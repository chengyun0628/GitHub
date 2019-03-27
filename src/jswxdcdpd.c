#include "su.h"
#include "segy.h"
#include "header.h"
#include <signal.h>

/*********************** self documentation **********************/
char *sdoc[] = {NULL};

 /* Trace header fields accessed:
 * Trace header fields modified: 
 */
/**************** end self doc ***********************************/
void wxd(int nt,int hnx,int fd,int hwid,int bts,int ets,int smpitnum,float dt,float iasmp,float **data,float *angx,float *minkj,float *maxkj);
segy tr;                        /* trace of input  						*/
segy tro;                       /* trace of output  					*/
int
main(int argc, char **argv)
{
 float **data;					/* ang_array of 2D real trace			*/
 float *datao;					/* ang_array of 1D real trace			*/
 float *angx;					/* array for scanned Dip angle			*/
 float *minkj;					/* array for the min aperture			*/
 float *maxkj;					/* array for the max aperture			*/
 float dt;						/* sample interval in time				*/
 float dx;						/* sample interval in space				*/
 float iasmp;					/* sample interval in angle				*/
 float kjl,kjr;					
 int mincdp;					/* cdp information on headers			*/
 int nt;                		/* number of points on input trace      */
 int ntr;						/* number of trace input				*/	
 int hntr;						/* half of the trace number				*/
 int bts,ets;					/* begain time and end time for scanning,(ms)			*/
 int smpitnum;					/* sample point number to express aperture. 3,5,7...31. */
 int widlth;					/* scanning window length in time		*/
 int hwid;						/* half of the scanning window length	*/
 int fd;						/* the dominant frequency of the data	*/
 int maxag;						/* the max Dip angle of angle gather	*/
 int it,ix; 					/* count number                         */
 int verbose;		            /* flag to get advisory messages	    */

/* file name */
 FILE *fp1;
 FILE *fp2;
 FILE *fp3;
 FILE *tracefp=NULL;	        /*  file to hold traces             	*/
 FILE *hfp=NULL;		        /* file to hold trace headers     		*/
 cwp_Bool seismic;	            /* is this seismic data?		        */
 cwp_Bool check_cdp=cwp_false;

/* Initialize */
 initargs(argc, argv);
 requestdoc(0);
 if (!getparint("verbose", &verbose))	verbose=0;

/* Get info from first trace */ 
 if (!gettr(&tr))  err("can't get first trace");
 nt = tr.ns;
 seismic = ISSEISMIC(tr.trid);		
 if (seismic) 
	{
	 if (verbose)	warn("input is seismic data, trid=%d",tr.trid);
	 dt = ((double) tr.dt)/1000000.0;
	}
 else 
	{
	 if (verbose)	warn("input is not seismic data, trid=%d",tr.trid);
	 dt = tr.d1;
	}
 if (!dt) 
	{
	 dt = .004;
	 if (verbose)	warn("dt or d1 not set, assumed to be .004");
	}

/* get optional parameters */
 if (!getparint("mincdp",&mincdp)) mincdp=tr.cdp; 
 if (!getparfloat("dx",&dx)) 
	{
	 dx=1.;
	 warn("dx not set,assumed to be 1.(m)");
	}
 if (!getparint("bts",&bts)) bts=0;
 if (!getparint("ets",&ets)) ets=nt*dt*1000;
 if (!getparint("smpitnum",&smpitnum)) smpitnum=13;
 if (!getparint("widlth",&widlth)) widlth = 7;
 if (!getparint("maxag",&maxag)) maxag=-tr.f2;
 if (!getparfloat("iasmp",&iasmp)) iasmp=tr.d2;
 if (!getparint("fd",&fd)) fd =20;
 bts=bts/(dt*1000);
 ets=ets/(dt*1000);
 smpitnum=1+0.5*(31-smpitnum);
 hwid=(widlth-1)*0.5;
 iasmp=1/iasmp;
/* Store traces in tmpfile while getting a count of number of traces */
 tracefp = etmpfile();
 hfp = etmpfile();
 ntr = 0;
 do 
	{
	 ++ntr;
	 /* read headers and data */
	 efwrite(&tr,HDRBYTES, 1, hfp);
	 efwrite(tr.data, FSIZE, nt, tracefp);
	} while (gettr(&tr));
 warn("ntr=%d",ntr);
 hntr =(ntr-1)*0.5;
/* rewind trace file pointer and header file pointer */
 erewind(tracefp);
 erewind(hfp);
/* Allocate space */
 data=ealloc2float(nt,ntr);
 datao=ealloc1float(nt);
 angx=ealloc1float(nt);
 minkj=ealloc1float(nt);
 maxkj=ealloc1float(nt);

/* Zero all arrays */
 memset((void *) data[0], 0,nt*ntr*FSIZE);
 memset((void *) datao, 0,nt*FSIZE);

/* load traces into the zero-offset array and close tmpfile */
 efread(*data, FSIZE, nt*ntr, tracefp);
 efclose(tracefp);
/* calculate the Fresnel zone*/
 wxd(nt,hntr,fd,hwid,bts,ets,smpitnum,dt,iasmp,data,angx,minkj,maxkj);
	for(it=0;it<nt;it++)
		{
		if(minkj[it]<-maxag)	minkj[it]=-maxag;
		if(maxkj[it]>maxag)	maxkj[it]=maxag;
		if(minkj[it]>maxkj[it])	minkj[it]=maxkj[it];
		kjl=minkj[it]*iasmp+hntr;
		kjr=maxkj[it]*iasmp+hntr;
		for(ix=kjl;ix<=kjr;ix++)
			datao[it]+=data[ix][it];
		}
/* output angle and aperture file in .txt */
 fp1=fopen("kjmin.txt","w");
 for(it=0;it<nt;it++)
	{
	fprintf(fp1,"%f ", it*dt);
	fprintf(fp1,"%f\n",minkj[it]);
	}
 fclose(fp1);
fp2=fopen("angx.txt","w");
for(it=0;it<nt;it++)
	{
	fprintf(fp2,"%f ", it*dt);
	fprintf(fp2,"%f\n",angx[it]);
	}
 fclose(fp2);
 fp3=fopen("kjmax.txt","w");
 for(it=0;it<nt;it++)
	{
	fprintf(fp3,"%f ", it*dt);
	fprintf(fp3,"%f\n",maxkj[it]);
	}
 fclose(fp3);
/* set header and output trace*/
 memset ((void *) &tro, (int) '\0', sizeof (tro));
 tro.trid = 1;
 tro.counit = 1;
 tro.f2=mincdp*dx;
 tro.d2 = dx;
 tro.ns = nt;
 tro.dt = dt*1000000;
 tro.cdp=mincdp;
 memcpy ((void *)tro.data, (const void *) datao, sizeof(float)*nt);
 puttr(&tro);
 /*free array*/
 free2float(data);
 free1float(datao);
 free1float(angx);
 free1float(minkj);
 free1float(maxkj);
 return(CWP_Exit());	
}
void wxd(int nt,int hnx,int fd,int hwid,int bts,int ets,int smpitnum,float dt,float iasmp,float **data,float *angx,float *minkj,float *maxkj)
{
 float *stk;				
 float *itint;
 float *itout;
 float *tmpx;
 float *tmpa;
 float *tmpl;
 float *tmpr;
 float *Atmpcut;
 float atmp0;
 float atmp1;
 float atmp2;
 float atmp3;
 float abg;
 float aed;
 float Tofd;
 float agtmp;
 int cys[15];
 int cyse[31];
 int cysz[255];
 int nter=10;
 int fdwid=30;
 int nx;
 int widmv;
 int widlth;
 int mnf;
 int itbg,ited;
 int iabg,iaed;
 int AG1,AG2;
 int FD1;
 int ix,inf,it,itx,ia,ifd,iter; 			/* count number */
 int j,k;
 for(iter=0;iter<15;iter++)
	cys[iter]=3+2*iter;
 widlth=hwid+hwid+1;
 nx=hnx+hnx+1;
 stk=ealloc1float(nt);
 itint=ealloc1float(nt);
 itout=ealloc1float(nt);
 tmpx=ealloc1float(nt);
 tmpa=ealloc1float(nt);
 tmpl=ealloc1float(nt);
 tmpr=ealloc1float(nt);
 Atmpcut=ealloc1float(widlth);
 memset((void *) stk, 0,nt*FSIZE);
 memset((void *) itint, 0,nt*FSIZE);
 memset((void *) tmpx, 0,nt*FSIZE);
 memset((void *) tmpa, 0,nt*FSIZE);
 memset((void *) tmpl, 0,nt*FSIZE);
 memset((void *) tmpr, 0,nt*FSIZE);
/* stacking data along angle */
 for(it=0;it<nt;it++)
    itout[it]=it*dt;
 for(ix=0;ix<nx;ix++)
	for(it=0;it<nt;it++)
		stk[it]+=data[ix][it];
 inf=0;
 for(iter=0;iter<15;iter++)
	{
	widmv=(int)(nt/cys[iter]+1);
 	for(it=0;it<nt;it+=widmv)
		{
		cysz[inf]=it;
		atmp0=fabs(stk[it]);
		for(itx=it+1;itx<it+widmv;itx++)
			{
			if(itx==nt-1) break;
			if(atmp0<fabs(stk[itx]))
				{
				cysz[inf]=itx;
				atmp0=fabs(stk[itx]);
				}
			}
		inf++;
		}
	}
 for(iter=0;iter<31;iter++)
	{
	cyse[iter]=1;
 	for(it=0;it<224;it++)
		{
		if(cysz[it]==cysz[224+iter])	cyse[iter]+=1;
		}
	}
 mnf=0;
 for(iter=0;iter<31;iter++)
	if(cyse[iter]>=smpitnum)
		{
		tmpx[mnf]=cysz[224+iter];
		mnf++;
		}
 warn("mnf=%d",mnf);
 k=1;
 for(inf=0;inf<mnf;inf++)
	{
	 FD1=fd;
	 it=tmpx[inf];
	 if(it<bts) continue;
	 if(it>ets) break;
	 itbg=it-hwid;
	 ited=it+hwid+1;
	for(iter=0;iter<nter;iter++)
		{
		Tofd=it*dt*FD1*2;
		atmp2=0;
	 	for(ix=0;ix<nx;ix++)
			{
		 	agtmp=(ix-hnx)*0.017453/iasmp;
		 	abg=(Tofd*sin(agtmp)-sqrt(Tofd+Tofd-1))/cos(agtmp)/(Tofd-1);
		 	aed=(Tofd*sin(agtmp)+sqrt(Tofd+Tofd-1))/cos(agtmp)/(Tofd-1);
	 	 	iabg=(int)(atan(abg)*57.29578*iasmp);
	 	 	iaed=(int)(atan(aed)*57.29578*iasmp);
		 	if(iabg>iaed)
				{
			 	j=iabg;
			 	iabg=iaed;
			 	iaed=j;
				}
	 	 	iabg=iabg>-hnx?iabg+hnx:0;
	 	 	iaed=iaed<hnx?iaed+hnx:nx-1;
		 	memset((void *) Atmpcut, 0,widlth*FSIZE);
		 	for(itx=itbg;itx<ited;itx++)
		 		for(ia=iabg;ia<=iaed;ia++)
			 		Atmpcut[itx-itbg]+=data[ia][itx];
		 	atmp1=0;
		 	for(itx=0;itx<widlth;itx++)
				atmp1+=fabs(Atmpcut[itx]);
		 	if(atmp1>atmp2)
				{
				atmp2=atmp1;
				tmpa[k]=ix;
				}
			}
		AG1=tmpa[k];
		atmp2=0;
		for(ifd=fd-fdwid;ifd<=fd+fdwid;ifd++)
			{
			Tofd=it*dt*ifd*2;
			agtmp=(AG1-hnx)*0.017453/iasmp;
		 	abg=(Tofd*sin(agtmp)-sqrt(Tofd+Tofd-1))/cos(agtmp)/(Tofd-1);
		 	aed=(Tofd*sin(agtmp)+sqrt(Tofd+Tofd-1))/cos(agtmp)/(Tofd-1);
	 	 	iabg=(int)(atan(abg)*57.29578*iasmp);
	 	 	iaed=(int)(atan(aed)*57.29578*iasmp);
		 	if(iabg>iaed)
				{
			 	j=iabg;
			 	iabg=iaed;
			 	iaed=j;
				}
	 	 	iabg=iabg>-hnx?iabg+hnx:0;
	 	 	iaed=iaed<hnx?iaed+hnx:nx-1;
			memset((void *) Atmpcut, 0,widlth*FSIZE);
		 	for(itx=itbg;itx<ited;itx++)
		 		for(ia=iabg;ia<=iaed;ia++)
			 		Atmpcut[itx-itbg]+=data[ia][itx];
		 	atmp1=0;
		 	for(itx=0;itx<widlth;itx++)
				atmp1+=fabs(Atmpcut[itx]);
			if(atmp1>atmp2)
				{
				tmpl[k]=iabg;
				tmpr[k]=iaed;
				atmp2=atmp1;
				FD1=ifd;
				}
			}
		if(iter==0)	
			{
			AG2=AG1;
			continue;
			}
		if(fabs(AG1-AG2)<2)	
			break;
		AG2=AG1;
		}
	iabg=tmpl[k];
	iaed=tmpr[k];
	atmp0=0;
	for(ia=iabg;ia<=iaed;ia++)
		atmp0+=data[ia][it];
	atmp1=atmp0;
	atmp2=fabs(atmp1);
	if(iaed==(nx-1)) iaed=nx-2;
	if(iabg==0)	iabg=1;
	if(atmp2<fabs(atmp0+data[iaed+1][it]))
		for(ia=iaed+1;ia<nx;ia++)
			{
			atmp1+=data[ia][it];
			atmp3=fabs(atmp1);
			if(atmp2<atmp3)
				{
				atmp2=atmp3;
				tmpr[k]=ia;
				}
			else break;
			} 
	else 
		for(ia=iaed;ia>hnx;ia--)
			{
			atmp1-=data[ia][it];
			atmp3=fabs(atmp1);
			if(atmp2<atmp3)
				{
				atmp2=atmp3;
				tmpr[k]=ia;
				}
			else break;
			}
	atmp1=atmp0;
	atmp2=fabs(atmp1);
	if(atmp2<fabs(atmp0+data[iabg-1][it]))
		for(ia=iabg-1;ia>=0;ia--)
			{
			atmp1+=data[ia][it];
			atmp3=fabs(atmp1);
			if(atmp2<atmp3)
				{
				atmp2=atmp3;
				tmpl[k]=ia;
				}
			else break;
			} 
	else 
		for(ia=iabg;ia<hnx;ia++)
			{
			atmp1-=data[ia][it];
			atmp3=fabs(atmp1);
			if(atmp2<atmp3)
				{
				atmp2=atmp3;
				tmpl[k]=ia;
				}
			else break;
			} 
	itint[k]=it*dt;
	k++;
	}
 itint[k]=nt*dt;
 tmpa[0]=tmpa[1];
 tmpa[k]=tmpa[k-1];
 tmpl[0]=tmpl[1];
 tmpl[k]=tmpl[k-1];
 tmpr[0]=tmpr[1];
 tmpr[k]=tmpr[k-1];
 k++;
 //intlin(k, itint, tmpa, tmpa[0], tmpa[k-1], nt, itout, angx);
 //intlin(k, itint, tmpl, tmpl[0], tmpl[k-1], nt, itout, minkj);
 //intlin(k, itint, tmpr, tmpr[0], tmpr[k-1], nt, itout, maxkj);
 float ydin[k][4];
 memset((void *) ydin[0], 0, k*4*FSIZE);
 cmonot(k,itint,tmpa,ydin);
 intcub(0,k,itint,ydin,nt,itout,angx);

 memset((void *) ydin[0], 0, k*4*FSIZE);
 cmonot(k,itint,tmpl,ydin);
 intcub(0,k,itint,ydin,nt,itout,minkj);

 memset((void *) ydin[0], 0, k*4*FSIZE);
 cmonot(k,itint,tmpr,ydin);
 intcub(0,k,itint,ydin,nt,itout,maxkj);
 for(it=0;it<nt;it++)
	{
	angx[it]=(angx[it]-hnx)/iasmp;
	maxkj[it]=(maxkj[it]-hnx)/iasmp;
	minkj[it]=(minkj[it]-hnx)/iasmp;
	}

/*free array*/
 free1float(stk);
 free1float(itint);
 free1float(itout);
 free1float(tmpx);
 free1float(tmpa);
 free1float(tmpl);
 free1float(tmpr);
 free1float(Atmpcut);
}
