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
void wxd(int nt,int nx,int fd,int hwid, int hnx,int smpitnum,float dt,float iasmp,float aratio,float arati,float **datact,float *angxx,float *minkjj,float *maxkjj);
segy tr;                        /* trace of input  						*/
segy tro;                       /* trace of output  					*/
int
main(int argc, char **argv)
{
 float **data;					/* ang_array of 2D real trace			*/
 float **datact;				/* ang_array of 2D real trace			*/
 float **datao;					/* ang_array of 2D real trace			*/
 float **minkj;
 float **maxkj;
 float **angx;
 float *minkjj;
 float *maxkjj;
 float *angxx;
 float dt;						/* sample interval						*/
 float dx;
 float iasmp;
 float aratio;
 float arati;
 float kjl,kjr;
 int smpitnum;					/* sample point number to express aperture. 3,5,7...31.*/
 int mincdp;
 int widlth;
 int hwid;
 int fd;
 int maxag;
 int nt;                		/* number of points on input trace      */
 int ntr;						/* number of trace input				*/
 int nx;
 int hnx;						
 int sntr;
 int hsntr;
 int itr,it,ix; 				/* count number                         */
 int verbose;		            /* flag to get advisory messages	    */

/* file name */
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
	 dx=6.;
	 warn("dx not set,assumed to be 6.(m)");
	}
 if (!getparint("smpitnum",&smpitnum)) smpitnum=13;
 if (!getparint("maxag",&maxag)) maxag=-tr.f2;
 if (!getparfloat("iasmp",&iasmp)) iasmp=tr.d2;
 if (!getparint("widlth",&widlth)) widlth = 21;
 if (!getparint("fd",&fd)) fd =20;
 if (!getparfloat("arati",&arati)) arati =0.95;
 if (!getparfloat("aratio",&aratio)) aratio =0.2;
 smpitnum=1+0.5*(31-smpitnum);
 hwid=(widlth-1)*0.5;
 iasmp=1/iasmp;
 nx=2*maxag*iasmp+1;
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
 sntr=ntr/nx;
 hsntr=(int)(0.5*sntr);
 hnx =(nx-1)*0.5;
/* rewind trace file pointer and header file pointer */
 erewind(tracefp);
 erewind(hfp);
/* Allocate space */
 data=ealloc2float(nt,ntr);
 datao=ealloc2float(nt,sntr);
 angx=ealloc2float(nt,sntr);
 minkj=ealloc2float(nt,sntr);
 maxkj=ealloc2float(nt,sntr);
 datact=ealloc2float(nt,nx);
 angxx=ealloc1float(nt);
 minkjj=ealloc1float(nt);
 maxkjj=ealloc1float(nt);

/* Zero all arrays */
 memset((void *) data[0], 0,nt*ntr*FSIZE);
 memset((void *) datao[0], 0,nt*sntr*FSIZE);
 memset((void *) angx[0], 0,nt*sntr*FSIZE);
 memset((void *) minkj[0], 0,nt*sntr*FSIZE);
 memset((void *) maxkj[0], 0,nt*sntr*FSIZE);
 memset((void *) datact[0], 0,nt*nx*FSIZE);
/* load traces into the zero-offset array and close tmpfile */
 efread(*data, FSIZE, nt*ntr, tracefp);
 efclose(tracefp);
/* loop */
 for(itr=hsntr+1;itr<sntr;itr++)
	{
	 for(ix=0;ix<nx;ix++)
		for(it=0;it<nt;it++)
			datact[ix][it]=data[itr*nx+ix][it];
	 memset((void *) angxx, 0,nt*FSIZE);
	 wxd(nt,nx,fd,hwid,hnx,smpitnum,dt,iasmp,aratio,arati,datact,angxx,minkjj,maxkjj);
	 for(it=0;it<nt;it++)
	 	{
		 angx[itr][it]=angxx[it];
		 minkj[itr][it]=minkjj[it];
		 maxkj[itr][it]=maxkjj[it];
		}
	}
 for(itr=hsntr;itr>=0;itr--)
	{
	 for(ix=0;ix<nx;ix++)
		for(it=0;it<nt;it++)
			datact[ix][it]=data[itr*nx+ix][it];
	 memset((void *) angxx, 0,nt*FSIZE);
	 wxd(nt,nx,fd,hwid,hnx,smpitnum,dt,iasmp,aratio,arati,datact,angxx,minkjj,maxkjj);
	 for(it=0;it<nt;it++)
		{
		angx[itr][it]=angxx[it];
		minkj[itr][it]=minkjj[it];
		maxkj[itr][it]=maxkjj[it];
		}
	}
 for(itr=0;itr<sntr;itr++)
	for(it=0;it<nt;it++)
		{
		if(minkj[itr][it]<-maxag)	minkj[itr][it]=-maxag;
		if(maxkj[itr][it]>maxag)	maxkj[itr][it]=maxag;
		if(minkj[itr][it]>maxkj[itr][it])	minkj[itr][it]=maxkj[itr][it];
		kjl=minkj[itr][it]*iasmp+hnx;
		kjr=maxkj[itr][it]*iasmp+hnx;
		for(ix=kjl;ix<=kjr;ix++)
			datao[itr][it]+=data[itr*nx+ix][it];
		}

/* set header */
 memset ((void *) &tro, (int) '\0', sizeof (tro));
 tro.trid = 1;
 tro.counit = 1;
 tro.f2=mincdp*dx;
 tro.d2 = dx;
 tro.ns = nt;
 tro.dt = dt*1000000;
for(itr=0;itr<sntr;itr++)
	{
	 tro.cdp=mincdp;
	 memcpy ((void *)tro.data, (const void *) datao[itr], sizeof(float)*nt);
	 puttr(&tro);
	 mincdp++;
	}
 /*free array*/
 free2float(data);
 free2float(datact);
 free2float(datao);
 free2float(minkj);
 free2float(maxkj);
 free1float(minkjj);
 free1float(maxkjj);
 free1float(angxx);
 return(CWP_Exit());	
}
void wxd(int nt,int nx,int fd,int hwid, int hnx,int smpitnum,float dt,float iasmp,float aratio,float arati,float **datact,float *angxx,float *minkjj,float *maxkjj)
{
 float *stk1;
 float *stk2;				
 float *itint;
 float *itout;
 float *tmpa;
 float *tmpx;
 float *tmpl;
 float *tmpr;
 float *Atmpcut;
 float atmpl;
 float atmpr;
 float atmp1;
 float atmp2;
 float atmp3;
 float tmpv;
 float tmpvx;
 float tmpvl;
 float tmpvr;
 float abg;
 float aed;
 float Tofd;
 float agtmp;
 float avg;
 int cys[15];
 int cyse[31];
 int cysz[255];
 int nter=10;
 int fdwid=30;
 int widmv;
 int widlthx;
 int itbg;
 int ited;
 int mnf;
 int iabg,iaed;
 int itr,inf,it,itx,ia,ifd,iter,AG1,FD1,AG2; 			/* count number */
 int j,k;
 FILE *fp1;
 FILE *fp2;
 FILE *fp4;
 for(iter=0;iter<15;iter++)
	cys[iter]=3+2*iter;
 widlthx=hwid+hwid+1;
 stk1=ealloc1float(nt);
 stk2=ealloc1float(nt);
 itint=ealloc1float(nt);
 itout=ealloc1float(nt);
 tmpl=ealloc1float(nt);
 tmpr=ealloc1float(nt);
 tmpx=ealloc1float(nt);
 tmpa=ealloc1float(nt);
 Atmpcut=ealloc1float(widlthx);
 memset((void *) stk1, 0,nt*FSIZE);
 memset((void *) tmpl, 0,nt*FSIZE);
 memset((void *) tmpr, 0,nt*FSIZE);
 memset((void *) tmpx, 0,nt*FSIZE);
 memset((void *) itint, 0,nt*FSIZE);
 memset((void *) tmpa, 0,nt*FSIZE);
/* stacking data along angle */
 for(it=0;it<nt;it++)
    itout[it]=it*dt;
 for(itr=0;itr<nx;itr++)
	for(it=0;it<nt;it++)
		stk1[it]+=datact[itr][it];
 ifd=0;
 for(iter=0;iter<15;iter++)
	{
	widmv=(int)(nt/cys[iter]+1);
 	for(it=0;it<nt;it+=widmv)
		{
		cysz[ifd]=it;
		atmpl=fabs(stk1[it]);
		for(itx=it+1;itx<it+widmv;itx++)
			{
			if(itx==nt-1) break;
			if(atmpl<fabs(stk1[itx]))
				{
				cysz[ifd]=itx;
				atmpl=fabs(stk1[itx]);
				}
			}
		ifd++;
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
	 itbg=it-hwid;
	 ited=it+hwid+1;
	for(iter=0;iter<nter;iter++)
		{
		Tofd=it*dt*FD1*2;
		atmp2=0;
	 	for(itr=0;itr<nx;itr++)
			{
		 	agtmp=(itr-hnx)*0.017453/iasmp;
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
		 	memset((void *) Atmpcut, 0,widlthx*FSIZE);
		 	for(itx=itbg;itx<ited;itx++)
		 		for(ia=iabg;ia<=iaed;ia++)
			 		Atmpcut[itx-itbg]+=datact[ia][itx];
		 	atmp1=0;
		 	for(itx=0;itx<widlthx;itx++)
				atmp1+=fabs(Atmpcut[itx]);
		 	if(atmp1>atmp2)
				{
				atmp2=atmp1;
				tmpa[k]=itr;
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
			memset((void *) Atmpcut, 0,widlthx*FSIZE);
		 	for(itx=itbg;itx<ited;itx++)
		 		for(ia=iabg;ia<=iaed;ia++)
			 		Atmpcut[itx-itbg]+=datact[ia][itx];
		 	atmp1=0;
		 	for(itx=0;itx<widlthx;itx++)
				atmp1+=fabs(Atmpcut[itx]);
			if(atmp1>atmp2)
				{
				atmp2=atmp1;
				FD1=ifd;
				}
			}
		if(iter==0)	
			{
			AG2=AG1;
			continue;
			}
		if(AG1-AG2<1||iter==nter-1)	
			{
			tmpl[k]=iabg;
			tmpr[k]=iaed;
			for(ia=iabg;ia<=iaed;ia++)
				tmpv+=datact[ia][it];
			break;
			}
		AG2=AG1;
		}
	if(iaed==(nx-1)) iaed=nx-2;
	if(iabg==0)	iabg=1;
	tmpvx=tmpv;
	tmpvl=fabs(tmpvx);
	if(tmpvl<fabs(tmpv+datact[iaed+1][it]))
		for(ia=iaed+1;ia<nx;ia++)
			{
			tmpvx+=datact[ia][it];
			tmpvr=fabs(tmpvx);
			if(tmpvl<tmpvr)
				{
				tmpvl=tmpvr;
				tmpr[k]=ia;
				}
			else break;
			} 
	/*else 
		for(ia=iaed;ia>=0;ia--)
			{
			tmpvx-=datact[ia][it];
			tmpvr=fabs(tmpvx);
			if(tmpvl<tmpvr)
				{
				tmpvl=tmpvr;
				tmpr[k]=ia;
				}
			else break;
			} */
	tmpvx=tmpv;
	tmpvl=fabs(tmpvx);
	if(tmpvl<fabs(tmpv+datact[iabg-1][it]))
		for(ia=iabg-1;ia>=0;ia--)
			{
			tmpvx+=datact[ia][it];
			tmpvr=fabs(tmpvx);
			if(tmpvl<tmpvr)
				{
				tmpvl=tmpvr;
				tmpl[k]=ia;
				}
			else break;
			} 
	/*else 
		for(ia=iabg;ia<nx;ia++)
			{
			tmpvx-=datact[ia][it];
			tmpvr=fabs(tmpvx);
			if(tmpvl<tmpvr)
				{
				tmpvl=tmpvr;
				tmpl[k]=ia;
				}
			else break;
			} 

	/*atmpl=0;
	for(itr=tmpa[k];itr>=0;itr--)
		atmpl+=datact[itr][it];
	atmpl=fabs(atmpl);
	atmpr=0;
	for(itr=tmpa[k]+1;itr<nx;itr++)
		atmpr+=datact[itr][it];
	atmpr=fabs(atmpr);
	atmp1=0;
	for(itr=tmpa[k];itr>=0;itr--)
		{
		atmp2=atmp1;
		atmp1+=datact[itr][it];
		atmp3=fabs(atmp1);
		if(atmp3>=atmpl*arati && fabs(atmp2)>=atmp3)
			{	
			 tmpl[k]=(itr-hnx)/iasmp;
		 	 break;
			}
		else	tmpl[k]=-hnx/iasmp;
		}
	atmp2=0;
	for(itr=tmpa[k]+1;itr<nx;itr++)
		{
		atmp1=atmp2;
		atmp2+=datact[itr][it];
		atmp3=fabs(atmp2);
		if(atmp3>=atmpr*arati && fabs(atmp1)>=atmp3)
			{	
			tmpr[k]=(itr-hnx)/iasmp;
		 	 break; 
			}
		else	tmpr[k]=hnx/iasmp;
		}*/
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
//intlin(k, itint, tmpa, tmpa[0], tmpa[k-1], nt, itout, angxx);
 //intlin(k, itint, tmpl, tmpl[0], tmpl[k-1], nt, itout, minkjj);
 //intlin(k, itint, tmpr, tmpr[0], tmpr[k-1], nt, itout, maxkjj);
 float ydin[k][4];
memset((void *) ydin[0], 0, k*4*FSIZE);
cmonot(k,itint,tmpa,ydin);
intcub(0,k,itint,ydin,nt,itout,angxx);

 memset((void *) ydin[0], 0, k*4*FSIZE);
 cmonot(k,itint,tmpl,ydin);
 intcub(0,k,itint,ydin,nt,itout,minkjj);

memset((void *) ydin[0], 0, k*4*FSIZE);
 cmonot(k,itint,tmpr,ydin);
 intcub(0,k,itint,ydin,nt,itout,maxkjj);
for(it=0;it<nt;it++)
	{
	angxx[it]=(angxx[it]-hnx)/iasmp;
	maxkjj[it]=(maxkjj[it]-hnx)/iasmp;
	minkjj[it]=(minkjj[it]-hnx)/iasmp;
	}
fp1=fopen("kjmin.txt","w");
for(it=0;it<nt;it++)
	{
	fprintf(fp1,"%f ", it*dt);
	fprintf(fp1,"%f\n",minkjj[it]);
	}
 fclose(fp1);
fp4=fopen("angx.txt","w");
for(it=0;it<nt;it++)
	{
	fprintf(fp4,"%f ", it*dt);
	fprintf(fp4,"%f\n",angxx[it]);
	}
 fclose(fp4);
fp2=fopen("kjmax.txt","w");
for(it=0;it<nt;it++)
	{
	fprintf(fp2,"%f ", it*dt);
	fprintf(fp2,"%f\n",maxkjj[it]);
	}
 fclose(fp2);

/*free array*/
free1float(stk1);
free1float(stk2);
free1float(itint);
free1float(itout);
free1float(tmpl);
free1float(tmpr);
free1float(tmpx);
free1float(tmpa);
free1float(Atmpcut);
}
