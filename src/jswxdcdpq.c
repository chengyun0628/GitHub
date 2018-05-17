#include "su.h"
#include "segy.h"
#include "header.h"
#include <signal.h>

/*********************** self documentation **********************/
char *sdoc[] = {NULL};

 /* Trace header fields accessed:
 * Trace header fields modified: 
 */
void wxd(int nt,int nx,int fd,int hwid, int hntr,int nstksm,float dt,float iasmp,float avg,float aratio,float *ang,float *stk1,float **datact,float *minkjj,float *maxkjj);
/**************** end self doc ***********************************/
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
 float *minkjj;
 float *maxkjj;
 float *stk1;
 float *ang;
 float avg=0;
 float dt;						/* sample interval						*/
 float dx;
 float iasmp;
 float aratio;
 int mincdp;
 int widlth;
 int hwid;
 int fd;
 int maxag;
 int hntr;
 int nt;                		/* number of points on input trace      */
 int ntr;						/* number of trace input				*/
 int sntr;
 int hsntr;
 int nx;
 int itr,it,ix; 				/* count number                         */
 int verbose;		            /* flag to get advisory messages	    */
 int ncdp;
 int nstksm;

/* file name */
 FILE *tracefp=NULL;	        /* angx file to hold traces             */
 FILE *hfp=NULL;		        /* angx file to hold trace headers      */
 FILE *fp1;
 FILE *fp2;
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
 if (!getparint("nstksm",&nstksm)) nstksm=2; 
 if (!getparfloat("dx",&dx)) 
	{
	 dx=20.;
	 warn("dx not set,assumed to be 20.(m)");
	}
 if (!getparint("maxag",&maxag)) maxag=-tr.f2;
 if (!getparfloat("iasmp",&iasmp)) iasmp=tr.d2;
 if (!getparint("widlth",&widlth)) widlth = 19;
 if (!getparint("fd",&fd)) fd =50;
 if (!getparfloat("aratio",&aratio)) aratio =0.2;
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
 hntr =(nx-1)*0.5;
/* rewind trace file pointer and header file pointer */
 erewind(tracefp);
 erewind(hfp);
 if (!getparint("ncdp", &ncdp))	ncdp=sntr;
/* Allocate space */
 data=ealloc2float(nt,ntr);
 datao=ealloc2float(nt,sntr);
 minkj=ealloc2float(nt,ncdp);
 maxkj=ealloc2float(nt,ncdp);
 datact=ealloc2float(nt,nx);
 minkjj=ealloc1float(nt);
 maxkjj=ealloc1float(nt);
 stk1=ealloc1float(nt);
 ang=ealloc1float(nt);

/* Zero all arrays */
 memset((void *) data[0], 0,nt*ntr*FSIZE);
 memset((void *) datao[0], 0,nt*sntr*FSIZE);
 memset((void *) minkj[0], 0,nt*ncdp*FSIZE);
 memset((void *) maxkj[0], 0,nt*ncdp*FSIZE);
 memset((void *) datact[0], 0,nt*nx*FSIZE);
 memset((void *) ang, 0,nt*FSIZE);
/* load traces into the zero-offset array and close tmpfile */
 efread(*data, FSIZE, nt*ntr, tracefp);
 efclose(tracefp);
 for(it=0;it<nt;it++)
	for(ix=0;ix<nx;ix++)
		stk1[it]+=data[hsntr*nx+ix][it];
/* Get the maxmum of every trace */
 for(it=1;it<nt;it++)
 	avg=avg>stk1[it]?avg:stk1[it];
/* loop */
 for(itr=hsntr;itr<sntr;itr++)
	{
	 for(ix=0;ix<nx;ix++)
		for(it=0;it<nt;it++)
			datact[ix][it]=data[itr*nx+ix][it];
	 memset((void *) stk1, 0,nt*FSIZE);
	 memset((void *) minkjj, 0,nt*FSIZE);
	 memset((void *) maxkjj, 0,nt*FSIZE);
	 wxd(nt,nx,fd,hwid,hntr,nstksm,dt,iasmp,avg,aratio,ang,stk1,datact,minkjj,maxkjj);
	 for(it=0;it<nt;it++)
	 	{
		 datao[itr][it]=stk1[it];
		 minkj[itr][it]=minkjj[it];
		 maxkj[itr][it]=maxkjj[it];
		}
	}
 for(itr=hsntr;itr>0;itr--)
	{
	 for(ix=0;ix<nx;ix++)
		for(it=0;it<nt;it++)
			datact[ix][it]=data[itr*nx+ix][it];
	 memset((void *) stk1, 0,nt*FSIZE);
	 memset((void *) minkjj, 0,nt*FSIZE);
	 memset((void *) maxkjj, 0,nt*FSIZE);
	 wxd(nt,nx,fd,hwid,hntr,nstksm,dt,iasmp,avg,aratio,ang,stk1,datact,minkjj,maxkjj);
	 for(it=0;it<nt;it++)
		{
	 	datao[itr][it]=stk1[it];
		minkj[itr][it]=minkjj[it];
		maxkj[itr][it]=maxkjj[it];
		}
	}
 fp1=fopen("kjmin","wb");
 fwrite(minkj[0],sizeof(float),nt*ncdp,fp1);
 fclose(fp1);
 fp2=fopen("kjmax","wb");
 fwrite(maxkj[0],sizeof(float),nt*ncdp,fp2);
 fclose(fp2);
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
 free1float(ang);
 free1float(stk1);
 return(CWP_Exit());	
}
void wxd(int nt,int nx,int fd,int hwid, int hntr,int nstksm,float dt,float iasmp,float avg,float aratio,float *ang,float *stk1,float **datact,float *minkjj,float *maxkjj)
{
 float *stk2;				
 float *itint;
 float *itout;
 float *angx;
 float *tmp;
 float *tmpx;
 float *Atmpcut;
 float atmp1;
 float atmp2;
 float atmp3;
 float abg;
 float aed;
 float T0;
 float dttmp;
 float Tofd;
 float agtmp;
 int widlthx;
 int ittmp;
 int itbg;
 int ited;
 int mnf;
 int iabg,iaed;
 int itr,inf,it,itx,ia; 			/* count number */
 int j;
 int k;
 widlthx=hwid*2+1;
 stk2=ealloc1float(nt);
 itint=ealloc1float(nt);
 itout=ealloc1float(nt);
 tmp=ealloc1float(nt);
 tmpx=ealloc1float(nt);
 angx=ealloc1float(nt);
 Atmpcut=ealloc1float(widlthx);
 memset((void *) tmp, 0,nt*FSIZE);
 memset((void *) tmpx, 0,nt*FSIZE);
 memset((void *) itint, 0,nt*FSIZE);
 memset((void *) itout, 0,nt*FSIZE);
 memset((void *) angx, 0,nt*FSIZE);

/* stacking data along angle */
for(itr=0;itr<nx;itr++)
	for(it=0;it<nt;it++)
		stk1[it]+=datact[itr][it];
 stk1[0]=0;
 for(it=0;it<nt;it++)
 	stk2[it]=stk1[it];
 for(it=1;it<nt;it++)
	if(stk1[it]<aratio*avg || stk1[it-1]>stk1[it])
		stk2[it]=0;
 mnf=0;
/*------------------------------*/
 for(it=1;it<nt;it++)
	if(stk2[it-1]>stk2[it])
		{
		 tmpx[mnf]=it-1;
		 mnf++;
		}
 if(mnf>=4)
	{
 	k=1;
 	for(inf=0;inf<mnf;inf++)
		{
	 	it=tmpx[inf];
	 	itbg=it-hwid;
	 	ited=it+hwid+1;
	 	T0=it*dt;
	 	Tofd=T0*fd*2;
	 	atmp2=0;
	 	for(itr=0;itr<nx;itr++)
			{
		 	agtmp=(itr-hntr)*0.017453/iasmp;
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
	 	 	iabg=iabg>-hntr?iabg+hntr:0;
	 	 	iaed=iaed<hntr?iaed+hntr:nx-1;
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
				tmp[k]=itr;
				}
			}
		itint[k]=it*dt;
		//warn("tmpk=%f",(tmp[k]-hntr)/iasmp);
		k++;
		}
 	itint[0]=0;
 	itint[k]=nt*dt;
 	tmp[k]=tmp[0]=hntr;
	k++;
 	float ydin[k][4];
 	memset((void *) ydin[0], 0, k*4*FSIZE);
 	for(it=0;it<nt;it++)
    	itout[it]=it*dt;
 	cakima(k,itint,tmp,ydin);
 	intcub(0,k,itint,ydin,nt,itout,angx);
 	//intlin(k, itint, tmp, tmp[0], tmp[k-1], nt, itout, angx);
	/*--------------------------------------*/
	for(inf=0;inf<mnf;inf++)
		{
	 	it=tmpx[inf];
		for(itx=it-hwid;itx<it+hwid;itx++)
			angx[itx]=tmp[inf+1];
		}
	/*-------------------------------------*/
	for(it=0;it<nt;it++)
 	 	ang[it]=angx[it];
	}
else
	{
	 for(it=0;it<nt;it++)
 	 	angx[it]=ang[it];
	}	
/* stack along Frenel */
 memset((void *) stk1, 0,nt*FSIZE);
 /*for(it=0;it<nt;it++)
	{
	 atmp2=0;
	 for(itr=angx[it];itr<nx;itr++)		
		{
		 atmp1=0;
		 for(itx=it-nstksm;itx<=it+nstksm;itx++)
			{
			 if(itx<0)	continue;
			 else if(itx>nt-1)	break;
			 atmp1+=fabs(datact[itr][itx]);
			}
		 atmp2+=atmp1;
		 if(atmp1/atmp2<0.03)
			{
			 maxkjj[it]=itr;
		 	 break;
			}
		}
	atmp3=0;
	 for(itr=angx[it];itr>0;itr--)
		{
		 atmp1=0;
		 for(itx=it-nstksm;itx<=it+nstksm;itx++)
			{
			 if(itx<0)	continue;
			 else if(itx>nt-1)	break;
			 atmp1+=fabs(datact[itr][itx]);
			}
		 atmp3+=atmp1;
		 if(atmp1/atmp3<0.03)
			{	
			 minkjj[it]=itr;
		 	 break;
			}
		}
	}
 for(it=0;it<nt;it++)
	{
	 for(itr=minkjj[it];itr<=maxkjj[it];itr++)
		stk1[it]+=datact[itr][it];
	 minkjj[it]=(minkjj[it]-hntr)/iasmp;
	 maxkjj[it]=(maxkjj[it]-hntr)/iasmp;
	}*/

 for(it=0;it<nt;it++)
	{
	 atmp1=0;
	 for(itr=angx[it];itr<nx;itr++)		
		{
		 atmp1+=datact[itr][it];
		 if(fabs(datact[itr][it]/atmp1)<0.03)
			{	
			 maxkjj[it]=(itr-hntr)/iasmp;
		 	 break;
			}
		}
	 atmp2=0;
	 for(itr=angx[it];itr>0;itr--)
		{
		 atmp2+=datact[itr][it];
		 if(fabs(datact[itr][it]/atmp2)<0.03)
			{	
			 minkjj[it]=(itr-hntr)/iasmp;
		 	 break;
			}
		}
	 stk1[it]=atmp1+atmp2;
	}

 /*for(it=0;it<nt;it++)
	{
	 atmp2=0;
	 for(itr=angx[it];itr<nx;itr++)		
		{
		 atmp1=atmp2;
		 atmp2+=datact[itr][it];
		 if(fabs(atmp1)>=fabs(atmp2))
			{	
			 maxkjj[it]=(itr-hntr)/iasmp;
		 	 break;
			}
		}
	 atmp3=0;
	 for(itr=angx[it];itr>0;itr--)
		{
		 atmp1=atmp3;
		 atmp3+=datact[itr][it];
		 if(fabs(atmp1)>=fabs(atmp3))
			{	
			 minkjj[it]=(itr-hntr)/iasmp;
		 	 break;
			}
		}
	 stk1[it]=atmp2+atmp3;
	}*/

 /*for(it=0;it<nt;it++)
	{
	 T0=it*dt;
	 Tofd=T0*fd*2;
	 if(Tofd<=1)	continue;
	 agtmp=(angx[it]-hntr)*0.017453/iasmp;
	 abg=(Tofd*sin(agtmp)-sqrt(Tofd+Tofd-1))/cos(agtmp)/(Tofd-1);
	 aed=(Tofd*sin(agtmp)+sqrt(Tofd+Tofd-1))/cos(agtmp)/(Tofd-1);
	 minkjj[it]=atan(abg)*57.29578;
	 maxkjj[it]=atan(aed)*57.29578;
	}
 for(it=0;it<nt;it++)
	{
	atmp1=minkjj[it]*iasmp+hntr;
	atmp2=maxkjj[it]*iasmp+hntr;
	atmp1=atmp1>0?atmp1:0;
	atmp2=atmp1<nx?atmp1:nx-1;
	for(itr=atmp1;itr<=atmp2;itr++)		
		stk1[it]+=datact[itr][it];
	}*/

/*free array*/
free1float(stk2);
free1float(itint);
free1float(itout);
free1float(tmp);
free1float(tmpx);
free1float(Atmpcut);
free1float(angx);
}
