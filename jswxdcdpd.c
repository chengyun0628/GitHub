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
segy tr;                        /* trace of input  						*/
int
main(int argc, char **argv)
{
 float **data;					/* angx_array of 2D real trace			*/
 float **shuchu;
 float *stk1;					/* angx_array of 1D real trace			*/
 float *stk2;					/* angx_array of 1D real trace			*/
 float *itint;
 float *itout;
 float *angx;
 float *tmp;
 float *Atmpcut;
 float atmp1;
 float atmp2;
 float aratio;
 float abg;
 float aed;
 float dt;						/* sample interval						*/
 float T0;
 float Tofd;
 float agtmp;
 float iasmp;
 int *tmpx;
 int fd;
 int hntr;
 int widlth,hwid;
 int itbg,ited;
 int iabg,iaed;
 int mnf;
 int nt;                		/* number of points on input trace      */
 int ntr;						/* number of trace input				*/
 int itr,inf,it,itx,ia; 		/* count number                         */
 int j;
 int k;
 float k1=0;
 int verbose;		            /* flag to get advisory messages	    */
 int nstksm;

/* file name */
 FILE *fp1;
 FILE *fp2;
 FILE *fp3;
 FILE *tracefp=NULL;	        /* angx file to hold traces             */
 FILE *hfp=NULL;		        /* angx file to hold trace headers      */
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
 if (!getparint("nstksm",&nstksm)) nstksm = 2;
 if (!getparint("widlth",&widlth)) widlth = 12;
 if (!getparint("fd",&fd)) fd =50;
 if (!getparint("itbg",&itbg)) itbg = 0.25/dt/fd+1;
 if (!getparfloat("aratio",&aratio)) aratio =0.2;
 if (!getparfloat("iasmp",&iasmp)) iasmp =tr.d2;
 hwid=widlth*0.5;
 if(itbg<hwid)	itbg=hwid;
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
 shuchu=ealloc2float(widlth,ntr);
 stk1=ealloc1float(nt);
 stk2=ealloc1float(nt);
 itint=ealloc1float(nt);
 itout=ealloc1float(nt);
 tmp=ealloc1float(nt);
 tmpx=ealloc1int(nt);
 angx=ealloc1float(nt);
 Atmpcut=ealloc1float(widlth);
/* Zero all arrays */
 memset((void *) data[0], 0,nt*ntr*FSIZE);
 memset((void *) tmp, 0,nt*FSIZE);
 memset((void *) tmpx, 0,nt*FSIZE);
 memset((void *) itint, 0,nt*FSIZE);
 memset((void *) itout, 0,nt*FSIZE);
 memset((void *) angx, 0,nt*FSIZE);

/* load traces into the zero-offset array and close tmpfile */
 efread(*data, FSIZE, nt*ntr, tracefp);
 efclose(tracefp);

/* stacking data along angle */
 for(itr=0;itr<ntr;itr++)
	for(it=1;it<nt;it++)
		stk1[it]+=data[itr][it];
 stk1[0]=0;
 for(it=0;it<nt;it++)
	{
	k1=k1>stk1[it]?k1:stk1[it];
 	stk2[it]=stk1[it];
	}
 for(it=1;it<nt;it++)
	if(stk1[it]<aratio*k1 || stk1[it-1]>stk1[it])
		stk2[it]=0;
 mnf=0;
 for(it=1;it<nt;it++)
	if(stk2[it-1]>stk2[it])
		{
		 tmpx[mnf]=it-1;
		 mnf++;
		}
 warn("mnf=%d",mnf);
 k=1;
 for(inf=0;inf<mnf;inf++)
	{
	 it=tmpx[inf];
	 itbg=it-hwid;
	 ited=it+hwid+1;
	 T0=it*dt;
	 Tofd=T0*fd*2;
	 memset((void *) shuchu[0], 0,widlth*ntr*FSIZE);
	 for(itr=0;itr<ntr;itr++)
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
	 	 iaed=iaed<hntr?iaed+hntr:ntr-1;
		 memset((void *) Atmpcut, 0,widlth*FSIZE);
		 for(itx=itbg;itx<ited;itx++)
		 	for(ia=iabg;ia<=iaed;ia++)
			 	Atmpcut[itx-itbg]+=data[ia][itx];
		 for(itx=0;itx<widlth;itx++)
			shuchu[itr][itx]=Atmpcut[itx];
		}
	atmp2=0;
	for(itr=0;itr<ntr;itr++)
		{
		 atmp1=0;
		for(itx=0;itx<widlth;itx++)
			atmp1+=fabs(shuchu[itr][itx]);
		if(atmp1>atmp2)
			{
			atmp2=atmp1;
			tmp[k]=itr;
			}
		}
	itint[k]=it*dt;
	warn("tmpk=%f,it=%d",(tmp[k]-hntr)/iasmp,it);
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
 //cmonot(k,itint,tmp,ydin);
 //csplin(k,itint,tmp,ydin);
 //chermite(k,itint,tmp,ydin);
 //intcub(0,k,itint,ydin,nt,itout,angx);
 intlin(k, itint, tmp, hntr, hntr, nt, itout, angx);
 fp1=fopen("angx.bin","wb");
 fwrite(angx,sizeof(float),nt,fp1);
 fclose(fp1);
 memset((void *) stk1, 0,nt*FSIZE);
 memset((void *) stk2, 0,nt*FSIZE);
 /*for(it=0;it<nt;it++)
	{
	 for(itr=angx[it];itr<ntr;itr++)		
		{
		 atmp1=stk1[it];
		 stk1[it]+=data[itr][it];
		 if(fabs(atmp1)>=fabs(stk1[it])) //if(fabs(data[itr][it]/stk1[it])<0.03)	break;	
			{
			 itint[it]=itr;
			 break;
			}
		}
		
	 for(itr=angx[it];itr>=0;itr--)
		{
		 atmp1=stk2[it];
		 stk2[it]+=data[itr][it];
		 if(fabs(atmp1)>=fabs(stk2[it]))
			{
			 itout[it]=itr;
			 break;
			}
		}
	 stk1[it]=stk1[it]+stk2[it];
	}*/
 for(it=0;it<nt;it++)
	{
	 atmp2=0;
	 for(itr=angx[it];itr<ntr;itr++)		
		{
		 atmp1=0;
		 for(itx=it-nstksm;itx<=it+nstksm;itx++)
			{
			 if(itx<0)	continue;
			 else if(itx>nt-1)	break;
			 atmp1+=fabs(data[itr][itx]);
			}
		 atmp2+=atmp1;
		 if(atmp1/atmp2<0.03)
			{	
			 itint[it]=itr;
		 	 break;
			}
		}
	stk1[it]+=atmp2;
	atmp2=0;
	 for(itr=angx[it];itr>0;itr--)
		{
		 atmp1=0;
		 for(itx=it-nstksm;itx<=it+nstksm;itx++)
			{
			 if(itx<0)	continue;
			 else if(itx>nt-1)	break;
			 atmp1+=fabs(data[itr][itx]);
			}
		 atmp2+=atmp1;
		 if(atmp1/atmp2<0.03)
			{	
			 itout[it]=itr;
		 	 break;
			}
		}
	 stk1[it]+=atmp2;
	}
 fp2=fopen("frenell.bin","wb");
 fwrite(itint,sizeof(float),nt,fp2);
 fclose(fp2);
 fp3=fopen("frenelr.bin","wb");
 fwrite(itout,sizeof(float),nt,fp3);
 fclose(fp3);
 /*fp2=fopen("stkwxd1.bin","wb");
 fwrite(stk1,sizeof(float),nt,fp2);
 fclose(fp2);*/
/*free array*/
efclose(hfp);
free2float(data);
free2float(shuchu);
free1float(stk1);
free1float(stk2);
free1float(itint);
free1float(itout);
free1float(tmp);
free1float(Atmpcut);
free1int(tmpx);
free1float(angx);
return(CWP_Exit());	
}

