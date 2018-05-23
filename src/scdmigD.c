#include"time.h"
#include "su.h"
#include "segy.h"
#include"header.h"

/*********************** self documentation **********************/
char *sdoc[] = {
" 									",
"                           		",
" 									",
" 					             	",
" 									",
" Required parameters:						",
"Refer to parfile 					",
" 									",
" Optional parameters:						",
"                                   ",
" 									",
" Notes:                                    ",
"Data must be sorted by offset 				",
" Caveats: 							    	",
" 									",
" Examples: 								",
" 									",
NULL};

/* Credit
 *
 * Note:
 *
 * Trace header fields accessed:
 * Trace header fields modified: 
 */
/**************** end self doc ***********************************/

static time_t t1,t2;

/* Prototype of function used internally */
void windtr(int nt,float *rtx,float ttt,float dt,float *firstt,float *datal);
void tanda(int ipx,int *anapx,int sx,int gx,float v,float vtt,float *ttt,float *qtmp);
void aperture(int it,int icdp,int mincdpout,int maxcdpout,int *bgc,int *edc,int anapxdx,float vt,float vtt,float h,float *angx1,float *angx2 );
void hammingFilter(int nf1,int nf2,int nf3,int nf4,int nf, float *filter);
#define LOOKFAC 2145
#define PFA_MAX 720720
segy tri;                        /* trace of input  						*/
segy tro;						/* trace of output 						*/
int
main(int argc, char **argv)
{
 register float *rt,*rtx;       /* real trace data and processed data	*/
 register complex *ct;          /* complex transformed trace            */
 register complex *hd;			/* half derivative						*/
 float *filter;              	/* filter array                         */
 float datal[8];				/* small cut of ununiform sample trace	*/ 
 float **data;					/* temp_array of 2D real trace			*/
 float **vel=NULL;              /* array for storing velocity           */
 float **kjmin=NULL;            /* array for storing the min aperture   */
 float **kjmax=NULL;            /* array for storing the max aperture   */
 float **mig=NULL;				/* array for storing migrated result	*/
 float ttt;						/* Travel time							*/
 float qtmp;					/* Amptitude of migration				*/	
 float va;    					/* ununiform sampled value				*/ 
 float firstt;					/* the first sample in datal[8]			*/
 float dt;						/* sample interval						*/
 float hdt;
 float T;						/* vertical time depth 					*/
 float p;						/* Boundary attenuation factor of amp	*/
 float df,dw;					/* freqency sample spacing				*/
 float tmax;					/* the max trace length					*/
 float v;						/* velocity								*/
 float vt;						/* vt=v*T								*/
 float vtt;						/* vtt=vt*vt							*/
 float ipxb;
 float f1,f2,f3,f4;				/* array of filter frequencies          */   
 int sx,gx;						/* coordinate of shot and geophone  	*/
 int offset,h;					/* offset and half of that				*/
 int *anapx;					/* coordinate of cdp for image space	*/
 int anapxmin;					/* min coordinate of imaging point		*/ 
 int anapxmax;					/* max coordinate of imaging point		*/
 int anapxdx;					/* spacing between imaging point		*/
 int startmt;					/* the start time for migration			*/
 int nstartmt;					/* nstartmt=startmt/dT					*/
 int endmt;						/* the end time for migration			*/
 int nendmt;					/* nendmt=endmt/D=dT					*/
 int tritvl;					/* trace interval						*/
 int itr,ix,ipx,it; 			/* count number                         */
 int icdp;						/* count number                         */
 int nf1,nf2,nf3,nf4;           /* nf1=(int)(f1/df)                     */
 int napmin;					/* napmin=(int)(anapxmin/anapxdx)		*/
 int oldcdp=0;	            	/* for temporary storage		        */
 int olddeltacdp=1;				/* for temporary storage                */
 int deltacdp;
 int bgc,edc;					/* begin and end of imaging	trace		*/		
 int mincdp,maxcdp;				/* mincdp and maxcdp of data input		*/
 int mincdpout;					/* the min cdp to output(image space)	*/
 int maxcdpout;					/* the min cdp to output(image space)	*/
 int mincdpo;					/* storing the mincdp for temporary to write header*/
 int firstcdp=0;	            /* first cdp in velocity file	    	*/
 int lastcdp=0;	                /* last cdp in velocity file	    	*/
 int ncdp;	                	/* number of cdps in the velocity file	*/ 
 int dcdp=0;	                /* number of cdps between consecutive traces */
 int nt;                		/* number of points on input trace      */
 int ntr;						/* number of trace input				*/
 int npx;						/* number of trace of imaging sapce		*/
 int nfft;                      /* number of points for fft trace       */
 int nf;                        /* number of frequencies (incl Nyq)     */
 int verbose;		            /* flag to get advisory messages	    */
 int nbjl,nbjr,nbj1,nbj2;
/* file name */
 char str[100],*path;
 char *vfile="";
 char *kjfile1="";
 char *kjfile2="";
 char *parfile="";
 FILE *fp;
 FILE *vfp=NULL;
 FILE *minbj=NULL;
 FILE *maxbj=NULL;
 FILE *tracefp=NULL;	        /* temp file to hold traces             */
 FILE *hfp=NULL;		        /* temp file to hold trace headers      */
 cwp_Bool seismic;	            /* is this seismic data?		        */
 cwp_Bool check_cdp=cwp_false;

/* Initialize */
 initargs(argc, argv);
 requestdoc(0);
 if (!getparint("verbose", &verbose))	verbose=0;

/* Get info from first trace */ 
 if (!gettr(&tri))  err("can't get first trace");
 nt = tri.ns;
 offset=tri.offset;
 seismic = ISSEISMIC(tri.trid);		
 if (seismic) 
	{
	 if (verbose)	warn("input is seismic data, trid=%d",tri.trid);
	 dt = ((double) tri.dt)/1000000.0;
	}
 else 
	{
	 if (verbose)	warn("input is not seismic data, trid=%d",tri.trid);
	 dt = tri.d1;
	}
 if (!dt) 
	{
	 dt = .004;
	 if (verbose)	warn("dt or d1 not set, assumed to be .004");
	}
 hdt=0.5*dt;
 tmax=(nt-1)*dt;
/* Get parameter*/
 if(!getparstring("path",&path))	err("path must be specified !");
 if(!getparstring("vfile",&vfile))	err("velocity file must be specified !");
 if(!getparstring("kjfile1",&kjfile1))	err("aperture file1 must be specified !");
 if(!getparstring("kjfile2",&kjfile2))	err("aperture file2 must be specified !");
 if(!getparstring("parfile",&parfile))	err("parameter file must be specified ! ");
 sprintf(str,"%s/par/%s",path,parfile);
 warn("str=%s",str);
 fp=fopen(str,"rb");
 fscanf(fp,"mincdp=%d\n",&mincdp);
 fscanf(fp,"maxcdp=%d\n",&maxcdp);
 fscanf(fp,"firstcdp=%d\n",&firstcdp);
 fscanf(fp,"lastcdp=%d\n",&lastcdp);
 fscanf(fp,"startmt=%d\n",&startmt);
 fscanf(fp,"endmt=%d\n",&endmt);
 fscanf(fp,"tritvl=%d\n",&tritvl);
 fscanf(fp,"f1=%f\n",&f1);
 fscanf(fp,"f2=%f\n",&f2);
 fscanf(fp,"f3=%f\n",&f3);
 fscanf(fp,"f4=%f\n",&f4);
 fscanf(fp,"anapxmin=%d\n",&anapxmin);
 fscanf(fp,"anapxmax=%d\n",&anapxmax);
 fscanf(fp,"anapxdx=%d\n",&anapxdx);
 efclose(fp);
 time(&t1);
 if (fmod(tritvl,anapxdx)!=0)	err("Check anapxdx value in parfile,cause that tritvl must be an integer number of anapxdx");

/*caculate image spacing*/
 napmin=(int)(anapxmin/anapxdx);
 npx=(int)(anapxmax/anapxdx)-napmin+1;
 anapx=ealloc1int(npx);
 anapx[0]=anapxmin;
 for(ipx=1;ipx<npx;ipx++)
 	anapx[ipx]=anapx[ipx-1]+anapxdx;
 mincdpout=mincdp-napmin;
 maxcdpout=maxcdp-napmin;
 nstartmt=ceil(0.001*startmt/dt);
 nendmt=floor(0.001*endmt/dt);
/* Store traces in tmpfile while getting a count of number of traces */
 tracefp = etmpfile();
 hfp = etmpfile();
 ntr = 0;
 do 
	{
	 ++ntr;
	 /* get new deltacdp value */
	 deltacdp=tri.cdp-oldcdp;

	 /* read headers and data */
	 efwrite(&tri,HDRBYTES, 1, hfp);
	 efwrite(tri.data, FSIZE, nt, tracefp);

	 /* error trappings. */
	 /* ...did cdp value interval change? */
	 if ((ntr>3) && (olddeltacdp!=deltacdp)) 
		{
		 if (verbose) 
			{
			 warn("cdp interval changed in data");	
			 warn("ntr=%d olddeltacdp=%d deltacdp=%d",ntr,olddeltacdp,deltacdp);
		 	 check_cdp=cwp_true;
			}
		}
		
	 /* save cdp and deltacdp values */
	 oldcdp=tri.cdp;
	 olddeltacdp=deltacdp;
	} while (gettr(&tri));
	warn("ntr=%d",ntr);
/* get last cdp  and dcdp */
 if (!getparint("dcdp",&dcdp))	dcdp=deltacdp;
 warn("dcp=%d",dcdp);
/* error trappings */
 if ( (firstcdp==lastcdp) 
	|| (dcdp==0)
	|| (check_cdp==cwp_true) )	warn("Check cdp values in data!");

/* rewind trace file pointer and header file pointer */
 erewind(tracefp);
 erewind(hfp);
/* total number of cdp's in data */
 ncdp=lastcdp-firstcdp+1;
 if(verbose)	warn("ncdp=%d",ncdp);

/* Set up FFT parameters */
 nfft=npfaro(nt,LOOKFAC*nt);
 if(nfft>=SU_NFLTS||nfft>=PFA_MAX)
    err("padded nt=%d--too big",nfft);
 nf=nfft/2+1;
 df=1.0/(nfft*dt);
 dw=2*PI*df;
  if (verbose)	warn("nf=%d,df=%f",nf,df);
/* Allocate space */
 rt=ealloc1float(nfft);
 rtx=ealloc1float(nfft);
 ct=ealloc1complex(nf);
 hd=ealloc1complex(nf);
 filter=ealloc1float(nf);
 data=ealloc2float(nt,ntr);
 vel=ealloc2float(nt,ncdp);
 kjmin=ealloc2float(nt,ncdp);
 kjmax=ealloc2float(nt,ncdp);
 mig=ealloc2float(nt,npx);

/* Zero all arrays */
 memset((void *) data[0], 0,nt*ntr*FSIZE);
 memset((void *) mig[0], 0,nt*npx*FSIZE);
 memset((void *) rt, 0, nfft*FSIZE);
 memset((void *) filter, 0, nf*FSIZE);
 
/* Read data from temporal array */
 for(itr=0;itr<ntr;++itr)
	{
	 efread(data[itr],FSIZE,nt,tracefp);
	}
	
/* read velocities */
 vfp=efopen(vfile,"r");
 efread(vel[0],FSIZE,nt*ncdp,vfp);
 efclose(vfp);
 minbj=efopen(kjfile1,"r");
 efread(kjmin[0],FSIZE,nt*ncdp,minbj);
 efclose(minbj);
 maxbj=efopen(kjfile2,"r");
 efread(kjmax[0],FSIZE,nt*ncdp,maxbj);
 efclose(maxbj);

/* Define half derivative*/
 for(ix=0;ix<nf;ix++)
	{
	 hd[ix].r=0.707107;
	 hd[ix].i=-0.707107;
	}

/* Get frequency parameter of filter */
 nf1=(int)(f1/df);
 nf2=(int)(f2/df);
 nf3=(int)(f3/df);
 nf4=(int)(f4/df);
 if (verbose)	warn("f1=%f,f2=%f,f3=%f,f4=%f,nf1=%d,nf2=%d,nf3=%d,nf4=%d",f1,f2,f3,f4,nf1,nf2,nf3,nf4);
 hammingFilter(nf1,nf2,nf3,nf4,nf,filter);
for(itr=0; itr<ntr; itr++)
	{
	 for(it=0; it<nt; ++it)
		 rt[it]=data[itr][it];

	/* zero array ct and rtx*/
	 memset((void *) rtx, 0, nfft*FSIZE);
	 memset((void *) ct, 0, nf*FSIZE);

	/* filtering the input trace and multiply by the half derivative*/
	 pfarc(1,nfft,rt,ct);
	 for(ix=0;ix<nf;ix++) 
		{	
		 if(ix>=nf1&&ix<nf4)	
			{
			ct[ix]=crmul(ct[ix],filter[ix-nf1]);
			ct[ix]=crmul(ct[ix],sqrt(ix*dw));
         	ct[ix]=cmul(ct[ix],hd[ix]);
			}
		 else	ct[ix].r=ct[ix].i=0.0;
		}
	 pfacr(-1,nfft,ct,rtx);
	 for(it=0;it<nt;it++)
		data[itr][it]=rtx[it]/nfft;
	}
/* Start the migration process */
/* Loop over input trace */
 warn("Starting migration process...\n");
 h=offset/2;
 for(ipx=mincdpout; ipx<=maxcdpout; ++ipx)
	{
	 T=nstartmt*hdt;
	 for(it=nstartmt;it<nendmt;it++)
		{
		 T+=hdt;
		 v=vel[ipx-firstcdp+napmin][it];
		 vt=v*T;
		 vtt=vt*vt;
		 aperture(it,ipx,mincdpout,maxcdpout,&bgc,&edc,anapxdx,vt,vtt,h,kjmin[ipx-firstcdp+napmin],kjmax[ipx-firstcdp+napmin]);
		 nbjl=(edc-bgc+1)*0.2;
		 nbj1=(bgc-nbjl)>mincdpout?nbjl:(bgc-mincdpout);
		 nbj2=(edc+nbjl)<maxcdpout?nbjl:(maxcdpout-edc);
		 if((bgc-nbj1-mincdpout)%dcdp!=0)	bgc=(bgc-nbj1-mincdpout)/dcdp+1;
		 else	bgc=(bgc-nbj1-mincdpout)/dcdp;
		 edc=(edc+nbj2-mincdpout)/dcdp;
		 nbj1=nbj1/dcdp;
		 nbj2=nbj2/dcdp;	
		 nbjl=edc-bgc-nbj1;
		 nbjr=edc-bgc-nbj2;
		 ipxb=(ipx-mincdpout)/dcdp;
		 ix=0;
		 /*------------------------------------------------------*/
		 if(ipxb>edc)		
			{
			sx=(mincdpout+napmin+(edc+1)*dcdp)*anapxdx-h;
		 	for(itr=edc;itr>=bgc;itr--) 
				{
				sx=sx-dcdp*anapxdx;
		 	 	gx=sx+offset;
     			tanda(ipx,anapx,sx,gx,v,vtt,&ttt,&qtmp);
             	if(ttt>=tmax)   break;
             	windtr(nt,data[itr],ttt,dt,&firstt,datal);
         	 	ints8r(8, dt, firstt, datal,
             	  	0.0, 0.0, 1, &ttt, &va);
			 	if(ix<nbj2)	p=sin(1.570796*ix/nbj2);
			 	else if(ix>nbjl)	p=cos(1.570796*(ix-nbjl)/nbj1);
			 	else p=1.0;
			 	ix++;
			 	mig[ipx][it]+=va*qtmp*p;
				}
			}
		 else if(ipxb<bgc)
			{
			sx=(mincdpout+napmin+(bgc-1)*dcdp)*anapxdx-h;
			for(itr=bgc;itr<=edc;itr++) 
				{
				sx=sx+dcdp*anapxdx;
		 	 	gx=sx+offset;
     			tanda(ipx,anapx,sx,gx,v,vtt,&ttt,&qtmp);
             	if(ttt>=tmax)   break;
             	windtr(nt,data[itr],ttt,dt,&firstt,datal);
         	 	ints8r(8, dt, firstt, datal,
             	  	0.0, 0.0, 1, &ttt, &va);
			 	if(ix<nbj1)	p=sin(1.570796*ix/nbj1);
			 	else if(ix>nbjr)	p=cos(1.570796*(ix-nbjr)/nbj2);
			 	else p=1.0;
			 	ix++;
			 	mig[ipx][it]+=va*qtmp*p;
				}
			}
		 else
			{
			sx=(mincdpout+napmin+(bgc-1)*dcdp)*anapxdx-h;
			for(itr=bgc;itr<=edc;itr++) 
				{
				sx=sx+dcdp*anapxdx;
		 	 	gx=sx+offset;
     			tanda(ipx,anapx,sx,gx,v,vtt,&ttt,&qtmp);
             	if(ttt>=tmax)   continue;
             	windtr(nt,data[itr],ttt,dt,&firstt,datal);
         	 	ints8r(8, dt, firstt, datal,
             	  	0.0, 0.0, 1, &ttt, &va);
			 	if(ix<nbj1)	p=sin(1.570796*ix/nbj1);
			 	else if(ix>nbjr)	p=cos(1.570796*(ix-nbjr)/nbj2);
			 	else p=1.0;
			 	ix++;
			 	mig[ipx][it]+=va*qtmp*p;
				}
			}
/*-----------------------------------------------------*/
		 /*sx=(mincdpout+napmin+(bgc-1)*dcdp)*anapxdx-h;
		 for(itr=bgc;itr<=edc;itr++)
			{
			 if(ix<nbj1)	p=sin(1.570796*ix/nbj1);
			 else if(ix>nbjr)	p=cos(1.570796*(ix-nbjr)/nbj2);
			 else p=1.0;
			 ix++;
			 sx=sx+dcdp*anapxdx;
		 	 gx=sx+offset;
			 tanda(ipx,anapx,sx,gx,v,vtt,&ttt,&qtmp);
             if(ttt>=tmax)   continue;
             windtr(nt,data[itr],ttt,dt,&firstt,datal);
         	 ints8r(8, dt, firstt, datal,
             	  0.0, 0.0, 1, &ttt, &va);
			 mig[ipx][it]+=va*qtmp*p;
			}*/
		}
	}
	memset ((void *) &tro, (int) '\0', sizeof (tro));
 	tro.trid = 1;
 	tro.counit = 1;
	tro.f2=mincdp*anapxdx;
	tro.d2 = anapxdx;
    tro.ns = nt;
    tro.dt = dt*1000000;
	/* Output migrated data */
	mincdpo=mincdp;
 	for(ipx=mincdpout; ipx<=maxcdpout; ++ipx)
    	{
     	tro.cdp = mincdpo;
     	memcpy ((void *) tro.data, (const void *)  mig[ipx],sizeof (float) * nt);
     	puttr(&tro);
     	mincdpo++;
    	}
 time(&t2);
 warn("time consuming in second = %f\n",difftime(t2,t1));
/*free array*/
 efclose(hfp);
 efclose(tracefp);
 free1float(rt);
 free1float(rtx);
 free1int(anapx);
 free2float(vel);
 free2float(kjmin);
 free2float(kjmax);
 free1float(filter);
 free1complex(ct);
 free1complex(hd);
 free2float(data);
 free2float(mig);
 return(CWP_Exit());
}

void hammingFilter(int nf1,int nf2,int nf3,int nf4,int nf, float *filter)
{
 int i;
 for(i=0;i<nf;i++)
	{
	 float tmp=0.0,tmpp=0.0;
	 if(i<(nf2-nf1))
		{
		 tmpp=PI*i/(nf2-nf1);
		 tmp=0.54-0.46*cos(tmpp);
		 filter[i]=tmp;
		}
	 else if(i>=(nf2-nf1)&&i<=(nf3-nf1))
		filter[i]=1.0;
	 else if(i>(nf3-nf1)&&i<=(nf4-nf1))
		{
		 tmpp=PI*(i+nf1+nf4-2*nf3+1)/(nf4-nf3+1);
		 tmp=0.54-0.46*cos(tmpp);
		 filter[i]=tmp;
		}
	 else filter[i]=0;
	}
}
        
void aperture(int it,int icdp,int mincdpout,int maxcdpout,int *bgc,int *edc,int anapxdx,float vt,float vtt,float h,float *angx1,float *angx2 )
{
 float x;
 float ang;
 float tanp;
 ang=angx1[it]*PI*0.005556;
 if(ang<-0.00001 || ang>0.00001)
	{
 	tanp=tan(ang)*tan(ang);
 	x=(vt*(tanp-1)+sqrt(vtt*(1+tanp)*(1+tanp)+4*h*h*tanp))*0.5/tan(ang);
 	*edc=MIN(maxcdpout,icdp - ceil(x/anapxdx));
	}
 else	*edc=icdp;
 ang=angx2[it]*PI*0.005556;
 if(ang<-0.00001 || ang>0.00001)
	{
 	tanp=tan(ang)*tan(ang);
 	x=(vt*(tanp-1)+sqrt(vtt*(1+tanp)*(1+tanp)+4*h*h*tanp))*0.5/tan(ang);
 	*bgc=MAX(mincdpout,icdp - ceil(x/anapxdx));
	}
 else	*bgc=icdp;
}

void tanda(int ipx,int *anapx,int sx,int gx,float v,float vtt,float *ttt,float *qtmp)
{
 float xxs,xxg;
 float ts,tg;
 xxs=anapx[ipx]-sx;
 xxg=anapx[ipx]-gx;
 ts=sqrt(xxs*xxs+vtt)/v;
 tg=sqrt(xxg*xxg+vtt)/v;
 *ttt=ts+tg;
 *qtmp=pow((ts/tg),1.5);
}

void windtr(int nt,float *rtx,float ttt,float dt,float *firstt,float *datal)
{
 int itt,itb,ite;
 itb=MAX(ceil(ttt/dt)-4,0);
 ite=MIN(itb+8,nt);
 *firstt=itb*dt;
 for(itt=itb;itt<ite;++itt)
	datal[itt-itb]=rtx[itt];
}

