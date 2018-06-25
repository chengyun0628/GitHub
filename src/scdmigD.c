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
void windtr(int nt,float *rtw,float ttt,float dt,float *firstt,float *datal);
void tanda(int ipx,int *anapx,float sx,float gx,float v,float vtt,float *ttt,float *qtmp);
void aperture(int it,int icdp,int mincdp,int maxcdp,int *bgc,int *edc,int anapxdx,float vt,float vtt,float h,float *angx1,float *angx2 );
void hammingFilter(int nf1,int nf2,int nf3,int nf4,int nf, float *filter);
#define LOOKFAC 2
#define PFA_MAX 720720
segy tri;                       /* trace of input  						*/
segy tro;						/* trace of output 						*/
int
main(int argc, char **argv)
{
 register float *rt;       		/* real trace data 						*/
 register complex *ct;          /* complex transformed trace            */
 register complex *hd;			/* half derivative						*/
 float datal[8];				/* small cut of ununiform sample trace	*/ 
 float *filter;              	/* filter array                         */
 float **data;					/* temp_array of 2D real trace			*/
 float **mig=NULL;				/* array for storing imaging gather		*/
 float **vel=NULL;              /* array for storing velocity           */
 float **kjmin=NULL;            /* array for storing the min aperture   */
 float **kjmax=NULL;            /* array for storing the max aperture   */

 float ttt;						/* Travel time							*/
 float qtmp;					/* Amptitude of migration				*/	
 float va;    					/* ununiform sampled value				*/ 
 float firstt;					/* the first sample in datal[8]			*/
 float dt,hdt;					/* sample interval & half of that		*/
 float T;						/* vertical time depth 					*/
 float p;						/* Boundary attenuation factor of amp	*/
 float df,dw;					/* freqency sample spacing				*/
 float tmax;					/* the max trace length					*/
 float v,vt,vtt;				/* velocity,vt=v*T,vtt=vt*vt			*/
 float f1,f2,f3,f4;				/* array of filter frequencies          */
 
 int *offx;  					/* store each co gather's offset		*/
 int *offarr;					/* store ~'s first trace num			*/ 
 int *mincdpx;					/* store ~'s min cdp num				*/
 int *maxcdpx;					/* store ~'s max cdp num				*/
 int *dcdp; 					/* store ~'s cdp interval				*/
 int *anapx;					/* coordinate of cdp for image space	*/ 

 int anapxmin,anapxmax;			/* min & max coordinate of image space	*/ 
 int anapxdx;					/* interval of image space				*/
 int napmin;					/* napmin=(int)(anapxmin/anapxdx)		*/
 int minoff,maxoff;				/* min & max offset of data input		*/
 int noff;						/* co group num for data input			*/ 
 int tritvl;					/* trace interval						*/
 int sx,gx;						/* coordinate of shot and geophone  	*/
 int oldcdp;					/* tmp value for cdp counter			*/
 int oldcdpt;					/* tmp value for cdp counter			*/
 int oldoffset;					/* tmp value for co group devide		*/
 int offset,h;					/* offset and half of that 				*/
 int startmt,endmt;				/* the start & end time for migration	*/
 int nf1,nf2,nf3,nf4;           /* nf1=(int)(f1/df)                     */
 int bgc,edc;					/* begin & end of imaging trace			*/		
 int mincdp,maxcdp;				/* min & max cdp of data input			*/
 int mincdpout,maxcdpout;		/* min & max cdp to output(image space) */
 int firstcdp,lastcdp;	        /* first & last cdp in velocity file	*/
 int ncdp;	                	/* number of cdps in the velocity file	*/ 
 int nkj;						/* number of trace in the aperture file	*/
 int nt;                		/* number of points on input trace      */
 int ntr;						/* number of trace input				*/
 int npx;						/* number of trace of imaging sapce		*/
 int nfft;                      /* number of points for fft trace       */
 int nf;                        /* number of frequencies (incl Nyq)     */
 int verbose;		            /* flag to get advisory messages	    */
 int icdp,itr,ipx,it,ix,jx; 	/* count number                         */
 int KG=1;						/* counter switch						*/
 int nbjl,nbjr,nbj1,nbj2;		/* acture range of migration			*/

/* file name */
 char str[100],*path;
 char *vfile="";
 char *kjfile1="";
 char *kjfile2="";
 char *parfile="";

 FILE *fp;						/* temp file to hold parfile      		*/
 FILE *vfp=NULL;				/* temp file to hold velocity      		*/
 FILE *minbj=NULL;				/* temp file to hold left aperture      */
 FILE *maxbj=NULL;				/* temp file to hold right aperture     */
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
 oldcdp=tri.cdp;
 oldoffset=tri.offset;
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
 fscanf(fp,"minoff=%d\n",&minoff);
 fscanf(fp,"maxoff=%d\n",&maxoff);
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
 if (fmod(tritvl,anapxdx)!=0)	err("Check anapxdx value in parfile,cause that tritvl must be an integer number of anapxdx");
 time(&t1);
 
/*caculate image spacing*/
 napmin=anapxmin/anapxdx;
 npx=anapxmax/anapxdx-napmin+1;
 anapx=ealloc1int(npx);
 anapx[0]=anapxmin;
 for(ipx=1;ipx<npx;ipx++)
 	anapx[ipx]=anapx[ipx-1]+anapxdx;

 mincdpout=mincdp-napmin;
 maxcdpout=maxcdp-napmin;

 startmt=0.001*startmt/dt;
 endmt = 0.001*endmt/dt;

 noff=(maxoff-minoff)/tritvl+2;
 ncdp=lastcdp-firstcdp+1;
 nkj=maxcdp-mincdp+1;
 
 /* Allocate space */
 offarr=ealloc1int(noff);
 offx=ealloc1int(noff);
 dcdp=ealloc1int(noff);
 mincdpx=ealloc1int(noff);
 maxcdpx=ealloc1int(noff);

 /* Zero all arrays */
 memset((void *) offarr, 0, noff*FSIZE);
 memset((void *) offx, 0, noff*FSIZE);
 memset((void *) dcdp, 0, noff*FSIZE);
 memset((void *) mincdpx, 0, noff*FSIZE);
 memset((void *) maxcdpx, 0, noff*FSIZE);

/* Store traces in tmpfile while getting a count of number of traces */
 tracefp = etmpfile();
 hfp = etmpfile();
 mincdpx[0]=oldcdp;
 ntr = 0;
 jx = 0;
 do 
	{
	 /* read headers and data */
	 efwrite(&tri,HDRBYTES, 1, hfp);
	 efwrite(tri.data, FSIZE, nt, tracefp);
	 if(ntr>0&&KG==1)	
		{
		 dcdp[jx]=tri.cdp-oldcdp;
		 offx[jx]=tri.offset;
		 KG=0;
		}

	 /* ...did offset value change? */ 
	 if ((ntr>0) && ( oldoffset!=tri.offset))
		{
		 mincdpx[jx+1]=tri.cdp;
		 maxcdpx[jx]=oldcdpt;
		 offarr[jx+1]=ntr;
		 oldcdp=tri.cdp;
		 KG=1;
		 ++jx;	
		}
	  
	  ++ntr;
	  oldcdpt=tri.cdp;
	  oldoffset=tri.offset;
	} while (gettr(&tri));

	maxcdpx[jx]=tri.cdp;
	jx++;
	offarr[jx]=ntr;
	warn("ntr=%d",ntr);

 /* coordinate transform from real sapce to image space */
 for(itr=0;itr<jx;itr++)
	 {
	 mincdpx[itr]-=napmin;
	 maxcdpx[itr]-=napmin;
	 }

/* rewind trace file pointer and header file pointer */
 erewind(tracefp);
 erewind(hfp);
 
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
 ct=ealloc1complex(nf);
 hd=ealloc1complex(nf);
 filter=ealloc1float(nf);
 data=ealloc2float(nt,ntr);
 vel=ealloc2float(nt,ncdp);
 kjmin=ealloc2float(nt,nkj);
 kjmax=ealloc2float(nt,nkj);
 mig=ealloc2float(nt,npx);

/* Zero all arrays */
 memset((void *) data[0], 0,nt*ntr*FSIZE);
 memset((void *) mig[0], 0,nt*npx*FSIZE);
 memset((void *) rt, 0, nfft*FSIZE);
 memset((void *) filter, 0, nf*FSIZE);
 
/* Define half derivative*/
 for(ix=0;ix<nf;ix++)
	{
	 hd[ix].r=0.707107;
	 hd[ix].i=-0.707107;
	}

/* Get frequency parameter of filter */
 nf1=f1/df;
 nf2=f2/df;
 nf3=f3/df;
 nf4=f4/df;
 hammingFilter(nf1,nf2,nf3,nf4,nf,filter);

/* Read data from temporal array */
 for(itr=0;itr<ntr;++itr)
	 efread(data[itr],FSIZE,nt,tracefp);
	
/* read velocities & aperture */
 vfp=efopen(vfile,"r");
 efread(vel[0],FSIZE,nt*ncdp,vfp);
 efclose(vfp);
 minbj=efopen(kjfile1,"r");
 efread(kjmin[0],FSIZE,nt*nkj,minbj);
 efclose(minbj);
 maxbj=efopen(kjfile2,"r");
 efread(kjmax[0],FSIZE,nt*nkj,maxbj);
 efclose(maxbj);

/* filter the input trace */
for(itr=0; itr<ntr; itr++)
	{
	 for(it=0; it<nt; ++it)
		 rt[it]=data[itr][it];

	/* zero array ct */
	 memset((void *) ct, 0, nf*FSIZE);

	/* FFT and multiply by the half derivative & filter*/
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
	 pfacr(-1,nfft,ct,rt);

	 for(it=0;it<nt;it++)
		data[itr][it]=rt[it]/nfft;

	/* zero array rt */
	 memset((void *) rt, 0, nfft*FSIZE);
	}

/* Start the migration process */
/* loop over each common offset gather*/
 warn("Starting migration process...\n");
 for(ix=0;ix<jx;ix++)
	{
	 offset=offx[ix];
	 h=offset*0.5;
	 for(ipx=mincdpout; ipx<=maxcdpout; ++ipx)
		{
	 	T=startmt*hdt;
	 	for(it=startmt;it<endmt;it++)
			{
		 	T+=hdt;
		 	v=vel[ipx-firstcdp+napmin][it];
		 	vt=v*T;
		 	vtt=vt*vt;
		 	aperture(it,ipx,mincdpx[ix],maxcdpx[ix],&bgc,&edc,anapxdx,vt,vtt,h,kjmin[ipx-mincdpout],kjmax[ipx-mincdpout]);
			if (bgc==999999)	continue;
			nbjl=(edc-bgc+1)*0.2;
		 	nbj1=(bgc-nbjl)>mincdpx[ix]?nbjl:(bgc-mincdpx[ix]);
		 	nbj2=(edc+nbjl)<maxcdpx[ix]?nbjl:(maxcdpx[ix]-edc);
		 	if((bgc-nbj1-mincdpx[ix])%dcdp[ix]!=0)	bgc=(bgc-nbj1-mincdpx[ix])/dcdp[ix]+1;
		 	else	bgc=(bgc-nbj1-mincdpx[ix])/dcdp[ix];
		 	edc=(edc+nbj2-mincdpx[ix])/dcdp[ix];
		 	nbj1=nbj1/dcdp[ix];
		 	nbj2=nbj2/dcdp[ix];	
		 	nbjl=edc-bgc-nbj1;
		 	nbjr=edc-bgc-nbj2;
			KG=0;
		 	if(ipx>=maxcdpx[ix])		
				{
				sx=(mincdpx[ix]+napmin+(edc+1)*dcdp[ix])*anapxdx-h;
		 		for(itr=offarr[ix]+edc;itr>=offarr[ix]+bgc;itr--) 
					{
					sx=sx-dcdp[ix]*anapxdx;
		 	 		gx=sx+offset;
     				tanda(ipx,anapx,sx,gx,v,vtt,&ttt,&qtmp);
             		if(ttt>=tmax)   break;
             		windtr(nt,data[itr],ttt,dt,&firstt,datal);
         	 		ints8r(8, dt, firstt, datal,
             	  		0.0, 0.0, 1, &ttt, &va);
			 		if(KG<nbj2)	p=sin(1.570796*KG/nbj2);
			 		else if(KG>nbjl)	p=cos(1.570796*(KG-nbjl)/nbj1);
			 		else p=1.0;
			 		KG++;
			 		mig[ipx][it]+=va*qtmp*p;
					}
				}
		 	else if(ipx<=mincdpx[ix])
				{
				sx=(mincdpx[ix]+napmin+(bgc-1)*dcdp[ix])*anapxdx-h;
				for(itr=offarr[ix]+bgc;itr<=offarr[ix]+edc;itr++) 
					{
					sx=sx+dcdp[ix]*anapxdx;
		 	 		gx=sx+offset;
     				tanda(ipx,anapx,sx,gx,v,vtt,&ttt,&qtmp);
             		if(ttt>=tmax)   break;
             		windtr(nt,data[itr],ttt,dt,&firstt,datal);
         	 		ints8r(8, dt, firstt, datal,
             	  		0.0, 0.0, 1, &ttt, &va);
			 		if(KG<nbj1)	p=sin(1.570796*KG/nbj1);
			 		else if(KG>nbjr)	p=cos(1.570796*(KG-nbjr)/nbj2);
			 		else p=1.0;
			 		KG++;
			 		mig[ipx][it]+=va*qtmp*p;
					}
				}
		 	else
				{
				sx=(mincdpx[ix]+napmin+(bgc-1)*dcdp[ix])*anapxdx-h;
				for(itr=offarr[ix]+bgc;itr<=offarr[ix]+edc;itr++) 
					{
					sx=sx+dcdp[ix]*anapxdx;
		 	 		gx=sx+offset;
     				tanda(ipx,anapx,sx,gx,v,vtt,&ttt,&qtmp);
             		if(ttt>=tmax)   continue;
             		windtr(nt,data[itr],ttt,dt,&firstt,datal);
         	 		ints8r(8, dt, firstt, datal,
             	  		0.0, 0.0, 1, &ttt, &va);
			 		if(KG<nbj1)	p=sin(1.570796*KG/nbj1);
			 		else if(KG>nbjr)	p=cos(1.570796*(KG-nbjr)/nbj2);
			 		else p=1.0;
			 		KG++;
			 		mig[ipx][it]+=va*qtmp*p;
					}
				}
			}
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
 for(ipx=mincdpout; ipx<=maxcdpout; ++ipx)
    {
     tro.cdp = mincdp;
     memcpy ((void *) tro.data, (const void *)  mig[ipx],sizeof (float) * nt);
     puttr(&tro);
     mincdp++;
    }
 time(&t2);
 warn("time consuming in second = %f\n",difftime(t2,t1));
/*free array*/
 efclose(hfp);
 efclose(tracefp);
 free1int(anapx);
 free1int(mincdpx);
 free1int(maxcdpx);
 free1int(dcdp);
 free1int(offx);
 free1int(offarr);
 free1float(rt);
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
        
void aperture(int it,int icdp,int mincdp,int maxcdp,int *bgc,int *edc,int anapxdx,float vt,float vtt,float h,float *angx1,float *angx2 )
{
 float x;
 float ang;
 float tanp;
 ang=angx1[it]*PI*0.005556;	
 if(ang<-0.00001 || ang>0.00001)
	{
 	tanp=tan(ang)*tan(ang);
 	x=(vt*(tanp-1)+sqrt(vtt*(1+tanp)*(1+tanp)+4*h*h*tanp))*0.5/tan(ang);
 	*edc=icdp - ceil(x/anapxdx);
	}
 else	*edc=icdp;

 ang=angx2[it]*PI*0.005556;
 if(ang<-0.00001 || ang>0.00001)
	{
 	tanp=tan(ang)*tan(ang);
 	x=(vt*(tanp-1)+sqrt(vtt*(1+tanp)*(1+tanp)+4*h*h*tanp))*0.5/tan(ang);
 	*bgc=icdp - ceil(x/anapxdx);
	}
 else	*bgc=icdp;

 if (*bgc>maxcdp || *edc<mincdp)	*bgc=999999;
 if (mincdp>*bgc)	*bgc=mincdp;
 if (*bgc!=999999 && *edc>maxcdp)	*edc=maxcdp;
}

void tanda(int ipx,int *anapx,float sx,float gx,float v,float vtt,float *ttt,float *qtmp)
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

void windtr(int nt,float *rtw,float ttt,float dt,float *firstt,float *datal)
{
 int itt,itb,ite;
 itb=MAX(ceil(ttt/dt)-4,0);
 ite=MIN(itb+8,nt);
 *firstt=itb*dt;
 for(itt=itb;itt<ite;++itt)
	datal[itt-itb]=rtw[itt];
}
