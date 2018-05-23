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
void tanda(int ipx,int *anapx,float sx,float gx,float v,float vtt,float *ttt,float *qtmp);
void aperture(int it,int icdp,int mincdp,int maxcdp,int *bgc,int *edc,int anapxdx,float vt,float vtt,float h,float *angx1,float *angx2 );
void hammingFilter(int nf1,int nf2,int nf3,int nf4,int nf, float *filter);
void calculateangle(int ximg,float angrange,float angdx,float angdxx,int h,float vt,float vtt,float *thit,int *nthit);
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
 float ***aglerst=NULL;			/* array for storing angle gather		*/
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
 float f1,f2,f3,f4;				/* array of filter frequencies          */
 float angrange;				/* angle range of angle gather			*/
 float angdx;					/* angle interval to image in 			*/
 float angdxx;					/* angdxx=1.0/angdx						*/
 float thit;   
 float ipxb;
 int ximg;
 int *anapx;					/* coordinate of cdp for image space	*/ 
 int *offarr;
 int *offx;  
 int *mincdpx;
 int *maxcdpx;
 int *dcdp;
 int anapxmin;					/* min coordinate of imaging point		*/ 
 int anapxmax;					/* max coordinate of imaging point		*/
 int anapxdx;					/* spacing between imaging point		*/
 int napmin;					/* napmin=(int)(anapxmin/anapxdx)		*/
 int minoff;
 int maxoff;
 int noff;
 int tritvl;					/* trace interval						*/
 int sx,gx;						/* coordinate of shot and geophone  	*/
 int oldcdp;
 int oldcdpt;
 int oldoffset;					/* offset for group devided	and output	*/
 int offset,h;
 int startmt;					/* the start time for migration			*/
 int nstartmt;					/* nstartmt=startmt/dT					*/
 int endmt;						/* the end time for migration			*/
 int nendmt;					/* nendmt=endmt/D=dT					*/
 int hmaxnthit;					/* hmaxnthit=ceil(angrange/angdx)		*/
 int maxnthit;					/* maxnthit=2*hmaxnthit 				*/
 int nthit;						/* angle trace number calculated		*/
 int iagle;
 int icdp;						/* count number                         */
 int itr,ix,ipx,it,jx; 			/* count number                         */
 int nf1,nf2,nf3,nf4;           /* nf1=(int)(f1/df)                     */
 int bgc,edc;					/* begin and end of imaging	trace		*/		
 int mincdp,maxcdp;				/* mincdp and maxcdp of data input		*/
 int mincdpout;					/* the min cdp to output(image space)	*/
 int maxcdpout;					/* the min cdp to output(image space)	*/
 int mincdpo;					/* storing the mincdp for temporary to write header*/
 int firstcdp=0;	            /* first cdp in velocity file	    	*/
 int lastcdp=0;	                /* last cdp in velocity file	    	*/
 int ncdp;	                	/* number of cdps in the velocity file	*/ 
 int nt;                		/* number of points on input trace      */
 int ntr;						/* number of trace input				*/
 int npx;						/* number of trace of imaging sapce		*/
 int nfft;                      /* number of points for fft trace       */
 int nf;                        /* number of frequencies (incl Nyq)     */
 int verbose;		            /* flag to get advisory messages	    */
 int KG=1;
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
 fscanf(fp,"angrange=%f\n",&angrange);
 fscanf(fp,"angdx=%f\n",&angdx);
 efclose(fp);
 time(&t1);
 
 if (fmod(tritvl,anapxdx)!=0)	err("Check anapxdx value in parfile,cause that tritvl must be an integer number of anapxdx");
 noff=(maxoff-minoff)/tritvl+2;
/*caculate image spacing*/
 napmin=anapxmin/anapxdx;
 npx=anapxmax/anapxdx-napmin+1;
 anapx=ealloc1int(npx);
 anapx[0]=anapxmin;
 for(ipx=1;ipx<npx;ipx++)
 	anapx[ipx]=anapx[ipx-1]+anapxdx;
 angdxx=1.0/angdx;
 hmaxnthit=ceil(angrange*angdxx);
 maxnthit=hmaxnthit+hmaxnthit+1;
 angrange=angrange+angdx*0.5;
 mincdpout=mincdp-napmin;
 maxcdpout=maxcdp-napmin;
 nstartmt=0.001*startmt/dt;
 nendmt = 0.001*endmt/dt;
 offarr=ealloc1int(noff);
 offx=ealloc1int(noff);
 dcdp=ealloc1int(noff);
 mincdpx=ealloc1int(noff);
 maxcdpx=ealloc1int(noff);
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
 jx=0;
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
	 /* error trappings. */
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
/* rewind trace file pointer and header file pointer */
 erewind(tracefp);
 erewind(hfp);
/* total number of cdp's in data */
 ncdp=lastcdp-firstcdp+1;
 if(verbose)	warn("ncdp=%d",ncdp);
 for(itr=0;itr<jx;itr++)
	 {
	 mincdpx[itr]-=napmin;
	 maxcdpx[itr]-=napmin;
	 }
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
 aglerst=ealloc3float(nt,maxnthit,npx);

/* Zero all arrays */
 memset((void *) data[0], 0,nt*ntr*FSIZE);
 memset((void *) aglerst[0][0], 0,nt*maxnthit*npx*FSIZE);
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
 nf1=f1/df;
 nf2=f2/df;
 nf3=f3/df;
 nf4=f4/df;
 if (verbose)	warn("f1=%f,f2=%f,f3=%f,f4=%f",f1,f2,f3,f4);
 hammingFilter(nf1,nf2,nf3,nf4,nf,filter);
 
/* Start the migration process */
/* Loop over input trace */
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
 warn("Starting migration process...\n");
for(ix=0;ix<jx;ix++)
	{
	 offset=offx[ix];
	 h=offset*0.5;
	 for(ipx=mincdpout; ipx<=maxcdpout; ++ipx)
		{
	 	T=nstartmt*hdt;
	 	for(it=nstartmt;it<nendmt;it++)
			{
		 	T+=hdt;
		 	v=vel[ipx-firstcdp+napmin][it];
		 	vt=v*T;
		 	vtt=vt*vt;
		 	aperture(it,ipx,mincdpx[ix],maxcdpx[ix],&bgc,&edc,anapxdx,vt,vtt,h,kjmin[ipx-firstcdp+napmin],kjmax[ipx-firstcdp+napmin]);
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
			ipxb=(ipx-mincdpx[ix])/dcdp[ix];
			KG=0;
			/*------------------------------------------------------*/
		 	if(ipxb>edc)		
				{
				sx=(mincdpout+napmin+(edc+1)*dcdp[ix])*anapxdx-h;
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
			 		icdp=(sx+gx)/2/anapxdx-napmin;
					ximg=(icdp-ipx)*anapxdx;
			 		if (ximg==0)	aglerst[ipx][hmaxnthit][it]+=va;
			 		else if (ximg<0)
						{
				 		ximg=-ximg;
			 	 		calculateangle(ximg,angrange,angdx,angdxx,h,vt,vtt,&thit,&nthit);
				 		if(nthit==1000000)	break;
				 		aglerst[ipx][hmaxnthit+nthit][it]+=va*qtmp*p;
						}
			 		else
						{
			 	 		calculateangle(ximg,angrange,angdx,angdxx,h,vt,vtt,&thit,&nthit);
				 		if(nthit==1000000)	break;
				 		aglerst[ipx][hmaxnthit-nthit][it]+=va*qtmp*p;
						}
					}
				}
		 	else if(ipxb<bgc)
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
			 		icdp=(sx+gx)/2/anapxdx-napmin;
					ximg=(icdp-ipx)*anapxdx;
			 		if (ximg==0)	aglerst[ipx][hmaxnthit][it]+=va;
			 		else if (ximg<0)
						{
				 		ximg=-ximg;
			 	 		calculateangle(ximg,angrange,angdx,angdxx,h,vt,vtt,&thit,&nthit);
				 		if(nthit==1000000)	break;
				 		aglerst[ipx][hmaxnthit+nthit][it]+=va*qtmp*p;
						}
			 		else
						{
			 	 		calculateangle(ximg,angrange,angdx,angdxx,h,vt,vtt,&thit,&nthit);
				 		if(nthit==1000000)	break;
				 		aglerst[ipx][hmaxnthit-nthit][it]+=va*qtmp*p;
						}
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
			 		icdp=(sx+gx)/2/anapxdx-napmin;
					ximg=(icdp-ipx)*anapxdx;
			 		if (ximg==0)	aglerst[ipx][hmaxnthit][it]+=va;
			 		else if (ximg<0)
						{
				 		ximg=-ximg;
			 	 		calculateangle(ximg,angrange,angdx,angdxx,h,vt,vtt,&thit,&nthit);
				 		if(nthit==1000000)	break;
				 		aglerst[ipx][hmaxnthit+nthit][it]+=va*qtmp*p;
						}
			 		else
						{
			 	 		calculateangle(ximg,angrange,angdx,angdxx,h,vt,vtt,&thit,&nthit);
				 		if(nthit==1000000)	break;
				 		aglerst[ipx][hmaxnthit-nthit][it]+=va*qtmp*p;
						}
					}
				}
			/*-----------------------------------------------------*/
			}
		}
	}
	memset ((void *) &tro, (int) '\0', sizeof (tro));
 	tro.trid = 1;
 	tro.counit = 1;
 	tro.f2 = angdx*0.5-angrange;
 	tro.d2 = angdx;
 	tro.ns = nt;
 	tro.dt = dt*1000000;
	/* Output migrated data */
	mincdpo=mincdp;
 	for(ipx=mincdpout; ipx<=maxcdpout; ++ipx)
    {
	for(iagle=0;iagle<maxnthit;iagle++)
    	{
     	 tro.cdp = mincdpo;
     	 tro.offset=tro.f2+iagle*angdx;
     	 memcpy ((void *) tro.data, (const void *)  aglerst[ipx][iagle],sizeof (float) * nt);
     	 puttr(&tro);
		}
     mincdpo++;
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
 free1float(rtx);
 free2float(vel);
 free2float(kjmin);
 free2float(kjmax);
 free1float(filter);
 free1complex(ct);
 free1complex(hd);
 free2float(data);
 free3float(aglerst);
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
 	*edc=MIN(maxcdp,icdp - ceil(x/anapxdx));
	}
 else	*edc=icdp;
 ang=angx2[it]*PI*0.005556;
 if(ang<-0.00001 || ang>0.00001)
	{
 	tanp=tan(ang)*tan(ang);
 	x=(vt*(tanp-1)+sqrt(vtt*(1+tanp)*(1+tanp)+4*h*h*tanp))*0.5/tan(ang);
 	*bgc=MAX(mincdp,icdp - ceil(x/anapxdx));
	}
 else	*bgc=icdp;
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

void windtr(int nt,float *rtx,float ttt,float dt,float *firstt,float *datal)
{
 int itt,itb,ite;
 itb=MAX(ceil(ttt/dt)-4,0);
 ite=MIN(itb+8,nt);
 *firstt=itb*dt;
 for(itt=itb;itt<ite;++itt)
	datal[itt-itb]=rtx[itt];
}
void calculateangle(int ximg,float angrange,float angdx,float angdxx,int h,float vt,float vtt,float *thit,int *nthit)
{
int nth;
float thitc,thitcc;
thitc=(ximg*ximg-h*h-vtt)*0.5/ximg/vt;
thitc=thitc+sqrt(thitc*thitc+1);
thitc=atan(thitc)*57.2957795;
*thit=thitc;
/*thit=28.6478897*(atan((ximg+h)/vt)+atan((ximg-h)/vt)); two method to caculate the angles,test shows,the effect of latter method is lower*/
if(thitc<angrange)
	{
	 thitc=thitc*angdxx;
	 nth=floor(thitc);
	 thitcc=thitc-nth;
	 if(thitcc>0.5)	nth=nth+1;
	}
 else	nth=1000000;
 *nthit=nth;
}

