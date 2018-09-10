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
void get_sx_gx_offset(float *sx, float *gx,float *offset);
void windtr(int nt,float *rtx,float ttt,float dt,float *firstt,float *datal);
void mutefct(float h,float hdt,float v1,float T1,float v2,float T2,float *ttt,float *qtmp);
void tanda(int ipx,float *anapx,float sx,float gx,float v1,float T1,float *ttt,float *qtmp);
void aperture(int it,int icdp,int mincdpout,int maxcdpout,int *bgc,int *edc,float anapxdx,float vt,float vtt,float h,float *angxb );
void angleintsmt(int nt,float dt,float tm1,float tm2,float tm3,float tm4,float ang1,float ang2,float ang3,float ang4,float *angx);
void hammingFilter(int nf1,int nf2,int nf3,int nf4,int nf, float *filter);
#define LOOKFAC 2
#define PFA_MAX 720720
segy tri;                        /* trace of input  					*/
segy tro;						/* trace of output 						*/
int
main(int argc, char **argv)
{
 register float *rt,*rtx;       /* real trace data and processed data	*/
 register complex *ct;          /* complex transformed trace            */
 register complex *hd;			/* half derivative						*/
 float *filter;              	/* filter array                         */
 float *angx;					/* array of given aperture angle 		*/
 float *angxb;					/* array of given aperture angle 		*/
 float *anapx;					/* coordinate of cdp for image space	*/
 float datal[8];				/* small cut of ununiform sample trace	*/ 
 float **data;					/* temp_array of 2D real trace			*/
 float **vel=NULL;              /* array for storing velocity           */
 float **mig=NULL;				/* array for storing migrated result	*/
 float ttt;						/* Travel time							*/
 float qtmp;					/* Amptitude of migration				*/	
 float va;    					/* ununiform sampled value				*/ 
 float firstt;					/* the first sample in datal[8]			*/
 float f1,f2,f3,f4;				/* array of filter frequencies          */   
 float sx,gx;					/* coordinate of shot and geophone  	*/
 float oldoffset;				/* offset for group devided	and output	*/
 float offset,h;				/* offset and half of that				*/
 float cdp;						/* coordinate of cdp					*/
 float tritvl;					/* trace interval						*/
 float ximg;					/* distance between cdp and image point	*/
 float dt;						/* sample interval						*/
 float hdt;						/* half of sample interval				*/
 float angtmp;					/* tmp variable for angle cacluation	*/
 float tanp;             		/* tanp=pow(tan(ang),2); 				*/
 float T1,T2;					/* vertical time depth 					*/
 float p;						/* Boundary attenuation factor of amp	*/
 float df,dw;					/* freqency sample spacing				*/
 float tmax=0;					/* the max trace length					*/
 float tmin=0;					/* the min time survive from mute		*/
 float v1;						/* velocity1							*/
 float v2;						/* velocity2							*/
 float vt;						/* vt=v*T								*/
 float vtt;						/* vtt=vt*vt							*/
 float ang1,ang2,ang3,ang4;		/* given migration aperture angle		*/
 float anapxmin;				/* min coordinate of imaging point		*/ 
 float anapxmax;				/* max coordinate of imaging point		*/
 float anapxdx;					/* spacing between imaging point		*/
 float tm1,tm2,tm3,tm4;			/* time corespond to apreture angle		*/
 float smute=0;					/* strech mute factor,smute=0 ->no mute */
 float osmute;					/* osmute=1/smute					 	*/
 float startmt;					/* the start time for migration			*/
 float endmt;					/* the end time for migration			*/
 int nstartmt;					/* nstartmt=startmt/dT					*/
 int nendmt;					/* nendmt=endmt/D=dT					*/
 int itr,ix,ipx,it; 			/* count number                         */
 int itmute;					/* the start migration time for mute	*/
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
 int mincdpo;					/* storing the mincdp for temporary to write header	*/
 int firstcdp=0;	            /* first cdp in velocity file	    	*/
 int lastcdp=0;	                /* last cdp in velocity file	    	*/
 int ncdp;	                	/* number of cdps in the velocity file	*/ 
 int dcdp=0;	                /* number of cdps between consecutive traces 		*/
 int nt;                		/* number of points on input trace      */
 int ntr;						/* number of trace input				*/
 int npx;						/* number of trace of imaging sapce		*/
 int nfft;                      /* number of points for fft trace       */
 int nf;                        /* number of frequencies (incl Nyq)     */
 int verbose;		            /* flag to get advisory messages	    */

/* file name */
 char str[100],*path;
 char *vfile="";
 char *parfile="";
 FILE *fp;
 FILE *fp1;
 FILE *vfp=NULL;
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
 if(!getparstring("parfile",&parfile))	err("parameter file must be specified ! ");
 sprintf(str,"%s/par/%s",path,parfile);
 warn("str=%s",str);
 fp=fopen(str,"rb");
 fscanf(fp,"mincdp=%d\n",&mincdp);
 fscanf(fp,"maxcdp=%d\n",&maxcdp);
 fscanf(fp,"firstcdp=%d\n",&firstcdp);
 fscanf(fp,"lastcdp=%d\n",&lastcdp);
 fscanf(fp,"startmt=%f\n",&startmt);
 fscanf(fp,"endmt=%f\n",&endmt);
 fscanf(fp,"tm1=%f\n",&tm1);
 fscanf(fp,"tm2=%f\n",&tm2);
 fscanf(fp,"tm3=%f\n",&tm3);
 fscanf(fp,"tm4=%f\n",&tm4);
 fscanf(fp,"ang1=%f\n",&ang1);
 fscanf(fp,"ang2=%f\n",&ang2);
 fscanf(fp,"ang3=%f\n",&ang3);
 fscanf(fp,"ang4=%f\n",&ang4);
 fscanf(fp,"tritvl=%f\n",&tritvl);
 fscanf(fp,"f1=%f\n",&f1);
 fscanf(fp,"f2=%f\n",&f2);
 fscanf(fp,"f3=%f\n",&f3);
 fscanf(fp,"f4=%f\n",&f4);
 fscanf(fp,"anapxmin=%f\n",&anapxmin);
 fscanf(fp,"anapxmax=%f\n",&anapxmax);
 fscanf(fp,"anapxdx=%f\n",&anapxdx);
 fscanf(fp,"smute=%f\n",&smute);
 efclose(fp);
 time(&t1);
 tm1=tm1/1000;
 tm2=tm2/1000;
 tm3=tm3/1000;
 tm4=tm4/1000;
 if (tm1>tmax||tm2>tmax||tm3>tmax||tm4>tmax)	err("Check tm values in data!");
 if (fmod(tritvl,anapxdx)!=0)	err("Check anapxdx value in parfile,cause that tritvl must be an integer number of anapxdx");
 if(smute!=0)	osmute=1/smute;

/*caculate image spacing*/
 napmin=(int)(anapxmin/anapxdx);
 npx=(int)(anapxmax/anapxdx)-napmin+1;
 anapx=ealloc1float(npx);
 anapx[0]=anapxmin;
 for(ipx=1;ipx<npx;ipx++)
 	anapx[ipx]=anapx[ipx-1]+anapxdx;
 mincdpout=mincdp-napmin;
 maxcdpout=maxcdp-napmin;
 nstartmt=itmute=ceil(0.001*startmt/dt)+1;
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
 if (!getparint("dcdp",&dcdp))	dcdp=deltacdp - 1;
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

/* Allocate space */
 angx=ealloc1float(nt);
 angxb=ealloc1float(nt);
 rt=ealloc1float(nfft);
 rtx=ealloc1float(nfft);
 ct=ealloc1complex(nf);
 hd=ealloc1complex(nf);
 filter=ealloc1float(nf);
 data=ealloc2float(nt,ntr);
 vel=ealloc2float(nt,ncdp);
 mig=ealloc2float(nt,npx);

/* Zero all arrays */
 memset((void *) data[0], 0,nt*ntr*FSIZE);
 memset((void *) mig[0], 0,nt*npx*FSIZE);
 memset((void *) angx, 0, nt*FSIZE);
 memset((void *) angxb, 0, nt*FSIZE);
 memset((void *) rt, 0, nfft*FSIZE);
 memset((void *) filter, 0, nf*FSIZE);
 
/* calculate aperture of migration */
 angleintsmt(nt,dt,tm1,tm2,tm3,tm4,ang1,ang2,ang3,ang4,angx);
 for(it=0;it<nt;it++)
	angxb[it]=angx[it]*1.25;

/* Read data from temporal array */
 for(itr=0;itr<ntr;++itr)
	{
	 efread(data[itr],FSIZE,nt,tracefp);
	}
	
/* read velocities */
 vfp=efopen(vfile,"r");
 efread(vel[0],FSIZE,nt*ncdp,vfp);
 efclose(vfp);

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
if(verbose)	
	{
	 fp1=fopen("filter.bin","wb");
	 fwrite(filter,sizeof(float),nf,fp1);
	 fclose(fp1);
	}
/* Start the migration process */
/* Loop over input trace */
 warn("Starting migration process...\n");
 for(itr=0; itr<ntr; ++itr)
	{
	 float perc;
	 perc=itr*100.0/(ntr-1);
	 if(fmod(itr*100.0,ntr-1)==0)
		warn("migrated %g\n",perc);
	 efread(&tri,HDRBYTES, 1, hfp);
	 /* caculate the coordinate of shot and geophone*/
	 get_sx_gx_offset(&sx,&gx,&offset);
	 if(oldoffset!=offset)
		{
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
     		 tro.offset=fabs(oldoffset);
     		 memcpy ((void *) tro.data, (const void *)  mig[ipx],sizeof (float) * nt);
     		 puttr(&tro);
     		 mincdpo++;
    		}
		 oldoffset=offset;
		 memset((void *) mig[0], 0,nt*npx*FSIZE);
		}

	 h=offset/2;	
     cdp=(sx+gx)/2;
     icdp=(int)(cdp/anapxdx)-napmin;

	 for(it=0; it<nt; ++it)
		 rt[it]=data[itr][it];

	/* zero array ct and rtx*/
	 memset((void *) rtx, 0, nfft*FSIZE);
	 memset((void *) ct, 0, nf*FSIZE);

	/* filtering the input trace*/
	 pfarc(1,nfft,rt,ct);
	 for(ix=0;ix<nf;ix++) 
		{	
		 if(ix>=nf1&&ix<nf4)	ct[ix]=crmul(ct[ix],filter[ix-nf1]);
		 else	ct[ix].r=ct[ix].i=0.0;
		}

	/* multiply by the half derivative*/ 
	 for(ix=0;ix<nf;ix++)
		{
		 ct[ix]=crmul(ct[ix],sqrt(ix*dw));
         ct[ix]=cmul(ct[ix],hd[ix]);
		}
	 pfacr(-1,nfft,ct,rtx);
	 for(it=0;it<nt;it++)
		rtx[it]=rtx[it]/nfft;

	/* determine index of first sample to survive mute */
 	 if(smute!=0)
		{
		tmin=0;
		itmute=nstartmt;
	 	T2=(nstartmt-2)*hdt;
	 	for(it=nstartmt;it<=nendmt;it++)
			{
		 	T2+=hdt;
		 	v2=vel[icdp-firstcdp+napmin][it-1];
		 	T1=T2+hdt;
		 	v1=vel[icdp-firstcdp+napmin][it];
		 	mutefct(h,hdt,v1,T1,v2,T2,&ttt,&qtmp);
		 	if	(qtmp<osmute)	tmin=ttt;
			}
	 	T1=(nendmt+1)*hdt;
	 	for(it=nendmt;it>=nstartmt;it--)
			{
		 	T1-=hdt;
		 	v1=vel[icdp-firstcdp+napmin][it];
 		 	ttt=2*hypotf(T1,h/v1);
		 	if(ttt<=tmin)
				{
				itmute=it;
				break;
				}	
			}
		}
	/* loop in the image space*/
	 T1=(itmute-1)*hdt;
	 for(it=itmute;it<=nendmt;it++)
		{
		 T1+=hdt;
		 v1=vel[icdp-firstcdp+napmin][it];
		 ttt=2*hypotf(T1,h/v1);
		 if(ttt>=tmax)   break;
		 vt=v1*T1;
		 vtt=vt*vt;
		 aperture(it,icdp,mincdpout,maxcdpout,&bgc,&edc,anapxdx,vt,vtt,h,angxb);
		 windtr(nt,rtx,ttt,dt,&firstt,datal);
		 /* sinc interpolate new data */
		 ints8r(8, dt, firstt, datal,
			 0.0, 0.0, 1, &ttt, &va);
		 mig[icdp][it]+=va;
		 p=1.0;
		 ximg=0;
		 for(ipx=icdp-1;ipx>bgc;ipx--)
		 	{
			 ximg=ximg+anapxdx;
			 v1=vel[ipx-firstcdp+napmin][it];
			 vt=v1*T1;
			 tanda(ipx,anapx,sx,gx,v1,T1,&ttt,&qtmp);
             if(ttt>=tmax)   break;
             windtr(nt,rtx,ttt,dt,&firstt,datal);
			 /* sinc interpolate new data */
         	 ints8r(8, dt, firstt, datal,
             	  0.0, 0.0, 1, &ttt, &va);
			 angtmp=(atan((ximg+h)/vt)+atan((ximg-h)/vt))*28.64789;
			 if(angtmp>angx[it])	p=cos(6.28318*(angtmp/angx[it]-1));
			 mig[ipx][it]+=va*qtmp*p;
			}
		 p=1.0;	
		 ximg=0;
		 for(ipx=icdp+1;ipx<edc;ipx++)
			{
			 ximg=ximg+anapxdx;
             v1=vel[ipx-firstcdp+napmin][it];
			 vt=v1*T1;
			 tanda(ipx,anapx,sx,gx,v1,T1,&ttt,&qtmp);
             if(ttt>=tmax)   break;
             windtr(nt,rtx,ttt,dt,&firstt,datal);
			 /* sinc interpolate new data */
              ints8r(8, dt, firstt, datal,
                  0.0, 0.0, 1, &ttt, &va);
			 angtmp=(atan((ximg+h)/vt)+atan((ximg-h)/vt))*28.64789;
			 if(angtmp>angx[it])	p=cos(6.28318*(angtmp/angx[it]-1));
			 mig[ipx][it]+=va*qtmp*p;
			}
		}
	if(itr==ntr-1)
		{
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
     		 tro.offset=fabs(oldoffset);
     		 memcpy ((void *) tro.data, (const void *)  mig[ipx],sizeof (float) * nt);
     		 puttr(&tro);
     		 mincdpo++;
    		}
		}
	}
 time(&t2);
 warn("time consuming in second = %f\n",difftime(t2,t1));
/*free array*/
 efclose(hfp);
 efclose(tracefp);
 free1float(rt);
 free1float(rtx);
 free1float(angx);
 free1float(angxb);
 free1float(anapx);
 free2float(vel);
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

void angleintsmt(int nt,float dt,float tm1,float tm2,float tm3,float tm4,float ang1,float ang2,float ang3,float ang4,float *angx)
{
 int it,ib;
 float *xin,*yin,*xout;
 float ydin[4][4];
 xin=ealloc1float(4);
 memset((void *) xin, 0, 4*FSIZE);
 yin=ealloc1float(4);
 memset((void *) yin, 0, 4*FSIZE);
 xout=ealloc1float(nt);
 memset((void *) xout, 0, nt*FSIZE);
 xin[0]=tm1;
 yin[0]=ang1;
 xin[1]=tm2;
 yin[1]=ang2;
 xin[2]=tm3;
 yin[2]=ang3;
 xin[3]=tm4;
 yin[3]=ang4;
 for(it=0;it<4;it++)
    for(ib=0;ib<4;ib++)
        ydin[it][ib]=0;
 for(it=0;it<nt;it++)
    xout[it]=dt*it;
 cakima(4,xin,yin,ydin);
// cmonot(4,xin,yin,ydin);
// csplin(4,xin,yin,ydin);
// chermite(4,xin,yin,ydin);
 intcub(0,4,xin,ydin,nt,xout,angx);
 free1float(xin);
 free1float(yin);
 free1float(xout);
}
        
void aperture(int it,int icdp,int mincdpout,int maxcdpout,int *bgc,int *edc,float anapxdx,float vt,float vtt,float h,float *angxb )
{
 float x;
 float ang;
 float tanp;
 ang=angxb[it]*PI*0.005556;
 tanp=tan(ang)*tan(ang);
 x=(vt*(tanp-1)+sqrt(vtt*(1+tanp)*(1+tanp)+4*h*h*tanp))/2/tan(ang);
 *bgc=MAX(mincdpout,icdp - ceil(x/anapxdx)); 
 *edc=MIN(maxcdpout,icdp + ceil(x/anapxdx));
}

void mutefct(float h,float hdt,float v1,float T1,float v2,float T2,float *ttt,float *qtmp)
{
 float ttmp;
 ttmp= hypotf(T2,h/v2);
 *ttt= hypotf(T1,h/v1);
 *qtmp=(*ttt-ttmp)/hdt;
 *ttt=*ttt+*ttt;
}

void tanda(int ipx,float *anapx,float sx,float gx,float v1,float T1,float *ttt,float *qtmp)
{
 float xxs,xxg;
 float ts,tg;
 xxs=anapx[ipx]-sx;
 xxg=anapx[ipx]-gx;
 ts = hypotf(T1,xxs/v1);
 tg = hypotf(T1,xxg/v1);
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

void get_sx_gx_offset(float *sx, float *gx,float *offset)
{ 
  /*****************************************************************************
 get_sx_gx_offset - get sx,gx and offset from headrs
  *****************************************************************************/
  float sy;		/* source coordinates */
  float gy;		/* geophone coordinates */
 if (tri.scalco) 
	{ 
	 /* if tri.scalco is set, apply value */
	 if (tri.scalco>0) 
		{	
		 *sx = (float) tri.sx*tri.scalco;
		 *gx = (float) tri.gx*tri.scalco;
		 sy = (float) tri.sy*tri.scalco;
		 gy = (float) tri.gy*tri.scalco;
		 *offset=(float) tri.offset*tri.scalco;
		}
	else 
		{ 
		 /* if tri.scalco is negative divide */
		 *sx = (float) tri.sx/ABS(tri.scalco);
		 *gx = (float) tri.gx/ABS(tri.scalco);
		 sy = (float) tri.sy/ABS(tri.scalco);
		 gy = (float) tri.gy/ABS(tri.scalco);
		 *offset=(float) tri.offset/ABS(tri.scalco);
		}
	} 
 else 
	{
	 *sx = (float) tri.sx;
	 *gx = (float) tri.gx;
	 sy = (float) tri.sy;
	 gy = (float) tri.gy;
	 *offset=(float) tri.offset;
	}
 if(tri.sy||tri.gy)
	{
	 /* use pythagorean theorem to remap radial direction to x-direction */
	 *sx = SGN(*sx-sy)*sqrt((*sx)*(*sx) + sy*sy);
	 *gx = SGN(*gx-gy)*sqrt((*gx)*(*gx) + gy*gy);
	}
 return;
}

