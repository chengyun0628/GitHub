#include"time.h"
#include "su.h"
#include"cwp.h"
#include "segy.h"
#include "math.h"
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
void tanda(int ipx,float *anapx,float sx,float gx,float v,float vtt,float *ttt,float *qtmp);
void aperture(int it,int icdp,int mincdp,int maxcdp,int *bgc,int *edc,float anapxdx,float vt,float vtt,float h,float *angxx );
void angleintsmt(int nt,float dt,float tm1,float tm2,float tm3,float tm4,float ang1,float ang2,float ang3,float ang4,float *angx);
void hammingFilter(int nf1,int nf2,int nf3,int nf4,int nf, float *filter);
void calculateangle(float ximg,float angrange,float angdx,float angdxx,float h,float vt,float vtt,float *thit,int *nthit);
#define LOOKFAC 2
#define PFA_MAX 720720
segy tr;                        /* trace of input  						*/
segy tro;						/* trace of output 						*/
int
main(int argc, char **argv)
{
 register float *rt,*rtx;       /* real trace data and processed data	*/
 register complex *ct;          /* complex transformed trace            */
 register complex *hd;			/* half derivative						*/
 float *filter;              	/* filter array                         */
 float *angx;					/* array of given aperture angle 		*/
 float *angxx;					/* array of given aperture angle 		*/
 float *anapx;					/* coordinate of cdp for image space	*/
 float datal[8];				/* small cut of ununiform sample trace	*/ 
 float **data;					/* temp_array of 2D real trace			*/
 float **vel=NULL;              /* array for storing velocity           */
 float ***aglerst=NULL;			/* array for storing angle gather		*/
 float ttt;						/* Travel time							*/
 float qtmp;					/* Amptitude of migration				*/	
 float va;    					/* ununiform sampled value				*/ 
 float firstt;					/* the first sample in datal[8]			*/
 float f1,f2,f3,f4;				/* array of filter frequencies          */   
 float sx,gx;					/* coordinate of shot and geophone  	*/
 float offset,h;				/* offset and half of that				*/
 float cdp;						/* coordinate of cdp					*/
 float tritvl;					/* trace interval						*/
 float dt;						/* sample interval						*/
 float hdt;
 float thit;
 float tanp;             		/* tanp=pow(tan(ang),2); 				*/
 float T;						/* vertical time depth 					*/
 float p;						/* Boundary attenuation factor of amp	*/
 float df,dw;					/* freqency sample spacing				*/
 float tmax;					/* the max trace length					*/
 float v;						/* velocity								*/
 float vt;						/* vt=v*T								*/
 float vtt;						/* vtt=vt*vt							*/
 float ang1,ang2,ang3,ang4;		/* given migration aperture angle		*/
 float tm1,tm2,tm3,tm4;			/* time corespond to apreture angle		*/
 float anapxmin;				/* min coordinate of imaging point		*/ 
 float anapxmax;				/* max coordinate of imaging point		*/
 float anapxdx;					/* spacing between imaging point		*/
 float startmt;					/* the start time for migration			*/
 float endmt;					/* the end time for migration			*/
 int nstartmt;					/* nstartmt=startmt/dT					*/
 int nendmt;					/* nendmt=endmt/D=dT					*/
 float ximg;					/* distance between cdp and image point	*/
 float angrange;				/* angle range of angle gather			*/
 float angdx;					/* angle interval to image in 			*/
 float angdxx;					/* angdxx=1.0/angdx						*/
 int nthit;						/* angle trace number calculated		*/
 int hmaxnthit;					/* hmaxnthit=ceil(angrange/angdx)		*/
 int maxnthit;					/* maxnthit=2*hmaxnthit 				*/
 int iagle;						/* count number for angle				*/
 int itr,ix,ipx,it; 			/* count number                         */
 int icdp;						/* count number                         */
 int nf1,nf2,nf3,nf4;           /* nf1=(int)(f1/df)                     */
 int napmin;					/* napmin=(int)(anapxmin/anapxdx)		*/
 int oldcdp=0;	            	/* for temporary storage		        */
 int olddeltacdp=1;				/* for temporary storage                */
 int deltacdp;
 int bgc,edc;					/* begin and end of imaging	trace		*/		
 int mincdp,maxcdp;				/* mincdp and maxcdp of data input		*/
 int mincdpx;					/* storing the mincdp for temporary to write header*/
 int firstcdp;	            	/* first cdp in velocity file	    	*/
 int lastcdp;	                /* last cdp in velocity file	    	*/
 int ncdp=0;	                /* number of cdps in the velocity file	*/ 
 int dcdp=0;	                /* number of cdps between consecutive traces */
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
 hdt=0.5*dt;
 tmax=(nt-1)*dt;
/* Get parameter*/
 if(!getparstring("path",&path))	err("path must be specified !");
 if(!getparstring("vfile",&vfile))	err("velocity file must be specified !");
 if(!getparstring("parfile",&parfile))	err("parameter file must be specified !");
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
 fscanf(fp,"angrange=%f\n",&angrange);
 fscanf(fp,"angdx=%f\n",&angdx);
 efclose(fp);
 time(&t1);
 tm1=tm1/1000;
 tm2=tm2/1000;
 tm3=tm3/1000;
 tm4=tm4/1000;
 if (tm1>tmax||tm2>tmax||tm3>tmax||tm4>tmax)	err("Check tm values in data!");
 if (fmod(tritvl,anapxdx)!=0)	err("Check anapxdx value in parfile,cause that tritvl must be an integer number of anapxdx");

/*caculate image spacing*/
 napmin=(int)(anapxmin/anapxdx);
 npx=(int)(anapxmax/anapxdx)-napmin+1;
 anapx=ealloc1float(npx);
 anapx[0]=anapxmin;
 for(ipx=1;ipx<npx;ipx++)
 	anapx[ipx]=anapx[ipx-1]+anapxdx;
 angdxx=1.0/angdx;
 hmaxnthit=ceil(angrange*angdxx);
 maxnthit=hmaxnthit+hmaxnthit+1;
 angrange=angrange+angdx*0.5;
 mincdp=mincdp-napmin;
 maxcdp=maxcdp-napmin;
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
	 deltacdp=tr.cdp-oldcdp;

	 /* read headers and data */
	 efwrite(&tr,HDRBYTES, 1, hfp);
	 efwrite(tr.data, FSIZE, nt, tracefp);

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
	 oldcdp=tr.cdp;
	 olddeltacdp=deltacdp;
	} while (gettr(&tr));
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
 if (verbose)	warn("nfft=%d,nf=%d,df=%f,dw=%f",nfft,nf,df,dw);

/* Allocate space */
 angx=ealloc1float(nt);
 angxx=ealloc1float(nt);
 rt=ealloc1float(nfft);
 rtx=ealloc1float(nfft);
 ct=ealloc1complex(nf);
 hd=ealloc1complex(nf);
 filter=ealloc1float(nf);
 data=ealloc2float(nt,ntr);
 vel=ealloc2float(nt,ncdp);
 aglerst=ealloc3float(nt,maxnthit,npx);

/* Zero all arrays */
 memset((void *) data[0], 0,nt*ntr*FSIZE);
 memset((void *) aglerst[0][0], 0,nt*maxnthit*npx*FSIZE);
 memset((void *) angx, 0, nt*FSIZE);
 memset((void *) angxx, 0, nt*FSIZE);
 memset((void *) rt, 0, nfft*FSIZE);
 memset((void *) filter, 0, nf*FSIZE);
 
/* calculate aperture of migration */
 angleintsmt(nt,dt,tm1,tm2,tm3,tm4,ang1,ang2,ang3,ang4,angx);
 for(it=0;it<nt;it++)
	angxx[it]=angx[it]*1.25;

/*if (verbose)	
	{
	 FILE *fp1;
	 fp1=fopen("angx.bin","wb");
	 fwrite(angx,sizeof(float),nt,fp1);
	 fclose(fp1);
	} */

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

/* Start the migration process */
/* Loop over input trace */
 if (verbose)	warn("Starting migration process...\n");
 for(itr=0; itr<ntr; ++itr)
	{
	 efread(&tr,HDRBYTES, 1, hfp);
	 float perc;
	 perc=itr*100.0/(ntr-1);
	 if(fmod(itr*100.0,ntr-1)==0 && verbose)
		warn("migrated %g\n",perc);
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
 
	/* caculate the coordinate of shot and geophone*/
	 get_sx_gx_offset(&sx,&gx,&offset);
	 h=offset/2;	
     cdp=(sx+gx)/2;
     icdp=(int)(cdp/anapxdx)-napmin;
	/* loop in the image space*/
	 T=nstartmt*hdt;
	 for(it=nstartmt;it<nendmt;it++)
		{
		 T+=hdt;
		 v=vel[icdp-firstcdp+napmin][it];
		 vt=v*T;
		 vtt=vt*vt;
		 aperture(it,icdp,mincdp,maxcdp,&bgc,&edc,anapxdx,vt,vtt,h,angxx);
		 tanda(icdp,anapx,sx,gx,v,vtt,&ttt,&qtmp);
		 if(ttt>=tmax)   break;
		 windtr(nt,rtx,ttt,dt,&firstt,datal);
		 /* sinc interpolate new data */
		 ints8r(8, dt, firstt, datal,
			 0.0, 0.0, 1, &ttt, &va);
		 aglerst[icdp][hmaxnthit][it]+=va;
		 p=1.0;
		 ximg=0;
		 for(ipx=icdp-1;ipx>=bgc;ipx--)
		 	{
			 ximg=ximg+anapxdx;
			 v=vel[ipx-firstcdp+napmin][it];
			 vt=v*T;
			 vtt=vt*vt;
			 calculateangle(ximg,angrange,angdx,angdxx,h,vt,vtt,&thit,&nthit);
			 if(nthit==1000000)	break;
			 tanda(ipx,anapx,sx,gx,v,vtt,&ttt,&qtmp);
             if(ttt>=tmax)   break;
             windtr(nt,rtx,ttt,dt,&firstt,datal);
			 /* sinc interpolate new data */
         	 ints8r(8, dt, firstt, datal,
             	 0.0, 0.0, 1, &ttt, &va);
			 if(thit>angx[it])	p=cos(6.28318*(thit/angx[it]-1));
			 aglerst[ipx][hmaxnthit-nthit][it]+=va*qtmp*p;
			}
		 p=1.0;	
		 ximg=0;
		 for(ipx=icdp+1;ipx<=edc;ipx++)
			{
			 ximg=ximg+anapxdx;
			 v=vel[ipx-firstcdp+napmin][it];
			 vt=v*T;
			 vtt=vt*vt;
			 calculateangle(ximg,angrange,angdx,angdxx,h,vt,vtt,&thit,&nthit);
			 if(nthit==1000000)    break;
			 tanda(ipx,anapx,sx,gx,v,vtt,&ttt,&qtmp);
             if(ttt>=tmax)   break;
             windtr(nt,rtx,ttt,dt,&firstt,datal);
			 /* sinc interpolate new data */
             ints8r(8, dt, firstt, datal,
                 0.0, 0.0, 1, &ttt, &va);
			 if(thit>angx[it])	p=cos(6.28318*(thit/angx[it]-1));
             aglerst[ipx][hmaxnthit+nthit][it]+=va*qtmp*p;
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
 mincdpx=mincdp+napmin;
 for(ipx=mincdp; ipx<=maxcdp; ++ipx)
	{
	 for(iagle=0;iagle<maxnthit;iagle++)
    	{
     	 tro.cdp = mincdpx;
		 tro.offset=tro.f2+iagle*angdx;
     	 memcpy ((void *) tro.data, (const void *)  aglerst[ipx][iagle],sizeof (float) * nt);
     	 puttr(&tro);
		}
     mincdpx++;
	}
 time(&t2);
 warn("time consuming in second = %f\n",difftime(t2,t1));
/*free array*/
 efclose(hfp);
 efclose(tracefp);
 free1float(rt);
 free1float(rtx);
 free1float(angx);
 free1float(angxx);
 free1float(anapx);
 free2float(vel);
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
        
void aperture(int it,int icdp,int mincdp,int maxcdp,int *bgc,int *edc,float anapxdx,float vt,float vtt,float h,float *angxx )
{
 float x;
 float ang;
 float tanp;
 ang=angxx[it]*PI*0.005556;
 tanp=tan(ang)*tan(ang);
 x=(vt*(tanp-1)+sqrt(vtt*(1+tanp)*(1+tanp)+4*h*h*tanp))/2/tan(ang);
 *bgc=MAX(mincdp,icdp - ceil(x/anapxdx)); 
 *edc=MIN(maxcdp,icdp + ceil(x/anapxdx));
}

void tanda(int ipx,float *anapx,float sx,float gx,float v,float vtt,float *ttt,float *qtmp)
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

void get_sx_gx_offset(float *sx, float *gx,float *offset)
{ 
  /*****************************************************************************
 get_sx_gx_offset - get sx,gx and offset from headrs
  *****************************************************************************/
  float sy;		/* source coordinates */
  float gy;		/* geophone coordinates */
 if (tr.scalco) 
	{ 
	 /* if tr.scalco is set, apply value */
	 if (tr.scalco>0) 
		{	
		 *sx = (float) tr.sx*tr.scalco;
		 *gx = (float) tr.gx*tr.scalco;
		 sy = (float) tr.sy*tr.scalco;
		 gy = (float) tr.gy*tr.scalco;
		 *offset=(float) tr.offset*tr.scalco;
		}
	else 
		{ 
		 /* if tr.scalco is negative divide */
		 *sx = (float) tr.sx/ABS(tr.scalco);
		 *gx = (float) tr.gx/ABS(tr.scalco);
		 sy = (float) tr.sy/ABS(tr.scalco);
		 gy = (float) tr.gy/ABS(tr.scalco);
		 *offset=(float) tr.offset/ABS(tr.scalco);
		}
	} 
 else 
	{
	 *sx = (float) tr.sx;
	 *gx = (float) tr.gx;
	 sy = (float) tr.sy;
	 gy = (float) tr.gy;
	 *offset=(float) tr.offset;
	}
 if(tr.sy||tr.gy)
	{
	 /* use pythagorean theorem to remap radial direction to x-direction */
	 *sx = SGN(*sx-sy)*sqrt((*sx)*(*sx) + sy*sy);
	 *gx = SGN(*gx-gy)*sqrt((*gx)*(*gx) + gy*gy);
	}
 return;
}

void calculateangle(float ximg,float angrange,float angdx,float angdxx,float h,float vt,float vtt,float *thit,int *nthit)
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
