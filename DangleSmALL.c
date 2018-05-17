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
void get_sx_gx_offset(int *sx, int *gx,int *offset);
void windtr(int nt,float *rtx,float ttt,float dt,float *firstt,float *datal);
void tanda(int ipx,int *anapx,int sx,int gx,float v,float vtt,float *ttt,float *qtmp);
void aperture(int it,int icdp,int mincdp,int maxcdp,int *bgc,int *edc,int anapxdx,float vt,float vtt,float h,float *angxx );
void angleintsmt(int nt,float dt,float tm1,float tm2,float tm3,float tm4,float ang1,float ang2,float ang3,float ang4,float *angx);
void hammingFilter(int nf1,int nf2,int nf3,int nf4,int nf, float *filter);
void calculateangle(float ximg,float angrange,float angdx,float angdxx,float h,float vt,float vtt,float *thit,int *nthit);
void wxd(int nt,int nx,int fd,int hwid, int hnx,float dt,float angdxx,float avg,float aratio,float *angarr,float *stk1,float **datact,float *minkjj,float *maxkjj);
#define LOOKFAC 2
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
 float *angx;					/* array of given aperture angle 		*/
 float *angxx;					/* array of given aperture angle 		*/
 int *anapx;					/* coordinate of cdp for image space	*/
 float datal[8];				/* small cut of ununiform sample trace	*/ 
 float **data;					/* temp_array of 2D real trace			*/
 float **vel=NULL;              /* array for storing velocity           */
 float ***aglerst=NULL;			/* array for storing angle gather		*/
 float **minkj;
 float **maxkj;
 float *minkjj;
 float *maxkjj;
 float *stk1;
 float *angarr;
 float ttt;						/* Travel time							*/
 float qtmp;					/* Amptitude of migration				*/	
 float va;    					/* ununiform sampled value				*/ 
 float firstt;					/* the first sample in datal[8]			*/
 float f1,f2,f3,f4;				/* array of filter frequencies          */   
 int sx,gx;					/* coordinate of shot and geophone  	*/
 int offset,h;				/* offset and half of that				*/
 float cdp;						/* coordinate of cdp					*/
 int tritvl;					/* trace interval						*/
 float dt;						/* sample interval						*/
 float hdt;
 float thit;
 float T;						/* vertical time depth 					*/
 float p;						/* Boundary attenuation factor of amp	*/
 float df,dw;					/* freqency sample spacing				*/
 float tmax;					/* the max trace length					*/
 float v;						/* velocity								*/
 float vt;						/* vt=v*T								*/
 float vtt;						/* vtt=vt*vt							*/
 float tm1,tm2,tm3,tm4;			/* time corespond to apreture angle		*/
 int ang1,ang2,ang3,ang4;		/* given migration aperture angle		*/
 int anapxmin;				/* min coordinate of imaging point		*/ 
 int anapxmax;				/* max coordinate of imaging point		*/
 int anapxdx;					/* spacing between imaging point		*/
 int startmt;					/* the start time for migration			*/
 int endmt;					/* the end time for migration			*/
 int nstartmt;					/* nstartmt=startmt/dT					*/
 int nendmt;					/* nendmt=endmt/D=dT					*/
 int ximg;					/* distance between cdp and image point	*/
 float angrange;				/* angle range of angle gather			*/
 float angdx;					/* angle interval to image in 			*/
 float angdxx;					/* angdxx=1.0/angdx						*/
 float aratio;
 float avg=0;
 int nthit;						/* angle trace number calculated		*/
 int itr,ix,ipx,it; 			/* count number                         */
 int icdp;						/* count number                         */
 int nf1,nf2,nf3,nf4;           /* nf1=(int)(f1/df)                     */
 int napmin;					/* napmin=(int)(anapxmin/anapxdx)		*/
 int oldcdp=0;	            	/* for temporary storage		        */
 int deltacdp;
 int bgc,edc;					/* begin and end of imaging	trace		*/		
 int mincdp,maxcdp;				/* mincdp and maxcdp of data input		*/
 int ncdpx=0;
 int mincdpo;					/* storing the mincdp for temporary to write header*/
 int firstcdp;	            	/* first cdp in velocity file	    	*/
 int lastcdp;	                /* last cdp in velocity file	    	*/
 int ncdp=0;	                /* number of cdps in the velocity file	*/ 
 int dcdp=0;	                /* number of cdps between consecutive traces */
 int nt;                		/* number of points on input trace      */
 int ntr;						/* number of trace input				*/
 int npx;						/* number of trace of imaging sapce		*/
 int nfft;                      /* number of points for fft trace       */
 int nf;                        /* number of frequencies (incl Nyq)     */
 int widlth;
 int fd;
 int hwid;
 int nx;
 int hnx;
 int hcdp;
 int verbose;		            /* flag to get advisory messages	    */

/* file name */
 char str[100],*path;
 char *vfile="";
 char *parfile="";
 FILE *fp;
 FILE *fp1;
 FILE *fp2;
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
 if(!getparstring("parfile",&parfile))	err("parameter file must be specified !");
 sprintf(str,"%s/par/%s",path,parfile);
 warn("str=%s",str);
 fp=fopen(str,"rb");
 fscanf(fp,"mincdp=%d\n",&mincdp);
 fscanf(fp,"maxcdp=%d\n",&maxcdp);
 fscanf(fp,"firstcdp=%d\n",&firstcdp);
 fscanf(fp,"lastcdp=%d\n",&lastcdp);
 fscanf(fp,"startmt=%d\n",&startmt);
 fscanf(fp,"endmt=%d\n",&endmt);
 fscanf(fp,"tm1=%f\n",&tm1);
 fscanf(fp,"tm2=%f\n",&tm2);
 fscanf(fp,"tm3=%f\n",&tm3);
 fscanf(fp,"tm4=%f\n",&tm4);
 fscanf(fp,"ang1=%d\n",&ang1);
 fscanf(fp,"ang2=%d\n",&ang2);
 fscanf(fp,"ang3=%d\n",&ang3);
 fscanf(fp,"ang4=%d\n",&ang4);
 fscanf(fp,"f1=%f\n",&f1);
 fscanf(fp,"f2=%f\n",&f2);
 fscanf(fp,"f3=%f\n",&f3);
 fscanf(fp,"f4=%f\n",&f4);
  fscanf(fp,"tritvl=%d\n",&tritvl);
 fscanf(fp,"anapxmin=%d\n",&anapxmin);
 fscanf(fp,"anapxmax=%d\n",&anapxmax);
 fscanf(fp,"anapxdx=%d\n",&anapxdx);
 fscanf(fp,"widlth=%d\n",&widlth);
 fscanf(fp,"fd=%d\n",&fd);
 fscanf(fp,"angrange=%f\n",&angrange);
 fscanf(fp,"angdx=%f\n",&angdx);
 fscanf(fp,"aratio=%f\n",&aratio);
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
 anapx=ealloc1int(npx);
 anapx[0]=anapxmin;
 for(ipx=1;ipx<npx;ipx++)
 	anapx[ipx]=anapx[ipx-1]+anapxdx;
 angdxx=1.0/angdx;
 hnx=ceil(angrange*angdxx);
 nx=hnx+hnx+1;
 angrange=angrange+angdx*0.5;
 mincdpo=mincdp;
 mincdp=mincdp-napmin;
 maxcdp=maxcdp-napmin;
 nstartmt=ceil(0.001*startmt/dt);
 nendmt=floor(0.001*endmt/dt);
 ncdpx=maxcdp-mincdp+1;
 hwid=(widlth-1)*0.5;
/* Store traces in tmpfile while getting a count of number of traces */
 tracefp = etmpfile();
 hfp = etmpfile();
 ntr = 0;
 do 
	{
	 ++ntr;

	 /* get new deltacdp value */
	 if(ntr==2)
	 	deltacdp=tri.cdp-oldcdp;

	 /* read headers and data */
	 efwrite(&tri,HDRBYTES, 1, hfp);
	 efwrite(tri.data, FSIZE, nt, tracefp);

	oldcdp=tri.cdp;
	} while (gettr(&tri));
	warn("ntr=%d",ntr);
/* get last cdp  and dcdp */
 if (!getparint("dcdp",&dcdp))	dcdp=deltacdp;
 if(verbose)	warn("dcdp=%d",dcdp);
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
 minkj=ealloc2float(nt,ncdpx);
 maxkj=ealloc2float(nt,ncdpx);
 minkjj=ealloc1float(nt);
 maxkjj=ealloc1float(nt);
 stk1=ealloc1float(nt);
 angarr=ealloc1float(nt);
 aglerst=ealloc3float(nt,nx,npx);

/* Zero all arrays */
 memset((void *) data[0], 0,nt*ntr*FSIZE);
 memset((void *) aglerst[0][0], 0,nt*nx*npx*FSIZE);
 memset((void *) minkj[0], 0,nt*ncdpx*FSIZE);
 memset((void *) maxkj[0], 0,nt*ncdpx*FSIZE);
 memset((void *) angarr, 0,nt*FSIZE);
 memset((void *) angx, 0, nt*FSIZE);
 memset((void *) angxx, 0, nt*FSIZE);
 memset((void *) rt, 0, nfft*FSIZE);
 memset((void *) filter, 0, nf*FSIZE);
 
/* calculate aperture of migration */
 angleintsmt(nt,dt,tm1,tm2,tm3,tm4,ang1,ang2,ang3,ang4,angx);
 for(it=0;it<nt;it++)
	angxx[it]=angx[it]*1.25;

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
	 efread(&tri,HDRBYTES, 1, hfp);
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
		 aglerst[icdp][hnx][it]+=va;
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
			 aglerst[ipx][hnx-nthit][it]+=va*qtmp*p;
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
             aglerst[ipx][hnx+nthit][it]+=va*qtmp*p;
			}
		}
	}
 hcdp=(int)(mincdp+maxcdp)*0.5;
 memset((void *) vel[0], 0,nt*ncdp*FSIZE);
 for(it=0;it<nt;it++)
	for(ix=0;ix<nx;ix++)
		stk1[it]+=aglerst[hcdp][ix][it];
/* Get the maxmum of every trace */
 for(it=1;it<nt;it++)
 	avg=avg>stk1[it]?avg:stk1[it];
 /* loop */
 for(itr=hcdp;itr<=maxcdp;itr++)
	{
	 memset((void *) stk1, 0,nt*FSIZE);
	 memset((void *) minkjj, 0,nt*FSIZE);
	 memset((void *) maxkjj, 0,nt*FSIZE);
	 wxd(nt,nx,fd,hwid,hnx,dt,angdxx,avg,aratio,angarr,stk1,aglerst[itr],minkjj,maxkjj);
	 for(it=0;it<nt;it++)
	 	{
		 vel[itr-mincdp][it]=stk1[it];
		 minkj[itr-mincdp][it]=minkjj[it];
		 maxkj[itr-mincdp][it]=maxkjj[it];
		}
	}
 for(itr=hcdp-1;itr>=mincdp;itr--)
	{
	 memset((void *) stk1, 0,nt*FSIZE);
	 memset((void *) minkjj, 0,nt*FSIZE);
	 memset((void *) maxkjj, 0,nt*FSIZE);
	 wxd(nt,nx,fd,hwid,hnx,dt,angdxx,avg,aratio,angarr,stk1,aglerst[itr],minkjj,maxkjj);
	 for(it=0;it<nt;it++)
		{
	 	vel[itr-mincdp][it]=stk1[it];
		minkj[itr-mincdp][it]=minkjj[it];
		maxkj[itr-mincdp][it]=maxkjj[it];
		}
	}
 fp1=fopen("kjmin","wb");
 fwrite(minkj[0],sizeof(float),nt*ncdpx,fp1);
 fclose(fp1);
 fp2=fopen("kjmax","wb");
 fwrite(maxkj[0],sizeof(float),nt*ncdpx,fp2);
 fclose(fp2);
/* set header */
 memset ((void *) &tro, (int) '\0', sizeof (tro));
 tro.trid = 1;
 tro.counit = 1;
 tro.f2=mincdpo*anapxdx;
 tro.d2 = anapxdx;
 tro.ns = nt;
 tro.dt = dt*1000000;
for(itr=0;itr<ncdpx;itr++)
	{
	 tro.cdp=mincdpo;
	 memcpy ((void *)tro.data, (const void *) vel[itr], sizeof(float)*nt);
	 puttr(&tro);
	 mincdpo+=dcdp;
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
 free1int(anapx);
 free2float(vel);
 free1float(filter);
 free1complex(ct);
 free1complex(hd);
 free2float(data);
 free2float(minkj);
 free2float(maxkj);
 free1float(minkjj);
 free1float(maxkjj);
 free1float(angarr);
 free1float(stk1);
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
        
void aperture(int it,int icdp,int mincdp,int maxcdp,int *bgc,int *edc,int anapxdx,float vt,float vtt,float h,float *angxx )
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

void get_sx_gx_offset(int *sx, int *gx,int *offset)
{ 
  /*****************************************************************************
 get_sx_gx_offset - get sx,gx and offset from headrs
  *****************************************************************************/
  int sy;		/* source coordinates */
  int gy;		/* geophone coordinates */
 if (tri.scalco) 
	{ 
	 /* if tri.scalco is set, apply value */
	 if (tri.scalco>0) 
		{	
		 *sx = (int) tri.sx*tri.scalco;
		 *gx = (int) tri.gx*tri.scalco;
		 sy = (int) tri.sy*tri.scalco;
		 gy = (int) tri.gy*tri.scalco;
		 *offset=(int) tri.offset*tri.scalco;
		}
	else 
		{ 
		 /* if tri.scalco is negative divide */
		 *sx = (int) tri.sx/ABS(tri.scalco);
		 *gx = (int) tri.gx/ABS(tri.scalco);
		 sy = (int) tri.sy/ABS(tri.scalco);
		 gy = (int) tri.gy/ABS(tri.scalco);
		 *offset=(int) tri.offset/ABS(tri.scalco);
		}
	} 
 else 
	{
	 *sx = (int) tri.sx;
	 *gx = (int) tri.gx;
	 sy = (int) tri.sy;
	 gy = (int) tri.gy;
	 *offset=(int) tri.offset;
	}
 if(tri.sy||tri.gy)
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

void wxd(int nt,int nx,int fd,int hwid, int hnx,float dt,float angdxx,float avg,float aratio,float *angarr,float *stk1,float **datact,float *minkjj,float *maxkjj)
{
 float *stk2;				
 float *itint;
 float *itout;
 float *angs;
 float *tmps;
 float *tmpx;
 float *Atmpcut;
 float atmp1;
 float atmp2;
 float abg;
 float aed;
 float T0;
 float Tofd;
 float agtmp;
 int widlthx;
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
 tmps=ealloc1float(nt);
 tmpx=ealloc1float(nt);
 angs=ealloc1float(nt);
 Atmpcut=ealloc1float(widlthx);
 memset((void *) tmps, 0,nt*FSIZE);
 memset((void *) tmpx, 0,nt*FSIZE);
 memset((void *) itint, 0,nt*FSIZE);
 memset((void *) itout, 0,nt*FSIZE);
 memset((void *) angs, 0,nt*FSIZE);

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
		 	agtmp=(itr-hnx)*0.017453/angdxx;
		 	abg=(Tofd*sin(agtmp)-sqrt(Tofd+Tofd-1))/cos(agtmp)/(Tofd-1);
		 	aed=(Tofd*sin(agtmp)+sqrt(Tofd+Tofd-1))/cos(agtmp)/(Tofd-1);
	 	 	iabg=(int)(atan(abg)*57.29578*angdxx);
	 	 	iaed=(int)(atan(aed)*57.29578*angdxx);
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
				tmps[k]=itr;
				}
			}
		itint[k]=it*dt;
		//warn("tmpk=%f",(tmps[k]-hnx)/angdxx);
		k++;
		}
 	itint[0]=0;
 	itint[k]=nt*dt;
 	tmps[k]=tmps[0]=hnx;
	k++;
 	float ydin[k][4];
 	memset((void *) ydin[0], 0, k*4*FSIZE);
 	for(it=0;it<nt;it++)
    	itout[it]=it*dt;
 	cakima(k,itint,tmps,ydin);
 	intcub(0,k,itint,ydin,nt,itout,angs);
 	//intlin(k, itint, tmps, tmps[0], tmps[k-1], nt, itout, angs);
	/*--------------------------------------*/
	for(inf=0;inf<mnf;inf++)
		{
	 	it=tmpx[inf];
		for(itx=it-hwid;itx<it+hwid;itx++)
			angs[itx]=tmps[inf+1];
		}
	/*-------------------------------------*/
	for(it=0;it<nt;it++)
 	 	angarr[it]=angs[it];
	}
else
	{
	 for(it=0;it<nt;it++)
 	 	angs[it]=angarr[it];
	}	
/* stack along Frenel */
 memset((void *) stk1, 0,nt*FSIZE);
 memset((void *) stk2, 0,nt*FSIZE);
 /*for(it=0;it<nt;it++)
	{
	 for(itr=angs[it];itr<nx;itr++)		
		{
		 stk1[it]+=datact[itr][it];
		 if(fabs(datact[itr][it]/stk1[it])<0.03)
			{	
			 maxkjj[it]=(itr-hnx)/angdxx;
		 	 break;
			}
		}
	 for(itr=angs[it];itr>0;itr--)
		{
		 stk2[it]+=datact[itr][it];
		 if(fabs(datact[itr][it]/stk2[it])<0.03)
			{	
			 minkjj[it]=(itr-hnx)/angdxx;
		 	 break;
			}
		}
	 stk1[it]=stk1[it]+stk2[it];
	}*/

 /*memset((void *) stk1, 0,nt*FSIZE);
 memset((void *) stk2, 0,nt*FSIZE);
 for(it=0;it<nt;it++)
	{
	 //maxkjj[it]=angs[it];
	 for(itr=angs[it];itr<nx;itr++)		
		{
		 atmp1=stk1[it];
		 stk1[it]+=datact[itr][it];
		 if(fabs(atmp1)>=fabs(stk1[it]))
			{	
			 maxkjj[it]=(itr-hnx)/angdxx;
		 	 break;
			}
		}
	 for(itr=angs[it];itr>0;itr--)
		{
		 atmp1=stk2[it];
		 stk2[it]+=datact[itr][it];
		 if(fabs(atmp1)>=fabs(stk2[it]))
			{	
			 minkjj[it]=(itr-hnx)/angdxx;
		 	 break;
			}
		}
	 stk1[it]=stk1[it]+stk2[it];
	}*/

for(it=0;it<nt;it++)
	{
	 T0=it*dt;
	 Tofd=T0*fd*2;
	 if(Tofd<=1)	continue;
	 agtmp=(angs[it]-hnx)*0.017453/angdxx;
	 abg=(Tofd*sin(agtmp)-sqrt(Tofd+Tofd-1))/cos(agtmp)/(Tofd-1);
	 aed=(Tofd*sin(agtmp)+sqrt(Tofd+Tofd-1))/cos(agtmp)/(Tofd-1);
	 minkjj[it]=atan(abg)*57.29578;
	 maxkjj[it]=atan(aed)*57.29578;
	}
 memset((void *) stk1, 0,nt*FSIZE);
 for(it=0;it<nt;it++)
	{
	atmp1=minkjj[it]*angdxx+hnx;
	atmp2=maxkjj[it]*angdxx+hnx;
	for(itr=atmp1;itr<=atmp2;itr++)		
		stk1[it]+=datact[itr][it];
	}

/*free array*/
free1float(stk2);
free1float(itint);
free1float(itout);
free1float(tmps);
free1float(tmpx);
free1float(Atmpcut);
free1float(angs);
}
