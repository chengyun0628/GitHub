/* CYFDMOD12: $Revision: 1.0 $ ; $Date: 2018/06/12 16:56:23 $        */

/* CYFDMOD12 - finite difference modeling by 12nd order explict method */

#include "par.h"
#include "su.h"
#include "segy.h"

/*********************** self documentation **********************/
char *sdoc[] = {
" 									",
" CYFDMOD12 - Finite-Difference MODeling (12nd order) for acoustic wave equation",
" Required Parameters:							",
" <vfile		file containing velocity[nx][nz]		",
" >wfile		file containing waves[nx][nz] for time steps	",
" nx=			number of x samples (2nd dimension)		",
" dx=			x sampling interval						",
" nz=			number of z samples (1st dimension)		",
" dz=			z sampling interval						",
" tmax=			maximum time					",
" 									",
" Optional Parameters:							",
" shot_num=1			total shot numbers	",  
" shot_interval=2		interval of each shot (grid sample)					",
" shot_startx=0			start position of shot in X direction (grid sample) ",
" shot_startz=0			start position of shot in Z direction (grid sample)	", 
" receiver_numx=nx		total receiver numbers in X direction 	", 
" receiver_numz=nz		total receiver numbers in Y direction 	",
" receiver_interval=1	interval of each receiver		",
" receiver_off=0		is receiver position fixed ?	",
" receiver_startx=0 	start position of receiver for HSD (grid sample)	",
" receiver_startz=0 	start position of receiver for VSP (grid sample)	",
"							",
" sstrength=1.0			strength of source				",
" mono=0				use ricker wavelet as source function 		",
" 		=1  			use single frequency src (freq=2*fpeak)	",
" 									",
" nt=1+tmax/dt			number of time samples (dt determined for stability)",
" fmax = vmin/(4.0*h)	maximum frequency in source 			",
" fpeak=0.5*fmax		peak frequency in ricker wavelet		",
" 									",
" pml_thick=100			PML boundary thickness in the up and down direction ",
" pml_width=100			PML boundary thickness in the left and right direction ",
" atten=50				PML boundary attenuation coefficient  ",
" movie=0				don't output wave field in time steps  ",
"       =1				output wave field in time steps  		",
" mt=2					number of time steps (dt) per output time step	",
" 									",
" hs1=			z coordinate of horizontal line of seismograms	",
" vs2=			x coordinate of vertical line of seismograms	",
" vsfile=		output file for vertical line of seismograms[nz][nt]",
" hsfile=		output file for horizontal line of seismograms[nx][nt]",
" ssfile=		output file for source point seismograms[nt]	",
" verbose=0		=1 for diagnostic messages, =2 for more		"
" Notes:								",
" 									",
" This program uses explicit 12 order differencing method",
" 									",
NULL};

/*
 * Authors:  NEPU:Cheng Yun
 *           

 * Trace header fields set: sx, gx, ns, delrt, tracl, tracr, offset, d1, d2,
 *                          sdepth, trid
 *
 * Technical reference:
 *
 * [1] 王维红, 柯璇, 裴江云. 完全匹配层吸收边界条件应用研究. 地球物理学进展, 2013, 28(5):2508-2514.
 * [2] sufdmod2.c by CWP:Dave Hale.
 *
 */
/**************** end self doc ********************************/

/* Prototypes for PML absorbing boundary conditions */
void mk_vpml(float **vpml,float **dvv,int du,int dd,int dl,int dr,int nx,int nz);
void mk_att(float **att1,int atten,int du,int dd,int dl,int dr,int nxx,int nzz);

/* Prototypes for finite differencing */
static float ricker (float t, float fpeak, int mono);
void ptsrc (float sstrength, int xs, int zs, float t, float fpeak, float tdelay, int mono, float **p);
void mk_coe(double *coe);
void extra (double *coe,int nxx, int nzz, float **att1,float **att2,float **att3,float **att4,float **pm, float **p, float **pp);

/* Globals for trace manipulation */
segy cubetr; 	/* data cube traces */
segy srctr;		/* source seismogram traces */
segy horiztr;	/* horizontal line seismogram traces */
segy verttr;	/* vertical line seismogram traces */

int
main(int argc, char **argv)
	{
	int ix,iz,it,is;		/* counters */
	int nx,nz,nt,mt;			/* x,z,t tsizes*/
	int du,dd,dl,dr;		/* PML boundary*/
	int pml_thick,pml_width;/* PML boundary thickness and width */
	int nxx,nzz;			/* nxx=nx+du+dd, nzz...*/
	int atten;				/* PML attenuation coefficent */

	int verbose;			/* is verbose? */
	int movie;				/* is movie? */
	int shot_num;			/* total shot numbers */ 
	int shot_startx;		/* start position of shot in X direction (grid sample) */
	int shot_startz;		/* start position of shot in Z direction (grid sample) */
	int shot_interval;		/* interval of each shot	*/
	int receiver_numx;		/* total receiver numbers in X direction */  
	int receiver_numz;		/* total receiver numbers in Y direction */  
	int receiver_startx; 	/* start position of receiver for HSD (grid sample) */
	int receiver_startz; 	/* start position of receiver for VSP (grid sample) */
	int receiver_endx;
	int receiver_endz;
	int receiver_interval;	/* interval of each receiver	*/
	int receiver_off;		/* is receiver position fixed	*/
	int vs2;				/* depth in samples of horiz rec line */
	int hs1;				/* horiz sample of vert rec line */

	float dx;				/* x sample interval */

	float dz;				/* z sample interval */
	float h;				/* minumum spatial sample interval */

	float dt;				/* time sample interval */
	float fmax;				/* maximum temporal frequency allowable */
	float fpeak;			/* peak frequency of ricker wavelet */
	float tdelay=0.0;		/* time delay of source beginning */

	float vmin;				/* minimum wavespeed in vfile */
	float vmax;				/* maximum wavespeed in vfile */

	float sstrength=0.0;	/* source strength */

	float tmax;				/* maximum time to compute */
	float t;				/* time */
	float **dvv=NULL;		/* array of velocity values from vfile */
	float **vpml=NULL;		/* array of velocity values for vpml */
	float **att1=NULL;		/* array of attenuation factor */
	float **att2=NULL;		/* temp array of attenuation factor */
	float **att3=NULL;		/* temp array of attenuation factor */
	float **att4=NULL;		/* temp array of attenuation factor */


	/* pressure field arrays */
	float **pm=NULL;		/* pressure field at t-dt */
	float **p=NULL;			/* pressure field at t */
	float **pp=NULL;		/* pressure field at t+dt */
	float **ptemp=NULL;		/* temp pressure array */

	/* output data arrays */
	float **ss=NULL;		/* source point seismogram array */
	float **hs=NULL;		/* seismograms from horiz receiver line */
	float **vs=NULL;		/* seismograms from vert receiver line */

	/* file names */
	char *vsfile="";		/* vert receiver seismogram line file  name */
	char *hsfile="";		/* horiz receiver seismogram line file name */
	char *ssfile="";		/* source point seismogram file name */

	/* input file pointers */
	FILE *velocityfp=stdin; /* pointer to input velocity data */

	/* output file pointers */
	FILE *hseisfp=NULL;		/* pointer to output horiz rec line file  */
	FILE *vseisfp=NULL;		/* pointer to output vert rec line file  */	
	FILE *sseisfp=NULL;		/* pointer to source point seis output file */

	/* SEGY fields */
	long tracl=0;			/* trace number within a line */
	long tracr=0;			/* trace number within a reel */

	int mono;				/* single fequency src flag */
	double *coe;
     	coe=alloc1double(7);

	/* hook up getpar to handle the parameters */
	initargs(argc,argv);
	requestdoc(0);
	
	/* get required parameters */
	if (!getparint("nx",&nx)) 	err("must specify nx!\n");
	if (!getparfloat("dx",&dx)) err("must specify dx!\n");
	if (!getparint("nz",&nz)) 	err("must specify nz!\n");
	if (!getparfloat("dz",&dz)) err("must specify dz!\n");
	if (!getparfloat("tmax",&tmax)) err("must specify tmax!\n");
	
	/* get optional parameters */
	if(!getparint("shot_num",&shot_num))	shot_num=1;
	if(!getparint("shot_interval",&shot_interval))	shot_interval=1;
	if(!getparint("shot_startx",&shot_startx))	shot_startx=0;
	if(!getparint("shot_startz",&shot_startz))	shot_startz=0;
	if(!getparint("receiver_numx",&receiver_numx))	receiver_numx=nx;
	if(!getparint("receiver_numz",&receiver_numz))	receiver_numz=nz;
	if(!getparint("receiver_interval",&receiver_interval))	receiver_interval=1;
	if(!getparint("receiver_off",&receiver_off))	receiver_off=0;
	if(!getparint("receiver_startx",&receiver_startx))	receiver_startx=0;
	if(!getparint("receiver_startz",&receiver_startz))	receiver_startz=0;
	if (!getparfloat("sstrength",&sstrength))	sstrength =1.0;
	if (!getparint("mono",&mono)) mono = 0;
	if (!getparint("nt",&nt)) nt = 0;
	if(!getparint("pml_thick",&pml_thick))	du=dd=100;
		else	du=dd=pml_thick;
	if(!getparint("pml_width",&pml_width))	dl=dr=100;
		else	dl=dr=pml_width;
	if(!getparint("atten",&atten))	atten=50;
	if(!getparint("movie",&movie))	movie=0;
	if (!getparint("mt",&mt)) mt = 2;
	if (!getparint("hs1",&hs1)) hs1 = 0;
	if (!getparint("vs2",&vs2)) vs2 = 0;
	getparstring("hsfile",&hsfile);
	getparstring("vsfile",&vsfile);
	getparstring("ssfile",&ssfile);

	if (!getparint("verbose",&verbose)) verbose = 0;

/*calculate the range of PML zone and in which the relative coordinate of shot && geophone */
	nxx=nx+dl+dr;
	nzz=nz+du+dd;
	shot_startx+=dl;
	shot_startz+=du;
	receiver_startx+=dl;
	receiver_startz+=du;
	hs1+=du;
	vs2+=dl;
/* print the input information */
	warn("shot_num=%d",shot_num);
	warn("shot_interval=%d",shot_interval);
	warn("shot_startx=%d",shot_startx-dl);
	warn("shot_startz=%d",shot_startz-du);
	warn("receiver_numx=%d",receiver_numx);
	warn("receiver_numz=%d",receiver_numz);
	warn("receiver_interval=%d",receiver_interval);
	warn("receiver_off=%d",receiver_off);
	warn("receiver_startx=%d",receiver_startx-dl);
	warn("receiver_startz=%d",receiver_startz-du);
	warn("pml_thick=%d",du);
	warn("pml_width=%d",dl);
	warn("pml_atten=%d",atten);
	
	/* allocate space */
	dvv = alloc2float(nz,nx);
	vpml= alloc2float(nzz,nxx);
	att1 = alloc2float(nzz,nxx);
	att2 = alloc2float(nzz,nxx);
	att3 = alloc2float(nzz,nxx);
	att4 = alloc2float(nzz,nxx);
	pm = alloc2float(nzz,nxx);
	p = alloc2float(nzz,nxx);
	pp = alloc2float(nzz,nxx);
	
	/* read velocities */
	fread(dvv[0],sizeof(float),nx*nz,velocityfp);
	
	/* determine minimum and maximum velocities */
	vmin = vmax = dvv[0][0];
	for (ix=0; ix<nx; ++ix) 
		for (iz=0; iz<nz; ++iz) 
			{
			vmin = MIN(vmin,dvv[ix][iz]);
			if (verbose==1 && dvv[ix][iz]==0) 
				warn("v=0 at (ix,iz)=(%i,%i)",ix,iz);
			vmax = MAX(vmax,dvv[ix][iz]);
			}
	
	/* determine mininum spatial sampling interval */
	h = MIN(ABS(dx),ABS(dz));
	
	/* determine time sampling interval to ensure stability */
	dt = h/(2.0*vmax);
	warn("stable dt=%g",dt);
	
	/* determine maximum temporal frequency to avoid dispersion */
	if (!getparfloat("fmax", &fmax))	fmax = vmin/(4.0*h);

	/* compute or set peak frequency for ricker wavelet */
	if (!getparfloat("fpeak", &fpeak))	fpeak = 0.5*fmax;

	/* determine number of time steps required to reach maximum time */
	if (nt==0) nt = 1+tmax/dt;

	/*------make shot & receiver coordinate-------*/
	int *source_x_cord;
	source_x_cord=alloc1int(shot_num);
	source_x_cord[0]=shot_startx;
	for(is=1;is<shot_num;is++)
		source_x_cord[is]=source_x_cord[is-1]+shot_interval;

	int *receiver_x_cord;
	receiver_x_cord=alloc1int(shot_num);
	receiver_x_cord[0]=receiver_startx;
	for(is=1;is<shot_num;is++)
		receiver_x_cord[is]=receiver_x_cord[is-1]+receiver_off*shot_interval;

	/* if requested, open file and allocate space for seismograms */
	if (*hsfile!='\0') 
		{
		if((hseisfp=fopen(hsfile,"w"))==NULL)
			err("cannot open hsfile=%s\n",hsfile);
		hs = alloc2float(nt,nx);
		} 
	else 
		hs = NULL;

	if (*vsfile!='\0') 
		{
		if((vseisfp=fopen(vsfile,"w"))==NULL)
			err("cannot open vsfile=%s\n",vsfile);
		vs = alloc2float(nt,nz);
		} 
	else
		vs = NULL;

	if (*ssfile!='\0') 
		{
		if((sseisfp=fopen(ssfile,"w"))==NULL)
			err("cannot open ssfile=%s\n",ssfile);
		ss = alloc2float(nt,shot_num);
		} 
	else
		ss = NULL;
	/* compute velocity^2 */	
	for (ix=0; ix<nx; ++ix) 
		for (iz=0; iz<nz; ++iz) 
			dvv[ix][iz] = dvv[ix][iz]*dvv[ix][iz];

	/* if verbose, print parameters */
	if (verbose) 
		{
		fprintf(stderr,"nx = %d\n",nx);
		fprintf(stderr,"dx = %g\n",dx);
		fprintf(stderr,"nz = %d\n",nz);
		fprintf(stderr,"dz = %g\n",dz);
		fprintf(stderr,"nt = %d\n",nt);
		fprintf(stderr,"dt = %g\n",dt);
		fprintf(stderr,"tmax = %g\n",tmax);
		fprintf(stderr,"fmax = %g\n",fmax);
		fprintf(stderr,"fpeak = %g\n",fpeak);
		fprintf(stderr,"vmin = %g\n",vmin);
		fprintf(stderr,"vmax = %g\n",vmax);
		fprintf(stderr,"mt = %d\n",mt);
		}

	/* make coefficient for PML and extrapolation */
	memset((void *) att1[0], 0,nxx*nzz*FSIZE);
	mk_vpml(vpml,dvv,du,dd,dl,dr,nx,nz);
	mk_att(att1,atten,du,dd,dl,dr,nxx,nzz); 
	mk_coe(coe);

	/* make coefficient att1,att2,att3,att4,to simply accelerate the extrapolation */
	for(ix=0;ix<nxx;ix++)
		for(iz=0;iz<nzz;iz++)
			{
			att1[ix][iz]*=dt;
			att2[ix][iz]=(2-att1[ix][iz]*att1[ix][iz])/(1+att1[ix][iz]);
			att3[ix][iz]=(1-att1[ix][iz])/(1+att1[ix][iz]);
			att4[ix][iz]=vpml[ix][iz]*dt*dt/(dx*dx*(1+att1[ix][iz]));
			att1[ix][iz]=att4[ix][iz]*dz*dz/dx/dx;
			}
	/* loop over shot and time steps */
	for(is=0;is<shot_num;is++)
		{
		memset((void *) pm[0], 0,nzz*nxx*FSIZE);
		memset((void *) p[0],  0,nzz*nxx*FSIZE);
		memset((void *) pp[0], 0,nzz*nxx*FSIZE);
		for (it=0,t=0.0; it<nt; ++it,t+=dt) 
			{
			/* update source function */
			ptsrc(sstrength,source_x_cord[is],shot_startz,t,fpeak,tdelay,mono,p);
			/* do one time step */
			extra(coe,nxx,nzz,att1,att2,att3,att4,pm,p,pp);

			/* write waves */
			/* if (it%mt==0) fwrite(pp[0],sizeof(float),nx*nz,stdout); */
			if (movie==1&&is==0&&it%mt==0) 
				{
				cubetr.sx = (shot_startx-dl)*dx;
				cubetr.sdepth = (shot_startz-du)*dz;
				cubetr.trid = 30 ;
				cubetr.ns = nz ;
				cubetr.d1 = dz ;
				cubetr.d2 = dx ;
				/* account for delay in source starting time */
				cubetr.delrt = - 1000.0 * tdelay;

				tracl = 0 ;

				for (ix=dl ; ix < nxx-dr ; ++ix) 
					{
					++tracl;
					++tracr;

					cubetr.offset = (ix-shot_startx)*dx;
					cubetr.gx = (ix-dl) * dx ;
					cubetr.tracl = (int) tracl;
					cubetr.tracr = (int) tracr;

					for (iz=du ; iz < nzz-dd ; ++iz) 
						{	
						cubetr.data[iz-du] = pp[ix][iz];
						}
					fputtr(stdout, &cubetr);
					}
				}

			/* if requested, save horizontal line of seismograms */
			if (hs!=NULL) 
				{
				receiver_endx=MIN(nx+dl,receiver_x_cord[is]+receiver_numx*receiver_interval);
				for (ix=receiver_x_cord[is]; ix<receiver_endx; ix+=receiver_interval)
					hs[ix-dl][it] = pp[ix][hs1];
				}
			/* if requested, save vertical line of seismograms */
			if (vs!=NULL) 
				{
				receiver_endz=MIN(nz+du,receiver_startz+receiver_numz*receiver_interval);
				for (iz=receiver_startz; iz<receiver_endz; iz+=receiver_interval)
					vs[iz-du][it] = pp[vs2][iz];
				}

			/* if requested, save seismograms at source locations */
			if (ss!=NULL) 
				ss[is][it] = pp[source_x_cord[is]][shot_startz];

			/* roll time slice pointers */
			ptemp=pm;
			pm = p;
			p = pp;
			pp=ptemp;
			}

		/* if requested, write horizontal line of seismograms */
		if (hs!=NULL) 
			{
			horiztr.sx = (source_x_cord[is]-dl)*dx;
			horiztr.sdepth = (shot_startz-du)*dz;
			horiztr.trid = 1;
			horiztr.ns = nt ;
			horiztr.dt = 1000000 * dt ;
			horiztr.d2 = dx ;

			/* account for delay in source starting time */
			horiztr.delrt = -1000.0 * tdelay ; 

			tracl = tracr = 0;

			for (ix=receiver_x_cord[is] ; ix < receiver_endx ; ix+=receiver_interval)
				{
				++tracl;
				++tracr;

				/* offset from first source location */
				horiztr.offset = (ix - source_x_cord[is])*dx;
				horiztr.gx = (ix-dl) * dx;

				horiztr.tracl = (int) tracl;
				horiztr.tracr = (int) tracr;

				for (it = 0 ; it < nt ; ++it)
					horiztr.data[it] = hs[ix-dl][it];
			
				fputtr(hseisfp , &horiztr);
				}

			}	
		
		/* if requested, write vertical line of seismograms */
		if (vs!=NULL) 
			{
			verttr.trid = 1;
			verttr.ns = nt ;
			verttr.sx = (source_x_cord[is]-dl)*dx;
			verttr.sdepth = (shot_startz-du)*dz;
			verttr.dt = 1000000 * dt ;
			verttr.d2 = dz ;
			/* account for delay source starting time */
			verttr.delrt = -1000.0 * tdelay ;

			tracl = tracr = 0;
			for (iz=receiver_startz; iz<receiver_endz; iz+=receiver_interval)
				{
				++tracl;
				++tracr;

				/* vertical line implies offset in z */
				verttr.offset = (iz-shot_startz)*dz;

				verttr.tracl = (int) tracl;
				verttr.tracr = (int) tracr;

				for (it = 0 ; it < nt ; ++it)
					verttr.data[it] = vs[iz-du][it];
			
				fputtr(vseisfp , &verttr);
				}	
			}

		/* if requested, write seismogram at source position */
		if (ss!=NULL) 
			{
			srctr.trid = 1;
			srctr.ns = nt ;
			srctr.dt = 1000000 * dt ;
			srctr.d2 = 1 ;
			srctr.delrt = -1000.0 * tdelay ;

			tracl = tracr = is;

			srctr.sx = (source_x_cord[is]-dl)*dx;
			srctr.sdepth = (shot_startz-du)*dz;
			srctr.tracl = (int) tracl;
			srctr.tracr = (int) tracr;

			for (it = 0 ; it < nt ; ++it)
				srctr.data[it] = ss[is][it];
			fputtr(sseisfp , &srctr);	
			}

		}
	if (*hsfile!='\0') fclose(hseisfp);	
	if (*vsfile!='\0') fclose(vseisfp);
	if (*ssfile!='\0') fclose(sseisfp);

	/* free space before returning */
	free1int(source_x_cord);
	free1int(receiver_x_cord);
	free2float(dvv);
	free2float(vpml);
	free2float(att1);
	free2float(att2);
	free2float(att3);
	free2float(att4);
	free2float(pm);
	free2float(p);
	free2float(pp);
	
	if (hs!=NULL) free2float(hs);
	if (vs!=NULL) free2float(vs);
	if (ss!=NULL) free2float(ss);
	
	return(CWP_Exit());
	}

void ptsrc (float sstrength, int xs, int zs, float t, float fpeak, float tdelay, int mono, float **p)
/*****************************************************************************
update source pressure function for a point source
******************************************************************************
Input:
sstrength 	the pressure strength
xs			x coordinate of point source
zs			z coordinate of point source
t			time at which to compute source function
fpeak		peak frequency

Output:
tdelay		time delay of beginning of source function
p			array[nx][nz] of source pressure at time t+dt
******************************************************************************
Author:  Dave Hale, Colorado School of Mines, 03/01/90
******************************************************************************/
	{
	int ix,iz;
	float ts;
	
	/* compute time-dependent part of source function */
	tdelay = 1.0/fpeak;
	if (t>2.0*tdelay) return;
	ts = ricker(t-tdelay,fpeak,mono);
	p[xs][zs]=sstrength*ts;
	}

static float ricker (float t, float fpeak, int mono)
/*****************************************************************************
Compute Ricker wavelet as a function of time
******************************************************************************
Input:
t		time at which to evaluate Ricker wavelet
fpeak		peak (dominant) frequency of wavelet
mono		=0 use ricker... =1 use single frequency (2*fpeak)  CLL 11/27/06
******************************************************************************
Notes:
The amplitude of the Ricker wavelet at a frequency of 2.5*fpeak is 
approximately 4 percent of that at the dominant frequency fpeak.
The Ricker wavelet effectively begins at time t = -1.0/fpeak.  Therefore,
for practical purposes, a causal wavelet may be obtained by a time delay
of 1.0/fpeak.
The Ricker wavelet has the shape of the second derivative of a Gaussian.
******************************************************************************
Author:  Dave Hale, Colorado School of Mines, 04/29/90
******************************************************************************/
	{
	float x,xx;
	
	x = PI*fpeak*t;
	xx = x*x;
	if (mono==0) 
		return exp(-xx)*(1.0-2.0*xx);
	else 
		return sin(x);
	}

/* make 12 order coefficient */
void mk_coe(double *coe)
    { 
	coe[0]=-3.12108522;
	coe[1]=1.83730507;
	coe[2]=-0.35408741;
	coe[3]=0.09988277;
	coe[4]=-0.02817135;
	coe[5]=0.006539;
	coe[6]=-0.00092547;
    }

/* 2D finite differencing subroutine */
void extra (double *coe,int nxx, int nzz, float **att1,float **att2,float **att3,float **att4,float **pm, float **p, float **pp)
/*****************************************************************************
One time step of FD solution (12nd order in space) to acoustic wave equation
******************************************************************************
Input:
nxx		number of x samples
nzz		number of z samples
pm		array[nxx][nzz] of pressure at time t-dt
p		array[nxx][nzz] of pressure at time t

Output:
pp		array[nxx][nzz] of pressure at time t+dt

******************************************************************************/
	{
	int ix,iz;
	
	/* do the finite-difference star */
	for (ix=6; ix<nxx-6; ++ix) 
		for (iz=6; iz<nzz-6; ++iz) 
			{
			 pp[ix][iz] = att2[ix][iz]*p[ix][iz]-att3[ix][iz]*pm[ix][iz] +
			 att4[ix][iz]*(
			 coe[0]*p[ix][iz]+
			 coe[1]*(p[ix-1][iz]+p[ix+1][iz])+
			 coe[2]*(p[ix-2][iz]+p[ix+2][iz])+
			 coe[3]*(p[ix-3][iz]+p[ix+3][iz])+
			 coe[4]*(p[ix-4][iz]+p[ix+4][iz])+
			 coe[5]*(p[ix-5][iz]+p[ix+5][iz])+
			 coe[6]*(p[ix-6][iz]+p[ix+6][iz]))+
			 att1[ix][iz]*(
			 coe[0]*p[ix][iz]+
			 coe[1]*(p[ix][iz-1]+p[ix][iz+1])+
			 coe[2]*(p[ix][iz-2]+p[ix][iz+2])+
			 coe[3]*(p[ix][iz-3]+p[ix][iz+3])+
			 coe[4]*(p[ix][iz-4]+p[ix][iz+4])+
			 coe[5]*(p[ix][iz-5]+p[ix][iz+5])+
			 coe[6]*(p[ix][iz-6]+p[ix][iz+6]));
			}

	}


void mk_vpml(float **vpml,float **dvv,int du,int dd,int dl,int dr,int nx,int nz)
{
int ix,iz;
/* four boundary zone */
//up
for(iz=0;iz<du;iz++)
	for(ix=0;ix<nx;ix++)
		vpml[ix+dl][iz]=dvv[ix][0];
//down
for(iz=nz+du;iz<nz+du+dd;iz++)
	for(ix=0;ix<nx;ix++)
		vpml[ix+dl][iz]=dvv[ix][nz-1];
//left
for(ix=0;ix<dl;ix++)
	for(iz=0;iz<nz;iz++)
		vpml[ix][iz+du]=dvv[0][iz];
//rigt
for(ix=nx+dl;ix<nx+dl+dr;ix++)
	for(iz=0;iz<nz;iz++)
		vpml[ix][iz+du]=dvv[nx-1][iz];
/* four corner zone */

//ul
for(ix=0;ix<dl;ix++)			
	for(iz=0;iz<du;iz++)
		vpml[ix][iz]=dvv[0][0];

//ur
for(ix=nx+dl;ix<nx+dl+dr;ix++)
	for(iz=0;iz<du;iz++)			
		vpml[ix][iz]=dvv[nx-1][0];

//dl
for(ix=0;ix<dl;ix++)
	for(iz=nz+du;iz<nz+du+dd;iz++)			
		vpml[ix][iz]=dvv[0][nz-1];

//dr
for(ix=nx+dl-1;ix<nx+dl+dr;ix++)
	for(iz=nz+du;iz<nz+du+dd;iz++)				
		vpml[ix][iz]=dvv[nx-1][nz-1];
/* middle zone */
for(ix=dl;ix<nx+dl;ix++)
	for(iz=du;iz<nz+du;iz++)
		vpml[ix][iz]=dvv[ix-dl][iz-du];
	
}

void mk_att(float **att1,int atten,int du,int dd,int dl,int dr,int nxx,int nzz) 
{
	int iz,ix;
	float *att_u;
	float *att_l;
	att_u=alloc1float(du);
	att_l=alloc1float(dl);
/* up && down */	
	for(iz=0;iz<du;iz++)
		att_u[iz]=atten*cos(0.5*PI*iz/du);
/* left && right */
	for(ix=0;ix<dl;ix++)
		att_l[ix]=atten*cos(0.5*PI*ix/dl);
	
//up
for(iz=0;iz<du;iz++)
	for(ix=0;ix<nxx;ix++)
		att1[ix][iz]=att_u[iz];
//down
for(iz=nzz-du;iz<nzz;iz++)
	for(ix=0;ix<nxx;ix++)
		att1[ix][iz]=att_u[nzz-iz-1];
//left
for(ix=0;ix<dl;ix++)
	for(iz=0;iz<nzz;iz++)
		att1[ix][iz]+=att_l[ix];
//right
for(ix=nxx-dl;ix<nxx;ix++)
	for(iz=0;iz<nzz;iz++)
		att1[ix][iz]+=att_l[nxx-ix-1];

free1float(att_u);
free1float(att_l);


}
