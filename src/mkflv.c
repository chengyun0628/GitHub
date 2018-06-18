/* MKFLV - Make the first layer's velocity for FBK modeling */

#include "su.h"

/*********************** self documentation **********************/
char *sdoc[] = {
" 									",
" mkflv - Make the first layer's velocity for FBK modeling  ",
" Required Parameters:							",
" nx=			number of x samples (2nd dimension)		",
" nz=			number of z samples (1st dimension)		",
" flv=			velocity for first layer				",
" 									",
" 									",
NULL};

/*
 * Authors:  NEPU:Cheng Yun
 *           
/**************** end self doc ********************************/

int
main(int argc, char **argv)
{
int ix,nx,iz,nz;
float flv;
float **vel;
FILE *vfp;

/* hook up getpar to handle the parameters */
initargs(argc,argv);
requestdoc(0);

if (!getparint("nx",&nx)) 	err("must specify nx!\n");
if (!getparint("nz",&nz)) 	err("must specify nz!\n");
if (!getparfloat("flv",&flv)) err("must specify flv!\n");
vel=alloc2float(nz,nx);

for(ix=0;ix<nx;ix++)
	for(iz=0;iz<nz;iz++)
		vel[ix][iz]=flv;
vfp=fopen("flvmdl","wb");
fwrite(vel,sizeof(float),nz*nx,vfp);

fclose(vfp);
free2float(vel);
return (0);
}
