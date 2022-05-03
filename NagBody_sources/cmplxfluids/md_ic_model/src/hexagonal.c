#include<stdio.h>
int main()
{
FILE *fpt;
	fpt=fopen("hexagonal.dat", "w");	
	int iix, iiy, iiz;
	int N;
	int ncx, ncy, ncz;
	int cx, cy, cz;
	double a, c;
	double xx, yy, zz;
	
	double r1x0, r1y0, r1z0;
	double r2x0, r2y0, r2z0;
	double r3x0, r3y0, r3z0;
	
	double r1x,r1y,r1z;
	double r2x,r2y,r2z;
	double r3x,r3y,r3z;

	
//numero de hexagonos

	ncx=8;
	ncy=6;
	ncz=1;
	
	a=1.0; //del hexagono
	c=1.0; //altuta en z		
	
	//contador de hexagonos

	cx=7;
	cy=5;
	cz=0;
	
	//posiciones iniciales
	
	r1x0=0.0;					 r1y0=0.86025;			r1z0=0.0;
	r2x0=0.0;		r2y0=0.28867;			r2z0=0.0;
	r3x0=a/2.0;		r3y0=0.0;			r3z0=0.0;
	
			
	xx=0.0; 	yy=0.0; 	zz=0.0;
	for(iiz=0; iiz<=cz; iiz++){
		yy=0.0;
		for(iiy=0; iiy<=cy; iiy++)	{
			xx=0.0;
			for(iix=0; iix<=cx; iix++)	{
			r1x=r1x0+a*xx;	r1y=r1y0+yy*1.73205;		r1z=r1z0+zz*c;
			r2x=r2x0+a*xx;	r2y=r2y0+yy*1.73205;		r2z=r2z0+zz*c;
			r3x=r3x0+a*xx;	r3y=r3y0+yy*1.73205;		r3z=r3z0+zz*c;
			

			fprintf(fpt, "Si	%.5f		%.5f		%.5f\n", r1x, r1y, r1z);
			fprintf(fpt, "Si	%.5f		%.5f		%.5f\n", r2x, r2y, r2z);
			fprintf(fpt, "Si	%.5f		%.5f		%.5f\n", r3x, r3y, r3z);
	
			xx=xx+1.0;
			}
			
		iix++;
		r1x=r1x0+a*xx;	r1y=r1y0+yy*1.73205;		r1z=r1z0+zz*c;
		r2x=r2x0+a*xx;	r2y=r2y0+yy*1.73205;		r2z=r2z0+zz*c;	
		r3x=r3x0+a*xx;	r3y=yy*1.73205;		r3z=r3z0+zz*c;	

		fprintf(fpt, "Si	%.5f		%.5f		%.5f\n", r1x, r1y, r1z);
		fprintf(fpt, "Si	%.5f		%.5f		%.5f\n", r2x, r2y, r2z);
		fprintf(fpt, "Si	%.5f		%.5f		%.5f\n", r3x, r3y, r3z);
			
				xx=0.0;
			for(iix=0; iix<=cx+1; iix++)	{
			r1x=r3x0+a*xx;	r1y=1.1547+yy*1.73205;		r1z=r1z0+zz*c;
			
			fprintf(fpt, "Si	%.5f		%.5f		%.5f\n", r1x, r1y, r1z);
				
			xx=xx+1.0;
		}
		yy=yy+1.0;
	
		}
		zz=zz+1.0;
	}
				
	fclose(fpt);
}
