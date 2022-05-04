#ifndef _NR_Q_H_
#define _NR_Q_H_


float qgaus(float (*func)(float), float a, float b);
float qgaus_01(float (*func)(float), float a, float b);
/* void qrdcmp(float **a, int n, float *c, float *d, int *sing); */
float qromb(float (*func)(float), float a, float b);
float qromo(float (*func)(float), float a, float b,
	float (*choose)(float (*)(float), float, float, int));
float qromo_01(float (*func)(float), float a, float b,
	float (*choose)(float (*)(float), float, float, int));
/* void qroot(float p[], int n, float *b, float *c, float eps); */
/* void qrsolv(float **a, int n, float c[], float d[], float b[]); */
/* void qrupdt(float **r, float **qt, int n, float u[], float v[]); */
float qsimp(float (*func)(float), float a, float b);
float qtrap(float (*func)(float), float a, float b);
float quad3d(float (*func)(float, float, float), float x1, float x2);
/*
void quadct(float x, float y, float xx[], float yy[], unsigned long nn,
	float *fa, float *fb, float *fc, float *fd);
void quadmx(float **a, int n);
void quadvl(float x, float y, float *fa, float *fb, float *fc, float *fd);
*/
void polint(float xa[], float ya[], int n, float x, float *y, float *dy);
void polint_01(float xa[], float ya[], int n, float x, float *y, float *dy);
float trapzd(float (*func)(float), float a, float b, int n);

float midexp(float (*funk)(float), float aa, float bb, int n);
float midexp_01(float (*funk)(float), float aa, float bb, int n);
float midinf(float (*funk)(float), float aa, float bb, int n);
float midinf_01(float (*funk)(float), float aa, float bb, int n);
float midpnt(float (*func)(float), float a, float b, int n);
float midpnt_01(float (*func)(float), float a, float b, int n);
float midsql(float (*funk)(float), float aa, float bb, int n);
float midsql_01(float (*funk)(float), float aa, float bb, int n);
float midsqu(float (*funk)(float), float aa, float bb, int n);
float midsqu_01(float (*funk)(float), float aa, float bb, int n);

#endif /* _NR_Q_H_ */
