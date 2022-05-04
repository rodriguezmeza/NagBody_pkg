

#include "stuff.h"

/*
 * 
 * Extracting displacements from the 2LPTic snapshots.
 * 
 * 
 * 
 */ 
  

void Displacements(void){
    
    dX=malloc(NumPart*sizeof(float));
    dY=malloc(NumPart*sizeof(float));
    dZ=malloc(NumPart*sizeof(float));
   
    obtainDisplacements(P,dX,dY,dZ);
    
}
 

float fix(float x,float q){ // Extract particle displacement from Eulerian and Lagrangian position. It takes care of periodic BC assuming box size is >> rms particle displacement (~15Mpc at z=0).
    q=putInBox(q);
    float s=x-q;
    
    if (-s>BoxSize/2.)s=x-(q-BoxSize);
    if (s>BoxSize/2.) s=(x-BoxSize)-q;
    
    return s;
    
}

void obtainDisplacements(struct particle_data *P,float *dX,float *dY,float *dZ){ // uses FORTRAN ordering

      int i,j,k,it;
      float ingrid=1./((float)NROW);
    
      float summ,f;
      summ=0;

      for (i=0;i<NROW;i++){
          for (j=0;j<NROW;j++){
              for (k=0;k<NROW;k++){
                it=i+j*NROW+k*NROW*NROW+1;

                f=fix(P[it].Pos[0],((float)i)*ingrid*BoxSize);
                WRtNROW(dX, i, j, k, f);
                
                summ+=f;
                f=fix(P[it].Pos[1],((float)j)*ingrid*BoxSize);
                WRtNROW(dY, i, j, k, f);

                f=fix(P[it].Pos[2],((float)k)*ingrid*BoxSize);
                WRtNROW(dZ, i, j, k, f);

      }}}
    // should be ->zero:
    //printf("==============AVERAGE 1D DISPL = %f\n",summ*powf(ingrid,3));
    meanM(dX,NROW*NROW*NROW,0.0);
    meanM(dY,NROW*NROW*NROW,0.0);
    meanM(dZ,NROW*NROW*NROW,0.0);
    
    
  }

/*
 * 
 * Change velocity variable from ds/dy to (ds/d\eta)/(aH(a))
 * 
 */ 


void velRSD(struct particle_data *P,float A){
        int i;
        float fac=A/Qfactor(A);
        for(i=1;i<=NumPart;i++){
            P[i].Vel[0]*=fac;
            P[i].Vel[1]*=fac;
            P[i].Vel[2]*=fac;
        }
}

/*
 * 
 *  Below follow the time dependent functions used in the code.
 *  Start with growth factors and derivatives.
 * 
 */



float growthD(float a){ // growth factor for LCDM
    return growthDtemp(a)/growthDtemp(1.0);
}


float growthDtemp(float a){
    // Decided to use the analytic expression for LCDM. More transparent if I change this to numerical integration?
    float x=-Om/(Om-1.0)/(a*a*a);

    
    
    float hyperP=0,hyperM=0;
    
    if (fabs(x-1.0) < 1.e-3) {
       // printf("mechka\n");
    hyperP=0.859596768064608 - 0.1016599912520404*(-1.0 + x) + 0.025791094277821357*pow(-1.0 + x,2) - 0.008194025861121475*pow(-1.0 + x,3) + 0.0029076305993447644*pow(-1.0 + x,4) - 0.0011025426387159761*pow(-1.0 + x,5) + 0.00043707304964624546*pow(-1.0 + x,6) - 0.0001788889964687831*pow(-1.0 + x,7);
    hyperM=1.1765206505266006 + 0.15846194123099624*(-1.0 + x) - 0.014200487494738975*pow(-1.0 + x,2) + 0.002801728034399257*pow(-1.0 + x,3) - 0.0007268267888593511*pow(-1.0 + x,4) + 0.00021801569226706922*pow(-1.0 + x,5) - 0.00007163321597397065*pow(-1.0 + x,6) +    0.000025063737576245116*pow(-1.0 + x,7);
    }
    else {
        if (x < 1.0) {
            hyperP=gsl_sf_hyperg_2F1(1.0/2.0,2.0/3.0,5.0/3.0,-x);
            hyperM=gsl_sf_hyperg_2F1(-1.0/2.0,2.0/3.0,5.0/3.0,-x);
        }
        x=1.0/x;
        if ((x < 1.0) && (x>1.0/30)) {
            
            hyperP=gsl_sf_hyperg_2F1(-1.0/6.0,0.5,5.0/6.0,-x);
            hyperP*=4*sqrt(x);
            hyperP+=-3.4494794123063873799*pow(x,2.0/3.0);
        
            hyperM=gsl_sf_hyperg_2F1(-7.0/6.0,-0.5,-1.0/6.0,-x);
            hyperM*=4.0/7.0/sqrt(x);
            hyperM+=pow(x,2.0/3.0)*(-1.4783483195598803057); //-(Gamma[-7/6]*Gamma[5/3])/(2*sqrt[Pi])
        }
        if (x<=1.0/30.0){
            hyperP=3.9999999999999996*sqrt(x) - 3.4494794123063865*pow(x,0.6666666666666666) + 0.3999999999999999*pow(x,1.5) -    0.13636363636363635*pow(x,2.5) + 0.07352941176470587*pow(x,3.5) - 0.04755434782608695*pow(x,4.5) +    0.033943965517241374*pow(x,5.5) - 0.02578125*pow(x,6.5) + 0.020436356707317072*pow(x,7.5) -    0.01671324384973404*pow(x,8.5) + 0.013997779702240564*pow(x,9.5) - 0.011945562847590041*pow(x,10.5) + 0.01035003662109375*pow(x,11.5) - 0.009080577904069926*pow(x,12.5);
            hyperM=0.5714285714285715/sqrt(x) + 2.000000000000001*sqrt(x) - 1.4783483195598794*pow(x,0.66666666666666666) +    0.10000000000000002*pow(x,1.5) - 0.022727272727272735*pow(x,2.5) + 0.009191176470588237*pow(x,3.5) -    0.004755434782608697*pow(x,4.5) + 0.002828663793103449*pow(x,5.5) - 0.0018415178571428578*pow(x,6.5) +    0.0012772722942073172*pow(x,7.5) - 0.0009285135472074472*pow(x,8.5) + 0.0006998889851120285*pow(x,9.5) -    0.0005429801294359111*pow(x,10.5) + 0.0004312515258789064*pow(x,11.5) - 0.00034925299631038194*pow(x,12.5);
        }
    }
    
    
    if (a>0.2) return sqrt(1.0 + (-1.0 + pow(a,-3))*Om)*(3.4494794123063873799*pow(-1.0 + 1.0/Om,0.666666666666666666666666666) + (hyperP*(4*pow(a,3)*(-1.0 + Om) - Om) - 7.0*pow(a,3)*hyperM*(-1.0 + Om))/(pow(a,5)*(-1.0+ Om) - pow(a,2)*Om));
    return (a*pow(1 - Om,1.5)*(1291467969*pow(a,12)*pow(-1 + Om,4) + 1956769650*pow(a,9)*pow(-1 + Om,3)*Om + 8000000000*pow(a,3)*(-1 + Om)*pow(Om,3) + 37490640625*pow(Om,4)))/(1.5625e10*pow(Om,5));
    
        }






float Qfactor(float a){ // Q\equiv a^3 H(a)/H0.
    return sqrt(Om/(a*a*a)+1.0-Om)*a*a*a;
}




float growthD2(float a){// Second order growth factor
    return growthD2temp(a)/growthD2temp(1.0);
}


float growthD2temp(float a){
    float d= growthD(a);
    float omega=Om/(Om+(1.0-Om)*a*a*a);
    return d*d*pow(omega,-1./143.);
}
 
float growthD2v(float a){ // explanation is in main()
    float d2= growthD2(a);
    float omega=Om/(Om+(1.0-Om)*a*a*a);
    return Qfactor(a)*(d2/a)*2.0*pow(omega,6./11.);
}

 

float decayD(float a){ // D_{-}, the decaying mode
    return sqrt(Om/(a*a*a)+1.0-Om);
}



double DprimeQ(double a,float nGrowth){ // returns Q*d(D_{+}^nGrowth*D_{-}^nDecay)/da, where Q=Qfactor(a)
        float nDecay=0.0;// not interested in decay modes in this code.
        float Nn=6.0*pow(1.0 - Om,1.5)/growthDtemp(1.0);
        return (pow(decayD(a),-1.0 + nDecay)*pow(growthD(a),-1.0 + nGrowth)*(nGrowth*Nn- (3.0*(nDecay + nGrowth)*Om*growthD(a))/(2.*a)));
        
    }
   


/*
 * 
 *  
 * Functions for our modified time-stepping (used when StdDA=0):
 * 
 * 
 */


double fun (double a, void * params) {
       
       double f;
       if (fullT==1) f = gpQ(a)/Qfactor(a); 
       else f = 1.0/Qfactor(a);
      
       return f;
     }
     
double   Sq(double ai,double af,double aRef)
     {
       gsl_integration_workspace * w 
         = gsl_integration_workspace_alloc (5000);
       
       double result, error;
        double alpha=0;
     
       gsl_function F;
       F.function = &fun;
       F.params = &alpha;
     
       gsl_integration_qag (&F, ai, af, 0, 1e-5, 5000,6,
                             w, &result, &error); 
     
       //printf ("result          = % .18f\n", result);
       //printf ("estimated error = % .18f\n", error);
      
     
       gsl_integration_workspace_free (w);
     
      
        //   result  = \int fun(a)da 
      
       if (fullT==1) return result/gpQ(aRef);
       else return result;
       
     }
     
     
double   Sphi(double ai,double af,double aRef)
     {
       double result;
       result=(gpQ(af)-gpQ(ai))*aRef/Qfactor(aRef)/DERgpQ(aRef);
     
     
       return result;
     }



double gpQ(double a){ 
    
        return pow(a,nLPT);
        
    }
double DERgpQ(double a){ // This must return d(gpQ)/da
        
        return nLPT*pow(a,nLPT-1);
        
    }
    




/*
 * 
 *  
 * Functions for Quinn et al time-stepping (used when StdDA=2):
 * 
 * 
 */ 


double funSqStd (double a, void * params) {
       
       double f = 1.0/Qfactor(a);
      
       return f;
     }
     
double   SqStd(double ai,double af)
     {
       gsl_integration_workspace * w 
         = gsl_integration_workspace_alloc (5000);
       
       double result, error;
        double alpha=0;
     
       gsl_function F;
       F.function = &funSqStd;
       F.params = &alpha;
     
       gsl_integration_qag (&F, ai, af, 0, 1e-5, 5000,6,
                             w, &result, &error); 
     
       //printf ("result          = % .18f\n", result);
       //printf ("estimated error = % .18f\n", error);
      
     
       gsl_integration_workspace_free (w);
     
      
       return result;
       
       
     }
     
     

double funSphiStd (double a, void * params) {
       
       double f = a/Qfactor(a);
      
       return f;
     }
     
double   SphiStd(double ai,double af)
     {
       gsl_integration_workspace * w 
         = gsl_integration_workspace_alloc (5000);
       
       double result, error;
        double alpha=0;
     
       gsl_function F;
       F.function = &funSphiStd;
       F.params = &alpha;
     
       gsl_integration_qag (&F, ai, af, 0, 1e-5, 5000,6,
                             w, &result, &error); 
     
       //printf ("result          = % .18f\n", result);
       //printf ("estimated error = % .18f\n", error);
      
     
       gsl_integration_workspace_free (w);
     
       return result;
       
       
     }





