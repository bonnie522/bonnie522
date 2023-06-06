/*   Waalkens, Liefting. June 2016 Appendix B.2 Bachelor Thesis "Time-Frequency Analysis of the Kuramoto Model"
	Calculates the wavelet transform of a single trajectory of the 3D system. The initial conditions are given by n1 and n2 which in a discrete way
   parametrize the energy surface. the wavelet transform is displayed by the program zljap in the same directory
*/

#include <math.h>
#include <stdio.h>

#define epsilon  0.04


/* parameters for the morlet wavelet */
#define sigma 2.0
#define lambda 1.0

#define ntmax 1000 /* 4000 */ //was 800000     // Total length signal

FILE *out,*out1,*out2,*out3,*in,*out_main,*out_freq;

double morlet_re(double);

double morlet_im(double);

double wave_trans(double, int, double *,double *, double *); 

double maxis(double,double);

double minis(double,double);

double theta(double,double,double);


int main(void)
{

  double eps=1.E-5;  // original: 1.E-14
  double h1=1.E-1,hmin=0,t1,t2,tmin = 0.,tmax = 440000.; /* 2200 */
  int nvar = 6,nok,nbad;
  int nt;
  
  double omega;
  
  /* for the wavelet transform*/
  
  double omega_x0,omega_y0,omega_z0;
  
  double wave_trans_x,wave_trans_y,wave_trans_z,wave_trans_0,wave_trans_sub,wave_trans_top;
  
  double omegamin = 0.2/(2.0*M_PI),omegamax = 0.45/(2.0*M_PI),domega; //omegamin = 0.01,omegamax = 1.0,domega;

  int nomega,nomegamax = 200;
  
  int nb,nbmin,nbmax;
  
    
   double x_array[ntmax];
   double px_array[ntmax];
   double t_array[ntmax];
  
   double omega_x0_array[ntmax];

   
  nbmin = 200;   // margins left & right
  nbmax = ntmax - nbmin;
  
  domega = (omegamax-omegamin)/nomegamax;
  
  printf("%f %f %f\n",omegamin,omegamax,domega);
  
  double xstart[nvar];

  printf("epsilon: %f\n",epsilon);

        /* reads trajectory from order.dat produced by orderpM.cpp */
          FILE*in_order = fopen("order.dat", "r"); // time - real part - imaginary part  

 	for(int k=0; k<ntmax; k++){
		fscanf(in_order, "%lf %lf %lf", &t_array[k], &x_array[k], &px_array[k] );  

printf("%lf %lf %lf \n", t_array[k], x_array[k], px_array[k]);
             } 

/* output for dat_2_ppm.C */ 
int l201=nbmin+1;
int l601=ntmax-2*nbmin+1;    

     out1 = fopen("out.dat","w");
    fprintf(out1, "%d %d \n", l201, l601);
     
     for(nb=nbmin;nb<=nbmax;nb++)
     {
              printf("nb=%d\n",nb);
    
      for(nomega=0;nomega<=nomegamax;nomega++)
      {
          
	 fprintf(out1,"%d %d ",nb,nomega);
	 
         
	 omega_x0 = omegamin +  nomega*(omegamax-omegamin)/nomegamax;
	 
          
         wave_trans_0     = wave_trans(omega_x0,nb,x_array,px_array,t_array);
	 fprintf(out1," %f\n",wave_trans_0); 
      }
      }
    
    fclose(out1); 
}


double morlet_re(double t)   /*  sigma=2.0, lambda=1.0  */
{   
   return 1./sigma/sqrt(2.*M_PI)*cos((2*M_PI*lambda*t))*exp(-t*t/(2.*sigma*sigma));
}

double morlet_im(double t)
{   
   return 1./sigma/sqrt(2.*M_PI)*sin((2*M_PI*lambda*t))*exp(-t*t/(2.*sigma*sigma));
}

double maxis(double x,double y)
{
  if(x>=y)
    return x;
  else
    return y;  
}

double minis(double x,double y)
{
  if(x<=y)
    return x;
  else
    return y;  
}

double wave_trans(double omega, int nb, double *array_signal, double *array_signal_H , double *array_t)
{
   double dt,a,b,tmin,tmax,mean = 0.,trans,trans_re = 0.,trans_im = 0.;
   double *array_si;
   int nt;
   
    a = (lambda + sqrt( lambda*lambda + 1./(2.*M_PI*M_PI*sigma*sigma) ))/(2.*omega);   
   
    b = array_t[nb];
   
    tmin = b - 3.*sqrt(2)*sigma*a;
    tmax = b + 3.*sqrt(2)*sigma*a;

  
      /* determine wavelet transform */		      
 
      dt = array_t[2] - array_t[1];
 
      for(nt=(int) maxis(0., tmin/dt -2.) ; nt<= (int) minis(ntmax,tmax/dt + 2.) ; nt++)
      {
            trans_re += (
                         array_signal[nt]   * morlet_re((array_t[nt]-b)/a) + 
	                 array_signal_H[nt] * morlet_im((array_t[nt]-b)/a)  
			)/sqrt(a) *dt ;
        			  
	    trans_im += (
                         array_signal_H[nt] * morlet_re((array_t[nt]-b)/a) -
			 array_signal[nt]   * morlet_im((array_t[nt]-b)/a)  
	                   
                        )/sqrt(a)* dt;	 
        	 
      }
      
      trans = sqrt(trans_re*trans_re + trans_im*trans_im); 
     
     return trans;
}


double theta(double q,double p,double h)
{
  double phi,thetais;
  
  phi = acos(q/pow(4.*h,0.25));
  
    thetais = M_PI/2.; /**ellf(phi,1./sqrt(2.))/ellK;*/
  
  if(q>=0 && p>=0)
    return 2*M_PI-thetais;
  else if(p>0 && q<=0)
     return 2*M_PI-(M_PI-thetais);
  else if(p<0 && q<=0)
     return 2*M_PI-(thetais+M_PI);
  else
     return 2*M_PI-(2.*M_PI-thetais);        

}
