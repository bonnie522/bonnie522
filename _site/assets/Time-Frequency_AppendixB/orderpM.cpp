/*   Waalkens, Liefting. June 2016 Appendix B.1 Bachelor Thesis "Time-Frequency Analysis of the Kuramoto Model" 
Computes the trajectory of the order parameter of the Kuramoto model.*/

#include <iostream>
#include <iomanip>
#include <cmath>
#include "nr.h"

#include <random>
#include <complex>


using namespace std;

	const int n=10;

	double w_array[n];

	double runtime=120;
	double K=5.0;


double phistart_array[n];



// Driver for routine odeint
DP dxsav;  // defining declarations
int kmax,kount;
Vec_DP *xp_p;
Mat_DP *yp_p;

int nrhs;   // counts function evaluations

void derivs(const DP t, Vec_I_DP &yvector, Vec_O_DP &dydx)
{	//double phi1dot, phi2dot, phi3dot;
	int i=0;
	int k=0;

        nrhs++;

for(i=0;i<n;i++){
dydx[i] = w_array[i];
	for(k=0;k<n;k++){
	dydx[i]+=K/n*sin(yvector[k]-yvector[i]);

	}
}

}

DP mod_2_PI(DP x)
{	while (x<=0) {
		x += 2.0*M_PI;
	}while (x>2*M_PI) {
		x -= 2.0*M_PI;
	}return x;
}


int main(void)
{	const int N=n, KMAX=100;
        int nbad,nok;
        DP eps=1.0e-7; //original: 1.0e-10
	DP h1=1E-10,hmin=0.0,x1=0.0,x2=0.0,dx = 0.03;
	DP x1max = runtime;
        Vec_DP ystart(N);
	Vec_DP ystart_old(N);
	int i=0;
	
ofstream outfile_w("omega.dat");
ofstream outfile_p("phases.dat");

  std::default_random_engine generator;
  std::normal_distribution<double> distribution(0.3,0.05); // (mean, standarddev)

	for(i=0;i<n;i++){

w_array[i]=distribution(generator);  // normal distribution: http://www.cplusplus.com/reference/random/normal_distribution/       
outfile_w << w_array[i]<<" "<<endl;
       
}
  std::default_random_engine generator2;
  std::uniform_real_distribution<double> distribution2(0,2*M_PI); // uniform distribution [a,b]

for(i=0;i<n;i++){	
ystart[i]=distribution2(generator2);
//ystart[i]=0;

outfile_p << ystart[i]<< " "<<endl;
        }
	outfile_w.close();
	outfile_p.close();

		
	nrhs = 0;
        dxsav=(x2-x1)/2.0;
        kmax=KMAX;
        xp_p=new Vec_DP(KMAX);
        yp_p=new Mat_DP(N,KMAX);

		

	ofstream outfile_s("trajectory.dat");
	ofstream outfile_o("order.dat");


	for(x1=0;x1<x1max;x1+=dx)
	{	cout << "time: " << x1 << endl;

        x2 = x1 +dx;
	NR::odeint(ystart,x1,x2,eps,h1,hmin,nok,nbad,derivs,NR::rkqs);
            
	for(i=0;i<n;i++){
	ystart[i]=mod_2_PI(ystart[i]);
	}

	outfile_s << x1;

	complex<double> im=-1;	
	im = sqrt(im);
for(i=0;i<n;i++){ 
	outfile_s<< " "<< ystart[i];
}
	outfile_s << endl;     			


	complex<double> order=exp(im*ystart[0]);
	for(i=1;i<n;i++){
	order += exp(im*ystart[i]);   //divide by n
	}

	double realorder=real(order)/n;
	double imagorder=imag(order)/n;

	outfile_o <<x1<<" "<< realorder <<" "<< imagorder<<" "<<endl;  

        }
	
	outfile_s.close();
	outfile_o.close();	
int length=runtime/dx+1;

cout << "length order.dat="<< length << endl;
		
	delete yp_p;
	delete xp_p;
return 0;
}
