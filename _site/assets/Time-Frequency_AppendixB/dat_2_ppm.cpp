/*   Waalkens, Liefting. June 2016 Appendix B.3 Bachelor Thesis "Time-Frequency Analysis of the Kuramoto Model"  */
#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>

using namespace std;

void hsv2rgb (double h, double s, double v, double *r, double *g, double *b)

{

	double q,p,f,t;
	
	int i;

   if (s==0){
      *r = *g = *b = v;
   }
   
	else{
	   h = h / 60;
	   i = (int)(h);
	   f = h - i;
	   p = v * (1 - s);
	   q = v * (1 - (s * f));
	   t = v * (1 - (s * (1 - f)));
	
	   switch(i){
	   
		  case 0:  *r = v; *g = t; *b = p;break;
	      case 1:  *r = q; *g = v; *b = p;break;
	      case 2:  *r = p; *g = v; *b = t;break;
	      case 3:  *r = p; *g = q; *b = v;break;
	      case 4:  *r = t; *g = p; *b = v;break;
	      case 5:  *r = v; *g = p; *b = q;break;
	   }
	   
	}
	
}

int main(void)

{ 	 
	 double max_fxy=-1.E-10,min_fxy=1.E10;
	 int Nx,Ny; 
	 double x,y,fxy;
	  
	 // find maximum and minimum values 	 	  
	 ifstream infile( "out.dat" ); 
	 if( infile )
	  {
			 // # grid points in x and y direction		
	     infile >> Nx >> Ny;
	 		 for( int ny=1 ; ny<=Ny ; ny++ )
			    for( int nx=1 ; nx<=Nx ; nx++ )
					 {
	 				   infile >> x >> y >> fxy;	
						 max_fxy =  fxy > max_fxy ? fxy : max_fxy;

						 min_fxy =  fxy < min_fxy ? fxy : min_fxy;
					 } 
	  }  
	 else 
    {
      cout << "file could not be opened";
      return (-1);
    }
	 infile.close(); 
	 cout << "Nx : " << Nx << " Ny : " << Ny << "\n";
	 cout << "Min : " << min_fxy << " Max : " << max_fxy << "\n";
	 ifstream infile_2( "out.dat" );
	 ofstream outfile( "out.ppm" ); 
	 outfile << "P3\n\n" << Nx << " " << Ny << "\n\n" << "255\n\n" ;  
	 
	 if( infile_2 )
	  {
			 // # grid points in x and y direction			 
	     infile_2 >> Nx >> Ny;	 
	     double r,g,b; 
	 		 for( int ny=1 ; ny<=Ny ; ny++ )
			    for( int nx=1 ; nx<=Nx ; nx++ )
					 {
	 				   infile_2 >> x >> y >> fxy;					
						 double r,g,b;
						 hsv2rgb( (max_fxy-fxy)/(max_fxy-min_fxy)*250,1.0,1.0,&r,&g,&b);
						 if(x*x+y*y>1)
						 {
						     hsv2rgb( (max_fxy-fxy)/(max_fxy-min_fxy)*250,1.0,1.0,&r,&g,&b);
						     
						     outfile << (int)( pow(fabs((fxy-min_fxy)/(max_fxy-min_fxy)),0.1)*  255*r) << " ";
						     
						     outfile << (int)( pow(fabs((fxy-min_fxy)/(max_fxy-min_fxy)),0.1)* 255*g) << " ";
						     
						     outfile << (int)( pow(fabs((fxy-min_fxy)/(max_fxy-min_fxy)),0.1)* 255*b) << "\n"; 
						 }
						 else 
						 {
						     hsv2rgb( (max_fxy-fxy)/(max_fxy-min_fxy)*250,0.4,1.0,&r,&g,&b);

						     outfile << (int)( pow(fabs((fxy-min_fxy)/(max_fxy-min_fxy)),0.1)*  255*r) << " ";
						     
						     outfile << (int)( pow(fabs((fxy-min_fxy)/(max_fxy-min_fxy)),0.1)* 255*g) << " ";
						     
						     outfile << (int)( pow(fabs((fxy-min_fxy)/(max_fxy-min_fxy)),0.1)* 255*b) << "\n"; 
						 }					 
					 } 
	  }  
	 else 
    {
      cout << "file could not be opened";
      return (-1);
    }
	 infile_2.close();
	 outfile.close(); 
	 return 0;
}








