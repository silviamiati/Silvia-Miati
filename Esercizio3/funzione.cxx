#include "funzione.h"
#include <cmath>
#define pi 3.14159

double erf(double x){
	double err;
	err = (2./sqrt(pi))*((x/sqrt(2)) - (pow((x/sqrt(2)),3)/3.) + (pow((x/sqrt(2)),5)/10.) -(pow((x/sqrt(2)),7)/42.) + (pow((x/sqrt(2)),9)/216.));
	return 0.5*(1. + err);
} 

double max(double x){
	double max1=0;
	if(x>0.){max1=x;
	}else{max1=0.;}
	return max1;
}
	
