#include "funzione.h"

void isto(int *NB, int bin, double m, double M, double R){
	double d = (M-m)/double(bin);
 	int i = floor(R/d) - floor(m/d);
	if(i< bin && i>=0) {
		NB[i]++;
	}
}

void null(int *NB, int bin){
	for(int j=0; j<bin; j++){
		NB[j]=0;
	}
}

