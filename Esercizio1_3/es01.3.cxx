#include "random.h"
#include "funzione.h"
#include <cmath>
#include <string>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <stdlib.h>
  

using namespace std;

int main(int argc, const char **argv) {

Random rnd;
   int seed[4];
   int p1, p2;
   ifstream Primes("Primes");
   if (Primes.is_open()){
      Primes >> p1 >> p2 ;
   } else cerr << "PROBLEM: Unable to open Primes" << endl;
   Primes.close();

   ifstream input("seed.in");
   string property;
   if (input.is_open()){
      while ( !input.eof() ){
         input >> property;
         if( property == "RANDOMSEED" ){
            input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
            rnd.SetRandom(seed,p1,p2);
         }
      }
      input.close();
   } else cerr << "PROBLEM: Unable to open seed.in" << endl;


///////////////////////////////////////////io/////////////////////////////////////

system("rm risultati.dat");


double x=0, theta=0;
double L=5/6.; //lunghezza ago
int d=1., block=100; // distanza tra sbarre, #blocchi
double pi=0, pi2=0;
const double Pi=3.141592;
int N=10000; // #lanci
int l = double(N)/(double)block; // # lanci per blocco
double S=0, S2=0, PI=0, PI2=0, dev_std=0, media=0;

ofstream uscita;

uscita.open("risultati.dat", ios::out);

for(int i=0; i<block; i++){	
	int n=0;
	x=0,theta=0;
	for(int j=0; j<l; j++){
		x= rnd.Rannyu(0,0.5);
		theta = rnd.Rannyu(0,Pi);
		if(x<=(L/2)*fabs(cos(theta))){
			n++;
		}
	}
	S += n/double(l) ;
	
	PI += (2*L)/((double)d*(n/double(l)));
	PI2 += pow((2*L)/((double)d*(n/double(l))),2);

	media= S/(double)(i+1);
	
	pi=(2*L)/(double(d)*media);
	pi2 = pow((2*L)/(double(d)*media),2);
	
	if(i==0){
		dev_std=0;
	}
	else{
	dev_std=pow(((PI2/(i+1)) - pi2)/i, 0.5);
	}

	uscita << pi << " " << dev_std << " " << i << endl;
}

uscita.close();

return 0;
}
