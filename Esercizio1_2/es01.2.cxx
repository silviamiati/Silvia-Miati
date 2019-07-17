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

system("rm isto_1.dat isto_2.dat isto_10.dat isto_100.dat");


int M= 10000; // #lanci
double r=0, rl=0, re=0;
int Nb=100; // #blocchi
double ms=-1., Ms=1.; // massimi e minimi delle funzioni random, esponenziale e lorentziana
double me=-4., Me=6.;
double ml=-15., Ml=15.;

double ds = (Ms-ms)/Nb; // passo per costruire l'istogramma
double de = (Me-me)/Nb;
double dl = (Ml-ml)/Nb;

int N[4] = {1,2,10,100};
int *vetr = new int[Nb];
int *vete = new int[Nb];
int *vetl = new int[Nb];

ofstream uscita;

string nome= "isto_";
string estensione= ".dat";

for(int i=0; i<4; i++){
	null(vetr,Nb);
	null(vete,Nb);
	null(vetl,Nb);
	for(int l=0; l<M; l++){ 
		r=0, re=0, rl=0;
			for(int k=0; k<N[i]; k++){
				r +=rnd.Rannyu(-1.,1.);
				re +=rnd.esponenziale(1.);
				rl +=rnd.lorenziana(1.,0.);
			}
		r= r/double(N[i]);
		re=re/double(N[i]);
		rl=rl/double(N[i]);

		
		isto(vetr, Nb, ms, Ms, r);
		isto(vete, Nb, me, Me, re);
		isto(vetl, Nb, ml, Ml, rl);
	}

	uscita.open(nome+to_string(N[i])+estensione);
	for(int j=0; j<Nb; j++){
		uscita<< ms + (double)j*ds << " " << vetr[j] << " " << me + (double)j*de << " " << vete[j] << " " <<ml + (double)j*dl << " " << vetl[j] << endl;
	}
	uscita.close();
}

delete[] vetr;
return 0;
}
