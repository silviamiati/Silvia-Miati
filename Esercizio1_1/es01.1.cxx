#include <cmath>
#include <string>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include "random.h"


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


system("rm risultati1.dat risultati2.dat");

int B=100; //numeri di blocchi
int M=100000; //numeri di lanci	
int L=M/B; //numeri random per ogni blocco 
double sum=0, sum_sigma=0, media=0, err=0; 
double S=0, S2=0, D=0, D2=0, dev_std=0, dev_std_err=0;
double r=0;

ofstream uscita1;
ofstream uscita2;

uscita1.open("risultati1.dat", ios::out);
uscita2.open("risultati2.dat", ios::out);


for(int i=0; i<B; i++){	
	sum=0;
	sum_sigma=0;
	for(int h=0; h<L; h++){
		r= rnd.Rannyu();
		sum += r;
		sum_sigma += pow(r-0.5,2);
		}
	S += sum/(L);
	S2 +=pow(sum/(double)(L),2);
	D += sum_sigma/(double)(L);
	D2 +=pow(sum_sigma/(double)(L),2);
	media= S/(double)(i+1);
	err=D/(double)(i+1);
	
	if(i==0){
		dev_std=0;
		dev_std_err=0;
	}
	else{
		dev_std=pow(((S2/(i+1) - media*media)/i), 0.5);
		dev_std_err=pow(((D2/(i+1) - err*err)/i), 0.5);
	}		
	uscita1 << media << " " << dev_std << " " << i <<endl;
	uscita2 << err << " " << dev_std_err << " " << i <<endl;
	}

uscita1.close();
uscita2.close();

return 0;
}
