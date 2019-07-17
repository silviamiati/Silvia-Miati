#include "math.h"
#include "random.h"
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#define pi 3.14

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

///////////////////////////////

system("rm dis_uni.dat dis_non.dat");

ofstream uscitauni,uscitanuni;
uscitauni.open("dis_uni.dat", ios::out);

double a=0.;
double b=1.;

double r=0, x=0, media_uni=0, dev_std_uni=0, S_uni=0, S_uni2=0, inte_uni=0, media_nuni=0, dev_std_nuni=0, S_nuni=0, S_nuni2=0, inte_nuni=0;
int M=100000, B=100; // #lanci, # blocchi
int L = M/B; // #lanci per blocco


for(int i=0; i<B; i++){
	inte_uni=0;
	x=rnd.Rannyu(); // sto utilizzando una variabile random!!!!!!!!!
	for(int h1=0; h1< L; h1++){
		inte_uni += (pi/2.)*cos((pi/2.)*x); // calcolo valore integrale 
		}
	S_uni += inte_uni/(double)L;
	S_uni2 += pow(inte_uni/(double)L,2);
	media_uni = S_uni/double(i+1);

	if(i==0){
		dev_std_uni=0;
	}
	else{
		dev_std_uni=pow((S_uni2/(i+1) - media_uni*media_uni)/i, 0.5);
	}
	uscitauni << media_uni << " " << dev_std_uni << " " << i <<endl;

}

uscitauni.close();
//////////////////////////////////////////////////////non uniforme//////////////////////

uscitanuni.open("dis_non.dat", ios::out);


for(int i=0; i<B; i++){
	inte_nuni=0;
	for(int h=0; h<L; h++){
		r= rnd.poli(); // sto utilizzando il metodo dell'important sampling -ho scelto un polinomio che approssimasse il cosino in [0,1]-
		inte_nuni += (pi/4.)*cos(pi*r*0.5)/(1-r);

		}
	S_nuni += inte_nuni/(double)L;
	S_nuni2 +=pow(inte_nuni/(double)L,2);
	media_nuni= S_nuni/double(i+1);
	
	if(i==0){
		dev_std_nuni=0;
	}
	else{
		dev_std_nuni=pow((S_nuni2/double(i+1) - media_nuni*media_nuni)/double(i), 0.5);
	}
	uscitanuni << media_nuni << " " << dev_std_nuni << " " << i <<endl;

	
	
	}
uscitanuni.close();

return 0;
}
