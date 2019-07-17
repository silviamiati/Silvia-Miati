#include "funzione.h"
#include "random.h"
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#define M_e 2.71828

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


system("rm S_directcall.dat S_directput.dat S_discretizedcall.dat S_discretizedput.dat");

double S0=100.; //asset price at t=0
double T=1.; //delivery time
double K=100.; // strike price
double r=0.1; // risk-free interest rate
double sigma=0.25; //volatility
double sigma2= sigma*sigma;

double d1= 1./(sigma*sqrt(T))*(log(S0/K) + (r+ (sigma2/2.))*T);
double d2 = d1 - sigma*sqrt(T);

double Ce=0;
double Cp=0; 

Ce= S0*erf(d1) - K*pow(M_e,-r*T)*erf(d2);
Cp= S0*(erf(d1)-1) -K*pow(M_e,-r*T)*(erf(d2)-1);

cout << Ce << " " << Cp << endl; 

//////////////////////////// S_direct //////////////////////
ofstream uscita1;
ofstream uscita2;

uscita1.open("S_directcall.dat", ios::out);
uscita2.open("S_directput.dat", ios::out);


int N= 10000, B=100, L = N/B;


double S=0, S2=0, S_direct=0, mediaC=0, mediaP=0, dev_stdC=0, dev_stdP=0;
double C_direct=0, P_direct=0, C=0, P=0, C2=0, P2=0;



for(int i=0; i<B; i++){		
	S_direct=0;
	C_direct=0;
	P_direct=0;
		for(int i=0; i< L; i++){
			S_direct = S0*pow(M_e, ((r-0.5*sigma2)*T) + sigma*rnd.Gauss(0,1)*sqrt(T));
			C_direct += pow(M_e, -r*T)* max(S_direct - K);
			P_direct += pow(M_e, -r*T)* max(-S_direct + K);
		} 
	C += C_direct/double(L);
	P += P_direct/double(L);
	
	
	C2 +=pow(C_direct/double(L),2);
	P2 +=pow(P_direct/double(L),2);

	
	mediaC = C/double(i+1);
	mediaP = P/double(i+1);

	if(i==0){
		dev_stdC=0;
		dev_stdP=0;
	}
	else{
		dev_stdC=pow((C2/double(i+1) - mediaC*mediaC)/double(i), 0.5);
		dev_stdP=pow((P2/double(i+1) - mediaP*mediaP)/double(i), 0.5);
	}
	uscita1 << mediaC << " " << dev_stdC << " " << i <<endl;
	uscita2 << mediaP << " " << dev_stdP << " " << i <<endl;
	}

uscita1.close();
uscita2.close();

///////////////////////////////////////////// S_discretized ///////////////////// 
ofstream uscita3;
ofstream uscita4;

uscita3.open("S_discretizedcall.dat", ios::out);
uscita4.open("S_discretizedput.dat", ios::out);

int tem=100;
double Sd=0, S2d=0, mediaCd=0, mediaPd=0, dev_stdCd=0, dev_stdPd=0;
double C_discre=0, P_discre=0, Cd=0, Pd=0, C2d=0, P2d=0;
double S_discre[B];
double t[B];
double delta = 0.01;

S_discre[0]=100.;
t[0] =0.;
for(int i=0; i<B; i++){	
	C_discre=0;
	P_discre=0;
		for(int j=0; j<L; j++){
			for(int i=1; i<=tem; i++){
        		t[i] = t[i-1] + delta;
        		S_discre[i] = S_discre[i-1]*exp((r-sigma2*0.5)*(t[i]-t[i-1]) + sigma*sqrt(t[i]-t[i-1])*rnd.Gauss(0,1));
			} 
		C_discre += pow(M_e, -r*T)* max(S_discre[B] - K);
		P_discre += pow(M_e, -r*T)* max(-S_discre[B] + K);
		
		}
		
	Cd += C_discre/double(L);
	Pd += P_discre/double(L);
	
	C2d +=pow(C_discre/double(L),2);
	P2d +=pow(P_discre/double(L),2);

	mediaCd = Cd/double(i+1);
	mediaPd = Pd/double(i+1);


	if(i==0){
		dev_stdCd=0;
		dev_stdPd=0;
	}
	else{
		dev_stdCd=pow((C2d/double(i+1) - mediaCd*mediaCd)/double(i), 0.5);
		dev_stdPd=pow((P2d/double(i+1) - mediaPd*mediaPd)/double(i), 0.5);
	}
	uscita3 << mediaCd << " " << dev_stdCd << " " << i <<endl;
	uscita4 << mediaPd << " " << dev_stdPd << " " << i <<endl;
	}

uscita3.close();
uscita4.close();

return 0;
}
