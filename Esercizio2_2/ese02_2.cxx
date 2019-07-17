#include "funzione.h"
#include "random.h"
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>

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

/////////////////////////////////////////////// RW RETICOLO ///////////////////////////////////////////////

system("rm rw.dat rwlibero.dat");

ofstream RWxyz;
ofstream RWtheta;

RWxyz.open("rw.dat", ios::out);  

double s=0,t=0;

vx[0]=0; // vettori per descrivere il reticolo
vy[0]=0;
vz[0]=0;


for(int j=0; j<B; j++){
	Reset(1);	
	for(int i=0; i<L; i++){
		Resetv(1);
		for(int k=0; k<nstep; k++){
			s=rnd.Rannyu(0.,3.);
			t=rnd.Rannyu(0.,2.);
			Movexyz(vx,vy,vz,k,1.,s,t); // movimento casuale nel reticolo (aggiorna la posizione)
			D[k] += Distanza(vx,vy,vz,k);
		}
	}
	for(int k=0; k<nstep; k++){ 
		stima_ar[k] = D[k]/(double)L;
		global_ar[k] += stima_ar[k];
		global_ar2[k]+=stima_ar[k]*stima_ar[k];
  }
}


for(int k=0; k<nstep; k++){
	RWxyz << sqrt(global_ar[k]/double(B))<< " " << sqrt(Error(global_ar[k], global_ar2[k], B)) << " " << k << endl;
}


RWxyz.close();

//////////////////////////////////////////////////////////// RW CONTINUO ///////////////////////////////////////////

RWtheta.open("rwlibero.dat", ios::out);

double r=0, z=0, a=1.;

vtheta[0]=0;
vphy[0]=0;
vrho[0]=0;

for(int j=0; j<B; j++){
	Reset(0);
	for(int l=0; l<L; l++){
		Resetv(0); 
		for(int k=0; k<nstep; k++){
			r=rnd.Dtheta();
			z=rnd.Rannyu(0.,2*pi);
			Movetheta(vtheta,vphy,vrho,k,a,r,z); //movimento casuale nel continuo (aggiorna la posizione)
			Dlibero[k] += Distanza(vtheta,vphy,vrho,k);
		}
	}
	for(int k=0; k<nstep; k++){ 
		stima_al[k] = Dlibero[k]/(double)L;
		global_al[k] += stima_al[k];
		global_al2[k] += stima_al[k]*stima_al[k];
  }
} 

for(int k=0; k<nstep; k++){
	RWtheta << sqrt(global_al[k]/double(B)) << " " << sqrt(Error(global_al[k], global_al2[k], B)) << " " << k << endl;
}

RWtheta.close();

return 0;
}


double Error(double sum, double sum2, int iblk){
    return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)iblk);
}


void Reset(bool libero){
	if(libero==1){ // libero=1 => reticolo se no continuo
		for (int i=0; i<nstep; i++){
			D[i]=0;
		}
	}else{
		for (int i=0; i<nstep; i++){
			Dlibero[i]=0;
		}
	}
}


void Resetv(bool libero){
		if(libero==1){
		for (int i=0; i<nstep; i++){
			vx[i]=0,vy[i]=0,vz[i]=0;
		}
	}else{
		for (int i=0; i<nstep; i++){
			vtheta[i]=0,vphy[i]=0,vrho[i]=0;
		}
	}
}
		

void Movexyz(double *vx, double *vy,double *vz, int i, double a, double r, double t){
	if(r<1){
		vy[i+1] =vy[i];
		vz[i+1] =vz[i];
		if(t<1){
			vx[i+1] = vx[i]-a;
		}else{vx[i+1] = vx[i]+a;}
	}
	if(r>=1 && r<2){
		vx[i+1] = vx[i];
		vz[i+1] =vz[i];
		if(t<1){
			vy[i+1]  =vy[i]-a;
		}else{vy[i+1] =vy[i]+a;}
		}
	if(r>=2){
		vy[i+1] =vy[i];
		vx[i+1] =vx[i];
		if(t<1){
			vz[i+1] = vz[i]-a;
		}else{vz[i+1] = vz[i]+a;}
	}
}


void Movetheta(double *vtheta, double *vphy, double *vrho, int i,double a, double r, double t){
	vtheta[i+1] = vtheta[i]+ a*sin(r)*cos(t); 
	vphy[i+1] = vphy[i] + a*sin(r)*sin(t);
	vrho[i+1] = vrho[i]+ a*cos(r);
}
	

double Distanza(double *vx, double *vy, double *vz, int i){
	double x;
	x = pow(vx[i],2) + pow(vy[i],2) + pow(vz[i],2);
	return x;
}

