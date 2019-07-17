#include "ese10.h"
#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <iomanip>
#include <string>
#include <cstdlib>
#include <stdlib.h>
#include <vector>
#include <algorithm>    // std::random_shuffle
//#include <ctime>        // std::time
//#include <valarray>
//#include <cstddef> 
//#include <unistd.h>
//#include <bits/stdc++.h> 
//#include <utility>


using namespace std;


int main(){

system("rm fitness.dat");

Input(); 
	 
Move(Pos,beta);

cout << "Il numero di interazioni è: " << tot << endl;

return 0;

}

void Input(void){

 ifstream ReadInput,ReadConfinal;

	int p1, p2; //Read seed for random numbers
	ifstream Primes("Primes");
  Primes >> p1 >> p2 ;
  Primes.close();

  ifstream input("seed.in");
  input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
  rnd.SetRandom(seed,p1,p2);
  input.close();

	ReadInput.open("input.dat");

	
	ReadInput >> tempe;
		beta = 1./tempe;
	ReadInput >> nstep;
	ReadInput >> first;
	ReadInput >> ncity;
	ReadInput >> rho; // raggio cerchio
	ReadInput >> restart;
	ReadInput >> percorso;

	

  ReadInput.close();

	if(first==1){
		ofstream Posi;
		double x=0.0,y=0.0;
		if(percorso==1){
			double theta=0.0;
			Posi.open("Posizioni_cerchio.dat", ios::app);
				for(int i=0; i<ncity; i++){
						theta=rnd.Rannyu(0,2.*pi);
						x=cos(theta);
						y=sin(theta);
						Posi << x << " " << y<< endl;
					}
		}else{
			Posi.open("Posizioni_quadrato.dat", ios::app);
				for(int i=0; i<ncity; i++){
					x = rnd.Rannyu();
					y = rnd.Rannyu();
					Posi << x << " " << y << endl;
				}
		}	
		Posi.close();
	}

	fstream pos;
	double a=0.0;
	if(percorso==1) pos.open("Posizioni_cerchio.dat");
	else pos.open("Posizioni_quadrato.dat"); // file da aprire per leggere le posizioni 

  for(int i=0; i<ncity; i++){
		pos >> a;
		xR[i]=a;
		pos >> a;
		yR[i]=a;  
  }
	for(int i=0; i<ncity; i++){
		Pos.push_back(make_pair(xR[i],yR[i])); // creazione del vettore di pair posizione
	}
	pos.close();

	if(restart){ 
   cout << "Read initial configuration from file config.final : " << endl << endl;
   ReadConfinal.open("config.final");
	 double x=0.0,y=0.0;
   for (int i=0; i<ncity; i++){
   	ReadConfinal >> x >> y;
		Pos[i].first=x;
		Pos[i].second=y;	
   }
   ReadConfinal.close();
	}

} // end Input()



double Fitness(vector<pair<double,double>> R){  
  
	L=0;
	for(int j=0; j<ncity; j++){
		if(j==ncity-1) L = L + sqrt(pow(R[ncity-1].first- R[1].first,2) +pow( R[ncity-1].second- R[1].second,2));
		else L = L + sqrt(pow(R[j+1].first - R[j].first,2) + pow(R[j+1].second - R[j].second,2));
	}
	return L;

}// end Fitness()
	


void Move(vector <pair<double,double> > &v , double beta){

ofstream Fitnes;
Fitnes.open("fitness.dat",ios::app); //attenzione il nome del file è sempre lo stesso sia per il cerchio che per il quadrato

	while(tempe>0.05){
		for(int istep=1; istep <= nstep; ++istep){
			double p=0, energy_old=0, energy_new=0;
			int rand = v.size()*rnd.Rannyu(); 
			int rand1 = v.size()*rnd.Rannyu(); 
	
			energy_old = Fitness(v);
		
			iter_swap(v.begin() + rand, v.begin() + rand1); 

			energy_new = Fitness(v);
	
			p = exp((1./tempe)*(energy_old-energy_new)); // M-H algoritmo
    	if(p < rnd.Rannyu()){ 
   			iter_swap(v.begin() + rand, v.begin() + rand1);
 	  	}
	  }
		tot++;
		Fitnes << L << endl;
		ConfFinal();
		tempe *= 0.98; // cooling rate
	}
Fitnes.close();
	
}  // end Move() 



void ConfFinal(void){

  ofstream WriteConf, WriteSeed;
  cout << "Print final configuration to file config.final " << endl << endl;
  WriteConf.open("config.final");//attenzione il nome del file è sempre lo stesso sia per il cerchio che per il quadrato
  for (int i=0; i<ncity; ++i){
    WriteConf << Pos[i].first << " " << Pos[i].second << endl;
  }
  WriteConf.close();
  rnd.SaveSeed();

} // end ConfFinal()
















		



