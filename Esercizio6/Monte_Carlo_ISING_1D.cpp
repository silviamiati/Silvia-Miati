/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <iomanip>
#include "Monte_Carlo_ISING_1D.h"

using namespace std;

int main(){
 
	//system("rm output.ene.0 output.Heat.0 output.Chi.0 output.Mag.0");


  Input(); //Inizialization
  
  for(int iblk=1; iblk <= nblk; ++iblk){ //Simulation
    Reset(iblk);   //Reset block averages
    for(int istep=1; istep <= nstep; ++istep){
      Move(metro);
      Measure();
      Accumulate(); //Update block averages
    }
    Averages(iblk);   //Print results for current block 
  }
  ConfFinal(); //Write final configuration

  return 0;

}



void Input(void){

  ifstream ReadInput,ReadConfinal;

  cout << "Classic 1D Ising model             " << endl;
  cout << "Monte Carlo simulation             " << endl << endl;
  cout << "Nearest neighbour interaction      " << endl << endl;
  cout << "Boltzmann weight exp(- beta * H ), beta = 1/T " << endl << endl;
  cout << "The program uses k_B=1 and mu_B=1 units " << endl;

//Read seed for random numbers
	int p1, p2;
  ifstream Primes("Primes");
  Primes >> p1 >> p2 ;
  Primes.close();

  ifstream input("seed.in");
  input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
  rnd.SetRandom(seed,p1,p2);
  input.close();
  
//Read input informations
  ReadInput.open("input.dat");

  ReadInput >> temp;
  beta = 1.0/temp;
  cout << "Temperature = " << temp << endl;

  ReadInput >> nspin;
  cout << "Number of spins = " << nspin << endl;

  ReadInput >> J;
  cout << "Exchange interaction = " << J << endl;

  ReadInput >> h;
  cout << "External field = " << h << endl << endl;
    
  ReadInput >> metro; // if=1 Metropolis else Gibbs

  ReadInput >> nblk;

  ReadInput >> nstep;
  ReadInput >> campo;

  ReadInput >> restart;

  if(metro==1) cout << "The program perform Metropolis moves" << endl;
  else cout << "The program perform Gibbs moves" << endl;
  cout << "Number of blocks = " << nblk << endl;
  cout << "Number of steps in one block = " << nstep << endl << endl;
  ReadInput.close();


//Prepare arrays for measurements
  iu = 0; //Energy
  ic = 1; //Heat capacity
  im = 2; //Magnetization
  ix = 3; //Magnetic susceptibility
  n_props = 4; //Number of observables

  if(restart ==0){
//initial configuration
  	for (int i=0; i<nspin; ++i){
  	  if(rnd.Rannyu() >= 0.5) s[i] = 1;
  	  else s[i] = -1;
 		}
  
//Evaluate energy etc. of the initial configuration
  Measure();

//Print initial values for the potential energy and virial
  cout << "Initial energy = " << walker[iu]/(double)nspin << endl; 
	}

	if(restart==1){ 
  	cout << "Read initial configuration from file config.final : " << endl << endl;
  	ReadConfinal.open("config.final");
  	for (int i=0; i<nspin; ++i){
    	ReadConfinal >> s[i];
  	}
  	ReadConfinal.close();
	}

}// end Input()



void Move(int metro){

  int o;
  double p;

  
  for(int i=0; i<nspin; ++i){
  //Select randomly a particle (for C++ syntax, 0 <= o <= nspin-1)
    o = (int)(rnd.Rannyu()*nspin);
    p=0;
    double r=rnd.Rannyu();
    
    if(metro==1){ //Metropolis
			p = s[o];
			s[o] = -s[o];
			metropolis(o,p,r,accepted);
    }
    else{ //Gibbs sampling
  		p = s[o];
			s[o] = - s[o];
	  	gibbs(o,p,r);	
    }
		attempted = attempted + 1;
  }

} // end Move()



double Boltzmann(int sm, int ip){

  double ene = -J * sm * ( s[Pbc(ip-1)] + s[Pbc(ip+1)] ) - h * sm;
  return ene;

} //end Boltzmann() 



void metropolis(int o, double p, double r,double &accepted){

	double exe= exp(2*beta*Boltzmann(p,o));
		if(exe>=1. or (exe<1. and r<=exe)) s[o] = s[o];
		else{
		 s[o] = p;
		 accepted = accepted + 1.0;
		}

}//end metropolis() - Algoritmo di M-H -
		


void gibbs(int o, double p, double r){

	double p_g = 1. / (1. + exp(Boltzmann(p,o)*2.*beta));
		if(r>= p_g) s[o] = s[o];
		else s[o] = -s[o];

} //end gibbs() - Algoritmo di Gibbs -



void Measure(){


 double u = 0.0, m = 0.0, chi=0.0;

//cycle over spins
 if(campo==1){
  for (int i=0; i<nspin; ++i){
       u += -J * s[i] * s[Pbc(i+1)] - 0.5 * h * (s[i] + s[Pbc(i+1)]);
			 u2 += (-J * s[i] * s[Pbc(i+1)] - 0.5 * h * (s[i] + s[Pbc(i+1)]) ) * (-J * s[i] * s[Pbc(i+1)] - 0.5 * h * (s[i] + s[Pbc(i+1)]));
       chi += s[i];
  }
  walker[iu] = u; 
	walker[ic] = u*u;
  walker[ix] = chi*chi*beta;
  }else{
    for (int i=0; i<nspin; ++i){ //così calcolo solo M.
			m += s[i];
    }
  walker[im] = m;
  }
 
} // end Measure()



void Reset(int iblk){ //Reset block averages
   
   if(iblk == 1){
   	for(int i=0; i<n_props; ++i){
      glob_av[i] = 0;
      glob_av2[i] = 0;
    }
   }

   for(int i=0; i<n_props; ++i){
     blk_av[i] = 0;
   }
   blk_norm = 0;
   attempted = 0;
   accepted = 0;

} //end Rese()



void Accumulate(void){//Update block averages

   for(int i=0; i<n_props; ++i){
     blk_av[i] = blk_av[i] + walker[i];
   }
   blk_norm = blk_norm + 1.0;

} // end Accumulate()



void Averages(int iblk){ //Print results for current block
  
   ofstream Ene, Heat, Mag, Chi;
  
   cout << "Block number " << iblk << endl;
   if(metro==1)cout << "Acceptance rate " << accepted/attempted << endl << endl;
    
 	 if(campo==1){
   	Ene.open("output.ene.0",ios::app);
    stima_u = blk_av[iu]/blk_norm/(double)nspin; //Energy
    glob_av[iu]  += stima_u;
    glob_av2[iu] += stima_u*stima_u;
    err_u=Error(glob_av[iu],glob_av2[iu],iblk);
    Ene << iblk << " "  << temp << " " << glob_av[iu]/(double)iblk << " " << err_u << endl;
    Ene.close();

    Heat.open("output.Heat.0",ios::app);
    stima_c = beta*beta*(blk_av[ic]/(blk_norm)/(double)(nspin) - (blk_av[iu]*blk_av[iu]/(blk_norm*blk_norm))/(double)nspin);//Heat
    glob_av[ic]  += stima_c;
    glob_av2[ic] += stima_c*stima_c;
    err_c=Error(glob_av[ic],glob_av2[ic],iblk);
    Heat << iblk << " " << temp << " " << glob_av[ic]/(double)iblk << " " << err_c << endl;
    Heat.close();

    Chi.open("output.Chi.0",ios::app);
    stima_x = blk_av[ix]/blk_norm/(double)nspin; 
    glob_av[ix]  += stima_x;
    glob_av2[ix] += stima_x*stima_x;
    err_x=Error(glob_av[ix],glob_av2[ix],iblk);
    Chi << iblk << " " << temp << " " << glob_av[ix]/(double)iblk << " " << err_x << endl; // Susceptibility
    Chi.close();
		}else{
    Mag.open("output.Mag.0",ios::app);
    stima_m = blk_av[im]/blk_norm/(double)nspin; 
    glob_av[im]  += stima_m;
    glob_av2[im] += stima_m*stima_m;
    err_m=Error(glob_av[im],glob_av2[im],iblk);
    Mag << iblk << " " << temp  << " " << glob_av[im]/(double)iblk << " " << err_m << endl; // Magnetization
    Mag.close();
	}

} // end Averages()


void ConfFinal(void){

  ofstream WriteConf;

  cout << "Print final configuration to file config.final " << endl << endl;
  WriteConf.open("config.final");
  for (int i=0; i<nspin; ++i){
    WriteConf << s[i] << endl;
  }
  WriteConf.close();

  rnd.SaveSeed();

} // end ConfFinal()



int Pbc(int i){  //Algorithm for periodic boundary conditions

    if(i >= nspin) i = i - nspin;
    else if(i < 0) i = i + nspin;
    return i;

} // end Pbc()



double Error(double sum, double sum2, int iblk){

    return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)iblk);

}// end Error()

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/	

