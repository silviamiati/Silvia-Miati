/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
/*///////////////////////////////////////////////////////////////////
AVVERTENZE PER L'UTENTE: 
In riga di compilazione scrivere:
- 1 x LIQUIDO;
- 2 x GAS;
- 3 X SOLIDO;
Scegliere poi quale elemento si vuole studiare:
- 1 x ARGON;
- 2 x KRYPTON;
*////////////////////////////////////////////////////////////////////

#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <iomanip>
#include "Monte_Carlo_NVT.h"
#include <string>

using namespace std;

int main(int argn,  char *argv[]){ 
	
	if(argn-1<2) {
		printf("numero inesatto di argomenti\n");
		exit(1);
	}

	F = atoi(argv[1]);
	E = atoi(argv[2]);


  Input(); //Inizialization
	Acceptance();
  int nconf = 1;
  for(int iblk=1; iblk <= nblk; ++iblk) { //Simulation
    Reset(iblk);   //Reset block averages
    for(int istep=1; istep <= nstep; ++istep) {
      Move();
      Measure();
      Accumulate(); //Update block averages
      if(istep%10 == 0){
       ConfXYZ(nconf);//Write actual configuration in XYZ format //Commented to avoid "filesystem full"! 
       nconf += 1;
      }
    }
    Averages(iblk);   //Print results for current block
  }
  ConfFinal(); //Write final configuration

  return 0;

}


void Input(void){

  ifstream ReadInput,ReadConf,ReadEle,ReadConfinal;
  
  cout << "Classic Lennard-Jones fluid        " << endl;
  cout << "Monte Carlo simulation             " << endl << endl;
  cout << "Interatomic potential v(r) = 4 * [(1/r)^12 - (1/r)^6]" << endl << endl;
  cout << "Boltzmann weight exp(- beta * sum_{i<j} v(r_ij) ), beta = 1/T " << endl << endl;
  cout << "The program uses Lennard-Jones units " << endl;

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
	if(F==1) ReadInput.open("inputL.dat");
	else if (F==2) ReadInput.open("inputG.dat");
	else if (F==3) ReadInput.open("inputS.dat");

  ReadInput >> temp;
  beta = 1.0/temp;
  cout << "Temperature = " << temp << endl;

  ReadInput >> npart;
  cout << "Number of particles = " << npart << endl;

  ReadInput >> rho;
  cout << "Density of particles = " << rho << endl;
  vol = (double)npart/rho;
  box = pow(vol,1.0/3.0);
  cout << "Volume of the simulation box = " << vol << endl;
  cout << "Edge of the simulation box = " << box << endl;

  ReadInput >> rcut;
  cout << "Cutoff of the interatomic potential = " << rcut << endl << endl;
    
  //Tail corrections for potential energy and pressure
  vtail = (8.0*pi*rho)/(9.0*pow(rcut,9)) - (8.0*pi*rho)/(3.0*pow(rcut,3));
  ptail = (32.0*pi*rho)/(9.0*pow(rcut,9)) - (16.0*pi*rho)/(3.0*pow(rcut,3));
  cout << "Tail correction for the potential energy = " << vtail << endl;
  cout << "Tail correction for the virial           = " << ptail << endl; 

  ReadInput >> nblk;

  ReadInput >> nstep;

  ReadInput >> Realsimulation;	

	ReadInput >> restart;
  
  cout << "The program perform Metropolis moves with uniform translations" << endl;
  cout << "Moves parameter = " << delta << endl;
  cout << "Number of blocks = " << nblk << endl;
  cout << "Number of steps in one block = " << nstep << endl << endl;
  ReadInput.close();

	if(E==1) ReadEle.open("Argon.dat");
	else if (E==2) ReadEle.open("Krypton.dat");
	ReadEle >> sigma;
		sigma *= 1E-9;
	ReadEle >> kb;
		kb *= 1E-23;
	ReadEle >> epsilon_kb;
	ReadEle.close();

  //Prepare arrays for measurements
  iv = 0; //Potential energy
  iw = 1; //Virial
 
  n_props = 2; //Number of observables

  igofr = 2; //measurement of g(r)
  nbins = 100;
  n_props = n_props + nbins;
  bin_size = (box/2.0)/(double)nbins;

	if(restart){ 
  cout << "Read initial configuration from file config.final : " << endl << endl;
  ReadConfinal.open("config.final");
  for (int i=0; i<npart; ++i){
    ReadConfinal >> x[i] >> y[i] >> z[i];
  	x[i] = Pbc( x[i] * box );
    y[i] = Pbc( y[i] * box );
    z[i] = Pbc( z[i] * box );
  }
  ReadConfinal.close();
	}else{	
  cout << "Read initial configuration from file config.0 " << endl << endl; //Read initial configuration
  ReadConf.open("config.0");
  	for (int i=0; i<npart; ++i){
    ReadConf >> x[i] >> y[i] >> z[i];
    x[i] = Pbc( x[i] * box );
    y[i] = Pbc( y[i] * box );
    z[i] = Pbc( z[i] * box );
  	}
  ReadConf.close();
	}
  
  Measure(); //Evaluate potential energy and virial of the initial configuration

  //Print initial values for the potential energy and virial
  cout << "Initial potential energy (with tail corrections) = " << walker[iv]/(double)npart + vtail << endl;
  cout << "Virial                   (with tail corrections) = " << walker[iw]/(double)npart + ptail << endl;
  cout << "Pressure                 (with tail corrections) = " << rho * temp + (walker[iw] + (double)npart * ptail) / vol << endl << endl;

} // end Input()



void Move(void){

  int o;
  double p, energy_old, energy_new;
  double xold, yold, zold, xnew, ynew, znew;

  for(int i=0; i<npart; ++i){
  //Select randomly a particle (for C++ syntax, 0 <= o <= npart-1)
    o = (int)(rnd.Rannyu()*npart);

  //Old
    xold = x[o];
    yold = y[o];
    zold = z[o];

    energy_old = Boltzmann(xold,yold,zold,o);

  //New
    xnew = Pbc( x[o] + delta*(rnd.Rannyu() - 0.5) );
    ynew = Pbc( y[o] + delta*(rnd.Rannyu() - 0.5) );
    znew = Pbc( z[o] + delta*(rnd.Rannyu() - 0.5) );

    energy_new = Boltzmann(xnew,ynew,znew,o);

  //Metropolis test
    p = exp(beta*(energy_old-energy_new));
    if(p >= rnd.Rannyu()){
    //Update
       x[o] = xnew;
       y[o] = ynew;
       z[o] = znew;
    
       accepted = accepted + 1.0;
    }
    attempted = attempted + 1.0; 
  }

} // end Move() - Algoritmo di M-H -



double Boltzmann(double xx, double yy, double zz, int ip){

  double ene=0.0;
  double dx, dy, dz, dr;

  for (int i=0; i<npart; ++i){
    if(i != ip) {
// distance ip-i in pbc
      dx = Pbc(xx - x[i]);
      dy = Pbc(yy - y[i]);
      dz = Pbc(zz - z[i]);

      dr = dx*dx + dy*dy + dz*dz;
      dr = sqrt(dr);

      if(dr < rcut)  ene += 1.0/pow(dr,12) - 1.0/pow(dr,6);   
    }
  }

  return 4.0*ene;

} // end Boltzmann()



void Measure(){

  int bin;
  double v = 0.0, w = 0.0;
  double rdr= 0.0, r= 0.0;
  double vij, wij;
  double dx, dy, dz, dr;

  ofstream Pres,Epot;  // Scommentare se si vuole vedere i valori istantanei

  Epot.open("gas_epot.dat",ios::app);
  Pres.open("gas_pres.dat",ios::app);
 

  for (int k=igofr; k<igofr+nbins; ++k) walker[k]=0.0; //reset the hystogram of g(r)

  for (int i=0; i<npart-1; ++i){ //cycle over pairs of particles
    for (int j=i+1; j<npart; ++j) {

     dx = Pbc(x[i] - x[j]); // distance i-j in pbc
     dy = Pbc(y[i] - y[j]);
     dz = Pbc(z[i] - z[j]);

     dr = dx*dx + dy*dy + dz*dz;
     dr = sqrt(dr);
     
     rdr= x[i]*x[i] + y[i]*y[i] + z[i]*z[i];
     rdr= sqrt(rdr);
     r= x[j]*x[j] + y[j]*y[j] + z[j]*z[j];
     r= sqrt(r);

     bin = dr/bin_size;

     walker[igofr+bin] += 2.; //update of the histogram of g(r) 

     if(dr < rcut){
       vij = 1.0/pow(dr,12) - 1.0/pow(dr,6);
       wij = 1.0/pow(dr,12) - 0.5/pow(dr,6);

       v += vij; // contribution to energy and virial
       w += wij;
     }
    }          
  }

  walker[iv] = 4.0 * v;
  walker[iw] = 48.0 * w / 3.0; 

  stimais_pot = 4.0 * v/(double)npart + vtail; 
  stimais_pres = rho * temp + ((48.0 * w/3.0)+ ptail * (double)npart) / vol;

  Pres << stimais_pres << endl;
  Epot << stimais_pot  << endl; 

  Epot.close();
  Pres.close();

} // Measure()



void Reset(int iblk){ //Reset block averages

   if(iblk == 1){
       for(int i=0; i<n_props; ++i){
           glob_av[i] = 0;
           glob_av2[i] = 0;
       }
   }

   for(int i=0; i<n_props; ++i) blk_av[i] = 0;
   
   blk_norm = 0;
   attempted = 0;
   accepted = 0;

} // end Reset()



void Accumulate(void){ //Update block averages

   for(int i=0; i<n_props; ++i) blk_av[i] = blk_av[i] + walker[i];
  
   blk_norm = blk_norm + 1.0;

} // end Accumulate()



void Averages(int iblk){ //Print results for current block

   double r, gdir;
   ofstream Gofr, Gave, Epot, Pres;
  
   cout << "Block number " << iblk << endl;
   cout << "Acceptance rate " << accepted/attempted << endl << endl;
    	
	 if(Realsimulation==0){

    Epot.open("output.epot.0",ios::app); Scommentare se si vuole vedere questi valori in unità riscalate
 	  Pres.open("output.pres.0",ios::app);
    Gofr.open("output.gofr.0",ios::app);
    Gave.open("output.gave.0",ios::app);
    
    stima_pot = blk_av[iv]/blk_norm/(double)npart + vtail; //Potential energy
    glob_av[iv] += stima_pot;
    glob_av2[iv] += stima_pot*stima_pot;
    err_pot=Error(glob_av[iv],glob_av2[iv],iblk);
    
    stima_pres = rho * temp + (blk_av[iw]/blk_norm + ptail * (double)npart) / vol; //Pressure
    glob_av[iw] += stima_pres;
    glob_av2[iw] += stima_pres*stima_pres;
    err_press=Error(glob_av[iw],glob_av2[iw],iblk);

		Epot <<iblk << " " << stima_pot << " " << glob_av[iv]/(double)iblk << " " << err_pot <<endl; //Potenziale
    Pres << iblk << " " << stima_pres << " " << glob_av[iw]/(double)iblk << " " << err_press <<endl;//Pressione

		for(int i=igofr; i<n_props; i++){ //sto calcolando la g(r) (Gofr), l'ultimo blocco lo stampo con Gave!
			r = double(i-igofr)*bin_size;
			double vr = (4./3.)*pi*(pow(r+bin_size,3) - pow(r,3));
			double norm = rho*npart*vr; // normalizzazione corretta
			gdir = blk_av[i]/(blk_norm*norm);
			glob_av[i] += gdir;
			glob_av2[i] += gdir*gdir;
			err_gdir= Error(glob_av[i],glob_av2[i],iblk);
			r = (double)(i-igofr)*bin_size + bin_size*0.5;
			Gofr << iblk << " " << i-igofr << " " << r << " " << gdir <<  " " << glob_av[i]/(double)iblk << " " << err_gdir << endl;
			if (iblk == nblk) Gave << i-igofr << " " << r << " " << gdir << " " << glob_av[i]/(double)iblk << " " << err_gdir << endl;
		}

		Epot.close();
    Pres.close();
    Gofr.close();
		Gofr.close();
		Gave.close();

	 }else if(Realsimulation==1){
    
		double pot,press,pot2,press2;

		ofstream is_Gdir, is_Epot, is_Press, is_Gdir;
		is_Gdir.open("ELE.gofr.dat" , ios::app);
		is_Epot.open("ELE.epot.dat",ios::app);
		is_Press.open("ELE.press.dat",ios::app);

		stima_pot = blk_av[iv]* epsilon_kb*kb/blk_norm/(double)npart + vtail; //Potential energy
    glob_av[iv] += stima_pot;
    glob_av2[iv] += stima_pot*stima_pot;
    err_pot=Error(glob_av[iv],glob_av2[iv],iblk);
    
    stima_pres = (rho * temp + (blk_av[iw]/blk_norm + ptail * (double)npart) / vol) * (epsilon_kb*kb/pow(sigma,3)); //Pressure
    glob_av[iw] += stima_pres;
    glob_av2[iw] += stima_pres*stima_pres;
    err_press=Error(glob_av[iw],glob_av2[iw],iblk);

		is_Epot << pot/(double)iblk << " " << Error(pot,pot2,iblk) << " " << iblk << endl;
		is_Press << press/(double)iblk << " " << Error(press,press2,iblk) << " " << iblk << endl;

		is_Epot.close();
		is_Press.close();


    for(int i=igofr; i<n_props; i++){ //sto calcolando la g(r) (Gofr), l'ultimo blocco lo stampo con Gave!
			r = double(i-igofr)*bin_size;
			double vr = (4./3.)*pi*(pow(r+bin_size,3) - pow(r,3));
			double norm = rho*npart*vr; // normalizzazione corretta
			gdir = blk_av[i]/(blk_norm*norm);
			glob_av[i] += gdir;
			glob_av2[i] += gdir*gdir;
			err_gdir= Error(glob_av[i],glob_av2[i],iblk);
			r = (double)(i-igofr)*bin_size + bin_size*0.5;
			if(iblk == nblk and Realsimulation == 1)is_Gdir << iblk << " " << i-igofr << " " << r*sigma << " " << glob_av[i]/(double)iblk << " " << err_gdir<< endl;
	  }
  
	  is_Gdir.close();	
		}

} // end Averages()



void ConfFinal(void){

  ofstream WriteConf, WriteSeed;

  cout << "Print final configuration to file config.final " << endl << endl;
  WriteConf.open("config.final");
  for (int i=0; i<npart; ++i){
    WriteConf << x[i]/box << "   " <<  y[i]/box << "   " << z[i]/box << endl;
  }
  WriteConf.close();

  rnd.SaveSeed();

} // end ConfFinal()



void ConfXYZ(int nconf){ //Write configuration in .xyz format
 
  ofstream WriteXYZ;

  WriteXYZ.open("frames/config_" + to_string(nconf) + ".xyz");
  WriteXYZ << npart << endl;
  WriteXYZ << "This is only a comment!" << endl;
  for (int i=0; i<npart; ++i){
    WriteXYZ << "LJ  " << Pbc(x[i]) << "   " <<  Pbc(y[i]) << "   " << Pbc(z[i]) << endl;
  }
  WriteXYZ.close();

} // end ConfXYZ()



double Pbc(double r){  //Algorithm for periodic boundary conditions with side L=box

    return r - box * rint(r/box);

} // end Pbc()



double Error(double sum, double sum2, int iblk){

    return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)iblk);

} // end Error()



double Acceptance(){
		 
		delta=0.01;
		double p=0;
		
		while(p<0.4 or p>0.65){	
			delta += 0.1;
			accepted=0;
			attempted=0;
			for(int i=0; i<1000; i++){
				Move();			
			}	
			p = accepted/attempted;
			cout << p << endl;	
		}
	return delta;

} // end Acceptance()


/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
