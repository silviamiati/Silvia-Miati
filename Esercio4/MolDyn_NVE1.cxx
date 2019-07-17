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
#include <stdlib.h>     // srand, rand: to generate random number
#include <iostream>     // cin, cout: Standard Input/Output Streams Library
#include <fstream>      // Stream class to both read and write from/to files.
#include <cmath>   
#include <string>     // rint, pow
#include "MolDyn_NVE.h"
#include "random.h"
#include <string>

using namespace std;


int main(int argn, const char *argv[]){ 

system("rm output_epot.dat output_ekin.dat output_temp.dat output_etot.dat output_pres.dat");
system("rm output.etot.0 output.ekin.0 output.epot.0 output.temp.0 output.pres.0 output.gofr.0 output.gave.0");

	if(argn-1<2) {
		printf("numero inesatto di argomenti\n");
		exit(1);
	}

	F = atoi(argv[1]);
	E = atoi(argv[2]);

	Input();             //Inizialization

  int nconf = 1;
  for(int istep=1; istep <= nstep; ++istep){
     Move();           //Move particles with Verlet algorithm
     if(istep%iprint == 0) cout << "Number of time-steps: " << istep << endl;
     if(istep%10 == 0){
        Measure();     //Properties measurement
//      ConfXYZ(nconf);//Write actual configuration in XYZ format //Commented to avoid "filesystem full"! 
      	Accumulate();
        if(istep%L==0) {Averages(block);
				block++;
				}
        nconf += 1;
     }
  }
  ConfFinal();         //Write final configuration to restart
 
 return 0;

}



void Input(void){ //Prepare all stuff for the simulation
   
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

  rnd.SaveSeed();
	
  ifstream ReadInput,ReadConf, ReadConfpre, ReadConfinal, ReadEle;
  double ep, ek, pr, et, vir;

  cout << "Classic Lennard-Jones fluid        " << endl;
  cout << "Molecular dynamics simulation in NVE ensemble  " << endl << endl;
  cout << "Interatomic potential v(r) = 4 * [(1/r)^12 - (1/r)^6]" << endl << endl;
  cout << "The program uses Lennard-Jones units " << endl;

  if(F==1) ReadInput.open("inputL.dat");
	else if (F==2) ReadInput.open("inputG.dat");
	else if (F==3) ReadInput.open("inputS.dat");

  ReadInput >> temp;

  ReadInput >> npart;
  cout << "Number of particles = " << npart << endl;

  ReadInput >> rho;
  cout << "Density of particles = " << rho << endl;
  vol = (double)npart/rho;
  cout << "Volume of the simulation box = " << vol << endl;
  box = pow(vol,1.0/3.0);
  cout << "Edge of the simulation box = " << box << endl;

  ReadInput >> rcut;
	ReadInput >> delta;
  ReadInput >> nblk;
  ReadInput >> nstep; 
		L = nstep/nblk;
		cout << L << " " << "L" << nblk << " " << nstep << endl;
  ReadInput >> iprint;
	ReadInput >> Realsimulation;
	ReadInput >> restart;

  cout << "The program integrates Newton equations with the Verlet method " << endl;
  cout << "Time step = " << delta << endl;
  cout << "Number of steps = " << nstep << endl << endl;
  ReadInput.close();

	if(E==1) ReadEle.open("Argon.dat");
	else if (E==2) ReadEle.open("Krypton.dat");
	ReadEle >> sigma;
		sigma *= 1E-9;
	ReadEle >> kb;
		kb *= 1E-23;
	ReadEle >> epsilon_kb;
	ReadEle.close();

//Prepare array for measurements
  iv = 0; //Potential energy
  ik = 1; //Kinetic energy
  ie = 2; //Total energy
  it = 3; //Temperature
	ip = 4; //Pression
  n_props = 5; //Number of observables

//prepare for g(r)
	igofr = 5;
  nbins = 100;
  n_props = n_props + nbins;
  bin_size = (box/2.0)/(double)nbins;

//Read initial configuration
if(restart==0){
  cout << "Read initial configuration from file config.0 : " << endl << endl;
  ReadConf.open("config.0");
  for (int i=0; i<npart; ++i){
    ReadConf >> x[i] >> y[i] >> z[i];
    x[i] = x[i] * box;
    y[i] = y[i] * box;
    z[i] = z[i] * box;
  }
  ReadConf.close();

	//Prepare initial velocities
	cout << "Prepare random velocities with center of mass velocity equal to zero " << endl << endl;
   double sumv[3] = {0.0, 0.0, 0.0};
   for (int i=0; i<npart; ++i){
     vx[i] = rnd.Rannyu() - 0.5;
     vy[i] = rnd.Rannyu() - 0.5;
     vz[i] = rnd.Rannyu() - 0.5;

     sumv[0] += vx[i];
     sumv[1] += vy[i];
     sumv[2] += vz[i];
   }
   for (int idim=0; idim<3; ++idim) sumv[idim] /= (double)npart;
   double sumv2 = 0.0, fs;
   for (int i=0; i<npart; ++i){
     vx[i] = vx[i] - sumv[0];
     vy[i] = vy[i] - sumv[1];
     vz[i] = vz[i] - sumv[2];

     sumv2 += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
   }
   sumv2 /= (double)npart;

   fs = sqrt(3 * temp / sumv2);   // fs = velocity scale factor 
   for (int i=0; i<npart; ++i){
     vx[i] *= fs;
     vy[i] *= fs;
     vz[i] *= fs;

     xold[i] = x[i] - vx[i] * delta;
     yold[i] = y[i] - vy[i] * delta;
     zold[i] = z[i] - vz[i] * delta; 
    
   } 

} else if(restart==1){
 
  cout << "Read initial configuration from file config.final : " << endl << endl;
  ReadConfinal.open("config.final");
  for (int i=0; i<npart; ++i){
    ReadConfinal >> x[i] >> y[i] >> z[i];
    x[i] = x[i] * box;
    y[i] = y[i] * box;
    z[i] = z[i] * box;
  }
  ReadConfinal.close();

  cout << "Read initial configuration from file config.prec : " << endl << endl;
  ReadConfpre.open("config.prec");
  for (int i=0; i<npart; ++i){
    ReadConfpre >> xold[i] >> yold[i] >> zold[i];
    xold[i] = xold[i] * box;
    yold[i] = yold[i] * box;
    zold[i] = zold[i] * box;
  }
  ReadConfpre.close();

	double xnew=0, ynew=0, znew=0, fxn[m_part], fyn[m_part], fzn[m_part];
  double sumv2_n = 0.0, fs_n;

  for(int i=0; i<npart; ++i){ //Force acting on particle i
    fxn[i] = Force(i,0);
    fyn[i] = Force(i,1);
    fzn[i] = Force(i,2);
  }

  for(int i=0; i<npart; ++i){ //Verlet integration scheme

    xnew = Pbc( 2.0 * x[i] - xold[i] + fxn[i] * pow(delta,2) );
    ynew = Pbc( 2.0 * y[i] - yold[i] + fyn[i] * pow(delta,2) );
    znew = Pbc( 2.0 * z[i] - zold[i] + fzn[i] * pow(delta,2) ); 

    vx_n[i] = Pbc(xnew - xold[i])/(2.0 * delta);
    vy_n[i] = Pbc(ynew - yold[i])/(2.0 * delta);
    vz_n[i] = Pbc(znew - zold[i])/(2.0 * delta); 
	
    sumv2_n += vx_n[i]*vx_n[i] + vy_n[i]*vy_n[i] + vz_n[i]*vz_n[i]; 

    x_n[i]=xnew;
    y_n[i]=ynew;
    z_n[i]=znew; 

    }

    sumv2_n /= double(npart);
    double T= sumv2_n/double(3.);
   
    fs_n = sqrt(temp/T);

		for(int i=0; i<npart; ++i){
			vx_n[i] *= fs_n;
			vy_n[i] *= fs_n;
			vz_n[i] *= fs_n;
		 
			xold[i] = Pbc(x_n[i] - 2*delta*vx_n[i]);
			yold[i] = Pbc(y_n[i] - 2*delta*vy_n[i]);
			zold[i] = Pbc(z_n[i] - 2*delta*vz_n[i]);
		}
}
   return;

} // end Input()



void Move(void){ //Move particles with Verlet algorithm

  double xnew, ynew, znew, fx[m_part], fy[m_part], fz[m_part];

  for(int i=0; i<npart; ++i){ //Force acting on particle i
    fx[i] = Force(i,0);
    fy[i] = Force(i,1);
    fz[i] = Force(i,2);
  }

  for(int i=0; i<npart; ++i){ //Verlet integration scheme

    xnew = Pbc( 2.0 * x[i] - xold[i] + fx[i] * pow(delta,2) );
    ynew = Pbc( 2.0 * y[i] - yold[i] + fy[i] * pow(delta,2) );
    znew = Pbc( 2.0 * z[i] - zold[i] + fz[i] * pow(delta,2) );

    vx[i] = Pbc(xnew - xold[i])/(2.0 * delta);
    vy[i] = Pbc(ynew - yold[i])/(2.0 * delta);
    vz[i] = Pbc(znew - zold[i])/(2.0 * delta);

    xold[i] = x[i];
    yold[i] = y[i];
    zold[i] = z[i];

    x[i] = xnew;
    y[i] = ynew;
    z[i] = znew;
  }
  return;

}// end Move()



double Force(int ip, int idir){ //Compute forces as -Grad_ip V(r)
  double f=0.0;
  double dvec[3], dr;

  for (int i=0; i<npart; ++i){
    if(i != ip){
      dvec[0] = Pbc( x[ip] - x[i] );  // distance ip-i in pbc
      dvec[1] = Pbc( y[ip] - y[i] );
      dvec[2] = Pbc( z[ip] - z[i] );

      dr = dvec[0]*dvec[0] + dvec[1]*dvec[1] + dvec[2]*dvec[2];
      dr = sqrt(dr);

      if(dr < rcut){
        f += dvec[idir] * (48.0/pow(dr,14) - 24.0/pow(dr,8)); // -Grad_ip V(r)
      }
    }
  }
  
  return f;

}//end Force()



void Measure(){ //Properties measurement
	
  int bin;
  double v, t, vij, pij, P;
  double dx, dy, dz, dr;
/* ofstream Epot, Ekin, Etot, Temp, Pres;

  Epot.open("output_epot.dat",ios::app);
  Ekin.open("output_ekin.dat",ios::app);
  Temp.open("output_temp.dat",ios::app);
  Etot.open("output_etot.dat",ios::app);
  Pres.open("output_pres.dat",ios::app);
*/
	for (int k=igofr; k<igofr+nbins; ++k) walker[k]=0.0;

  v = 0.0; //reset observables
  t = 0.0;
  P = 0.0;

//cycle over pairs of particles
  for (int i=0; i<npart-1; ++i){
    for (int j=i+1; j<npart; ++j){

     dx = Pbc( x[i] - x[j] );
     dy = Pbc( y[i] - y[j] );
     dz = Pbc( z[i] - z[j] );

     dr = dx*dx + dy*dy + dz*dz;
     dr = sqrt(dr);

		 bin = dr/(bin_size);
	
		 //update of the histogram of g(r) 
     walker[igofr+bin] += 2.;

     if(dr < rcut){
       vij = 4.0/pow(dr,12) - 4.0/pow(dr,6);
			 pij = (48./3.)/pow(dr,12) - (24.0/3.)/pow(dr,6);

//Pressione
       P += pij;

//Potential energy
       v += vij;
     }
    }          
  }

//Kinetic energy
  for (int i=0; i<npart; ++i) t += 0.5 * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);
   
    stima_pot = v/(double)npart; //Potential energy
    stima_kin = t/(double)npart; //Kinetic energy
    stima_temp = (double)2.0/(double)3.0 * t/(double)npart; //Temperature
    stima_etot = (t+v)/(double)npart; //Total enery
    

		walker[iv] = stima_pot;
		walker[ik] = stima_kin;
		walker[ie] = stima_etot;
		walker[it] = stima_temp;
		
		stima_pres = rho*walker[it] + P/(vol); //Pression
		walker[ip] = stima_pres;
/*
    Epot << stima_pot  << endl;
    Ekin << stima_kin  << endl;
    Temp << stima_temp << endl;
    Etot << stima_etot << endl;
    Pres << stima_pres << endl;

   
 		Etot.clear();
    Epot.close();
    Ekin.close();
    Temp.close();
    Etot.close();
    Pres.close();
*/
    return;
} // end Measure()



void Accumulate(void){ //Update block averages

   for(int i=0; i<n_props; ++i){
     blk_av[i] = blk_av[i] + walker[i];
   }
   blk_norm = blk_norm + 1.0;

} // end Accumulate()



void ConfFinal(void){ //Write final configuration

  ofstream WriteConf;

  cout << "Print final configuration to file config.final " << endl << endl;
  WriteConf.open("config.final");

  for (int i=0; i<npart; ++i){
    WriteConf << x[i]/box << "   " <<  y[i]/box << "   " << z[i]/box << endl;
  }
  WriteConf.close();

	ofstream WriteConfprec;
	cout << "Print final configuration to file config.prec " << endl << endl;
	WriteConfprec.open("config.prec");

	 for (int i=0; i<npart; ++i){
   		 WriteConfprec << xold[i]/box << "   " <<  yold[i]/box << "   " << zold[i]/box << endl;
 	 }
  	WriteConfprec.close(); 	
  
  return;

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



void Averages(int iblk){ //Print results for current block
		
    double r, gdir;

if(Realsimulation==0){ 
	
    ofstream Gofr, Gave, Epot, Pres, Etot, Kin, Tem; 
    cout << "Block number " << iblk << endl;

		Etot.open("output.etot.0",ios::app);
		Kin.open("output.ekin.0",ios::app);
    Epot.open("output.epot.0",ios::app);
		Tem.open("output.temp.0",ios::app);
    Pres.open("output.pres.0",ios::app);
    Gofr.open("output.gofr.0",ios::app);
    Gave.open("output.gave.0",ios::app);
    
    stima_k = blk_av[ik] /blk_norm; //Kinetical energy
    glob_av[ik] += stima_k ; 
    glob_av2[ik] += stima_k*stima_k;
    err_k=Error(glob_av[ik],glob_av2[ik],iblk);
    
    stima_et = blk_av[ie]/blk_norm; // Energia totale
    glob_av[ie] += stima_et; 
    glob_av2[ie] += stima_et*stima_et;
    err_et=Error(glob_av[ie],glob_av2[ie],iblk);

		stima_p = blk_av[ip]/blk_norm; // Energia potenziale
		glob_av[iv] += stima_p; 
    glob_av2[iv] += stima_p*stima_p;
    err_p=Error(glob_av[iv],glob_av2[iv],iblk);

		stima_t = blk_av[it]/blk_norm; // temperatura
		glob_av[it] += stima_t;
    glob_av2[it] += stima_t*stima_t;
    err_t=Error(glob_av[it],glob_av2[it],iblk);

		stima_pr = blk_av[ip]/blk_norm; // Pressione
		glob_av[ip] += stima_pr; 
    glob_av2[ip] += stima_pr*stima_pr;
    err_pr=Error(glob_av[ip],glob_av2[ip],iblk);

    Epot << " "  << iblk << " " << stima_p << " "  << glob_av[iv]/(double)iblk << " "  << err_p << endl;
    Tem << " " << iblk <<  " " << stima_t << " "  << glob_av[it]/(double)iblk << " " << err_t << endl;
		Pres << " " << iblk <<  " " << stima_pr << " "  << glob_av[ip]/(double)iblk << " " << err_pr<< endl;
		Kin << " " << iblk <<  " " << stima_k << " "  << glob_av[ik]/(double)iblk << " " << err_k << endl;
		Etot << " " << iblk <<  " " << stima_et << " "  << glob_av[ie]/(double)iblk << " " << err_et << endl;

		Epot.close();
		Tem.close();
		Pres.close();
		Kin.close();
		Etot.close();
		
		for(int i=igofr; i<n_props; i++){
			r = double(i-igofr)*bin_size;
			double vr = (4./3.)*pi*(pow(r+bin_size,3) - pow(r,3));
			double norm = rho*npart*vr;
			gdir = blk_av[i]/(blk_norm*norm);
			glob_av[i] += gdir;
			glob_av2[i] += gdir*gdir;
			err_gdir= Error(glob_av[i],glob_av2[i],iblk);
			r = (double)(i-igofr)*bin_size + bin_size*0.5;
			Gofr << iblk << " " << i-igofr << " " << r << " " << gdir <<  " " << glob_av[i]/(double)iblk << " " << err_gdir <<  endl;
			if (iblk == nblk) Gave << i-igofr << " " << r*sigma << " " << gdir << " " << glob_av[i]/(double)iblk << " " << err_gdir << endl;
		}

		Gofr.close();
	
	}else if(Realsimulation==1){ // SIMULAZIONE
		ofstream is_Gdir,is_Epot, is_Pres, is_Etot, is_Kin, is_Tem;
		
		is_Epot.open("ELE.epot.dat",ios::app);
		is_Pres.open("ELE.press.dat",ios::app);
		is_Etot.open("ELE.etot.dat",ios::app);
		is_Kin.open("ELE.kin.dat",ios::app);
		is_Tem.open("ELE.tem.dat",ios::app);
		is_Gdir.open("Realoutput.gofr.dat" , ios::app);

		stima_k = blk_av[ik]* epsilon_kb*kb /blk_norm; //Kinetical energy
    glob_av[ik] += stima_k ; // cinetica in unità SI
    glob_av2[ik] += stima_k*stima_k;
    err_k=Error(glob_av[ik],glob_av2[ik],iblk);
    
    stima_et = blk_av[ie] * epsilon_kb*kb/blk_norm; // Energia totale
    glob_av[ie] += stima_et; //totale in unità SI
    glob_av2[ie] += stima_et*stima_et;
    err_et=Error(glob_av[ie],glob_av2[ie],iblk);

		stima_p = blk_av[ip] * epsilon_kb*kb /blk_norm; // Energia potenziale
		glob_av[iv] += stima_p; // potenziale in unità SI
    glob_av2[iv] += stima_p*stima_p;
    err_p=Error(glob_av[iv],glob_av2[iv],iblk);

		stima_t = blk_av[it] * epsilon_kb /blk_norm; // temperatura
		glob_av[it] += stima_t; // temperatura in unità SI
    glob_av2[it] += stima_t*stima_t;
    err_t=Error(glob_av[it],glob_av2[it],iblk);

		stima_pr = blk_av[ip]* (epsilon_kb*kb/pow(sigma,3))/blk_norm; // Pressione
		glob_av[ip] += stima_pr; // pressione in unità SI
    glob_av2[ip] += stima_pr*stima_pr;
    err_pr=Error(glob_av[ip],glob_av2[ip],iblk);

    is_Epot << " "  << iblk << " " << stima_p << " "  << glob_av[iv]/(double)iblk << " "  << err_p << endl;
    is_Tem << " " << iblk <<  " " << stima_t << " "  << glob_av[it]/(double)iblk << " " << err_t << endl;
		is_Pres << " " << iblk <<  " " << stima_pr << " "  << glob_av[ip]/(double)iblk << " " << err_pr<< endl;
		is_Kin << " " << iblk <<  " " << stima_k << " "  << glob_av[ik]/(double)iblk << " " << err_k << endl;
		is_Etot << " " << iblk <<  " " << stima_et << " "  << glob_av[ie]/(double)iblk << " " << err_et << endl;

		is_Epot.close();
		is_Pres.close();
		is_Etot.close();
		is_Kin.close();
		is_Tem.close();

  	for(int i=igofr; i<n_props; i++){
			r = double(i-igofr)*bin_size;
			double vr = (4./3.)*pi*(pow(r+bin_size,3) - pow(r,3));
			double norm = rho*npart*vr;
			gdir = blk_av[i]/(blk_norm*norm);
			glob_av[i] += gdir;
			glob_av2[i] += gdir*gdir;
			err_gdir= Error(glob_av[i],glob_av2[i],iblk);
			r = (double)(i-igofr)*bin_size + bin_size*0.5;
			if(iblk == nblk && Realsimulation == 1) is_Gdir << iblk << " " << i-igofr << " " << r*sigma << " " << glob_av[i]/(double)iblk << " " << err_gdir <<endl;
		}
	
  	is_Gdir.close();
	}

	Reset(iblk);

} // end Averages()



double Error(double sum, double sum2, int iblk){

    return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)iblk);

} // end Error()



void Reset(int iblk){ //Reset block averages

   if(iblk == 1) {
       for(int i=0; i<n_props; ++i) {
           walker[i] = 0;
           walker[i] = 0;
       }
   }
   for(int i=0; i<n_props; ++i){
     blk_av[i] = 0;
   }
   blk_norm = 0;

} // end Reset()



/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
