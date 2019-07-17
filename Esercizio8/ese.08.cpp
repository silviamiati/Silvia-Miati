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
#include "ese.08.h"

using namespace std;

int main(){ 

system("rm output.histo.0 output.e.0 ");

  Input(); //Inizialization

  Acceptance();
	
  for(int iblk=1; iblk <= nblk; ++iblk){ //Simulation
   Reset(iblk);   //Reset block averages
    for(int istep=1; istep <= nstep; ++istep){
      Move();
			Measure();
      Accumulate(); //Update block averages
    }	
   Averages(iblk);   //Print results for current block
	 Eblock(iblk); //scommentare se si vuole l'energia per blocco
  }
  ConfFinal(); //Write final configuration
	Histo(); //scommentare se si vuole l'istogramma
	//EUltimo();
	
  return 0;

}



void Input(void){

	 double x;
	 ifstream ReadInput,ReadConf,ReadConfinal;
 
   int p1, p2;
   ifstream Primes("Primes");
   Primes >> p1 >> p2 ;
   Primes.close();

   ifstream input("seed.in");
   input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
   rnd.SetRandom(seed,p1,p2);
   input.close();
  
  ReadInput.open("input.dat");

  ReadInput >> nblk;
  ReadInput >> nstep; 
	ReadInput >> mu;
		cout << "la mu è :" << mu << endl; 
  ReadInput >> sigma2;
		cout << "la sigma quadro è :" << sigma2 << endl;
	ReadInput >> restart;
	ReadInput >> length;
  
  ReadInput.close();

	ie=0;

	n_props = 1;

	nbins = 100;
	igofr = 1;
  n_props = n_props + igofr;

  bin_size = length/(double)nbins;

	if(restart==1){ 
   cout << "Read initial configuration from file config.final : " << endl << endl;
   ReadConfinal.open("config.final");
   for (int i=0; i<npart; ++i) ReadConfinal >> xold;   
   ReadConfinal.close();
 }else{
	 cout << "Read initial configuration from file config.0 " << endl << endl;
   ReadConf.open("config.0");
   ReadConf >> x ; 
   xold = x; 
   ReadConf.close();
 }
   rnd.SaveSeed();

	for (int k=0; k<nbins; k++){
	 histo[k]=0.0;
	 glob_is[k]=0.0;
	 glob_is2[k]=0.0;
	}

}// end Input()



void Move(void){

  double p, wave_old, wave_new;
  double xnew;

    wave_old = Wave2(xold,mu,sigma2);
		
    xnew = xold + delta*(rnd.Rannyu(-0.5, 0.5));
		
    wave_new = Wave2(xnew,mu,sigma2);
	
    p = wave_new/wave_old;

    if(p >= rnd.Rannyu()){
       xold = xnew; 
			 accepted = accepted + 1.0;
    }
    attempted = attempted + 1.0; 

} // end Move(void)



double Wave2(double xx, double mu, double sigma2){

  double w=0.0;
  w = exp(-pow(xx-mu,2.)/sigma2) + exp(-pow(xx+mu,2.)/sigma2) + 2*exp(-(xx*xx + mu*mu)/sigma2); 
	return w;

}// end Wave2() -modulo quadro della funzione d'onda-



double Wave(double xx, double mu, double sigma2){

	double w=0.0;
	w = exp(-pow(xx-mu,2.)/(2.*sigma2)) + exp(-pow(xx+mu,2.)/(2.*sigma2));
	return w;
	
} // end Wave() -funzione d'onda-



void Measure(){

	ofstream Eis;
	int bin;
	//Eis.open("output.eis.0",ios::app);
  double e = 0.0;
	e = (+1./(2.*sigma2))*(exp(-(xold-mu)*(xold-mu)/(2.*sigma2)) + exp(-(xold+mu)*(xold+mu)/(2.*sigma2))) - (1./(2.*sigma2*sigma2))*((xold-mu)*(xold-mu)*exp(-(xold-mu)*(xold-mu)/(2.*sigma2)) + (xold+mu)*(xold+mu)*exp(-(xold+mu)*(xold+mu)/(2.*sigma2))) + (pow(xold,4) -2.5*xold*xold)*Wave(xold,mu,sigma2); //derivate seconde
	e = e/Wave(xold,mu,sigma2);
	walker[ie] = e;
	//Eis << e << endl;
	//Eis.close();  


	if(xold<3. && xold> -3.){
		bin = floor((xold + 3)/6.*nbins); //sto traslando per fare andare l'istogramma da 0 a 6
		histo[int(bin)] += 1;
	}
		
} // end Measure() - calcola energia e riempe l'istogramma per |(\phy)^2|



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

}// end Reset()



void Accumulate(void){ //Update block averages

   for(int i=0; i<n_props; ++i){
     blk_av[i] = blk_av[i] + walker[i];
   }		 
   blk_norm = blk_norm + 1.0;

} // end Accumulate() 



void Averages(int iblk){ //Print results for current block  
 
    cout << "Block number " << iblk << endl;
    cout << "Acceptance rate " << accepted/attempted << endl << endl;
    
    stima_e = blk_av[ie]/blk_norm; //Potential energy
    glob_av[ie] += stima_e;
    glob_av2[ie] += stima_e*stima_e;
    err_e=Error(glob_av[ie],glob_av2[ie],iblk);
    
		for(int i=0; i<nbins; i++){
			histo[i] /= nstep;
			glob_is[i] += histo[i];
			glob_is2[i] += histo[i]*histo[i];
			histo[i]=0;		
		}
	
} // end Averages() -stampa valori dell'energia -
 


void ConfFinal(void){

  ofstream WriteConf, WriteSeed;
	
  cout << "Print final configuration to file config.final " << endl << endl;
  WriteConf.open("config.final");

  WriteConf << xold << endl;
  
  WriteConf.close();

  rnd.SaveSeed();

} // end ConfFinal()



void ConfXYZ(void){ //Write configuration in .xyz format

  ofstream WriteXYZ;
  WriteXYZ.open("confi.x",ios::app);
    WriteXYZ << xold <<endl;
  
  WriteXYZ.close();

} // end ConfXYZ()



double Error(double sum, double sum2, int iblk){

    return sqrt(abs((sum2/(double)iblk - pow(sum/(double)iblk,2)))/(double)iblk);

}// end Error()



double Acceptance(){
		 
		delta=0.1;
		double p=0;
		
		while(p<0.4 or p>0.6){	
			accepted=0;
			attempted=0;
			for(int i=0; i<1000; i++){
				Move();	
			}	
			p = (double)accepted/attempted;
			delta += 0.5;
		}
	return delta;

} //end Acceptance() -calcola il passo tale per cui l'accettanza è circa il 50%-
			


void Histo(){	
	
	double norma=0.0;
	ofstream Histo;
	Histo.open("output.histo.0", ios::app);

	for(int i=0; i< nbins; i++) norma += glob_is[i]/nblk;

	norma *= bin_size;

	for(int i=0; i< nbins; i++){
	 Histo << bin_size*i - 3. << " " << glob_is[i]/(double)nblk/norma  << " " << Error(glob_is[i],glob_is2[i+1],nblk)/norma <<endl;	
  }

	Histo.close();

} // end Histo() 



void EUltimo(){
	
	ofstream Eultimo;
	Eultimo.open("output.e.ultimo",ios::app);  // mi faccio stampare l'ultimo valore dell'energia associato ad una certa mu e sigma
	Eultimo << glob_av[ie]/(double)nblk << " " << mu << " " << sigma2 << endl; 
	Eultimo.close();		
               
} // end EUltimo()  



void Eblock(int iblk){ 
	
	ofstream E;   
	E.open("output.e.0",ios::app);  
	E << iblk << " " << stima_e << " " << glob_av[ie]/(double)iblk << " " << err_e << endl;
	E.close(); 

} // end Eblock() 

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
