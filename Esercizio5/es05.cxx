#include "random.h"
#include "es05.h"
#include <cmath>
#include <string>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <iostream>
#include <fstream>


using namespace std;

int main(int argc, const char **argv) { 
system ("rm R_m_gauss.out");

Input();
Acceptance(); //mi trova il valore della delta in modo tale che l'accettanza sia al 50%

int L=nstep/nblock;


for(int i=0; i<nstep; i++){ 
	Move();
	Measure();
	//ConfXYZ();
	if (i%L==0 and i!=0){
		Block(block);
		block++;
	}	
	
}
ConfFinal(); 

return 0;
}



void Input(void){

	 ifstream ReadInput,ReadConf,ReadConfinal;

   int p1, p2;
   ifstream Primes("Primes");
   Primes >> p1 >> p2 ;
   Primes.close();

   ifstream input("seed.in");
   input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
   rnd.SetRandom(seed,p1,p2);
   input.close();
   rnd.SaveSeed();

	if(restart==1){ 
  	cout << "Read initial configuration from file config.final : " << endl << endl;
 		ReadConfinal.open("config.final");
  	ReadConfinal >> xin >> yin >> zin; 
		xpre = xin;
		ypre = yin;
		zpre = zin;
  	ReadConfinal.close();
	}else{
		cout << "Read initial configuration from file config.0 " << endl << endl;
  	ReadConf.open("config.0");
		 
    ReadConf >> xin >> yin >> zin; 
		xpre=xin;
		ypre=yin;
		zpre=zin;
    ReadConf.close();
	}
  
	 ReadInput.open("input.dat"); //Read input

   ReadInput >> passo;
   ReadInput >> gr;
   ReadInput >> nstep;
   ReadInput >> nblock;
	 ReadInput >> restart;


   ReadInput.close(); 

	 ir = 0; //raggio
   n_props = 1; //Number of observables

} // end Input()
 


double Funz(double x, double y, double z){

	double p;
	if(gr==1){  // sto calcolando la funzione 100
	double r = sqrt(pow(x,2) + pow(y,2) + pow(z,2));
	p = exp(-2.*r);
	}
	else{ // sto calcolando la funzione 210
	double r = sqrt(pow(x,2) + pow(y,2) + pow(z,2));
	double ctheta = fabs(z/r);
	p = exp(-r)*pow(r,2)*pow(ctheta,2);

	}
	return p;

} // end Funz() - calcola il modulo quadro della funzione d'onda -



void Measure(void){

	ave_raggio += sqrt(pow(xpre,2) + pow(ypre,2) + pow(zpre,2));
	blk_norm = blk_norm + 1;

} // end Measure() - calcola il raggio -



void Block(int iblk){

  ofstream rmedio;
	rmedio.open("prova.out",ios::app);

	cout << "accettanza : " << accepted/attempted << endl << endl;
  
  media = ave_raggio/blk_norm;
	glob_av[ir] += media;
	glob_av2[ir] += media*media;
	err_r=Error(glob_av[ir],glob_av2[ir],iblk);
 
  rmedio << glob_av[ir]/(double)iblk << " " << err_r << " " << iblk << endl;

  rmedio.close();

	Reset();

} // end Block() 



void Reset(){ //Reset block averages  

   ave_raggio=0.0;
	 blk_norm=0.0;
   attempted = 0;
   accepted = 0;

} // end Reset()



double Error(double sum, double sum2, int iblk){

    return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)iblk);

} // end Error()



void Move(void){

	double xnew,ynew,znew,x,y,z,alpha;

	if(passo==1){
		x = delta*(rnd.Rannyu() - 0.5) ;
		y = delta*(rnd.Rannyu() - 0.5) ;
		z = delta*(rnd.Rannyu() - 0.5) ;
	
	}else{
		x = rnd.Gauss(0.,delta/2.);
		y = rnd.Gauss(0.,delta/2.);
		z = rnd.Gauss(0.,delta/2.);
	}

	xnew = x+xpre;
	ynew = y+ypre;
	znew = z+zpre;

  alpha = Funz(xnew,ynew,znew)/Funz(xpre,ypre,zpre);

  if(alpha >= rnd.Rannyu()){
    //Update
     xpre = xnew;
     ypre = ynew;
     zpre = znew;
    
		 accepted = accepted + 1.0;
    }
	attempted = attempted + 1.0; 

} // end Move() - Implementa l'algoritmo di M-H -


double Acceptance(){
		 
		delta=0.0;
		double p=0;
		
		while(p<0.4 or p>0.6){	
			accepted=0;
			attempted=0;
			for(int i=0; i<1000; i++){
				Move();	
				
			}	
			p = accepted/attempted;
			delta += 0.5;
		
		}

	return delta;

} // end Acceptance()



void ConfFinal(void){

  ofstream WriteConf, WriteSeed;
	
  cout << "Print final configuration to file config.final " << endl << endl;
  WriteConf.open("config.final");

  WriteConf << xpre << " " << ypre << " " << zpre <<  endl;
  
  WriteConf.close();

}// end ConfFinal()



void ConfXYZ(void){ //Write configuration in .xyz format

  ofstream WriteXYZ;
  WriteXYZ.open("confi.x",ios::app);
    WriteXYZ << xpre << " " << ypre << " " << zpre << endl;
  
  WriteXYZ.close();

} // end ConfXYZ()









