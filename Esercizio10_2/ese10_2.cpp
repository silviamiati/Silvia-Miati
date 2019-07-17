#include "ese10_2.h"
#include "random.h"
#include <iostream>
#include <fstream>
#include <ostream>
#include <iomanip>
#include <vector>
#include <algorithm>    // std::random_shuffle 
#include <mpi.h>

using namespace std;


int main(int argc, char *argv[]){

MPI_Init(&argc,&argv);

MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
MPI_Comm_rank(MPI_COMM_WORLD, &myid);

if(numprocs > 4){cout<<"Hai scelto troppi processi"<<endl;
return 1;}

	Input(myid); 
	Move(Pos,beta);
	double irecv;

//informazione da passare (&L) al processore 0. Mi restituisce il minimo dei valori
MPI_Reduce(&L, &irecv, 1, MPI_DOUBLE, MPI_MIN, 0, MPI::COMM_WORLD); 

if(myid==0) cout<< "Fitness migliore per i 4 processi: " << irecv << endl;


MPI_Finalize();

return 0;
}


void Input(int myid){

ifstream ReadInput,ReadConfinal;
 
//Read seed for random numbers
   int pa, ps;
	  
   ifstream Primes("Primes");
	 for (int j=0; j<numprocs; j++){
   Primes >> pa;
	 p1[j]= pa;
	 Primes >> ps;
	 p2[j] = ps;
	 }
   Primes.close();

	 ifstream input("seed.in");
   input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
   rnd.SetRandom(seed,p1[myid],p2[myid]);
   input.close();
	 
rnd.SaveSeed();

ReadInput.open("input.dat");

	ReadInput >> tempe;
		beta = 1./tempe;
	ReadInput >> nstep;
	ReadInput >> first;
	ReadInput >> ncity;
	ReadInput >> rho; 
	ReadInput >> fine;
	ReadInput >> restart;

  ReadInput.close();


if(first==1){
ofstream Posi;
double x=0.0,y=0.0,theta=0.0;

Posi.open("Posizioni1.dat", ios::app);

		for(int i=0; i<ncity; i++){
			theta=rnd.Rannyu(0,2.*pi);
			x=cos(theta);
			y=sin(theta);
			Posi << x << " " << y<< endl;
			}
		/*for(int i=0; i<ncity; i++){
			x = rnd.Rannyu();
			y = rnd.Rannyu();
			Posi << x << " " << y << endl;
			}*/
Posi.close();
}

if(restart==1){ 
   cout << "Read initial configuration from file config.final : " << endl << endl;
  ReadConfinal.open("config.final");
  for (int i=0; i<ncity; ++i){
    ReadConfinal >> Pos[i].first >> Pos[i].second;
  }
  ReadConfinal.close();
	}


/////////////////////posizioni////////////////

fstream pos;
double a=0.0;

pos.open("Posizioni1.dat");


for(int i=0; i<30; i++){
	pos >> a;
	xR[i]=a;
	pos >> a;
	yR[i]=a;  
}

for(int i=0; i<ncity; i++){
	Pos.push_back(make_pair(xR[i],yR[i]));
}
pos.close();


}


double Fitness(vector<pair<double,double> > R){    
	L=0;
	for(int j=0; j<ncity; j++){
		if(j==ncity-1) L = L + sqrt(pow(R[ncity-1].first- R[1].first,2) +pow( R[ncity-1].second- R[1].second,2));
		else L = L + sqrt(pow(R[j+1].first - R[j].first,2) + pow(R[j+1].second - R[j].second,2));
	}
	return L;
}
	

void Move(vector <pair<double,double> > &v , double beta){

	while(tempe>0.1){
		for(int istep=1; istep <= nstep; ++istep){
			double p=0, energy_old=0, energy_new=0;
			int rand = v.size()*rnd.Rannyu(); 
			int rand1 = v.size()*rnd.Rannyu(); 
	
			energy_old = Fitness(v);
		
			iter_swap(v.begin() + rand, v.begin() + rand1); 

			energy_new = Fitness(v);
	
			p = exp((1./tempe)*(energy_old-energy_new));
    	if(p < rnd.Rannyu()){ 
   			iter_swap(v.begin() + rand, v.begin() + rand1);
 	  	}
	  }
		tempe *= 0.8;
	}

	
}   
















		



