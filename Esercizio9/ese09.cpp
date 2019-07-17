#include "ese09.h"
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
#include <valarray>
#include <cstddef> 
#include <unistd.h>
#include <bits/stdc++.h> 


using namespace std;


int main(int argc, const char **argv){

system("rm popolazione.dat Fitness.dat fine_path.dat EvoFitne.dat");

	Input(); 
	ofstream Route, Best_path;
	Route.open("popolazione.dat", ios::app);
	for(int j=0; j<npopulation; j++){
		Population();  // Vettore R (pair di vector (x,y))
	}	
// scrivo su file popolazione.dat la mia popolazione iniziale	
	for(unsigned int i=0; i<R.size(); i++) Route<< R[i].first << " " << R[i].second<< endl; 
	Route.close();
// mi calcolo la fitness iniziale
	Fitness();	

	for(int i=0; i<generation; i++){
		if(i != 0)  Fitness();	
		while(populi.size() < N){ // finchè non ho raggiunto la stessa dimensione della popolazione
			Selector(); 
			Crossover();	
		}
	R.assign(populi.begin(),populi.end());

//Riassegno a populi (Vettore di pair di appoggio contenente la popolazione) la dimensione 0
	populi.resize(0);

// Quando raggiungo la fine della mia simulazione scrivo la popolazione finale in fine.dat!!+
	if(i == generation - 1 ){
		int pos = (npopulation-2)*ncity; // sto stampando il penultimo percorso
		Best_path.open("fine_path.dat");
		for(unsigned int i=0; i<ncity; i++) Best_path << R[i+pos].first << " " << R[i+pos].second << endl;
		Best_path.close();

	}
	Reset();
}


return 0;
}



void Input(void){

ifstream ReadInput;
 
//Read seed for random numbers
	int p1, p2;
	ifstream Primes("Primes");
  Primes >> p1 >> p2 ;
  Primes.close();

  ifstream input("seed.in");
  input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
  rnd.SetRandom(seed,p1,p2);
  input.close();

  ReadInput.open("input.dat");

	ReadInput >> first;
	ReadInput >> ncity;
	ReadInput >> rho; 
	ReadInput >> percorso;

  ReadInput.close();

	if(first==1){
		ofstream Posi;
	  double x=0.0,y=0.0;
	  if(percorso==1){ // Cerchio
	    Posi.open("Posizioni_cerchio.dat", ios::app);
		  double theta=0.0;
		    for(int i=0; i<ncity; i++){
			    theta=rnd.Rannyu(0,2.*pi);
			    x=cos(theta);
			    y=sin(theta);
			    Posi << x << " " << y<< endl;
			  }
		}else{ // Quadrato
			 Posi.open("Posizioni_quadrato.dat", ios::app);
		   for(int i=0; i<ncity; i++){
			   x = rnd.Rannyu();
			   y = rnd.Rannyu();
			   Posi << x << " " << y << endl;
			 }
		}
		Posi.close();
	}

	rnd.SaveSeed();


	fstream pos;
	double a=0.0;
	if(percorso==1)pos.open("Posizioni_cerchio.dat");
	else pos.open("Posizioni_quadrato.dat");

	for(int i=0; i<ncity; i++){
		pos >> a;
		xR[i]=a;
		pos >> a;
		yR[i]=a;  
	}

	for(int i=0; i<ncity; i++) Pos.push_back(make_pair(xR[i],yR[i]));

}//end Input()



void Population(){ // data la posizione delle 30 città genera una popolazione shuffolando (casualmente)

	shuffle(Pos.begin(), Pos.end(), default_random_engine(0));
  for(vector<pair<double,double>>::const_iterator it=Pos.begin(); it!=Pos.end(); ++it)
	R.push_back(make_pair(it->first, it->second));
	
} // end Population() - Crea popolazione-



void Fitness(){  //Calcola la fitness e scrive i valori in appoggio su Fitness.dat (F= 1/L)
  
	int i=0,l=0,p=0;
	double fit_medio=0; 

	ofstream Fitne, Fitnend, Fitmedio;

	Fitne.open("Fitness.dat",ios::app); // File in cui ci sono tutte le fitness
	Fitnend.open("EvoFitne.dat", ios::app); // File in cui viene stampata fitness migliore
	Fitmedio.open("Fit_medio.dat", ios::app); // File in cui c'è la media su metà popolazione della fitness
	while(i<npopulation*ncity){
		for(int j=i; j<ncity+i; j++){
			if(j==ncity-1) L[l] = L[l] + sqrt(pow(R[ncity-1].first- R[1].first,2) + pow(R[ncity-1].second- R[1].second,2));
			else { L[l] = L[l] + sqrt(pow(R[j+1].first - R[j].first,2) + pow(R[j+1].second - R[j].second,2));
			}
			p++;
		}
	l++;
	i = i + 30;
	}
	for(int i=0; i<npopulation; i++){
		Fitne << L[i] << endl; 

			if(i == npopulation - 2) Fitnend << L[i] << endl;
				
			if(i>npopulation/2) fit_medio += L[i];

		L[i] = 1./L[i];
		F.push_back(make_pair(L[i],i)); // vector di pair con in prima colonna le fitness e in seconda la posizione
	}
	
	Fitmedio << fit_medio/(double)npopulation * 2. << endl;

	Fitne.close();
	Fitnend.close(); 
	Fitmedio.close();	

} // end Fitness()



void Selector(){ // scelta genitori tramite il metodo della roulette truccata 

	double somma=0.0;
	double sum[npopulation];

	for(int i=0; i<npopulation; i++) sum[i]=0.0;
	
	for(int i=0; i<mapa; i++) parents[i]=0;

	for(int l=0; l<npopulation; l++){  //mi ricavo la normalizzazione per le probabilità
		somma = somma + L[l];
	}

	sort(F.begin(), F.end()); // ordino le Fitness secondo le fitness
	
	for(int p=0; p<npopulation-1; p++){ //probabilità normalizzata a 1
		if(p==0) sum[0] = F[0].first/somma;    
		sum[p+1] =  F[p+1].first/somma + sum[p];
	}
	
	for(int i=0; i<mapa; i++){      //scelta dei genitori
		double p=rnd.Rannyu();
			for(int k=0; k<npopulation; k++){
				if (p<sum[k]){
					parents[i] = F[k].second;  //tiene memoria della posizione dei genitori
					//cout << " E' stato scelto il genitore: " << F[k].second << " " << F[k].first << endl;
				}
			}	
	}

// Carico i genitori! Alla fine del crossover verranno posti a 0 (anche lunghezza)	
	Padre.insert(Padre.begin(),R.begin() + parents[0]*ncity, R.begin() + (parents[0]+1)*ncity); 
	Madre.insert(Madre.begin(),R.begin() + parents[1]*ncity, R.begin() + (parents[1]+1)*ncity); 
	
}// end Selector()



double Crossover(){ // Generazione figli	

	int taglio = rnd.Rannyu(0,ncity-1); // taglio casuale nei genitori
	if(rnd.Rannyu() > 0.30){

		for(int j=0; j<taglio; j++) son_1.push_back(make_pair(Madre[j].first,Madre[j].second)); // carico in son_1 Madre

// carico in appo_1 la parte d Madre da cercare in Padre	  
		for(int k=taglio; k<ncity; k++) appo_1.push_back(make_pair(Madre[k].first,Madre[k].second)); 
			 	
	  for(int j=0; j<taglio; j++) son_2.push_back(make_pair(Padre[j].first,Padre[j].second)); // carico in son_2 Padre 

// carico in appo_2 la parte d Padre da cercare in Madre		
		for(int k=taglio; k< ncity; k++) appo_2.push_back(make_pair(Padre[k].first,Padre[k].second));	
		
		for(int l=0; l<ncity-taglio; l++){
			auto it = find(Madre.begin(), Madre.end(), appo_2[l]); //sto cercando la parte "scartata" di son_2 in tutta la Madre
			int result1 = distance(Madre.begin(), it);	
			index1.push_back(result1); // carico le posizioni delle città mancanti in index
		}		
	
		sort(index1.begin(), index1.end()); // ordino le posizioni 

// completo son_2 con le città che gli mancavano in ordine di come sono messe in Madre
		for(int j=0; j<ncity-taglio; j++) son_2.push_back(make_pair(Madre[index1[j]].first,Madre[index1[j]].second)); 
				
// STESSA COSA PER son_1
	  for(int l=0; l<ncity-taglio; l++){
			auto it =find(Padre.begin(), Padre.end(), appo_1[l]); //sto cercando la parte scartata di padre in tutta la madre
			int result2= distance(Padre.begin(), it);
			index2.push_back(result2);
		}

		sort(index2.begin(), index2.end());

	  for(int j=0; j<ncity-taglio; j++) son_1.push_back(make_pair(Padre[index2[j]].first,Padre[index2[j]].second));		

// Applico le mutazioni!!			
		Mutation(son_1);
		Mutation(son_2);

// Aggiorno la mia popolazione!!!!!!!	
		for(int j=0; j<ncity; j++) populi.push_back(make_pair(son_2[j].first,son_2[j].second));
		for(int j=0; j<ncity; j++) populi.push_back(make_pair(son_1[j].first,son_1[j].second));
	
		son_1.resize(0),son_2.resize(0),Madre.resize(0),Padre.resize(0),appo_1.resize(0),appo_2.resize(0);

		index1.resize(0),index2.resize(0);
}	

	return populi.size();

} //end Crossover() - Continuo finchè non ho aggiornato tutta la mia vecchia popolazione-
 


void Mutation(vector <pair<double,double> > &v){

	double x=rnd.Rannyu(), y=rnd.Rannyu(), z=rnd.Rannyu(), w=rnd.Rannyu(), f=rnd.Rannyu();
	if(x<=0.1){
		int rand = v.size()*rnd.Rannyu();
		int rand1 = v.size()*rnd.Rannyu(); 
		iter_swap(v.begin() + rand, v.begin() + rand1);
	} else if(y<=0.1){
		int rot=6;
		int rand = rnd.Rannyu(0,ncity/2.);
		rotate(v.begin() + rand, v.begin()+rot + rand, v.end());
	}else if(z<=0.1){
		int rot=2 ;
		int rand =rnd.Rannyu(ncity/2, ncity-rot);
		rotate(v.begin() + rand, v.begin()+rot +rand, v.end());
	} else if(w<=0.1){
		int rand = rnd.Rannyu(0,ncity/2.);
		int rand1 = rnd.Rannyu(ncity/2., ncity); 
	  reverse(v.begin()+rand,v.begin()+rand1);	
	} 
		
} // end Mutation()
		


void Reset(){

	for(int k=0; k<npopulation; k++) L[k] = 0;
	F.resize(0);

}
























		



