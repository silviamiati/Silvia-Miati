#ifndef __Ga_
#define __Ga_

#include "random.h"
#include <vector>
#include <utility>


int seed[4];
Random rnd;
bool first,percorso;

const int n_city=30, npopulation=300, mapa=2, generation=600; // #città, #individui, #genitori, #generazioni
unsigned int N = npopulation*n_city;
int ncity,rho,dim;
int taglio;

double prob[npopulation];
int parents[mapa];


double xR[n_city],yR[n_city];
std::vector <double> L(npopulation);
std::vector <std::pair<double,int>> F;
std::vector <int> index1,index2;
std::vector <std::pair<double,double>> Pos;
std::vector <std::pair<double,double>> R,Madre,Padre,son_1,son_2,appo_1,appo_2,populi;


const double pi=3.14159;


void Input(void);
void Population();
void Fitness();
void Selector();
double Crossover();
void Mutation(std::vector<std::pair<double,double>> &);
void Reset();
#endif
