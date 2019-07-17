#ifndef __Ga_
#define __Ga_

#include "random.h"
#include <vector>
#include <utility>


int seed[4];
Random rnd;
bool first,fine,restart;

const int n_city=30, n_primes=384;
int ncity,rho,dim,generation;
int nstep,myid,numprocs;
int p1[n_primes],p2[n_primes];

double  L, beta,accepted, attempted, tempe;

double xR[n_city],yR[n_city];
std::vector <std::pair<double,double> > Pos;



const double pi=3.14159;


void Input(int);
double Fitness(std::vector<std::pair<double,double > >);
void Move(std::vector<std::pair<double,double > > &, double);
void ConfFinal(int);
#endif
