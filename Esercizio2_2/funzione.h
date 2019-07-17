#include <cmath>
#include <string>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include "random.h"


int nstep=100;  //numero di passi
int B=100;  //numero di blocchi
int N=10000; // numero simulazioni
int L=N/B; //numero di simulazioni per blocco

bool libero;

const int m_N=10000, m_nstep=100, m_L=100;
const double pi=3.14159;

double vx[m_L],vy[m_L],vz[m_L],vtheta[m_L],vphy[m_L],vrho[m_L];
double D[m_nstep],Dlibero[m_nstep];
double global_ar[m_nstep]={0},global_ar2[m_nstep]={0},stima_ar[m_nstep],global_al[m_nstep],global_al2[m_nstep],stima_al[m_nstep];



void Reset(bool);
void Resetv(bool);
void Movexyz(double *,double *,double *, int ,double ,double ,double t);
void Movetheta(double *,double *,double *, int ,double , double , double);
double Distanza(double *, double *, double *, int);
double Error(double, double, int);

