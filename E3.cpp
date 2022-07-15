#include <iostream>
#include "E3.h"
using namespace std;

double dxdt(double, double); //This class can integrate I order differential equations
double dVdT(double, double, double); //This class can integrate II order differential equations


int main(){

  //----Integration of I order differential equations----//
  
  RungeKutta a(2, 4);
  a.RungeK(&dxdt, 1000, 0.01);

  Euler b(2, 4);
  b.EulerB(&dxdt, 1000, 0.01);

   //----Integration of II order differential equations----//

  EulerCromer c(2, 4, 0);
  c.EulerC(&dVdT, 1000, 0.01);

  MidPoint d(2, 4, 0);
  d.MidP(&dVdT, 1000, 0.01);

  LeapFrog e(2, 4, 0);
  e.LeapF(&dVdT, 1000, 0.01);

  Verlet f(2, 4, 0);
  f.Verl(&dVdT, 1000, 0.01);

  //----Integration of II order differential equations using Euler and RungeKutta----//

  Euler2 g(2, 4, 0);
  g.Euler2B(&dVdT, 1000, 0.01);

  RungeKutta2 h(2, 4, 0);
  h.RungeK2(&dVdT, 1000, 0.01);

  Exponential i(3, 4);
  double o = i.Expo(0.1);
  cout << o << endl;

  return 0;
}


//Here we write the I order differential equation we want to integrate
//in the form dx/dt = ax + b
double dxdt(double t, double x){
  return x*t + 3;
}

//Here we write the II order differential equation we want to integrate
//in the form dV/dT = aX + bV + c
double dVdT(double T, double X, double V){
  return - 0.5*X;
}


