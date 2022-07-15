#include <iostream>
#include "F.h"
using namespace std;


int main(){
 
  Function f1;
  RungeKutta* r1 = new RungeKutta(1, 2);
  Euler* e1 = new Euler(1, 2);
  
  r1 -> Integrator(1000, 0.01, &f1, &Function::dxdt);
  e1 -> Integrator(1000, 0.01, &f1, &Function::dxdt);

  /////////Integrazione del II ordine////////////
  
  Function2 f2;
  EulerCromer* eu = new EulerCromer(0, 0, 0);
  MidPoint* mp = new MidPoint(0, 0, 0);
  LeapFrog* lp = new LeapFrog(0, 0, 0);
  Verlet* v = new Verlet(0, 0, 0);
  Euler2* e2 = new Euler2(0, 0, 0);
  RungeKutta2* r2 = new RungeKutta2(0, 0, 0);

  eu -> Integrator(1000, 0.01, &f2, &Function2::dVdT);
  mp -> Integrator(1000, 0.01, &f2, &Function2::dVdT);
  lp -> Integrator(1000, 0.01, &f2, &Function2::dVdT);
  v -> Integrator(1000, 0.01, &f2, &Function2::dVdT);
  e2 -> Integrator(1000, 0.01, &f2, &Function2::dVdT);
  r2 -> Integrator(1000, 0.01, &f2, &Function2::dVdT);
  
  delete r1;
  delete e1;
  delete eu;
  delete mp;
  delete lp;
  delete v;
  delete e2;
  delete r2;
  
  return 0;
}




 

