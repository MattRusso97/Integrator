#include <iostream>
#include <cmath>
#include <fstream>
#include "E3.h"
using namespace std;


//----RungeKutta class----//

RungeKutta::RungeKutta(double t0, double x0){
  t0_ = t0;
  x0_ = x0;
}

RungeKutta::~RungeKutta(){
  cout << "RungeKutta is destroyed!" << endl;
}

double RungeKutta::Get_t0() const{
  return t0_;
}

double RungeKutta::Get_x0() const{
  return x0_;
}

void RungeKutta::Set_t0(double t0){
  t0_ = t0;
}

void RungeKutta::Set_x0(double x0){
  x0_ = x0;
}

void RungeKutta::print() const{
  cout << "t0 = " << Get_t0() << endl
       << "x0 = " << Get_x0() << endl;
}

void RungeKutta::RungeK(double (*f) (double t, double x), double N, double dt){
 
  double K1 = 0, K2 =0;

  double t0 = Get_t0();  
  double x0 = Get_x0(); 
  
  ofstream myfile("dataR.txt");
  
  for(int i = 0; i < N; i++){

    myfile << x0 << "\t" << t0 << endl;  

    K1 = f(t0, x0);
    K2 = f(t0 + dt/2., x0 + K1*dt/2.);

    x0 = x0 + K2*dt;  
    t0 += dt;         
  }
  myfile.close();
}


//----Euler class----//

Euler::Euler(double t0, double x0){
  t0_ = t0;
  x0_ = x0;
}

Euler::~Euler(){
  cout << "Euler is destroyed!" << endl;
}

double Euler::Get_t0() const{
  return t0_;
}

double Euler::Get_x0() const{
  return x0_;
}

void Euler::Set_t0(double t0){
  t0_ = t0;
}

void Euler::Set_x0(double x0){
  x0_ = x0;
}

void Euler::print() const{
  cout << "t0 = " << Get_t0() << endl
       << "x0 = " << Get_x0() << endl;
}

void Euler::EulerB(double (*f) (double t, double x), double N, double dt){
  
  double t0 = Get_t0();  
  double x0 = Get_x0(); 

  ofstream myfile("dataE.txt");
  
  for(int i = 0; i < N; i++){
    
    myfile << x0 << "\t" << t0 << endl; 

    x0 = x0 + f(t0, x0)*dt;  
    t0 += dt;                
  }
  myfile.close();
}


//----Euler-Cromer class----//

EulerCromer::EulerCromer(double T0, double X0, double V0){
  T0_ = T0;
  X0_ = X0;
  V0_ = V0;
}

EulerCromer::~EulerCromer(){
  cout << "EulerCromer is destroyed!" << endl;
}

double EulerCromer::Get_T0() const{
  return T0_;
}

double EulerCromer::Get_X0() const{
  return X0_;
}

double EulerCromer::Get_V0() const{
  return V0_;
}

void EulerCromer::Set_T0(double T0){
  T0_ = T0;
}

void EulerCromer::Set_X0(double X0){
  X0_ = X0;
}

void EulerCromer::Set_V0(double V0){
  V0_ = V0;
}

void EulerCromer::print() const{
  cout << "t0 = " << Get_T0() << endl
       << "x0 = " << Get_X0() << endl
       << "v0 = " << Get_V0() << endl;
}

void EulerCromer::EulerC(double (*f) (double T, double X, double V), double N, double dt){
  
  double T0 = Get_T0();  
  double X0 = Get_X0();  
  double V0 = Get_V0(); 

  ofstream myfile("dataC.txt");
  
  for(int i = 0; i < N; i++){
    
    myfile << X0 << "\t" << V0 << "\t" << T0 << endl;  
    
    X0 = X0 + V0*dt;            
    V0 = V0 + f(T0, X0, V0)*dt;  
    T0 += dt;                    
  }
  myfile.close();
}


//----MidPoint class----//

MidPoint::MidPoint(double T0, double X0, double V0){
  T0_ = T0;
  X0_ = X0;
  V0_ = V0;
}

MidPoint::~MidPoint(){
  cout << "MidPoint is destroyed!" << endl;
}

double MidPoint::Get_T0() const{
  return T0_;
}

double MidPoint::Get_X0() const{
  return X0_;
}

double MidPoint::Get_V0() const{
  return V0_;
}

void MidPoint::Set_T0(double T0){
  T0_ = T0;
}

void MidPoint::Set_X0(double X0){
  X0_ = X0;
}

void MidPoint::Set_V0(double V0){
  V0_ = V0;
}

void MidPoint::print() const{
  cout << "t0 = " << Get_T0() << endl
       << "x0 = " << Get_X0() << endl
       << "v0 = " << Get_V0() << endl;
}

void MidPoint::MidP(double (*f) (double T, double X, double V), double N, double dt){

  double V1 = 0;
  
  double T0 = Get_T0();  
  double X0 = Get_X0(); 
  double V0 = Get_V0();  

  ofstream myfile("dataMP.txt");
  
  for(int i = 0; i < N; i++){
    
    myfile << X0 << "\t" << V0 << "\t" << T0 << endl;  
    
    V1 = V0 + f(T0, X0, V0)*dt;  
    X0 = X0 + (V0 + V1)/2.*dt;  

    V0 = V1;                     
    T0 += dt;                   
  }
  myfile.close();
}


//----LeapFrog class----//

LeapFrog::LeapFrog(double T0, double X0, double V0){
  T0_ = T0;
  X0_ = X0;
  V0_ = V0;
}

LeapFrog::~LeapFrog(){
  cout << "LeapFrog is destroyed!" << endl;
}

double LeapFrog::Get_T0() const{
  return T0_;
}

double LeapFrog::Get_X0() const{
  return X0_;
}

double LeapFrog::Get_V0() const{
  return V0_;
}

void LeapFrog::Set_T0(double T0){
  T0_ = T0;
}

void LeapFrog::Set_X0(double X0){
  X0_ = X0;
}

void LeapFrog::Set_V0(double V0){
  V0_ = V0;
}

void LeapFrog::print() const{
  cout << "t0 = " << Get_T0() << endl
       << "x0 = " << Get_X0() << endl
       << "v0 = " << Get_V0() << endl;
}

void LeapFrog::LeapF(double (*f) (double T, double X, double V), double N, double dt){

  double V12 = 0;
  
  double T0 = Get_T0(); 
  double X0 = Get_X0(); 
  double V0 = Get_V0(); 

  ofstream myfile("dataLF.txt");

  V12 = V0 + f(T0, X0, V0)*dt/2.;  
  
  for(int i = 0; i < N; i++){
    
    myfile << X0 << "\t" << V0 << "\t" << T0 << endl;  
    
    X0 = X0 + V12*dt;                
    V0 = V12 + f(T0, X0, V0)*dt;     
     
    V12 = V0;                        
    T0 += dt;                       
  }
  myfile.close();
}


//----Verlet class----//

Verlet::Verlet(double T0, double X0, double V0){
  T0_ = T0;
  X0_ = X0;
  V0_ = V0;
}

Verlet::~Verlet(){
  cout << "Verlet is destroyed!" << endl;
}

double Verlet::Get_T0() const{
  return T0_;
}

double Verlet::Get_X0() const{
  return X0_;
}

double Verlet::Get_V0() const{
  return V0_;
}

void Verlet::Set_T0(double T0){
  T0_ = T0;
}

void Verlet::Set_X0(double X0){
  X0_ = X0;
}

void Verlet::Set_V0(double V0){
  V0_ = V0;
}

void Verlet::print() const{
  cout << "t0 = " << Get_T0() << endl
       << "x0 = " << Get_X0() << endl
       << "v0 = " << Get_V0() << endl;
}

void Verlet::Verl(double (*f) (double T, double X, double V), double N, double dt){

  double V12 = 0;
  
  double T0 = Get_T0();  
  double X0 = Get_X0();  
  double V0 = Get_V0();  

  ofstream myfile("dataV.txt");
  
  for(int i = 0; i < N; i++){
    
    myfile << X0 << "\t" << V0 << "\t" << T0 << endl;  
    
    V12 = V0 + f(T0, X0, V0)*dt*0.5;  
    X0 = X0 + V12*dt;                 
    V0 = V12 + f(T0, X0, V0)*dt*0.5;  
     
    T0 += dt;                        
  }
  myfile.close();
}


//----Euler2 class----//

Euler2::Euler2(double t0, double x0, double v0){
  t0_ = t0;
  x0_ = x0;
  v0_ = v0;
}

Euler2::~Euler2(){
  cout << "Euler2 is destroyed!" << endl;
}

double Euler2::Get_t0() const{
  return t0_;
}

double Euler2::Get_x0() const{
  return x0_;
}

double Euler2::Get_v0() const{
  return v0_;
}

void Euler2::Set_t0(double t0){
  t0_ = t0;
}

void Euler2::Set_x0(double x0){
  x0_ = x0;
}

void Euler2::Set_v0(double v0){
  v0_ = v0;
}

void Euler2::print() const{
  cout << "t0 = " << Get_t0() << endl
       << "x0 = " << Get_x0() << endl
       << "v0 = " << Get_v0() << endl;
}

void Euler2::Euler2B(double (*f) (double t, double x, double v), double N, double dt){
  
  double t0 = Get_t0(); 
  double x0 = Get_x0(); 
  double v0 = Get_v0(); 
  
  ofstream myfile("dataE2.txt");
  
  for(int i = 0; i < N; i++){
    
    myfile << x0 << "\t" << v0 << "\t" << t0 << endl; 

    x0 = x0 + v0*dt;            
    v0 = v0 + f(t0, x0, v0)*dt;  
    t0 += dt;                    
  }
  myfile.close();
}


//----RungeKutta2 class----//

RungeKutta2::RungeKutta2(double t0, double x0, double v0){
  t0_ = t0;
  x0_ = x0;
}

RungeKutta2::~RungeKutta2(){
  cout << "RungeKutta2 is destroyed!" << endl;
}

double RungeKutta2::Get_t0() const{
  return t0_;
}

double RungeKutta2::Get_x0() const{
  return x0_;
}

double RungeKutta2::Get_v0() const{
  return v0_;
}

void RungeKutta2::Set_t0(double t0){
  t0_ = t0;
}

void RungeKutta2::Set_x0(double x0){
  x0_ = x0;
}

void RungeKutta2::Set_v0(double v0){
  v0_ = v0;
}

void RungeKutta2::print() const{
  cout << "t0 = " << Get_t0() << endl
       << "x0 = " << Get_x0() << endl
       << "v0 = " << Get_x0() << endl;
}

void RungeKutta2::RungeK2(double (*f) (double t, double x, double v), double N, double dt){
 
  double K1x = 0, K2x =0, K1v = 0, K2v = 0;

  double t0 = Get_t0();  
  double x0 = Get_x0(); 
  double v0 = Get_v0();  
  
  ofstream myfile("dataR2.txt");
  
  for(int i = 0; i < N; i++){

    myfile << x0 << "\t" << v0 << "\t" << t0 << endl;  

    K1x = f(t0, x0, v0);
    K1v = v0;
    
    K2x = f(t0 + dt/2., x0 + K1x*dt/2., v0);
    K2v = v0 + K1v*dt/2.;

    x0 = x0 + K2x*dt; 
    v0 = v0 + K2v*dt;  
    t0 += dt;         
  }
  myfile.close();
}



