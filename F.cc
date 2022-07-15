#include <iostream>
#include <cmath>
#include <fstream>
#include "F.h"
using namespace std;


//----Function class----//
//I write here the I order differential equation I want to be integrated

double Function::dxdt(double t0, double x0){
  return t0 + 3; 
}


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

void RungeKutta::Set_t0(double t0) {
  t0_ = t0;
}

void RungeKutta::Set_x0(double x0) {
  x0_ = x0;
}

void RungeKutta::print() const{
  cout << "t0 = " << Get_t0() << endl
       << "x0 = " << Get_x0() << endl;
}

void RungeKutta::Integrator(int N, double dt, Function* obj, double(Function::*fp)(double, double)) const{
  
  double K1 = 0, K2 = 0;
  double t0 = Get_t0(), x0 = Get_x0();
  
  ofstream myfile("dataR.txt");
  
  for(int i = 0; i < N; i++){
    
    myfile << t0 << "\t" << x0 << endl;  
    
    K1 = (obj->*fp)(t0,x0);
    K2 = (obj->*fp)(t0 + dt/2., x0 + K1*dt/2.);
    
    x0 = x0 + K2*dt;  
    t0 += dt;       
  }
  myfile.close();
}


//----Euler class----//

Euler::Euler(double t0, double x0){
  t0_ = t0;
  x0_ = x0;
};

Euler::~Euler(){
  cout << "Euler is destroyed!" << endl;
}

double Euler::Get_t0() const{
  return t0_;
}

double Euler::Get_x0() const{
  return x0_;
}

void Euler::Set_t0(double t0) {
  t0_ = t0;
}

void Euler::Set_x0(double x0) {
  x0_ = x0;
}

void Euler::print() const{
  cout << "t0 = " << Get_t0() << endl
       << "x0 = " << Get_x0() << endl;
}


void Euler::Integrator(int N, double dt, Function* obj, double(Function::*fp)(double, double)) const {

  double t0 = Get_t0(), x0 = Get_x0();
  
  ofstream myfile("dataE.txt");
  
  for(int i = 0; i < N; i++){

    myfile << t0 << "\t" << x0  << endl; 
   
    x0 = x0 + (obj->*fp)(t0,x0)*dt;     
    t0 += dt;                           
  }
  myfile.close();
}



////////////////////////////////////////////


//----Function2 class----//
//I write here the II order differential equation I want to be integrated

double Function2::dVdT(double t0, double x0, double v0){
  return t0 + 3;  
}



//----EulerCromer class----//

EulerCromer::EulerCromer(double t0, double x0, double v0){
  t0_ = t0;
  x0_ = x0;
  v0_ = v0;
}

EulerCromer::~EulerCromer(){
  cout << "EulerCromer is destroyed!" << endl;
}

double EulerCromer::Get_t0() const{
  return t0_;
}

double EulerCromer::Get_x0() const{
  return x0_;
}

double EulerCromer::Get_v0() const{
  return v0_;
}

void EulerCromer::Set_t0(double t0){
  t0_ = t0;
}

void EulerCromer::Set_x0(double x0){
  x0_ = x0;
}

void EulerCromer::Set_v0(double v0){
  v0_ = v0;
}

void EulerCromer::print() const{
  cout << "t0 = " << Get_t0() << endl
       << "x0 = " << Get_x0() << endl
       << "v0 = " << Get_v0() << endl;
}

void EulerCromer::Integrator(int N, double dt, Function2* obj, double(Function2::*fp)(double t, double x, double v)) const{
  
  double t0 = Get_t0();  
  double x0 = Get_x0();  
  double v0 = Get_v0();  

  ofstream myfile("dataEC.txt");
  
  for(int i = 0; i < N; i++){
    
    myfile << x0 << "\t" << v0 << "\t" << t0 << endl;  
    
    x0 = x0 + v0*dt;                     
    v0 = v0 + (obj->*fp)(t0, x0, v0)*dt; 
    t0 += dt;                              
  }
  myfile.close();
}


//---- MidPoint class----//

MidPoint::MidPoint(double t0, double x0, double v0){
  t0_ = t0;
  x0_ = x0;
  v0_ = v0;
}

MidPoint::~MidPoint(){
  cout << "MidPoint is destroyed!" << endl;
}

double MidPoint::Get_t0() const{
  return t0_;
}

double MidPoint::Get_x0() const{
  return x0_;
}

double MidPoint::Get_v0() const{
  return v0_;
}

void MidPoint::Set_t0(double t0){
  t0_ = t0;
}

void MidPoint::Set_x0(double x0){
  x0_ = x0;
}

void MidPoint::Set_v0(double v0){
  v0_ = v0;
}

void MidPoint::print() const{
  cout << "t0 = " << Get_t0() << endl
       << "x0 = " << Get_x0() << endl
       << "v0 = " << Get_v0() << endl;
}

void MidPoint::Integrator(int N, double dt, Function2* obj, double(Function2::*fp)(double t, double x , double v)) const{

  double v1 = 0;
  
  double t0 = Get_t0();  
  double x0 = Get_x0();  
  double v0 = Get_v0();  

  ofstream myfile("dataMP.txt");
  
  for(int i = 0; i < N; i++){
    
    myfile << x0 << "\t" << v0 << "\t" << t0 << endl;  
    
    v1 = v0 + (obj->*fp)(t0, x0, v0)*dt; 
    x0 = x0 + (v0 + v1)/2.*dt;            

    v0 = v1;                          
    t0 += dt;                            
  }
  myfile.close();
}


//----LeapFrog class----//

LeapFrog::LeapFrog(double t0, double x0, double v0){
  t0_ = t0;
  x0_ = x0;
  v0_ = v0;
}

LeapFrog::~LeapFrog(){
  cout << "LeapFrog is destroyed!" << endl;
}

double LeapFrog::Get_t0() const{
  return t0_;
}

double LeapFrog::Get_x0() const{
  return x0_;
}

double LeapFrog::Get_v0() const{
  return v0_;
}

void LeapFrog::Set_t0(double t0){
  t0_ = t0;
}

void LeapFrog::Set_x0(double x0){
  x0_ = x0;
}

void LeapFrog::Set_v0(double v0){
  v0_ = v0;
}

void LeapFrog::print() const{
  cout << "t0 = " << Get_t0() << endl
       << "x0 = " << Get_x0() << endl
       << "v0 = " << Get_v0() << endl;
}

void LeapFrog::Integrator(int N, double dt, Function2* obj, double(Function2::*fp)(double t, double x , double v)) const{

  double v12 = 0;
  
  double t0 = Get_t0();  
  double x0 = Get_x0();  
  double v0 = Get_v0();  

  ofstream myfile("dataLF.txt");

  v12 = v0 + (obj->*fp)(t0, x0, v0)*dt/2.; 
  
  for(int i = 0; i < N; i++){
    
    myfile << x0 << "\t" << v0 << "\t" << t0 << endl;  
    
    x0 = x0 + v12*dt;                        
    v0 = v12 + (obj->*fp)(t0, x0, v0)*dt;    
     
    v12 = v0;                                
    t0 += dt;                               
  }
  myfile.close();
}


//----Verlet class----//

Verlet::Verlet(double t0, double x0, double v0){
  t0_ = t0;
  x0_ = x0;
  v0_ = v0;
}

Verlet::~Verlet(){
  cout << "Verlet is destroyed!" << endl;
}

double Verlet::Get_t0() const{
  return t0_;
}

double Verlet::Get_x0() const{
  return x0_;
}

double Verlet::Get_v0() const{
  return v0_;
}

void Verlet::Set_t0(double t0){
  t0_ = t0;
}

void Verlet::Set_x0(double x0){
  x0_ = x0;
}

void Verlet::Set_v0(double v0){
  v0_ = v0;
}

void Verlet::print() const{
  cout << "t0 = " << Get_t0() << endl
       << "x0 = " << Get_x0() << endl
       << "v0 = " << Get_v0() << endl;
}

void Verlet::Integrator(int N, double dt, Function2* obj, double(Function2::*fp)(double t, double x , double v)) const{

  double v12 = 0;
  
  double t0 = Get_t0();  
  double x0 = Get_x0();  
  double v0 = Get_v0();  

  ofstream myfile("dataV.txt");
  
  for(int i = 0; i < N; i++){
    
    myfile << x0 << "\t" << v0 << "\t" << t0 << endl;  
    
    v12 = v0 + (obj->*fp)(t0, x0, v0)*dt*0.5;  
    x0 = x0 + v12*dt;                          
    v0 = v12 + (obj->*fp)(t0, x0, v0)*dt*0.5; 
     
    t0 += dt;                                 
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

void Euler2::Integrator(int N, double dt, Function2* obj, double(Function2::*fp)(double t, double x , double v)) const{
  
  double t0 = Get_t0();  
  double x0 = Get_x0();  
  double v0 = Get_v0(); 
  
  ofstream myfile("dataE2.txt");
  
  for(int i = 0; i < N; i++){
    
    myfile << x0 << "\t" << v0 << "\t" << t0 << endl;  

    x0 = x0 + v0*dt;                     
    v0 = v0 + (obj->*fp)(t0, x0, v0)*dt;  
    t0 += dt;                            
  }
  myfile.close();
}


//----RungeKutta2 class ----//

RungeKutta2::RungeKutta2(double t0, double x0, double v0){
  t0_ = t0;
  x0_ = x0;
  v0_ = v0;
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

void RungeKutta2::Integrator(int N, double dt, Function2* obj, double(Function2::*fp)(double t, double x , double v)) const{
 
  double K1x = 0, K2x = 0, K1v = 0, K2v = 0;

  double t0 = Get_t0();  
  double x0 = Get_x0();  
  double v0 = Get_v0();  

  cout << v0 << endl;
  
  ofstream myfile("dataR2.txt");
  
  for(int i = 0; i < N; i++){

    myfile << x0 << "\t" << v0 << "\t" << t0 << endl;  

    K1x = dt*v0;
    K1v = dt*(obj->*fp)(t0, x0, v0);

    K2x = dt*v0;
    K2v = dt*(obj->*fp)(t0 + dt/2., x0, v0 + K1v*dt/2.);
    
    x0 = x0 + K1x;  
    v0 = v0 + K2v;  
    t0 += dt;          
  }
  myfile.close();
  }
		   


