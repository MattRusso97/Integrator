#include <iostream>
using namespace std;


///////INTEGRAZIONE EQ. DIFF. I ORDINE///////////

//----Classe per l'equazione differenziale da integrare----//

class Function{
 public:
  double dxdt(double t0, double x0);
};


//----Classe di tipo Abstract per il metodo di integrazione scelto----//

class Method{
 public: 
  virtual void Integrator(int N, double dt, Function* obj, double(Function::*fp)(double, double)) const = 0;
};


//----Classe per l'integrazione con RungeKutta----//

class RungeKutta: public Method{
 public:
  RungeKutta(double t0, double x0);   //Constructor
  virtual ~RungeKutta();              //Destructor

  virtual void Integrator(int N, double dt, Function* obj, double(Function::*fp)(double, double)) const;

  double Get_t0() const;    //Getter: restituisce il tempo iniziale t0
  double Get_x0() const;    //Getter: restituisce la posizione iniziale x0
  
  void Set_t0(double t0);   //Setter: permette di settare il tempo iniziale
  void Set_x0(double x0);   //Setter: permette di settare la posizione iniziale

  void print() const;       //Stampa le condizioni iniziali

 private:
  double t0_;               //Tempo iniziale (in [s])
  double x0_;               //Posizione iniziale (in [m])
};


//----Implementazione della classe Euler----//

class Euler: public Method{
 public:

  Euler(double t0, double x0);    //Constructor
  virtual ~Euler();               //Destructor
  
  virtual void Integrator(int N, double dt, Function* obj, double(Function::*fp)(double, double)) const;

  double Get_t0() const;   //Getter: restituisce il tempo iniziale t0
  double Get_x0() const;   //Getter: restituisce la posizione iniziale x0
  
  void Set_t0(double t0);  //Setter: permette di settare il tempo iniziale
  void Set_x0(double x0);  //Setter: permette di settare la posizione iniziale

  void print() const;      //Stampa le condizioni iniziali
  
 private:
 double t0_;               //Tempo iniziale (in [s])
 double x0_;               //Posizione iniziale (in [m])
};


//////////INTEGRAZIONE EQ. DIFF. DEL II ORDINE//////////

class Function2{
 public:
  double dVdT(double t0, double x0, double v0);
};


//----Classe di tipo Abstract per il metodo di integrazione scelto----//

class Method2{
 public: 
  virtual void Integrator(int N, double dt, Function2* obj, double(Function2::*fp)(double t, double x, double v)) const = 0;
};


//----Metodo di integrazione Eulero-Cromer----//

class EulerCromer: public Method2{
 public:
  EulerCromer(double t0, double x0, double v0);
  virtual ~EulerCromer();

  virtual void Integrator(int N, double dt, Function2* obj, double(Function2::*fp)(double t, double x , double v)) const;

  double Get_t0() const;
  double Get_x0() const;
  double Get_v0() const;

  void Set_t0(double t0);
  void Set_x0(double x0);
  void Set_v0(double v0);

  void print() const;

private:
  double t0_;
  double x0_;
  double v0_;
};


class MidPoint{
 public:
  MidPoint(double t0, double x0, double v0);
  virtual ~MidPoint();

   virtual void Integrator(int N, double dt, Function2* obj, double(Function2::*fp)(double t, double x , double v)) const;
  
  double Get_t0() const;
  double Get_x0() const;
  double Get_v0() const;

  void Set_t0(double t0);
  void Set_x0(double x0);
  void Set_v0(double v0);

  void print() const;

private:
  double t0_;
  double x0_;
  double v0_;
};


class LeapFrog{
 public:
  LeapFrog(double t0, double x0, double v0);
  virtual ~LeapFrog();

  void Integrator(int N, double dt, Function2* obj, double(Function2::*fp)(double t, double x , double v)) const;

  double Get_t0() const;
  double Get_x0() const;
  double Get_v0() const;

  void Set_t0(double t0);
  void Set_x0(double x0);
  void Set_v0(double v0);

  void print() const;

private:
  double t0_;
  double x0_;
  double v0_;
};


class Verlet{
 public:
  Verlet(double t0, double x0, double v0);
  virtual ~Verlet();

  void Integrator(int N, double dt, Function2* obj, double(Function2::*fp)(double t, double x , double v)) const;

  double Get_t0() const;
  double Get_x0() const;
  double Get_v0() const;

  void Set_t0(double t0);
  void Set_x0(double x0);
  void Set_v0(double v0);

  void print() const;

private:
  double t0_;
  double x0_;
  double v0_;
};


class Euler2{
 public:
  Euler2(double t0, double x0, double v0);
  virtual ~Euler2();

  void Integrator(int N, double dt, Function2* obj, double(Function2::*fp)(double t, double x , double v)) const;

  double Get_t0() const;
  double Get_x0() const;
  double Get_v0() const;

  void Set_t0(double t0);
  void Set_x0(double x0);
  void Set_v0(double v0);

  void print() const;

private:
  double t0_;
  double x0_;
  double v0_;
};


class RungeKutta2{
 public:
  RungeKutta2(double t0, double x0, double v0);
  ~RungeKutta2();

   void Integrator(int N, double dt, Function2* obj, double(Function2::*fp)(double t, double x , double v)) const;

  double Get_t0() const;
  double Get_x0() const;
  double Get_v0() const;

  void Set_t0(double t0);
  void Set_x0(double x0);
  void Set_v0(double v0);

  void print() const;

private:
  double t0_;
  double x0_;
  double v0_;
};


///Stampareeeee/////////

 

