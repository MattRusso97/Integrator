#include <iostream>

class RungeKutta{
 public:
  RungeKutta(double t0, double x0);
  ~RungeKutta();

  void RungeK(double (*f) (double t, double x), double N, double dt);

  double Get_t0() const;
  double Get_x0() const;

  void Set_t0(double t0);
  void Set_x0(double x0);

  void print() const;

private:
  double t0_;
  double x0_;
};

class Euler{
 public:
  Euler(double t0, double x0);
  ~Euler();

  void EulerB(double (*f) (double t, double x), double N, double dt);

  double Get_t0() const;
  double Get_x0() const;

  void Set_t0(double t0);
  void Set_x0(double x0);

  void print() const;

private:
  double t0_;
  double x0_;
};


class EulerCromer{
 public:
  EulerCromer(double T0, double X0, double V0);
  ~EulerCromer();

  void EulerC(double (*f) (double T, double X, double V), double N, double dt);

  double Get_T0() const;
  double Get_X0() const;
  double Get_V0() const;

  void Set_T0(double T0);
  void Set_X0(double X0);
  void Set_V0(double V0);

  void print() const;

private:
  double T0_;
  double X0_;
  double V0_;
};


class MidPoint{
 public:
  MidPoint(double T0, double X0, double V0);
  ~MidPoint();

  void MidP(double (*f) (double T, double X, double V), double N, double dt);

  double Get_T0() const;
  double Get_X0() const;
  double Get_V0() const;

  void Set_T0(double T0);
  void Set_X0(double X0);
  void Set_V0(double V0);

  void print() const;

private:
  double T0_;
  double X0_;
  double V0_;
};


class LeapFrog{
 public:
  LeapFrog(double T0, double X0, double V0);
  ~LeapFrog();

  void LeapF(double (*f) (double T, double X, double V), double N, double dt);

  double Get_T0() const;
  double Get_X0() const;
  double Get_V0() const;

  void Set_T0(double T0);
  void Set_X0(double X0);
  void Set_V0(double V0);

  void print() const;

private:
  double T0_;
  double X0_;
  double V0_;
};


class Verlet{
 public:
  Verlet(double T0, double X0, double V0);
  ~Verlet();

  void Verl(double (*f) (double T, double X, double V), double N, double dt);

  double Get_T0() const;
  double Get_X0() const;
  double Get_V0() const;

  void Set_T0(double T0);
  void Set_X0(double X0);
  void Set_V0(double V0);

  void print() const;

private:
  double T0_;
  double X0_;
  double V0_;
};


class Euler2{
 public:
  Euler2(double t0, double x0, double v0);
  ~Euler2();

  void Euler2B(double (*f) (double t, double x, double v), double N, double dt);

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

  void RungeK2(double (*f) (double t, double x, double v), double N, double dt);

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

 
