// compile with
// g++ ode_demo.cpp bla/calcinverse.cpp bla/exception.cpp bla/localheap.cpp -std=c++11 -o ode_demo -fmax-errors=2 -fdiagnostics-color=auto
//


// a simple linear algebra library
#include "bla/bla.hpp"
using namespace ngbla;

#include "ode.hpp"
#include "RK_orig.hpp"


/* ******************* The ODEs I want to solve  ******************** */


//  f(t,y) = a y 
class My_First_ODE_Function : public ODE_Function
{
  double a;

public:
  My_First_ODE_Function (double aa) { a = aa; }

  virtual void Eval (double t, const Vector<> & y, Vector<> & f) const
  {
    f = a * y;
  }
};





class MassSpring_ODE_Function : public ODE_Function
{
  double m;
  double k;
public:
  MassSpring_ODE_Function (double am, double ak) { m = am; k = ak; }

  // position, momentum
  virtual void Eval (double t, const Vector<> & y, Vector<> & f) const
  {
    f(0) = 1.0/m * y(1);
    f(1) = -k * y(0);
  }
};

class Pendulum_ODE_Function : public ODE_Function
{
  double length;
  double mass;
  double gravity=9.81;
 
public:
  Pendulum_ODE_Function(double _length, double _mass) : length(_length), mass(_mass){}
 
  virtual void Eval (double time, const Vector<> & y, Vector<> & f) const
  {
    f(0)=y(1);
    f(1)=-sin(y(0))*gravity/length;
  }
};


int main ()
{
cout <<"hello1" <<endl;
  ImplicitMidPoint imp_mp;
  ImplicitGauss imp_gauss;

  Vector<> y0(2);
  Pendulum_ODE_Function func(1.0,1.0);

  ofstream out("pendel_mp.data");
  y0(0)=M_PI/4;
  y0(1)=0.0;
  ODESolver (func, imp_mp, 0, y0, 10000, 0.0001, out);

  ofstream out2("pendel_gauss.data");
  y0(0)=M_PI/4;
  y0(1)=0.0;
  ODESolver (func, imp_gauss, 0, y0, 10000, 0.0001, out2);

  return 0;
}
