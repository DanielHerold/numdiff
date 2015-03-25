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





class B13_ODE_Function : public ODE_Function
{
public:


  B13_ODE_Function (int i) {}

  virtual void Eval (double t, const Vector<> & y, Vector<> & f) const
  {
    f(0) =0.001*exp((sin(t)-y(0))/0.2)-0.001-y(0);
  }
};







int main ()
{
  ExplicitEuler expl_euler;
  ImprovedEuler impr_euler;

  ofstream out2("B13.out");
  B13_ODE_Function func(1);
  Vector<> y0(1);
  y0(0) = 1.0;
  ODESolver (func, impr_euler, 0, y0, 30, 0.01, out2);

  return 0;
}
