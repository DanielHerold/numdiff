// a simple ODE - solver library
// Joachim Schoeberl


// the base class for the right-hand-side f(t,y)
class ODE_Function
{
public:
  // must be overloaded by derived class
  virtual void Eval (double t, const Vector<> & y, Vector<> & f) const = 0;

  virtual void EvalDfDy (double t, const Vector<> & y, Matrix<> & dfdy) const
  {
    // numerical differentiation
    int n = y.Size();
    Vector<> yr(n), yl(n), fr(n), fl(n);
    double eps = 1e-6;
    for (int i = 0; i < n; i++)
      {
        yl = y;  yl(i) -= eps;
        yr = y;  yr(i) += eps;
        Eval (t, yl, fl);
        Eval (t, yr, fr);

	dfdy.Col(i) = 1.0/(2*eps) * (fr - fl);
      }
  }
};



// base class for the single-step method
class SSM
{
public:
  // do the step
  virtual void Step (double t, double h, const ODE_Function & func, 
                     const Vector<> & yold, Vector<> & ynew) const = 0;
};





// the time integration loop
void ODESolver (const ODE_Function & func, const SSM & ssm, 
		double t0, Vector<> & y0, double tend, double h,
		ostream & out)
{
  double t = t0;
  int n = y0.Size();
  double temp=0.0;
  int temp2=-1;

  Vector<> yold(n), ynew(n);

  yold = y0;
  while (t < tend)
    {
     if(t-temp>0.1)
     {
      out << t;
      for (int i = 0; i < n; i++)
        out << " " << yold(i);
      out << endl;
      temp=t;
     }

      if(t-temp2>1){printf("\r Bereits %f Sekunden, bzw %f Prozent...", t, t/(tend-t0)*100);temp2=t;}

      ssm.Step (t, h, func, yold, ynew);
      yold = ynew;
      t += h;
    }
}








/* *************** Here are the specific single-step methods *************** */



class ExplicitEuler : public SSM
{
public:
  virtual void Step (double t, double h, const ODE_Function & func,
                     const Vector<> & yold, Vector<> & ynew) const
  {
    Vector<> f(yold.Size());

    func.Eval (t, yold, f);
    ynew = yold + h * f;
  }
};


class ImprovedEuler : public SSM
{
public:
  virtual void Step (double t, double h, const ODE_Function & func,
                     const Vector<> & yold, Vector<> & ynew) const
  {
    Vector<> f(yold.Size());

    func.Eval (t, yold, f);
    ynew = yold + h/2.0 * f;

    func.Eval (t+h/2.0, ynew, f);
    ynew = yold + h * f;
  }
};

/*
class ImplicitEuler : public SSM
{
public:
  virtual void Step (double t, double h, const ODE_Function & func,
                     const Vector<> & yold, Vector<> & ynew) const
  {
    Matrix<> DF(yold.Size());
    Identity Id(yold.Size());
    Matrix<> funcdfdy(yold.Size());
    Matrix<> InvDF(yold.Size());

    Vector<> f(yold.Size());
    Vector<> update(yold.Size());

    double err=1;
    double epsilon = 1e-8;

    ynew=yold;

    int cnt=0;

    while ((err > epsilon) && (cnt <20))
      {
 	func.EvalDfDy(t+h,ynew,funcdfdy);
      
	DF = Id - h * funcdfdy;
	CalcInverse(DF,InvDF);
      
	func.Eval(t+h,ynew,f);

	update = InvDF * (ynew-yold-h*f);
	
	err = L2Norm(update);
	
	ynew -= update;
	cnt++;
      }

  }
};
*/
