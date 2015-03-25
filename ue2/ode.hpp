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
  int write_old=-1;

  Vector<> yold(n), ynew(n);

  yold = y0;
  while (t < tend)
    {
      if((int)(4*t)>(write_old)){
      write_old++;
      out << t;
      for (int i = 0; i < n; i++)
        out << " " << yold(i);
      out << endl;
printf("\r Bereits %f Sekunden, bzw. %f ...",t , t/(tend-t0)*100 );
//cout<<"now"<<endl;
      }


      ssm.Step (t, h, func, yold, ynew);
      yold = ynew;
      t += h;
    }
}


/* Implicit Euler */

class Impl_Euler_ODE_Function
{

  double h;

  Vector<> yold;
  int size;

public:
  Impl_Euler_ODE_Function (double h,  const Vector<> yold) : h(h), yold(yold) {size=yold.Size();}

  virtual void Eval(double time, const Vector<> &y, const ODE_Function & f, Vector<> & F) const
  {
   Vector<> f_eval(size);
   f.Eval(time, y, f_eval);
   for (int i=0; i<size; i++)
	{
	 F(i)=yold(i)-y(i)+h*f_eval(i);
/////////////////////////////
//////////	 cout<<"Auswertung F("<<i<<"= "<<F(i)<< endl;

	};
  }

  virtual void EvalDfDy (double t, const Vector<> & y, const ODE_Function & f,Matrix<> & dfdy) const
  {
    // numerical differentiation
    int n = y.Size();
    Vector<> yr(n), yl(n), fr(n), fl(n);
    double eps = 1e-6;
    for (int i = 0; i < n; i++)
      {
        yl = y;  yl(i) -= eps;
        yr = y;  yr(i) += eps;
        Eval (t, yl, f,fl);
        Eval (t, yr, f,fr);

	dfdy.Col(i) = 1.0/(2*eps) * (fr - fl);
      }
  }
};



void newton(const Impl_Euler_ODE_Function & func, const ODE_Function & ode_f, double time, const Vector<> &x, double tol, Vector<> &zero)
{
  int size = x.Size();
  Vector<> delta_x(size);
  delta_x=0.0;
  Vector<> delta_x_old(size);
  delta_x_old=0.0;
  Matrix<> df(size);
  Vector<> eval_f(size);
  double q=0.0;

  func.EvalDfDy(time, x,ode_f, df);
////////////////////7
////////cout <<"df= "<<df<<endl;
if(((double)L2Norm(df))<=tol){zero=x; return;};
  CalcInverse(df,df);
///////////////////
//cout <<"1/df= "<<df<<endl;

  func.Eval(time, x, ode_f, eval_f);

/////////////////////
/////////cout <<"f= "<<eval_f<<endl;

  delta_x_old=df*eval_f;
  zero=x-delta_x_old;
///////////////////
//cout << delta_x_old<< endl;

  do{
//     func.EvalDfDy(time, zero, ode_f,df);
//     CalcInverse(df,df);
     func.Eval(time, zero, ode_f, eval_f);

     delta_x=df*eval_f;
     zero=zero-delta_x;

//cout<< delta_x<<endl;

     q=L2Norm(delta_x)/L2Norm(delta_x_old);
     if (q>=1) { cout << "Newtonverfahren konvergiert nicht." <<delta_x<<"|"<<delta_x_old<< endl; break;}
     delta_x_old=delta_x;
////////////////////
/////////cout<<delta_x<<endl;


}while(q/(1-q)*L2Norm(delta_x)>tol);

////////////////////
//cout<<"Newtonverfahren konv bei "<< zero<<endl;

}; 







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

class ImplicitEuler : public SSM
{
public:
  virtual void Step (double t, double h, const ODE_Function & func,
                     const Vector<> & yold, Vector<> & ynew) const
  {
    Impl_Euler_ODE_Function F(h,yold);

///////////////////

    newton(F,func,t+h, yold, 1e-2, ynew);
////////////////
/////////cout << "Impl Euler mit startwerten, schritten, zeit und endwerten: "<< yold<<h<<t<<ynew<<endl;
  }
};





 
