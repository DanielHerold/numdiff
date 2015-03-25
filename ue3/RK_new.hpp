// Runge Kutta methods 
// Joachim Schoeberl

class RungeKuttaMethod : public SSM
{
protected:
  int stages;
  Matrix<> A;
  Vector<> b,c;
public:
  RungeKuttaMethod (int s)
    : stages(s), A(s,s), b(s), c(s) 
  { ; }
 
  void SetAbc (const Matrix<> & aa,
	       const Vector<> & ab,
	       const Vector<> & ac)
  {
    A = aa;
    b = ab;
    c = ac;
  }
};



class ExplicitRKMethod : public RungeKuttaMethod
{
public:
  ExplicitRKMethod (int as) : RungeKuttaMethod (as){ ; }

  
  virtual void Step (double t, double h, const ODE_Function & func,
                     const Vector<> & yold, Vector<> & ynew) const
  {
    int n = yold.Size();

    Vector<> yi(n), ki(n);
    Vector<> all_ki(stages*n);

    for (int i = 0; i < stages; i++)
      {
	double ti = t + h * c(i);
	yi = yold;
	for (int j = 0; j < i; j++)
	  yi += h * A(i,j) * all_ki.Range (j*n, (j+1)*n);
	func.Eval (ti, yi, ki);
	all_ki.Range(i*n, (i+1)*n) = ki;
      }


    ynew = yold;
    for (int i = 0; i < stages; i++)
      ynew += h * b(i) * all_ki.Range (i*n, (i+1)*n);
  }
};


class NewtonF
{
  double t;
  Vector<> yold;
  Vector<> c;
  Matrix<> A;
  double h;

  public:
    NewtonF(double tt, double hh, 
            const Vector<> & yyold, const Matrix<> & AA, const Vector<> & cc)
	    : t(tt), yold(yyold), c(cc), A(AA), h(hh){};

    void Eval (const Vector<> & all_ki_old, const ODE_Function & func, Vector<> & all_ki_new) const
      {
         int s=c.Size();
         int n=yold.Size();
         Vector<> sum(n);
	 Vector<> evalfunc(n);
         for (int i=0; i<s; i++)
         {
            sum = 0.0;
            for (int l=0; l<s; l++)
	    {sum += A(i,l)*all_ki_old.Range(i*n,(i+1)*n);}
            func.Eval(t+c(i)*h, yold+h*sum, evalfunc);
            all_ki_new.Range(i*n, (i+1)*n) = all_ki_new.Range(i*n,(i+1)*n)-evalfunc;
         }

      }
     
    void EvalDfDy (const Vector<> & y, const ODE_Function & func, Matrix<> & dfdy) const
    {
    int n = y.Size();
    Vector<> yr(n), yl(n), fr(n), fl(n);
    double eps = 1e-6;
    for (int i = 0; i < n; i++)
      {
        yl = y;  yl(i) -= eps;
        yr = y;  yr(i) += eps;
        Eval (yl, func, fl);
        Eval (yr, func, fr);

        dfdy.Col(i) = 1.0/(2*eps) * (fr - fl);
      }

    }

};





void Newton(const NewtonF & F, const ODE_Function & func, Vector<> & zero)
{
  double tol=1e-6;


  int sn= zero.Size();
  Vector<> delta(sn);
  Matrix<> dF(sn);
  Vector<> evalF(sn);
  Vector<> delta_old(sn);
  double q=0;

  F.EvalDfDy(zero, func, dF);


if(((double)L2Norm(dF))<= tol){zero=0.0; return;}
  CalcInverse(dF, dF);
  F.Eval(zero, func, evalF);


  delta=-dF*evalF;
  zero+=delta;
  delta_old=delta;

int test;

  do{
    F.Eval(zero, func, evalF);
    F.EvalDfDy(zero, func, dF);
    CalcInverse(dF, dF);
    delta=-dF*evalF;
    zero+=delta;

cout << "delta="<<delta<<endl;

    q=L2Norm(delta)/L2Norm(delta_old);
    delta_old=delta;
    if (q>=1) {cout << "Newton konvergiert nicht bei t="<<F.t<<endl;}
  }while(q/(1-q)*L2Norm(delta)>tol);

};



  

class ImplicitRKMethod : public RungeKuttaMethod
{
public:
  ImplicitRKMethod (int as) : RungeKuttaMethod (as){ ; }

  
  virtual void Step (double t, double h, const ODE_Function & func,
                     const Vector<> & yold, Vector<> & ynew) const
  {
    int n = yold.Size();

    Vector<> all_ki(stages*n);
    all_ki=0.0;

    NewtonF F(t, h, yold, A, c);
    Newton(F, func, all_ki);

    ynew = yold;
    for (int i = 0; i < stages; i++)
      ynew += h * b(i) * all_ki.Range (i*n, (i+1)*n);
  }
};

//////////////////////////////////////////////////////////////////////  

class ClassicalRK : public ExplicitRKMethod
{
public:
  ClassicalRK() : ExplicitRKMethod (4)
  {
    Matrix<> A(4,4);
    Vector<> b(4), c(4);

	c = { 0, 0.5, 0.5, 1 };
    
	b = { 1.0 / 6, 2.0 / 6, 2.0 / 6, 1.0 / 6 };

    A = 0.0;
    A(1,0) = 0.5;
    A(2,1) = 0.5;    
    A(3,2) = 1;
    
    SetAbc (A, b, c);
  }
};


class ImprovedEulerRK : public ExplicitRKMethod
{
public:
  ImprovedEulerRK() : ExplicitRKMethod (2)
  {
    Matrix<> A(2,2);
    Vector<> b(2), c(2);

    c = { 0, 0.5 };
    b = { 0.0, 1.0 };

    A = 0.0;
    A(1,0) = 0.5;
    
    SetAbc (A, b, c);
  }
};


class ImplicitMidPoint : public ImplicitRKMethod
{
public:
  ImplicitMidPoint() : ImplicitRKMethod (1)
  {
    Matrix<> A(1,1);
    Vector<> b(1),c(1);
    A=0.5;
    b=1;
    c=0.5;

    SetAbc (A,b,c);
  }
};


class ImplicitGauss : public ImplicitRKMethod
{
public:
  ImplicitGauss() : ImplicitRKMethod (2)
  {
    Matrix<> A(2,2);
    Vector<> b(2),c(2);

    A=0.25;
    A(1,0)=0.25+sqrt(3)/6;
    A(0,1)=0.25-sqrt(3)/6;
    b=0.5;
    c={0.5-sqrt(3)/6, 0.5+sqrt(3)/6};

    SetAbc (A,b,c);
  }
};
