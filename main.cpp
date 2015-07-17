#include "rk4.hpp"
#include <iostream>
#include <cmath>
using namespace std;

class F
{
public:
    double operator()(double t, double y) {
	return y-t*t+1;
    }
};

class Y
{
public:
    double operator()(double t) {
	return t*t+2*t+1-0.5*exp(t);
    }
};

void analysis(double t0, double t1, double y0, double h, double tol)
{
    RK4<F> rk4;
    F func;
    Y soln;
    
    rk4.setInitialValue(t0,y0);
    rk4.setTimeStep(h);
    rk4.setFunction(func);

    cout<<"h = "<<h<<"\n";
    cout<<"t Exact Numerical error\n";
    for(double t=t0+h; t<t1+h/2; t+=h) {
	rk4.next();
	double y = rk4.getSolution();
	double yexact = soln(t);
	cout<<t<<" "<<yexact<<" "<<y<<" "<<fabs(yexact-y)<<"\n";
    }
    cout<<"\n\n\n";

    RKF4<F> rkf4;
    rkf4.setInitialValue(t0,y0);
    rkf4.setTimeStep(h);
    rkf4.setFunction(func);
    rkf4.setTolerance(tol);
    
    cout<<"h = "<<h<<", tol = "<<tol<<"\n";
    cout<<"h t Exact Numerical error\n";
    double t = rkf4.getCurrentTime();
    h = rkf4.getTimeStep();
    while(t < t1) {
	bool res = rkf4.next();
	if(!res) continue;
	t = rkf4.getCurrentTime();
	h = rkf4.getTimeStep();
	double y = rkf4.getSolution();
	double yexact = soln(t);
	cout<<h<<" "<<t<<" "<<yexact<<" "<<y<<" "<<fabs(yexact-y)<<"\n";
    }
    cout<<"\n\n\n";

}

int main()
{
    double t0 = 0;
    double t1 = 2;
    double y0 = 0.5;
    double h = 0.5;
    double tol = 1e-6;

    analysis(t0,t1,y0,h,tol);

    h = 0.2;
    analysis(t0,t1,y0,h,tol);

    h = 0.05;
    analysis(t0,t1,y0,h,tol);

    return 0;
}
