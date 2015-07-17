#include <cmath>

template<typename Func>
class RK4
{
public:

    // constructors
    RK4():t(0), h(1), y(0), f(0) {}
    ~RK4() {}

    // set 
    void setInitialValue(double t0, double y0) {
	t = t0;
	y = y0;
    }
    void setTimeStep(double dt) {h = dt;}
    void setFunction(Func& fu) {f = &fu;}

    // get
    double getSolution() const {return y;}
    double getCurrentTime() const {return t;}

    // time step
    void next() {
	double k1 = h*(*f)(t,y);
	double k2 = h*(*f)(t+h/2.0,y+k1/2.0);
	double k3 = h*(*f)(t+h/2.0,y+k2/2.0);
	double k4 = h*(*f)(t+h,y+k3);

	y += (k1+2*k2+2*k3+k4)/6.0;
	t += h;
    }

    
private:
    double t,h,y;
    Func* f;
};

template<typename Func>
class RKF4
{
public:

    // constructors
    RKF4():t(0), h(1), y(0), f(0), tol(1e-6) {}
    ~RKF4() {}

    // set 
    void setInitialValue(double t0, double y0) {
	t = t0;
	y = y0;
    }
    void setTimeStep(double dt) {h = dt;}
    void setFunction(Func& fu) {f = &fu;}
    void setTolerance(double t) {tol = t;}

    // get
    double getSolution() const {return y;}
    double getCurrentTime() const {return t;}
    double getTimeStep() const {return h;}

    // time step
    bool next() {
	double k1 = h*(*f)(t,y);
	double k2 = h*(*f)(t+h/4.0,y+k1/4.0);
	double k3 = h*(*f)(t+3*h/8.0,y+3*k1/32.0+9*k2/32.0);
	double k4 = h*(*f)(t+12*h/13.0,y+1932*k1/2197.0-7200*k2/2197.0+7296*k3/2197.0);
	double k5 = h*(*f)(t+h,y+439*k1/216.0-8*k2+3680*k3/513.0-845*k4/4104.0);
	double k6 = h*(*f)(t+h/2.0,y-8*k1/27.0+2*k2-3544*k3/2565.0+1859*k4/4104.0-11*k5/40.0);
	double dy = 25*k1/216.0+1408*k3/2565.0+2197*k4/4104.0-k5/5.0;
	double dy1 = 16*k1/135.0+6656*k3/12825.0+28561*k4/56430.0-9*k5/50.0+2*k6/55.0;
	double R = fabs(dy-dy1)/h;
	double delta = 0.84*pow(tol/R,0.25);

	if(R <= tol) {
	    y += dy;
	    t += h;
	    h *= delta;
	    return true;
	} else {
	    h *= delta;
	    return false;
	}

	return true;
    }

    
private:
    double t,h,y,tol;
    Func* f;
};
