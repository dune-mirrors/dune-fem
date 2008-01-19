 #define USERK 0
int order = 2;

#include <cmath>
#include <cassert>


#ifdef NOFEL
const int size = NOFEL;
#else
#warning SETTING SIZE to 100
const int size = 1;
#endif

/*
double Len = 2.*M_PI;
const double MAXT = 2.;
const double savestep = 0.05;
*/
bool usePAdapt = false;
double* h; 

#include <config.h>
#include <iostream>
#include <fem/operator/common/spaceoperatorif.hh>
#include <fem/solver/rungekutta.hh>
#include <fem/solver/multistepAdams.hh>

using namespace Dune;
using namespace DuneODE;
using namespace std;

#include "grid.hh"
#include "bucklev.hh"

template <int N>
struct Dest : public FieldVector<double,N> {
  typedef FieldVector<double,N> BaseType;
  typedef double DomainFieldType;
  typedef double RangeFieldType;
  typedef SpaceDummy DiscreteFunctionSpaceType;
  typedef Dest<N> ScalarDestinationType;
  Dest(string,const SpaceDummy&) {}
  Dest() {}
  void clear() {
    BaseType::operator=(0.);
  }
  void assign(const Dest& other) {
    BaseType::operator=(other);
  }
  void addScaled(const Dest& other,RangeFieldType l) {
    (*this).axpy(l,other);
  }
  double operator()(int i) const {
    if (i<0 || i>this->size-1) {
      std::cout << "ERROR: Accessing element " << i << std::endl;
    }
    return (*this)[i];
  }
};

Dest<size> x;

struct Heat : SpaceOperatorInterface<Dest<size> > {
  double w0;
  double visc(double x) const {
    return 1.;
    return 1.-log(1.-x*w0);
    return ((atan(-(x-M_PI))+1.6)/3.)*0.5*(sin(x)+1.25);
  }
  double dvisc(double x) const {
    return 0.;
    return w0/(1.-x*w0);
  }
  double lsg(double t,double x) const {
    // return 1.+0.5*t*t;
    // return cos(x) + t*t*t; // sin(t*M_PI);
    return cos(x)*exp(-t*visc(x)) + sin(t*M_PI);
  }
  TimeProvider* tp_;
  SpaceType space_;
  mutable int eval;
  Heat() : eval(0), tp_(0), w0(-0.855/M_PI) {}
  const double length() const {
    return 2.*M_PI;
  }
  const double maxTime() const {
    return 10.;
  }
  const double saveStep() const {
    return 0.1;
  }
  const SpaceType& space() const {return space_;}
  double mass(double x,double gridh) {
    return gridh*gridh/visc(x);
  }
  double cfl(double dt,double h,double gridh) {
    return dt/h*length();
  }
  double flux(int i,int j,const DestinationType& u) const {
    ++eval;
    // return 0;
    // return u[j] - u[i];
    // int i0  = (i<j)?i:j;
    // int i1  = i0+1;
    // int i2  = (i1==size-1)i1:i1+1;
    // int i_1 = (i0==0)?i0:i0-1;
    // int i2  = (i1==size-1)?0:i1+1;
    // int i_1 = (i0==0)?size-1:i0-1;
    int i0  = i;
    int i1  = j;
    int i2  = (i1==size-1)?0:i1+1;
    int i_1 = (i0==0)?size-1:i0-1; 
    // Higher order
    double ret = (-u(i2)+u(i_1)+15.*u(i1)-15.*u(i0))/12.;
    // Second order
    // double ret = u[i1]-u[i0];
    return ret;
  }
  void applyIntersection (const  
	       SpaceType::IntersectionIteratorType& iter,
	       double time,
	       const DestinationType& x, 
	       double& ret) const {
    SpaceType::IteratorType::EntityType in = iter.inside();
    SpaceType::IteratorType::EntityType out = iter.outside();
    if (iter.numberInSelf()==1)
      ret = flux(in,out,x);
    else {
      ret = flux(out,in,x);
      ret *= -1;
    }
  }
  double elementFlux (int i,
		      double time,
		      const DestinationType& u) const {
    ++eval;
    // return 0.;
    // return time;
    return M_PI*cos(M_PI*time);
    // return -visc(x[i])*u[i]; 
    //        +time+visc(x[i])*0.5*time*time; 
    double f= -2.*visc(x[i])*time*
      dvisc(x[i])*exp(-time*visc(x[i]))*sin(x[i])+
      visc(x[i])*time*dvisc(x[i])*dvisc(x[i])*u[i]*(1.-time);
    if (i==0)
      f -= (u[0]-
	    lsg(time,2.*x[0]-x[1]))/h[0];
    if (i==size-1)
      f += (lsg(tp_->time(),2.*x[size-1]-x[size-2])-
	    u[size-1])/h[size-1];
    return f;
  }
  void applyElement (const 
               SpaceType::IteratorType::EntityType& en,
	       double time,
	       const DestinationType& u,
	       double& ret) const {
    ret = elementFlux(en,time,u);
  }
		     
  void operator()(const DestinationType& u,
                  DestinationType& y) const {
    for (int i=0;i<size;i++) {
      y[i]   = elementFlux(i,tp_->time(),u);
    }
    for (int i=0;i<size-1;i++) {
      double f = flux(i,i+1,u);
      y[i]   += f/h[i];
      y[i+1] -= f/h[i+1];
    }
    double f = flux(size-1,0,u);
    y[size-1] += f/h[size-1];
    y[0]      -= f/h[0];
  }
  void timeProvider(TimeProvider* tp) {
    tp_ = tp;
  }
};
struct Advect : SpaceOperatorInterface<Dest<size> > {
  double a;
  TimeProvider* tp_;
  SpaceType space_;
  mutable int eval;
  double lsg(double t,double x) const {
    double s = x-a*t;
    while (s>length()) s-=length();
    while (s<0) s+=length();
    return cos(4.*M_PI*s);
    if (s>0.1 && s<0.3)
      return 1.;
    else
      return 0.;
  }
  Advect() : eval(0), tp_(0), a(2.) {
  }
  const double length() const {
    return 1.;
  }
  const double maxTime() const {
    return 5.;
  }
  const double saveStep() const {
    return 0.1;
  }
  const SpaceType& space() const {return space_;}
  double mass(double x,double gridh) {
    return gridh/a;
  }
  double cfl(double dt,double h,double gridh) {
    return dt/h;
  }
  double flux(double ui,double uj) const {
    ++eval;
    double ret;
    ret = -ui;
    return ret;
  }
  double minmod(double a,double b) const {
    double ret=0;
    // return ret;
    if (a*b>0) {
      if (fabs(a)<fabs(b))
	ret=a;
      else 
	ret=b;
    } else ret=0;
    return ret;
  }
  void applyIntersection (const  
	       SpaceType::IntersectionIteratorType& iter,
	       double time,
	       const DestinationType& u, 
	       double& ret,
         int& localOrder) const {
    localOrder=2;
    SpaceType::IteratorType::EntityType in = iter.inside();
    SpaceType::IteratorType::EntityType out = iter.outside();
    typedef SpaceType::IntersectionIteratorType 
      IntersectionIteratorType;
    double aIn,aOut;
    {
      double W0,W[2];
      W0=u[in];
      IntersectionIteratorType nbend=space().iend(in);
      for (IntersectionIteratorType nbiter=space().ibegin(in);
	   !(nbiter==nbend);++nbiter) {
	W[nbiter.numberInSelf()] = u[nbiter.outside()];
      }
      aIn = minmod(W0-W[0],W[1]-W0);
    }
    {
      double W0,W[2];
      W0=u[out];
      IntersectionIteratorType nbend=space().iend(out);
      for (IntersectionIteratorType nbiter=space().ibegin(out);
	   !(nbiter==nbend);++nbiter) {
	W[nbiter.numberInSelf()] = u[nbiter.outside()];
      }
      aOut = minmod(W0-W[0],W[1]-W0);
    }
    /*
    if (iter.numberInSelf()==1)
      ret = flux(u[in],u[out]);
    else {
      ret = flux(u[out],u[in]);
      ret *= -1;
    }
    */
    if (iter.numberInSelf()==1)
      ret = flux(u[in]+0.5*aIn,u[out]-0.5*aOut);
    else {
      ret = flux(u[out]+0.5*aOut,u[in]-0.5*aIn);
      ret *= -1;
    }
  }
  double elementFlux (int i,
		      double time,
		      const DestinationType& u) const {
    ++eval;
    double f = 0.;
    return f;
  }
  void applyElement (const 
               SpaceType::IteratorType::EntityType& en,
	       double time,
	       const DestinationType& u,
	       double& ret) const {
    ret = elementFlux(en,time,u);
  }
		     
  void operator()(const DestinationType& u,
                  DestinationType& y) const {
    y.clear();
    SpaceType::IteratorType end=space().end();
    for (SpaceType::IteratorType 
	   iter=this->space().begin();
	 !(iter==end);++iter) {
      double value;
      applyElement(*iter,tp_->time(),u,value);
      y[*iter] += value;
      SpaceType::IntersectionIteratorType nbend=
	space().iend(*iter);
      for (SpaceType::IntersectionIteratorType nbiter=
	     space().ibegin(*iter);
	   !(nbiter==nbend);++nbiter) {
	if (nbiter.inside()<nbiter.outside()) {
	  double flux;
    int order;
	  applyIntersection(nbiter,tp_->time(),u,flux,order);
	  y[nbiter.inside()]  += flux/h[nbiter.inside()];
	  y[nbiter.outside()] -= flux/h[nbiter.outside()];
	}
      }
    }
  }
  void timeProvider(TimeProvider* tp) {
    tp_ = tp;
  }
};
struct BuckleyLevert : SpaceOperatorInterface<Dest<size> > {
  // Size=400 dt=0.001
  U0BuckLev sol;
  TimeProvider* tp_;
  SpaceType space_;
  mutable int eval;
  double lsg(double t,double x) const {
    double u;
    sol.evaluate(t,x,u);
    return u;
  }
  BuckleyLevert() : eval(0), tp_(0) {
  }
  const double length() const {
    return 1.;
  }
  const double maxTime() const {
    return 0.5;
  }
  const double saveStep() const {
    return 0.05;
  }
  double maxTime() {
    return 0.25;
  }
  const SpaceType& space() const {return space_;}
  double mass(double x,double gridh) {
    return gridh;
  }
  double cfl(double dt,double h,double gridh) {
    return dt/h*2.1;
  }
  double flux(double ui,double uj) const {
    ++eval;
    double ret;
    ret = -sol.f(ui);
    return ret;
  }
  double minmod(double a,double b) const {
    double ret=0;
    // return ret;
    if (a*b>0) {
      if (fabs(a)<fabs(b))
	ret=a;
      else 
	ret=b;
    } else ret=0;
    return ret;
  }
  void applyIntersection (const  
	       SpaceType::IntersectionIteratorType& iter,
	       double time,
	       const DestinationType& u, 
	       double& ret,
	       int& localOrder) const {
    SpaceType::IteratorType::EntityType in = iter.inside();
    SpaceType::IteratorType::EntityType out = iter.outside();
    typedef SpaceType::IntersectionIteratorType 
      IntersectionIteratorType;
    double aIn=0,aOut=0;
    if (!usePAdapt || fabs(u[in]-u[out])<0.1) {
      {
	double W0,W[2];
	W0=u[in];
	IntersectionIteratorType nbend=space().iend(in);
	for (IntersectionIteratorType 
	       nbiter=space().ibegin(in);
	     !(nbiter==nbend);++nbiter) {
	  W[nbiter.numberInSelf()] = u[nbiter.outside()];
	}
	aIn = minmod(W0-W[0],W[1]-W0);
      }
      {
	double W0,W[2];
	W0=u[out];
	IntersectionIteratorType nbend=space().iend(out);
	for (IntersectionIteratorType 
	       nbiter=space().ibegin(out);
	     !(nbiter==nbend);++nbiter) {
	  W[nbiter.numberInSelf()] = u[nbiter.outside()];
	}
	aOut = minmod(W0-W[0],W[1]-W0);
	localOrder = 2;
      }
    } else {
      localOrder = 1;
    }
    if (iter.numberInSelf()==1)
      ret = flux(u[in]+0.5*aIn,u[out]-0.5*aOut);
    else {
      ret = flux(u[out]+0.5*aOut,u[in]-0.5*aIn);
      ret *= -1;
    }
  }
  double elementFlux (int i,
		      double time,
		      const DestinationType& u) const {
    ++eval;
    double f = 0.;
    return f;
  }
  void applyElement (const 
               SpaceType::IteratorType::EntityType& en,
	       double time,
	       const DestinationType& u,
	       double& ret) const {
    ret = elementFlux(en,time,u);
  }
		     
  void operator()(const DestinationType& u,
                  DestinationType& y) const {
    y.clear();
    SpaceType::IteratorType end=space().end();
    for (SpaceType::IteratorType 
	   iter=this->space().begin();
	 !(iter==end);++iter) {
      double value;
      applyElement(*iter,tp_->time(),u,value);
      y[*iter] += value;
      SpaceType::IntersectionIteratorType nbend=
	space().iend(*iter);
      for (SpaceType::IntersectionIteratorType nbiter=
	     space().ibegin(*iter);
	   !(nbiter==nbend);++nbiter) {
	if (nbiter.inside()<nbiter.outside()) {
	  double flux;
	  int localOrder;
	  applyIntersection(nbiter,tp_->time(),u,flux,
			    localOrder);
	  y[nbiter.inside()]  += flux/h[nbiter.inside()];
	  y[nbiter.outside()] -= flux/h[nbiter.outside()];
	}
      }
    }
  }
  void timeProvider(TimeProvider* tp) {
    tp_ = tp;
  }
};

int main (int argc, char **argv) {
  if(argc != 2)
  {
    fprintf(stderr,"usage: %s dt \n",argv[0]);
    exit(1);
  }
  // Heat rhs;
  Advect rhs;
  // BuckleyLevert rhs;
  TimeProvider tp(0.,1.);
  rhs.timeProvider(&tp);
  tp.setDeltaT(atof(argv[1]));
  double cfl = 1.;
  #if USERK
  ExplicitRungeKuttaSolver<Dest<size> > rk(rhs,tp,order,true);
  // ExplicitMultiStepSolver<Dest<size> > rk(rhs,tp,order,true); cfl /= double(order);
  // cfl /= 16.;
  #else
  ExplicitMultiStepLTSolver<Dest<size> > rk(rhs,tp,order,true);
  // cfl /= 16.;
  cfl /= double(order);
  #endif
  h = new double[size];
  double minh = 10.;
  tp.setCfl(cfl);
  double meshWidth = rhs.length()/double(size);
  Dest<size> u;
  for (int i=0;i<size;i++) {
    x[i] = rhs.length()*(double(i)+0.5)/double(size);
    u[i] = rhs.lsg(0,x[i]);
    h[i] = rhs.mass(x[i],meshWidth);
    if (h[i]<minh) minh = h[i];
  }
  /* two odes 
  u[0] = 0.6;
  u[1] = 0.4;
  h[0] = 0.1;
  h[1] = 0.5;
  */
  /* one ode 
  u[0] = 1.;
  h[0] = 1.;
  */
  std::cout << "#  " << minh << " " << " " << tp.deltaT() 
	    << " "
	    << tp.deltaT()/minh*rhs.length() << std::endl; 
  rk.initialize(u);
  int steps = 0,oldfluxes = 0;
  std::cout.precision(16);
  std::cout << scientific;
  double savetime = 0.;
  double mxerror=0;
  double l1error=0;
  double mass=0;
  while (1) {
    if (tp.time()>=savetime || size<=4) {
      l1error=0;
      mxerror=0;
      mass=0;
      std::cout <<  " ############# Step : " << steps << std::endl;
      std::cerr <<  " ############# Step : " << steps << std::endl;
      for (int i=0;i<size;i++) {
#if USERK
	double ldt  = tp.deltaT();
#else
	double ldt  = (*rk.dt_[rk.dt_.size()-1])[i];
#endif
	double lerr = fabs(u[i] - rhs.lsg(tp.time(),x[i]));
	mass += u[i]*meshWidth;
	l1error += lerr*meshWidth;
	mxerror = (mxerror<lerr)?lerr:mxerror;
	std::cout << x[i] 
		  << " " << u[i] 
		  << " " << rhs.lsg(tp.time(),x[i])
		  << " " << tp.time()
	          << " " << tp.deltaT()
		  << " " << ldt
		  << " " << rhs.cfl(ldt,h[i],meshWidth)
		  << " " << steps
		  << std::endl;
      }
      /* two odes
      std::cout << tp.time() << " " 
		<< u[0] << " " << u[1] 
		<< std::endl;
      */
    }
    if (tp.time()>=savetime) {
      std::cout << "# " << tp.time() << " " << tp.deltaT() 
		<< " "
		<< rhs.eval << " " << steps
		<< " # output" << std::endl;
      std::cout << "# " << tp.time() << " " 
		<< " Linf-Error: " << mxerror << " "
		<< " L1-Error: " << l1error << " ";
      std::cout << "  Mass: " << mass << std::endl;
    }
    if (tp.time()>=savetime) {
      std::cout << std::endl << std::endl;
      savetime += rhs.saveStep();
    }
    if (oldfluxes == rhs.eval) {
      std::cout << "ERROR: EMPTY UPDATE CYCLE!" << std::endl;
      // return 0;
    }
    oldfluxes = rhs.eval; 
    rk.solve(u);
#if USERK
    tp.augmentTime();
#endif
    steps++;
    if (USERK) {
       steps++;
       // if (steps==4 && USERK) 
       // tp.setDeltaT(tp.deltaT()*0.25*0.25);
    }
    if (tp.time()>rhs.maxTime()) 
      break;
  }
  mxerror=0;
  l1error=0;
  mass=0;
  for (int i=0;i<size;i++) {
#if USERK
    double ldt  = tp.deltaT();
#else
    double ldt  = (*rk.dt_[0])[i];
#endif
    double lerr = fabs(u[i] - rhs.lsg(tp.time(),x[i]));
    l1error += lerr*meshWidth;
    mass += u[i]*meshWidth;
    mxerror = (mxerror<lerr)?lerr:mxerror;
    std::cout << x[i] 
	      << " " << u[i] 
	      << " " << rhs.lsg(tp.time(),x[i])
	      << " " << tp.time()
	      << " " << tp.deltaT()
	      << " " << ldt
	      << " " << rhs.cfl(ldt,h[i],meshWidth)
	      << std::endl;
  }
  /* two odes
     std::cout << tp.time() << " " 
     << u[0] << " " << u[1] 
     << std::endl;
  */
  std::cout << "# " << tp.time() << " " << tp.deltaT() 
	    << " "
	    << rhs.eval << " " << steps
	    << " # output" << std::endl;
  std::cout << "# Linf-Error: " << mxerror << " "
	    << "# L1-Error: " << l1error << " ";
  std::cout << "  Mass: " << mass << std::endl;
}
