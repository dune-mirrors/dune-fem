#include <config.h>
#include <iostream>
#include <fem/operator/common/spaceoperatorif.hh>
#include <fem/solver/rungekutta.hh>
#include <fem/solver/multistepAdams.hh>

using namespace Dune;
using namespace DuneODE;
using namespace std;

struct SpaceDummy {
};
struct Dest {
  typedef double DomainFieldType;
  typedef double RangeFieldType;
  typedef SpaceDummy DiscreteFunctionSpaceType;
  Dest(string,const SpaceDummy&) {}
  Dest() {}
  RangeFieldType operator[](int i) const {
    return d_;
  }
  RangeFieldType& operator[](int i) {
    return d_;
  }
  void assign(const Dest& other) {
    d_ = other.d_;
  }
  void addScaled(const Dest& other,RangeFieldType l) {
    d_ += other.d_*l;
  }
  Dest& operator*=(RangeFieldType l) {
    d_ *= l;
    return *this;
  }
  Dest& operator+=(const Dest& other) {
    d_ += other.d_;
    return *this;
  }
  Dest& operator-=(const Dest& other) {
    d_ += other.d_;
    return *this;
  }
  RangeFieldType d_;
};
struct RHS : SpaceOperatorInterface<Dest> {
  TimeProvider* tp_;
  SpaceType space_;
  mutable int eval;
  RHS() : eval(0) {}
  const SpaceType& space() const {return space_;}
  void operator()(const DestinationType& x,
                  DestinationType& y) const {
    y[0]=x[0];
    ++eval;
  }
  void timeProvider(TimeProvider* tp) {
    tp_ = tp;
  }
};

int main() {
  RHS rhs;
  TimeProvider tp(0.,1.);
  tp.setDeltaT(0.1);
  double cfl = 1.;
  ExplicitMultiStepSolver<Dest> rk(rhs,tp,2,true);
  cfl /= 2.;
  // ExplicitRungeKuttaSolver<Dest> rk(rhs,tp,2,true);

  Dest u;
  u[0] = 1.;
  rk.initialize(u);
  tp.setCfl(cfl);
  std::cout << tp.time() << " " << rhs.eval 
            << " " << u[0] << " " 
            << 1.0 << " " << tp.deltaT() << std::endl;
  while (tp.time()<1.) {
    double dt = tp.deltaT();
    rk.solve(u);
    tp.augmentTime();
    
    tp.setDeltaT(dt*(1.+0.5*(2.*double(random())/double(RAND_MAX)-1.0)));
    std::cout << tp.time() << " " << rhs.eval 
              << " " << u[0] << " " 
              << tp.deltaT()/dt << " " << tp.deltaT() << std::endl;
  }
}
