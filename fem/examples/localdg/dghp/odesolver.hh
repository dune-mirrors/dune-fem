#ifndef RUNGEKUTTA_ODE_SOLVER_HH
#define RUNGEKUTTA_ODE_SOLVER_HH

#include <iostream>
#include <cmath>
#include <dune/fem/misc/timeutility.hh>
namespace DuneODE {
  using namespace Dune;
  using namespace std;
#include "localfunction.hh"
template<class Operator,int STEPS>
class ExplRungeKutta : public TimeProvider {
 public:
  enum {maxord=10,ord_=STEPS};
  typedef typename Operator::SpaceType SpaceType;
  typedef typename Operator::SpaceType FunctionSpaceType;
  typedef typename Operator::DestinationType DestinationType;
  typedef typename DestinationType::LocalFunctionType LocalFunctionType; 
 private:
  double cfl_;
  double **a;
  double *b;
  double *c;
  double dt;
  // int ord_;
  const SpaceType& spc_;
  DestinationType* U0;
  LocalFunctionType* LU0;
  std::vector<DestinationType*> Upd;
  std::vector<LocalFunctionType*> LUpd;
  mutable LocalFuncHelper<LocalFunctionType> tmp_;
public:
  typedef typename DestinationType::RangeType RangeType;
  typedef typename DestinationType::DomainType DomainType;
  typedef typename DestinationType::JacobianRangeType JacobianRangeType;
  ExplRungeKutta(Operator& op,int pord,double cfl) :
    op_(op),
    spc_(op.space()),
    tmp_(op.space()),
    cfl_(cfl), 
    // ord_(pord), 
    Upd(0),
    savetime_(0.0), savestep_(1)
  {
    op.timeProvider(this);
    assert(ord_>0);
    a=new double*[ord_];
    for (int i=0;i<ord_;i++)
      a[i]=new double[ord_];
    b=new double [ord_];
    c=new double [ord_];
    switch (ord_) {
    case 4 :
      a[0][0]=0.;     a[0][1]=0.;     a[0][2]=0.;    a[0][3]=0.;
      a[1][0]=1.0;    a[1][1]=0.;     a[1][2]=0.;    a[1][3]=0.;
      a[2][0]=0.25;   a[2][1]=0.25;   a[2][2]=0.;    a[2][3]=0.;
      a[3][0]=1./6.;  a[3][1]=1./6.;  a[3][2]=2./3.; a[3][3]=0.;
      b[0]=1./6.;     b[1]=1./6.;     b[2]=2./3.;    b[3]=0.;
      c[0]=0.;        c[1]=1.0;       c[2]=0.5;      c[3]=1.0;
      break;
    case 3 :
      a[0][0]=0.;     a[0][1]=0.;     a[0][2]=0.;
      a[1][0]=1.0;    a[1][1]=0.;     a[1][2]=0.;
      a[2][0]=0.25;   a[2][1]=0.25;   a[2][2]=0.;
      b[0]=1./6.;     b[1]=1./6.;     b[2]=2./3.;
      c[0]=0.;        c[1]=1;         c[2]=0.5;
      break;
    case 2 :
      a[0][0]=0.;     a[0][1]=0.;
      a[1][0]=1.0;    a[1][1]=0.;
      b[0]=0.5;       b[1]=0.5;
      c[0]=0;         c[1]=1;
      break;
    case 1:
      a[0][0]=0.;
      b[0]=1.;
      c[0]=0.;
      break;
    default : std::cerr << "Runge-Kutta method of this order not implemented" 
			<< std::endl;
              abort();
    }
    U0 = new DestinationType("Start",op_.space());
    LU0 = new LocalFunctionType(*U0);
    for (int i=0;i<ord_;i++) {
      Upd.push_back(new DestinationType("URK",op_.space()) );
      LUpd.push_back(new LocalFunctionType(*(Upd[i])));
    }
    Upd.push_back(new DestinationType("Ustep",op_.space()) );
  }
  void timecoeff(double s,double *ret) const {
    assert(0<=s && s<=1.);
    switch (ord_) {
    case 4:
      ret[0]=2./3.*s*s*s-3./2.*s*s+s;
      ret[1]=-1./3.*s*s*s+1./2.*s*s;
      ret[2]=-4./3.*s*s*s+2.*s*s;
      ret[3]=s*s*s-s*s;
      break;
    case 3 :
      ret[0]=(6.*c[0]-3.)*b[0]*s*s+(4.-6.*c[0])*b[0]*s;   
      ret[1]=(6.*c[1]-3.)*b[1]*s*s+(4.-6.*c[1])*b[1]*s;   
      ret[2]=(6.*c[2]-3.)*b[2]*s*s+(4.-6.*c[2])*b[2]*s;   
      break;
    case 2:
      ret[0]=(b[0]-1.)*s*s+s;
      ret[1]=b[1]*s*s;
      break;
    case 1:
      ret[0]=s;
      break;
    }
  }
  void dtimecoeff(double s,double *ret) const {
    assert(0<=s && s<=1.);
    switch (ord_) {
    case 4:
      ret[0]=2.*s*s-3.*s+1.;
      ret[1]=-s*s+s;
      ret[2]=-4.*s*s+4.*s;
      ret[3]=3.*s*s-2.*s;
    case 3 :
      ret[0]=(6.*c[0]-3.)*b[0]*2.*s+(4.-6.*c[0])*b[0];   
      ret[1]=(6.*c[1]-3.)*b[1]*2.*s+(4.-6.*c[1])*b[1];   
      ret[2]=(6.*c[2]-3.)*b[2]*2.*s+(4.-6.*c[2])*b[2];   
      break;
    case 2:
      ret[0]=(b[0]-1.)*2.*s+1.;
      ret[1]=b[1]*2.*s;
      break;
    case 1:
      ret[0]=1.;
      break;
    }
  }
  /*****************************************************************/
  template <class EntityType>
  void setEntity(const EntityType& en) {
    LU0->init(en);
    for (int i=0;i<ord_;i++) {
      LUpd[i]->init(en);
    }
    tmp_.init(en);
  }
  int numDofs() {
    return LU0->numDofs();
  }
  const SpaceType& space() {
    return spc_;
  }
  /*
  template <class EntityType,class QuadratureType>
  void ucomponents(const EntityType& en,
		   const QuadratureType& quad,int l,double t,
		   std::vector<RangeType>& comp) const {
    assert(0<=t && t<=1.);
    RangeType ret(0);
    double coeff[ord_];
    timecoeff(t,coeff);
    tmp_.assign(*LU0);
    for (int j=0;j<ord_;j++) {
      tmp_.addscaled((*LUpd[j]),dt*coeff[j]);
    }
    tmp_.components(en,quad,l,comp);
  }
  */
  template <class EntityType,class QuadratureType>
  RangeType uval(const EntityType& en,
		 const QuadratureType& quad,int l,double t,
		 int maxp) const {
    assert(0<=t && t<=1.);
    RangeType ret(0);
    double coeff[ord_];
    timecoeff(t,coeff);
    tmp_.assign(*LU0);
    for (int j=0;j<ord_;j++) {
      tmp_.addscaled((*LUpd[j]),dt*coeff[j]);
    }
    tmp_.evaluateLocal(en,quad,l,maxp,ret);
    return ret;
  }
  template <class EntityType,class QuadratureType>
  JacobianRangeType dxuval(const EntityType& en,
			   const QuadratureType& quad,int l,double t,
         int maxp) const {
    assert(0<=t && t<=1.);
    JacobianRangeType ret(0);
    double coeff[ord_];
    timecoeff(t,coeff);
    tmp_.assign(*LU0);
    for (int j=0;j<ord_;j++) {
      tmp_.addscaled(*(LUpd[j]),dt*coeff[j]);
    }
    tmp_.jacobianLocal(en,quad,l,maxp,ret);
    return ret;
  }
  template <class EntityType,class QuadratureType>
  RangeType dtuval(const EntityType& en,
		   const QuadratureType& quad,int l,double t,
       int maxp) const {
    assert(0<=t && t<=1.);
    RangeType ret(0);
    double coeff[ord_];
    dtimecoeff(t,coeff);
    tmp_.assign(0.);
    for (int j=0;j<ord_;j++) {
      tmp_.addscaled(*(LUpd[j]),dt*coeff[j]);     
    }
    tmp_.evaluateLocal(en,quad,l,maxp,ret);
    ret /= dt;
    return ret;
  }
  double solve(typename Operator::DestinationType& U) {
    resetTimeStepEstimate();
    double t=time();
    U0->assign(U);
    // Compute Steps
    op_(U,*(Upd[0]));
    dt=cfl_*timeStepEstimate();
    for (int i=1;i<ord_;i++) {
      (Upd[ord_])->assign(U);
      for (int j=0;j<i;j++) 
	(Upd[ord_])->addScaled(*(Upd[j]),(a[i][j]*dt));
      setTime(t+c[i]*dt);
      op_(*(Upd[ord_]),*(Upd[i]));
      double ldt=cfl_*timeStepEstimate();
    }
    // Perform Update
    for (int j=0;j<ord_;j++) {
      U.addScaled(*(Upd[j]),(b[j]*dt));
    }
    setTime(t+dt);
    return time();
  }
  template <class Projection>
  void project(typename Operator::DestinationType& U,
	       Projection& proj) {
    typedef typename Operator::DestinationType DestType;
    typedef typename DestType::LocalFunctionType LFuncType;
    typedef typename DestType::FunctionSpaceType::IteratorType IteratorType;
    typedef typename DestType::FunctionSpaceType::BaseFunctionSetType BaseFunctionSetType;
    IteratorType endit = U.space().end();
    for(IteratorType it = U.space().begin(); 
	it != endit ; ++it) {
      int deg = proj.usePolDeg(*it);
      LFuncType lU = U.localFunction(*it);
      Dune::DofConversionUtility< PointBased > dofConversion(RangeType::size);
      for (int i = 0; i < lU.numDofs(); ++i) {
	int phiord = dofConversion.containedDof(i);
	if (phiord>=numPol[deg]) {
	  lU[i] = 0.;
	}
      }
    }
  }
  template <class Projection>
  double solve(typename Operator::DestinationType& U,
	       typename Operator::DestinationType& V,
	       Projection& proj) {
    resetTimeStepEstimate();
    double t=time();
    U0->assign(U);
    project(*U0,proj);
    // Compute Steps
    op_(*U0,*(Upd[0]));
    dt=cfl_*timeStepEstimate();
    for (int i=1;i<ord_;i++) {
      (Upd[ord_])->assign(U);
      for (int j=0;j<i;j++) 
	(Upd[ord_])->addScaled(*(Upd[j]),(a[i][j]*dt));
      setTime(t+c[i]*dt);
      project(*(Upd[ord_]),proj);
      op_(*(Upd[ord_]),*(Upd[i]));
      double ldt=cfl_*timeStepEstimate();
    }
    // Perform Update
    V.assign(*U0);
    for (int j=0;j<ord_;j++) {
      V.addScaled(*(Upd[j]),(b[j]*dt));
    }
    project(V,proj);
    setTime(t+dt);
    return time();
  }
  void printGrid(int nr, 
		 const typename Operator::DestinationType& U) {
    if (time()>=savetime_) {
      printSGrid(time(),savestep_*10+nr,op_.space(),U);
      ++savestep_;
      savetime_+=0.001;
    }
  }
  void printmyInfo(string filename) const {
    std::ostringstream filestream;
    filestream << filename;
    std::ofstream ofs(filestream.str().c_str(), std::ios::app);
    ofs << "ExplRungeKutta, steps: " << ord_ << "\n\n";
    ofs << "                cfl: " << cfl_ << "\\\\\n\n";
    ofs.close();
    op_.printmyInfo(filename);
  }
 private:
  const Operator& op_;
  int savestep_;
  double savetime_;
};

}
#endif
