

#include "modeldefault.hh"
#include "mhd_eqns.hh"
#include "rotator.hh"
namespace EULERCHORIN {
#include "chorjo.hh"
}
namespace EULERNUMFLUX {
#include "euler_flux.hpp"
}

template <int dimDomain>
class EulerFlux {
  enum { e = dimDomain+1};
 public:
  typedef FieldVector<double,dimDomain+2> RangeType;
  typedef FieldMatrix<double,dimDomain+2,dimDomain> FluxRangeType;
  typedef FieldVector<double, dimDomain> DomainType;
  inline void analyticalFlux(const double gamma,
			     const RangeType& u,
			     FluxRangeType& f) const;
  inline  void jacobian(const double gamma,
			const RangeType& u, 
			const FluxRangeType& du,
			RangeType& A) const; /* {
    A = 0.;
    DomainType velocity(u[0]);
    du.umv(velocity,A);
    }*/
  inline double pressure(const double gamma,const RangeType& u) const;
  inline double maxSpeed(const double gamma,
			 const FieldVector<double,dimDomain>& n,
			 const FieldVector<double,dimDomain+2>& u) const;
};
// ************************************************
template <class Model>
class DWNumFlux;
template <class Model>
class HLLNumFlux;
// ************************************************
template <int dimDomain>
class ConsVec : public FieldVector<double,dimDomain+2> {
 public:
  explicit ConsVec (const double& t) : FieldVector<double,dimDomain+2>(t) {}
  ConsVec () : FieldVector<double,dimDomain+2>(0) {}
};
template <class GridPart,int dimRange2,
	  int dimRange1=dimRange2*GridPart::GridType::dimensionworld>
class EulerModelTraits {
 public:
  typedef GridPart GridPartType;
  typedef typename GridPart::GridType GridType;
  enum { dimDomain = GridType::dimensionworld };
  enum { dimRange = dimRange2, dimGradRange = dimRange1 };
  typedef FieldVector<double, dimDomain> DomainType;
  typedef FieldVector<double, dimDomain-1> FaceDomainType;
  typedef FieldVector<double,dimRange> RangeType;
  // typedef ConsVec<dimDomain> RangeType;
  typedef FieldVector<double,dimGradRange> GradientType;
  typedef FieldMatrix<double,dimRange,dimDomain> FluxRangeType;
  typedef FieldMatrix<double,dimGradRange,dimDomain> DiffusionRangeType;
  typedef typename GridPart::IntersectionIteratorType IntersectionIterator;
  typedef typename GridType::template Codim<0>::Entity EntityType;
};
// ************************************************
template <class GridPartType,class ProblemType>
class EulerModel {
 public:
  enum { dimDomain = GridPartType::GridType::dimensionworld };
  enum { dimRange = dimDomain+2};
  typedef EulerModelTraits<GridPartType,dimRange> Traits;
  typedef typename Traits::GridType GridType;
  typedef typename Traits::RangeType RangeType;
  typedef typename Traits::DomainType DomainType;
  typedef typename Traits::FluxRangeType FluxRangeType;
  //typedef typename Traits::GradientType GradientType;
  //typedef typename Traits::DiffusionRangeType DiffusionRangeType;
 public:
  EulerModel(GridType& grid,const ProblemType& problem) :
    gamma_(problem.gamma),
    problem_(problem)
  {}

  // has flux 
  bool hasFlux() const { return true; }
  
  // has flux 
  bool hasSource() const { return false; }

  // needed for limitation
  inline  void velocity(
             const typename Traits::EntityType& en,
             double time,
             const DomainType& x,
             const RangeType& u,
             DomainType& velocity) const
  {
    for(int i=0; i<dimDomain; ++i) 
    {
      // check again 
      velocity[i] = u[i+1];
    }
  }

  inline  void analyticalFlux(typename Traits::EntityType& en,
			      double time,  
			      const typename Traits::DomainType& x,
			      const RangeType& u, 
			      FluxRangeType& f) const 
  {
    EulerFlux<dimDomain>().analyticalFlux(gamma_,u,f);
  }
  inline  void jacobian(const typename Traits::EntityType& en,
			double time,  
			const typename Traits::DomainType& x,
			const RangeType& u, 
			const FluxRangeType& du,
			RangeType& A) const {
    EulerFlux<dimDomain>().jacobian(gamma_,u,du,A);
  }
  inline bool hasBoundaryValue(typename Traits::IntersectionIterator& it,
			       double time, 
			       const typename Traits::FaceDomainType& x) const {
    return (abs(it.boundaryId()) != 4);
  }
  inline double boundaryFlux(typename Traits::IntersectionIterator& it,
			     double time, 
			     const typename Traits::FaceDomainType& x,
			     const RangeType& uLeft, 
			     RangeType& gLeft) const  {
    DomainType xgl=it.intersectionGlobal().global(x);
    const typename Traits::DomainType normal = it.integrationOuterNormal(x); 
    double p = EulerFlux<dimDomain>().pressure(gamma_,uLeft);
    gLeft = 0;
    for (int i=0;i<dimDomain; ++i) {
      gLeft[i+1] = normal[i]*p;
    }
    return 0.;
  }
  inline  void boundaryValue(typename Traits::IntersectionIterator& it,
			     double time, 
			     const typename Traits::FaceDomainType& x,
			     const RangeType& uLeft, 
			     RangeType& uRight) const {
    uRight=uLeft;
    if (abs(it.boundaryId()) == 1) {
      const typename Traits::DomainType normal = it.integrationOuterNormal(x);  
      double len2 = normal.two_norm2();
      double y = 0.;
      for (int i=0;i<dimDomain;i++) 
        y += uLeft[i+1]*normal[i];
      y *= 2./len2;
      for (int i=0;i<dimDomain;i++) 
        uRight[i+1] -= y*normal[i];
    } else if (abs(it.boundaryId()) == 2) {
      DomainType xgl=it.intersectionGlobal().global(x);
      problem_.evaluate(time,xgl,uRight);
    } else if (abs(it.boundaryId()) == 3) {
    } else {
      std::cerr << "Wrong Boundary ID " << it.boundaryId() << "\n\n";
    }
  }
  inline void maxSpeed(const typename Traits::DomainType& normal,
		       double time,  
		       const typename Traits::DomainType& x,
		       const RangeType& u,
		       double& advspeed,double& totalspeed) const {
    advspeed=EulerFlux<dimDomain>().maxSpeed(gamma_,normal,u);
    totalspeed=advspeed;
  }
  inline const ProblemType& problem() const {
    return problem_;
  }
 protected:
  double gamma_;
  double tstep_eps;
  ProblemType problem_;
  friend class DWNumFlux<EulerModel<GridPartType,ProblemType> >;
  friend class HLLNumFlux<EulerModel<GridPartType,ProblemType> >;
};
// ***********************
template <class GridPartType,class ProblemType>
class DWNumFlux<EulerModel<GridPartType,ProblemType> > {
 public:
  typedef typename GridPartType::GridType GridType;
  enum { dimDomain = GridType::dimensionworld };
  typedef Mhd::MhdSolver SolverType;
  typedef EulerModel<GridPartType,ProblemType> Model;
  typedef typename Model::Traits Traits;
  enum { dimRange = Model::dimRange };
  typedef typename Model::RangeType RangeType;
  typedef typename Model::FluxRangeType FluxRangeType;
  DWNumFlux(const Model& mod) : 
    model_(mod),
    eos(SolverType::Eosmode::me_ideal),
    numFlux_(eos,mod.gamma_),
    rot_(1) 
  {
    for(int i=0; i<9; ++i) 
    {
      ulmhd[i] = 0.0;
      urmhd[i] = 0.0;
      retmhd[i] = 0.0;
    }
  }
  double numericalFlux(typename Traits::IntersectionIterator& it,
		       double time,
		       const typename Traits::FaceDomainType& x,
		       const RangeType& uLeft,
		       const RangeType& uRight,
		       RangeType& gLeft,
		       RangeType& gRight) const;
  const Model& model() const {return model_;}
 private:
  Adi::FieldRotator<Model> rot_;
  const Model& model_;
  const typename SolverType::Eosmode::meos_t eos;
  mutable SolverType numFlux_;
  mutable SolverType::Vec9 ulmhd,urmhd,retmhd;
};
template <class GridPartType,class ProblemType>
double
DWNumFlux<EulerModel<GridPartType,ProblemType> > :: 
numericalFlux(typename DWNumFlux<EulerModel<GridPartType,ProblemType> >::Traits::
	      IntersectionIterator& it,
	      double time,
	      const typename 
	      DWNumFlux<EulerModel<GridPartType,ProblemType> >::Traits::
	      FaceDomainType& x,
	      const typename 
	      DWNumFlux<EulerModel<GridPartType,ProblemType> >:: RangeType& uLeft,
	      const typename
	      DWNumFlux<EulerModel<GridPartType,ProblemType> > :: RangeType& uRight,
	      typename
	      DWNumFlux<EulerModel<GridPartType,ProblemType> > :: RangeType& gLeft,
	      typename
	      DWNumFlux<EulerModel<GridPartType,ProblemType> > :: RangeType& gRight)
  const {
    typename Traits::DomainType normal = it.integrationOuterNormal(x);
    // double len = normal.two_norm();
    double len = it.intersectionGlobal().integrationElement(x);
    normal *= 1./len;
    RangeType& ul = const_cast<RangeType&>(uLeft);
    RangeType& ur = const_cast<RangeType&>(uRight);
    // RangeType ul(uLeft);
    // RangeType ur(uRight);
    rot_.rotateForth(ul, normal);
    rot_.rotateForth(ur, normal);
    /*
    ulmhd[0] = ulmhd[1] = ulmhd[2] = ulmhd[3] = ulmhd[4] =
      ulmhd[5] = ulmhd[6] = ulmhd[7] = ulmhd[8] = 0.;
    urmhd[0] = urmhd[1] = urmhd[2] = urmhd[3] = urmhd[4] =
      urmhd[5] = urmhd[6] = urmhd[7] = urmhd[8] = 0.;
    */
    ulmhd[0] = ul[0];
    urmhd[0] = ur[0];
    for (int i=0;i<dimDomain;++i) {
      ulmhd[i+1] = ul[i+1];
      urmhd[i+1] = ur[i+1];
    }
    ulmhd[7] = ul[1+dimDomain];
    urmhd[7] = ur[1+dimDomain];
    double p[3];
    double ldt=numFlux_(ulmhd,urmhd,p,retmhd);
    
    gLeft[0] = retmhd[0];
    for (int i=0;i<dimDomain;++i)
      gLeft[i+1] = retmhd[i+1];
    
    gLeft[1+dimDomain] = retmhd[7];
    
    rot_.rotateBack(gLeft,normal);
    gLeft *= len;
    
    gRight = gLeft;
    return ldt*len;
}
template <class GridPartType,class ProblemType>
class HLLNumFlux<EulerModel<GridPartType,ProblemType> > {
 public:
  typedef typename GridPartType::GridType GridType;
  enum { dimDomain = GridType::dimensionworld };
  typedef EulerModel<GridPartType,ProblemType> Model;
  typedef typename Model::Traits Traits;
  enum { dimRange = Model::dimRange };
  typedef typename Model::RangeType RangeType;
  typedef typename Model::FluxRangeType FluxRangeType;
  HLLNumFlux(const Model& mod) : 
    model_(mod),
    numFlux_(mod.gamma_)
  {}
  double numericalFlux(typename Traits::IntersectionIterator& it,
		       double time,
		       const typename Traits::FaceDomainType& x,
		       const RangeType& uLeft,
		       const RangeType& uRight,
		       RangeType& gLeft,
		       RangeType& gRight) const;
  const Model& model() const {return model_;}
 private:
  const Model& model_;
  EULERNUMFLUX::EulerFlux<dimDomain,EULERNUMFLUX::HLL> numFlux_;
};
template <class GridPartType,class ProblemType>
double
HLLNumFlux<EulerModel<GridPartType,ProblemType> > :: 
numericalFlux(typename HLLNumFlux<EulerModel<GridPartType,ProblemType> >::Traits::
	      IntersectionIterator& it,
	      double time,
	      const typename 
	      HLLNumFlux<EulerModel<GridPartType,ProblemType> >::Traits::
	      FaceDomainType& x,
	      const typename 
	      HLLNumFlux<EulerModel<GridPartType,ProblemType> >:: RangeType& uLeft,
	      const typename
	      HLLNumFlux<EulerModel<GridPartType,ProblemType> > :: RangeType& uRight,
	      typename
	      HLLNumFlux<EulerModel<GridPartType,ProblemType> > :: RangeType& gLeft,
	      typename
	      HLLNumFlux<EulerModel<GridPartType,ProblemType> > :: RangeType& gRight)
  const {
  typename Traits::DomainType normal = it.integrationOuterNormal(x);
  double len = it.intersectionGlobal().integrationElement(x);
  normal *= 1./len;
  double ldt = numFlux_.num_flux((&(uLeft[0])),
				 (&(uRight[0])),
				 (&(normal[0])),
				 (&(gLeft[0])));
  gLeft *= len;
  gRight = gLeft;
  return ldt*len;
}
// ***********************
template <>
inline
void EulerFlux<1>::analyticalFlux(const double gamma,
				  const FieldVector<double,1+2>& u,
				  FieldMatrix<double,1+2,1>& f) const {
  assert(u[0]>1e-10);
  double rhoeps = u[e]-0.5*u[1]/u[0]*u[1];
  assert(rhoeps>1e-10);
  double p = (gamma-1)*rhoeps;
  f[0][0] = u[1];
  f[1][0] = u[1]/u[0]*u[1]+p;
  f[e][0] = u[1]/u[0]*(u[e]+p);
}
template <>
inline
double EulerFlux<1>::pressure(const double gamma,const RangeType& u) const {
  assert(u[0]>1e-10);
  double rhoeps = u[e]-0.5*(u[1]*u[1])/u[0];
  /*
  if (rhoeps<1e-10) 
    cerr << "negative internal energy density " << rhoeps 
	 << " in analyticalFlux: "
	 << u << endl;
   */
  assert(rhoeps>1e-10);
  return (gamma-1)*rhoeps;
}
template <>
inline
double EulerFlux<2>::pressure(const double gamma,const RangeType& u) const {
  assert(u[0]>1e-10);
  double rhoeps = u[e]-0.5*(u[1]*u[1]+u[2]*u[2])/u[0];
  /*
  if (rhoeps<1e-10) 
    cerr << "negative internal energy density " << rhoeps 
	 << " in analyticalFlux: "
	 << u << endl;
   */
  assert(rhoeps>1e-10);
  return (gamma-1)*rhoeps;
}
template <>
inline
void EulerFlux<2>::analyticalFlux(const double gamma,
				  const FieldVector<double,2+2>& u,
				  FieldMatrix<double,2+2,2>& f) const {
  assert(u[0]>1e-10);
  const double rhoeps = u[e]-0.5*(u[1]*u[1]+u[2]*u[2])/u[0];
  /*
  if (rhoeps<1e-10) 
    cerr << "negative internal energy density " << rhoeps 
	 << " in analyticalFlux: "
	 << u << endl;
  */
  assert(rhoeps>1e-10);
  double v[2] = {u[1]/u[0],u[2]/u[0]};
  double p = (gamma-1)*rhoeps;
  f[0][0] = u[1];            f[0][1] = u[2];
  f[1][0] = v[0]*u[1]+p;     f[1][1] = v[1]*u[1];
  f[2][0] = v[0]*u[2];       f[2][1] = v[1]*u[2]+p;
  f[e][0] = v[0]*(u[e]+p);   f[e][1] = v[1]*(u[e]+p);
}

template <>
inline
double EulerFlux<3>::pressure(const double gamma,const RangeType& u) const {
  assert(u[0]>1e-10);
  const double rhoeps = u[e]-0.5*(u[1]*u[1] + u[2]*u[2] + u[3]*u[3])/u[0];
  /*
  if (rhoeps<1e-10) 
    cerr << "negative internal energy density " << rhoeps 
	 << " in analyticalFlux: "
	 << u << endl;
   */
  assert(rhoeps>1e-10);
  return (gamma-1)*rhoeps;
}
template <>
inline
void EulerFlux<3>::analyticalFlux(const double gamma,
				  const FieldVector<double,3+2>& u,
				  FieldMatrix<double,3+2,3>& f) const {
  assert(u[0]>1e-10);
  double rhoeps = u[e]-0.5*(u[1]*u[1]+u[2]*u[2]+u[3]*u[3])/u[0];
  assert(rhoeps>1e-10);
  double p = (gamma-1)*rhoeps;
  f[0][0]=u[1];              f[0][1]=u[2];             f[0][2]=u[3];
  f[1][0]=u[1]/u[0]*u[1]+p;  f[1][1]=u[2]/u[0]*u[1];   f[1][2]=u[3]/u[0]*u[1];
  f[2][0]=u[1]/u[0]*u[2];    f[2][1]=u[2]/u[0]*u[2]+p; f[2][2]=u[3]/u[0]*u[2];
  f[3][0]=u[1]/u[0]*u[3];    f[3][1]=u[2]/u[0]*u[2];   f[3][2]=u[3]/u[0]*u[3]+p;
  f[e][0]=u[1]/u[0]*(u[e]+p);f[e][1]=u[2]/u[0]*(u[e]+p);f[e][2]=u[3]/u[0]*(u[e]+p);
}
template <>
inline
void EulerFlux<2>::jacobian(const double gamma,
			const EulerFlux<2>::RangeType& u, 
			const EulerFlux<2>::FluxRangeType& du,
			EulerFlux<2>::RangeType& A) const {
  assert(u[0]>1e-10);
  /*
  double rhoeps = u[e]-0.5*(u[1]*u[1]+u[2]*u[2])/u[0];
  if (rhoeps<1e-10) 
    cerr << "negative internal energy density " << rhoeps 
	 << " in analyticalFlux: "
	 << u << endl;
  assert(rhoeps>1e-10);
  */
  double v[2] = {u[1]/u[0],u[2]/u[0]};
  
  A[0] = du[1][0] + du[1][1];
  A[1] = du[0][0]*((gamma-3.0)/2.0*v[0]*v[0] +(gamma-1.0)/2.0*v[1]*v[1])-
    du[0][1]*v[0]*v[1];
  A[1] += du[1][0]*(3.0-gamma)*v[0] + du[1][1]*v[1];
  A[1] += du[2][0]*(1.0-gamma)*v[1] + du[2][1]*v[0];
  A[1] += du[3][0]*(gamma-1.0);
  A[2] = du[0][1]*((gamma-3.0)/2.0*v[1]*v[1] +(gamma-1.0)/2.0*v[0]*v[0])-
    du[0][0]*v[0]*v[1];
  A[2] += du[1][1]*(1.0-gamma)*v[0] + du[1][0]*v[1];
  A[2] += du[2][1]*(3.0-gamma)*v[1] + du[2][0]*v[0];
  A[2] += du[3][1]*(gamma-1.0);
  A[3] = du[0][0]*(-gamma*u[3]*v[0]/u[0]+(gamma-1.0)*v[0]*(v[0]*v[0]+v[1]*v[1]))+
    du[0][1]*(-gamma*u[3]*v[1]/u[0]+(gamma-1.0)*v[1]*(v[0]*v[0]+v[1]*v[1]));
  A[3] += du[1][0]*(gamma*u[3]/u[0]-(gamma-1.0)/2.0*(3.0*v[0]*v[0]+v[1]*v[1]))-
    du[1][1]*(gamma-1.0)*v[0]*v[1];
  A[3] += -du[2][0]*(gamma-1.0)*v[0]*v[1]+
    du[2][1]*(gamma*u[3]/u[0]-(gamma-1.0)/2.0*(v[0]*v[0]+3.0*v[1]*v[1]));
  A[3] += du[3][0]*gamma*v[0] + du[3][1]*gamma*v[1];
}
template <>
inline
double EulerFlux<1>::maxSpeed(const double gamma,
			      const FieldVector<double,1>& n,
			      const FieldVector<double,1+2>& u) const {
  assert(u[0]>1e-10);
  double u_normal = u[1]*n[0] / u[0];
  double rhoeps = u[e]-0.5*u[1]*u[1]/u[0];
  assert(rhoeps>1e-10);
  double p = (gamma-1)*rhoeps;
  double c2 = gamma*p/u[0]*n.two_norm2();
  assert(c2>1e-10);
  if (0>c2) {
    printf("Negative sound speed!\n");
    // abort();
  }
  return fabs(u_normal) + sqrt(c2);
}
template <>
inline
double EulerFlux<2>::maxSpeed(const double gamma,
			      const FieldVector<double,2>& n,
			      const FieldVector<double,2+2>& u) const {
  assert(u[0]>1e-10);
  double u_normal = (u[1]*n[0]+u[2]*n[1]) / u[0];
  double rhoeps = u[e]-0.5*(u[1]*u[1]+u[2]*u[2])/u[0];
  assert(rhoeps>1e-10);
  double p = (gamma-1)*rhoeps;
  double c2 = gamma*p/u[0]*n.two_norm2();
  assert(c2>1e-10);
  if (0>c2) {
    // printf("Negative sound speed!\n");
    // abort();
    // c2 = 0;
  }
  return fabs(u_normal) + sqrt(c2);
}
template <>
inline
double EulerFlux<3>::maxSpeed(const double gamma,
			      const FieldVector<double,3>& n,
			      const FieldVector<double,3+2>& u) const {
  assert(u[0]>1e-10);
  double u_normal = (u[1]*n[0]+u[2]*n[1]+u[3]*n[2]) / u[0];
  double rhoeps = u[e]-0.5*(u[1]*u[1]+u[2]*u[2]+u[3]*u[3])/u[0];
  assert(rhoeps>1e-10);
  double p = (gamma-1)*rhoeps;
  double c2 = gamma*p/u[0]*n.two_norm2();
  assert(c2>1e-10);
  if (0>c2) {
    printf("Negative sound speed!\n");
    abort();
  }
  return fabs(u_normal) + sqrt(c2);
}

/*****************************************************************/
// Initial Data
class U0Smooth1D {
public:
  U0Smooth1D(double,bool diff_timestep=true) : gamma(1.4),myName("Advection") {}

  void printmyInfo(string filename) {}
  double endtime() {
    return 0.2;
  }
  double gamma;
  string myName;

  template <class DomainType, class RangeType>
  void evaluate(const DomainType& arg, RangeType& res) const {
    evaluate(0,arg,res);
  }
  template <class DomainType, class RangeType>
  void evaluate(const DomainType& arg,double t, RangeType& res) const 
  {
    evaluate(t,arg,res);
  }
  template <class DomainType, class RangeType>
  void evaluate(double t,const DomainType& arg, RangeType& res) const {
    DomainType c(0.0);
    DomainType x = arg;
    x -= c;
    double r2 = 0.25*0.25;
    res=0;
    if (x*x<r2) {
      res[0] = 1.0;
      res[3] = 1.0/(gamma-1.);
    } else {
      res[0] = 0.125;
      res[3] = 0.1/(gamma-1.);
    }
    return;
    res[1] = cos(0.2*M_PI);
    res[2] = sin(0.2*M_PI);
    x[0] -= t*res[1];
    x[1] -= t*res[2];
    if (x[0]*x[0]+x[1]*x[1]<r2) {
      res[0] = cos((x[0]*x[0]+x[1]*x[1])/r2*M_PI)+1.;
      res[0] *= res[0];
      res[0] /= 4.;
    } else {
      res[0] = 0.;
    }
    res[0] += 0.5;
    /*
    if (arg[0]*arg[0] < 0.25) {
      // res[0] = cos(arg[0]*M_PI*2.)+2;
      res[0] = -8.*(arg[0]*arg[0]-0.25)+1.;
    }
    else {
      res[0] = 1.0;
    }
    */
    res[1] = cos(0.2*M_PI);
    res[2] = sin(0.2*M_PI);
    res[3] = 0.3/(gamma-1.0);
    res[1] *= res[0];
    res[2] *= res[0];
    res[3] += 0.5*(res[1]*res[1]+res[2]*res[2])/res[0];
  }
};
class U0RotatingCone {
public:
  U0RotatingCone() : gamma(1.4) {}
  template <class DomainType, class RangeType>
  void evaluate(const DomainType& arg, RangeType& res) const {
    evaluate(0,arg,res);
  }
  template <class DomainType, class RangeType>
  void evaluate(double t,const DomainType& arg, RangeType& res) const {
    res*=0.;
    DomainType c(0.5);
    DomainType x=arg;
    x-=c;
    double r2=0.04;
    if (x*x < r2) {
      res[0] =cos(x*x/r2*M_PI)+2;
    }
    else {
      res[0] = 1.0;
    }
    x=arg;
    x-=DomainType(1.0);
    if (DomainType::size>1) {
      res[1] = x[1]*res[0];
      res[2] = -x[0]*res[0];
    } else {
      res[1] = -1.*res[0];
    }
    res[DomainType::size+1] = 2.;
    if (arg.size>1) {
      res[DomainType::size+1] += 
	0.5*(res[1]*res[1]+res[2]*res[2])/res[0];
    } else {
      res[DomainType::size+1] += 
	0.5*(res[1]*res[1])/res[0];
    }
  }
  double gamma;
};
class U0VW {
public:
  U0VW() : gamma(1.4) {}
  template <class DomainType, class RangeType>
  void evaluate(const DomainType& arg, RangeType& res) const {
    evaluate(0,arg,res);
  }
  template <class DomainType, class RangeType>
  void evaluate(double t,const DomainType& arg, RangeType& res) const {
    if (arg[0]<0.25) {
      res[0]=1.;
      res[1]=-1.;
      res[1]=0.;
      res[2]=0.;
      res[3]=1./(1.4-1.0);
    } else {
      res[0]=1.; 
      res[1]=1.;
      res[2]=0.;
      res[3]=1.0/(1.4-1.0);
    }
    res[1] *= res[0];
    res[2] *= res[0];
    res[3] += 0.5*res[1]*res[1]/res[0];
  }
  double gamma;
};
class U0Sod {
  double ql,qr,ul,ur,pl,pr;
public:
  U0Sod() : gamma(1.4) {
    myName = "RP-Sod";
    ql = 1.0;
    qr = 0.125;
    ul = 0.;
    ur = 0.;
    pl = 1.0;
    pr = 0.1;
  }
  U0Sod(double eps,int flag,bool diff_timestep=true) 
    : gamma(1.4) {
    myName = "RP-Sod";
    ql = 1.0;
    qr = 0.125;
    ul = 0.;
    ur = 0.;
    pl = 1.0;
    pr = 0.1;
  }
  double endtime() {
    return 0.4;
  }
  double saveinterval() {
    return 0.01;
  }
  template <class DomainType, class RangeType>
  void evaluate(const DomainType& arg, RangeType& res) const {
    evaluate(0.,arg,res);
  }
  template <class DomainType, class RangeType>
  void evaluate(double t,const DomainType& arg, 
		RangeType& res) const {
    double x = arg[0] - 0.5;
    if (t>1e-8)
      chorin(t,x,res[0],res[1],res[2]);
    else {
      if (x<0.) {
	res[0]=ql; // 0.5;
	res[1]=ul;
	// res[2]=0.;
	res[2]=pl;
      } else {
	res[0]=qr;
	res[1]=ur;
	// res[2]=0.;
	res[2]=pr;
      }
    }
    res[1] *= res[0];
    // res[2] *= res[0];
    res[2] = res[2]/(gamma-1.0)+0.5*(res[1]*res[1])/res[0];
  }
  
  void chorin(double t,double x,
	      double& q_erg,double& u_erg,double& p_erg) const
  {
    EULERCHORIN::
      lsg(x,t,&q_erg,&u_erg,&p_erg,ql,qr,ul,ur,pl,pr,gamma);
  }

  void printmyInfo(std::string filename)
  {
    std::ostringstream filestream;
    filestream << filename;

    std::ofstream ofs(filestream.str().c_str(), std::ios::app);
	
    ofs << "Problem: " << myName << "\n\n"
	<< "gamma = " << gamma << "\n\n";
    ofs	<< "\n\n";
	
    ofs.close();
			
  }
  double gamma;
  std::string myName;
};
class FFS {
public:
  FFS(double eps,int flag,bool diff_timestep=true) : gamma(1.4) {myName = "FFS";}
  double endtime() {
    return 3.0;
  }
  double saveinterval() {
    return 0.1;
  }
  template <class DomainType, class RangeType>
  void evaluate(const DomainType& arg, RangeType& res) const {
    evaluate(0,arg,res);
  }
  template <class DomainType, class RangeType>
  void evaluate(double t,const DomainType& arg, RangeType& res) const {
    res[0]=1.4;
    res[1]=3.0*1.4;
    res[2]=0.;
    res[3]=8.8;
  }
  void printmyInfo(std::string filename)
  {
    std::ostringstream filestream;
    filestream << filename;

    std::ofstream ofs(filestream.str().c_str(), std::ios::app);
	
    ofs << "Problem: " << myName << "\n\n"
	<< "gamma = " << gamma << "\n\n";
    ofs	<< "\n\n";
	
    ofs.close();
			
  }
  double gamma;
  std::string myName;
};
