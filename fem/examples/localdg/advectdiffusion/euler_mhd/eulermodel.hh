

#include "modeldefault.hh"
#include "mhd_eqns.hh"
#include "rotator.hh"

template <int dimDomain>
class EulerFlux {
  enum { e = dimDomain+1};
 public:
  inline void analyticalFlux(const double gamma,
			     const FieldVector<double,dimDomain+2>& u,
			     FieldMatrix<double,dimDomain+2,dimDomain>& f) const;
  inline double maxSpeed(const double gamma,
			 const FieldVector<double,dimDomain>& n,
			 const FieldVector<double,dimDomain+2>& u) const;
};
// ************************************************
template <class Model>
class DWNumFlux;
// ************************************************
template <int dimDomain>
class ConsVec : public FieldVector<double,dimDomain+2> {
 public:
  explicit ConsVec (const double& t) : FieldVector<double,dimDomain+2>(t) {}
  ConsVec () : FieldVector<double,dimDomain+2>(0) {}
};
template <class Grid,int dimRange2,
	  int dimRange1=dimRange2*Grid::dimensionworld>
class EulerModelTraits {
 public:
  typedef Grid GridType;
  enum { dimDomain = GridType::dimensionworld };
  enum { dimRange = dimRange2, dimGradRange = dimRange1 };
  typedef FieldVector<double, dimDomain> DomainType;
  typedef FieldVector<double, dimDomain-1> FaceDomainType;
  typedef FieldVector<double,dimRange> RangeType;
  // typedef ConsVec<dimDomain> RangeType;
  typedef FieldVector<double,dimGradRange> GradientType;
  typedef FieldMatrix<double,dimRange,dimDomain> FluxRangeType;
  typedef FieldMatrix<double,dimGradRange,dimDomain> DiffusionRangeType;
  typedef typename GridType::Traits::IntersectionIterator IntersectionIterator;
  typedef typename GridType::template Codim<0>::Entity EntityType;
};
// ************************************************
template <class GridType,class ProblemType>
class EulerModel {
 public:
  enum { dimDomain = GridType::dimensionworld };
  enum { dimRange = dimDomain+2};
  typedef EulerModelTraits<GridType,dimRange> Traits;
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
  inline  void analyticalFlux(typename Traits::EntityType& en,
			      double time,  
			      const typename Traits::DomainType& x,
			      const RangeType& u, 
			      FluxRangeType& f) const {
    EulerFlux<dimDomain>().analyticalFlux(gamma_,u,f);
  }
  inline bool hasBoundaryValue(typename Traits::IntersectionIterator& it,
			       double time, 
			       const typename Traits::FaceDomainType& x) const {
    return true;
  }
  inline double boundaryFlux(typename Traits::IntersectionIterator& it,
			     double time, 
			     const typename Traits::FaceDomainType& x,
			     const RangeType& uLeft, 
			     RangeType& gLeft) const  {
    DomainType xgl=it.intersectionGlobal().global(x);
    const typename Traits::DomainType normal = it.integrationOuterNormal(x);  
    RangeType u;
    problem_.evaluate(time,xgl,u);
    FluxRangeType f;
    EulerFlux<dimDomain>().analyticalFlux(gamma_,u,f);
    gLeft*=0;
    f.umv(normal,gLeft);
    return 0.;
  }
  inline  void boundaryValue(typename Traits::IntersectionIterator& it,
			     double time, 
			     const typename Traits::FaceDomainType& x,
			     const RangeType& uLeft, 
			     RangeType& uRight) const {
    // uRight=uLeft;
    DomainType xgl=it.intersectionGlobal().global(x);
    problem_.evaluate(time,xgl,uRight);
  }
  inline void maxSpeed(const typename Traits::DomainType& normal,
		       double time,  
		       const typename Traits::DomainType& x,
		       const RangeType& u,
		       double& advspeed,double& totalspeed) const {
    advspeed=EulerFlux<dimDomain>().maxSpeed(gamma_,normal,u);
    totalspeed=advspeed;
  }
 protected:
  double gamma_;
  double tstep_eps;
  ProblemType problem_;
  friend class DWNumFlux<EulerModel<GridType,ProblemType> >;
};
// ***********************
template <class GridType,class ProblemType>
class DWNumFlux<EulerModel<GridType,ProblemType> > {
 public:
  enum { dimDomain = GridType::dimensionworld };
  typedef Mhd::MhdSolver SolverType;
  typedef EulerModel<GridType,ProblemType> Model;
  typedef typename Model::Traits Traits;
  enum { dimRange = Model::dimRange };
  typedef typename Model::RangeType RangeType;
  typedef typename Model::FluxRangeType FluxRangeType;
  DWNumFlux(const Model& mod) : 
    model_(mod),
    eos(SolverType::Eosmode::me_ideal),
    numFlux_(eos,mod.gamma_),
    rot_(1) {
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
};
template <class GridType,class ProblemType>
double
DWNumFlux<EulerModel<GridType,ProblemType> > :: 
numericalFlux(typename DWNumFlux<EulerModel<GridType,ProblemType> >::Traits::
	      IntersectionIterator& it,
	      double time,
	      const typename 
	      DWNumFlux<EulerModel<GridType,ProblemType> >::Traits::
	      FaceDomainType& x,
	      const typename 
	      DWNumFlux<EulerModel<GridType,ProblemType> >:: RangeType& uLeft,
	      const typename
	      DWNumFlux<EulerModel<GridType,ProblemType> > :: RangeType& uRight,
	      typename
	      DWNumFlux<EulerModel<GridType,ProblemType> > :: RangeType& gLeft,
	      typename
	      DWNumFlux<EulerModel<GridType,ProblemType> > :: RangeType& gRight)
  const {
    typename Traits::DomainType normal = it.integrationOuterNormal(x);
    double len = normal.two_norm();
    normal /= len;

    RangeType ul,ur;
    ul = uLeft;
    ur = uRight;
    rot_.rotateForth(ul, normal);
    rot_.rotateForth(ur, normal);
    SolverType::Vec9 ulmhd,urmhd,retmhd;
    ulmhd[0] = ulmhd[1] = ulmhd[2] = ulmhd[3] = ulmhd[4] =
      ulmhd[5] = ulmhd[6] = ulmhd[7] = ulmhd[8] = 0.;
    urmhd[0] = urmhd[1] = urmhd[2] = urmhd[3] = urmhd[4] =
      urmhd[5] = urmhd[6] = urmhd[7] = urmhd[8] = 0.;
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
void EulerFlux<2>::analyticalFlux(const double gamma,
				  const FieldVector<double,2+2>& u,
				  FieldMatrix<double,2+2,2>& f) const {
  assert(u[0]>1e-10);
  double rhoeps = u[e]-0.5*(u[1]*u[1]+u[2]*u[2])/u[0];
  if (rhoeps<1e-10) 
      cerr << "negative internal energy density in analyticalFlux: "
	   << u << " " << rhoeps << endl;
  assert(rhoeps>1e-10);
  double p = (gamma-1)*rhoeps;
  f[0][0] = u[1];                 f[0][1] = u[2];
  f[1][0] = u[1]/u[0]*u[1]+p;     f[1][1] = u[2]/u[0]*u[1];
  f[2][0] = u[1]/u[0]*u[2];       f[2][1] = u[2]/u[0]*u[2]+p;
  f[e][0] = u[1]/u[0]*(u[e]+p);   f[e][1] = u[2]/u[0]*(u[e]+p);
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
  return fabs(u_normal) + sqrt(c2);
}

/*****************************************************************/
// Initial Data
class U0Smooth1D {
public:
  U0Smooth1D() : gamma(1.4) {}
  template <class DomainType, class RangeType>
  void evaluate(const DomainType& arg, RangeType& res) const {
    evaluate(0,arg,res);
  }
  template <class DomainType, class RangeType>
  void evaluate(double t,const DomainType& arg, RangeType& res) const {
    if (arg[0]*arg[0] < 0.25) {
      // res[0] = cos(arg[0]*M_PI*2.)+2;
      res[0] = -8.*(arg[0]*arg[0]-0.25)+1.;
    }
    else {
      res[0] = 1.0;
    }
    res[1] = 1.5;
    res[2] = 0.;
    res[3] = 10.;
    res[1] *= res[0];
    res[2] *= res[0];
    res[3] += 0.5*res[1]*res[1]/res[0];
  }
  double gamma;
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
    DomainType c(1.25);
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
public:
  U0Sod() : gamma(1.4) {}
  template <class DomainType, class RangeType>
  void evaluate(const DomainType& arg, RangeType& res) const {
    evaluate(0,arg,res);
  }
  template <class DomainType, class RangeType>
  void evaluate(double t,const DomainType& arg, RangeType& res) const {
    if (arg[0]<0.) {
      res[0]=1.;
      res[1]=0.;
      res[1]=0.;
      res[2]=0.;
      res[3]=1./(1.4-1.0);
    } else {
      res[0]=0.125;
      res[1]=0.;
      res[2]=0.;
      res[3]=0.1/(1.4-1.0);
    }
    res[1] *= res[0];
    res[2] *= res[0];
    res[3] += 0.5*res[1]*res[1]/res[0];
  }
  double gamma;
};
