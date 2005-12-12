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
template <class GridType>
class EulerModel {
 public:
  enum { dimDomain = GridType::dimensionworld };
  enum { dimRange = dimDomain+2};
  typedef ModelTraits<GridType,dimRange> Traits;
  typedef typename Traits::RangeType RangeType;
  typedef typename Traits::DomainType DomainType;
  typedef typename Traits::FluxRangeType FluxRangeType;
  //typedef typename Traits::GradientType GradientType;
  //typedef typename Traits::DiffusionRangeType DiffusionRangeType;
 public:
  EulerModel(double gamma,bool difftimestep=true) : gamma_(gamma)
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
    return 0.;
  }
  inline  void boundaryValue(typename Traits::IntersectionIterator& it,
			     double time, 
			     const typename Traits::FaceDomainType& x,
			     const RangeType& uLeft, 
			     RangeType& uRight) const {
    uRight=uLeft;
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
  friend class DWNumFlux<EulerModel<GridType> >;
};
// ***********************
template <class GridType>
class DWNumFlux<EulerModel<GridType> > {
 public:
  enum { dimDomain = GridType::dimensionworld };
  typedef Mhd::MhdSolver SolverType;
  typedef EulerModel<GridType> Model;
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
  inline  double numericalFlux(typename Traits::IntersectionIterator& it,
                                     double time,
                                     const typename Traits::FaceDomainType& x,
                                     const RangeType& uLeft,
                                     const RangeType& uRight,
                                     RangeType& gLeft,
                                     RangeType& gRight) const {
    typename Traits::DomainType normal = it.integrationOuterNormal(x);
    double len = normal.two_norm();
    normal /= len;
    if (std::fabs(normal.two_norm() - 1.0) > 1.0e-8)
      cerr << normal << " " << len << endl;

    RangeType ul,ur,ret;
    rot_.rotateForth(uLeft,  ul, normal);
    rot_.rotateForth(uRight, ur, normal);
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
    ret[0] = retmhd[0];
    for (int i=0;i<dimDomain;++i)
      ret[i+1] = retmhd[i+1];
    ret[1+dimDomain] = retmhd[7];
    rot_.rotateBack(ret,gLeft,normal);
    gLeft *= len;
    gRight = gLeft;
    return ldt*len;
  }  
  const Model& model() const {return model_;}
 private:
  Adi::FieldRotator<Model> rot_;
  const Model& model_;
  const typename SolverType::Eosmode::meos_t eos;
  mutable SolverType numFlux_;
};

// ***********************
template <>
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
void EulerFlux<2>::analyticalFlux(const double gamma,
				  const FieldVector<double,2+2>& u,
				  FieldMatrix<double,2+2,2>& f) const {
  assert(u[0]>1e-10);
  double rhoeps = u[e]-0.5*(u[1]*u[1]+u[2]*u[2])/u[0];
  if (rhoeps<1e-10) 
      cerr << u << " " << rhoeps << endl;
  assert(rhoeps>1e-10);
  double p = (gamma-1)*rhoeps;
  f[0][0] = u[1];                 f[0][1] = u[2];
  f[1][0] = u[1]/u[0]*u[1]+p;     f[1][1] = u[2]/u[0]*u[1];
  f[2][0] = u[1]/u[0]*u[2];       f[2][1] = u[2]/u[0]*u[2]+p;
  f[e][0] = u[1]/u[0]*(u[e]+p);   f[e][1] = u[2]/u[0]*(u[e]+p);
}
template <>
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
double EulerFlux<1>::maxSpeed(const double gamma,
			      const FieldVector<double,1>& n,
			      const FieldVector<double,1+2>& u) const {
  assert(u[0]>1e-10);
  double u_normal = u[1]*n[0] / u[0];
  double rhoeps = u[e]-0.5*u[1]*u[1]/u[0];
  assert(rhoeps>1e-10);
  double p = (gamma-1)*rhoeps;
  double c2 = gamma*p/u[0];
  assert(c2>1e-10);
  return fabs(u_normal) + sqrt(c2);
}
template <>
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
double EulerFlux<3>::maxSpeed(const double gamma,
			      const FieldVector<double,3>& n,
			      const FieldVector<double,3+2>& u) const {
  assert(u[0]>1e-10);
  double u_normal = (u[1]*n[0]+u[2]*n[1]+u[3]*n[2]) / u[0];
  double rhoeps = u[e]-0.5*(u[1]*u[1]+u[2]*u[2]+u[3]*u[3])/u[0];
  assert(rhoeps>1e-10);
  double p = (gamma-1)*rhoeps;
  double c2 = gamma*p/u[0];
  assert(c2>1e-10);
  return fabs(u_normal) + sqrt(c2);
}

