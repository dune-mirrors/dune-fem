#include "modeldefault.hh"

template <class GridType>
class BurgersModel {
 public:
  enum { dimDomain = GridType::dimensionworld };
  enum { dimRange = 1};
  typedef ModelTraits<GridType,dimRange> Traits;
  typedef typename Traits::RangeType RangeType;
  typedef typename Traits::DomainType DomainType;
  typedef typename Traits::FluxRangeType FluxRangeType;
  typedef typename Traits::GradientType GradientType;
  typedef typename Traits::DiffusionRangeType DiffusionRangeType;
 public:
  BurgersModel(double eps,bool difftimestep=true) : 
    epsilon(eps), tstep_eps((difftimestep)?eps:0) {}
  inline  void analyticalFlux(typename Traits::EntityType& en,
				    double time,  
				    const typename Traits::DomainType& x,
				    const RangeType& u, 
				    FluxRangeType& f) const {
    f *= 0;
    f[0] = u*u*0.5;
  }
  inline  void diffusion(typename Traits::EntityType& en,
			       double time, 
			       const DomainType& x,
			       const RangeType& u, 
			       DiffusionRangeType& a) const {
    a*=0;
    for (int i=0;i<dimDomain;i++)
      a[i][i]=u;
  }
  inline double diffusion(typename Traits::EntityType& en,
			  double time, 
			  const DomainType& x,
			  const RangeType& u, 
			  const GradientType& v,
			  FluxRangeType& A) const {
    A *= 0;
    A[0] = epsilon*v[0];
    return tstep_eps;
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
			     const GradientType& vLeft, 
			     RangeType& gLeft) const  {
    return 0.;
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
    advspeed=std::abs(normal[0]*u);
    totalspeed=advspeed+tstep_eps;
  }
 protected:
  double epsilon;
  double tstep_eps;
};
// ***********************
template <class Model>
class UpwindFlux;
template <class GridType>
class AdvectionDiffusionModel {
 public:
  enum { dimDomain = GridType::dimensionworld };
  enum { dimRange = 1};
  typedef ModelTraits<GridType,dimRange,dimRange*dimDomain> Traits;
  typedef typename Traits::DomainType DomainType;
  typedef typename Traits::RangeType RangeType;
  typedef typename Traits::GradientType GradientType;
  typedef typename Traits::FluxRangeType FluxRangeType;
  typedef typename Traits::DiffusionRangeType DiffusionRangeType;
 public:
  AdvectionDiffusionModel(DomainType& velo,double eps,bool diff_timestep=true) :
    velocity(velo), epsilon(eps), tstep_eps((diff_timestep)?eps:0) {}
  inline  void analyticalFlux(typename Traits::EntityType& en,
			      double time,  
			      const typename Traits::DomainType& x,
				    const RangeType& u, 
				    FluxRangeType& f) const {
    f[0] = velocity;
    f *= u;
  }
  inline  void diffusion(typename Traits::EntityType& en,
			       double time, 
			       const DomainType& x,
			       const RangeType& u, 
			       DiffusionRangeType& a) const {
    a*=0;
    for (int i=0;i<dimDomain;i++)
      a[i][i]=u;
  }
  inline double diffusion(typename Traits::EntityType& en,
			 double time, 
			 const DomainType& x,
			 const RangeType& u, 
			 const GradientType& v,
			 FluxRangeType& A) const {
    
    A[0] = v;
    A *= epsilon;
    return tstep_eps;
  }
  inline bool hasBoundaryValue(typename Traits::IntersectionIterator& it,
			       double time, 
			       const typename Traits::FaceDomainType& x) const {
    const DomainType normal = it.integrationOuterNormal(x);
    int boundaryId = (std::abs(normal[0])>1e-10)? 1:2;
    return (boundaryId==1);
  }
  inline double boundaryFlux(typename Traits::IntersectionIterator& it,
			     double time, 
			     const typename Traits::FaceDomainType& x,
			     const RangeType& uLeft, 
			     const GradientType& vLeft, 
			     RangeType& gLeft) const  {
    gLeft*=0.;
    return 0.;
  }
  inline double boundaryFlux(typename Traits::IntersectionIterator& it,
			     double time, 
			     const typename Traits::FaceDomainType& x,
			     const RangeType& uLeft, 
			     RangeType& gLeft) const  {
    gLeft*=0.;
    return 0.;
  }
  inline  void boundaryValue(typename Traits::IntersectionIterator& it,
			     double time, 
			     const typename Traits::FaceDomainType& x,
			     const RangeType& uLeft, 
			     RangeType& uRight) const {
    
    uRight*=0; //uLeft;
  }
  inline  double maxSpeed(const typename Traits::DomainType& normal,
			  double time,  
			  const typename Traits::DomainType& x,
			  const RangeType& u,
			  double& advspeed,double& totalspeed) const {
    advspeed=std::abs(normal*velocity);
    totalspeed=advspeed+tstep_eps;
  }
 protected:
  DomainType velocity;
  double epsilon;
  double tstep_eps;
  friend class UpwindFlux<AdvectionDiffusionModel<GridType> >;
};
// Numerical Upwind-Flux
template <class GridType>
class UpwindFlux<AdvectionDiffusionModel<GridType> > {
 public:
  typedef AdvectionDiffusionModel<GridType> Model;
  typedef typename Model::Traits Traits;
  enum { dimRange = Model::dimRange };
  typedef typename Model::RangeType RangeType;
  typedef typename Model::FluxRangeType FluxRangeType;
  typedef typename Model::DiffusionRangeType DiffusionRangeType;
 public:
  UpwindFlux(const Model& mod) : model_(mod) {}
  const Model& model() const {return model_;}
  inline  double numericalFlux(typename Traits::IntersectionIterator& it,
				     double time, 
				     const typename Traits::FaceDomainType& x,
				     const RangeType& uLeft, 
				     const RangeType& uRight,
				     RangeType& gLeft,
				     RangeType& gRight) const {
    const typename Traits::DomainType normal = it.integrationOuterNormal(x);    
    double upwind = normal*model_.velocity;
    if (upwind>0)
      gLeft = uLeft;
    else
      gLeft = uRight;
    gLeft *= upwind;
    gRight = gLeft;
    return std::abs(upwind)+model_.tstep_eps;
  }
 private:
  const Model& model_;
};
