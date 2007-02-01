#ifndef DUNE_SCALARMODELS_HH
#define DUNE_SCALARMODELS_HH

#include "modeldefault.hh"

template <class GridPartType,class ProblemType>
class BurgersModel {
 public:
  typedef typename GridPartType	:: GridType GridType;
  enum { dimDomain = GridType::dimensionworld };
  enum { dimRange = 1};
  typedef ModelTraits<GridPartType,dimRange> Traits;
  typedef typename Traits::RangeType RangeType;
  typedef typename Traits::DomainType DomainType;
  typedef typename Traits::FluxRangeType FluxRangeType;
  typedef typename Traits::GradientType GradientType;
  typedef typename Traits::DiffusionRangeType DiffusionRangeType;
 public:
  BurgersModel(GridType& grid,
			  const ProblemType& problem) :
    problem_(problem),
    epsilon(problem.epsilon), 
    tstep_eps((problem.diff_tstep)?problem.epsilon:0) {}
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
  inline double diffusionTimeStep() const {
    return 2.*tstep_eps;
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
    // Outflow
    uRight=uLeft;
    // Dirichlet
    // DomainType xgl=it.intersectionGlobal().global(x);
    // problem_.evaluate(time,xgl,uRight);
  }
  inline void maxSpeed(const typename Traits::DomainType& normal,
		       double time,  
		       const typename Traits::DomainType& x,
		       const RangeType& u,
		       double& advspeed,double& totalspeed) const {
    advspeed=std::abs(normal[0]*u);
    totalspeed=advspeed; // +tstep_eps;
  }
 protected:
  const ProblemType& problem_;
  double epsilon;
  double tstep_eps;
};
// ***********************
template <class Model>
class UpwindFlux;
template <class GridPartType,class ProblemType>
class AdvectionDiffusionModel {
 public:
  typedef typename GridPartType :: GridType GridType;	 
  enum { dimDomain = GridType::dimensionworld };
  enum { dimRange = 1};
  typedef ModelTraits<GridPartType,dimRange,dimRange*dimDomain> Traits;
  typedef typename Traits::DomainType DomainType;
  typedef typename Traits::RangeType RangeType;
  typedef typename Traits::GradientType GradientType;
  typedef typename Traits::FluxRangeType FluxRangeType;
  typedef typename Traits::DiffusionRangeType DiffusionRangeType;
 public:
  AdvectionDiffusionModel(GridType& grid,
			  const ProblemType& problem) :
    problem_(problem),
    velocity(problem.velocity), epsilon(problem.epsilon), 
    tstep_eps((problem.diff_tstep)?problem.epsilon:0) {}
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
    return true;
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
    
    // Dirichlet
    DomainType xgl=it.intersectionGlobal().global(x);
    problem_.evaluate(time,xgl,uRight);
  }
	inline double diffusionTimeStep() const {
		  return 2.*tstep_eps;
	}
  inline void maxSpeed(const typename Traits::DomainType& normal,
		       double time,  
		       const typename Traits::DomainType& x,
		       const RangeType& u,
		       double& advspeed,double& totalspeed) const {
    advspeed=std::abs(normal*velocity);
    totalspeed=advspeed; // +tstep_eps;
  }
 protected:
  const ProblemType& problem_;
  DomainType velocity;
  double epsilon;
  double tstep_eps;
  friend class UpwindFlux<AdvectionDiffusionModel<GridType,ProblemType> >;
};
// Numerical Upwind-Flux
template <class GridPartType,class ProblemType>
class UpwindFlux<AdvectionDiffusionModel<GridPartType,ProblemType> > {
 public:
  typedef AdvectionDiffusionModel<GridPartType,ProblemType> Model;
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
    return std::abs(upwind); // +model_.tstep_eps;
  }
 private:
  const Model& model_;
};

/*********************************************************************
template <class GridType>
class U0 {
public:
  enum { dimDomain = GridType::dimensionworld };  
  typedef FieldVector<double,dimDomain> DomainType;
  typedef FieldVector<double,1> RangeType;
  U0(double eps,bool diff_timestep=true) :
    velocity(0), epsilon(eps), diff_tstep(diff_timestep) {
      velocity[0]=0.8;
    }
  void evaluate(const DomainType& arg, RangeType& res) const {
    evaluate(0,arg,res);
  }
  void evaluate(double t,const DomainType& arg, RangeType& res) const {
    if (arg[0]*arg[0] < 0.25) {
      res = cos(arg[0]*M_PI*2.)+1;
    }
    else {
      res = 0.0;
    }
    double diffusion_ = 0.01;
    double t0 = 1.0;
    // res = 1./sqrt(4.*M_PI*diffusion_*t0)*
    //       exp(-(arg[0]+0.8-0.8*t0)*(arg[0]+0.8-0.8*t0)/(4.*diffusion_*t0));
    // res = 1.0;
  }
  DomainType velocity;
  double epsilon;
  bool diff_tstep;
};
*/
#endif
