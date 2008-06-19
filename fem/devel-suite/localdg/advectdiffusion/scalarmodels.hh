#ifndef DUNE_SCALARMODELS_HH
#define DUNE_SCALARMODELS_HH

#include "modeldefault.hh"

template <class GridPartType,class ProblemType>
class BurgersModel {
 public:
  enum { ConstantVelocity = false };
  typedef typename GridPartType :: GridType GridType;
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
    epsilon(problem.epsilon) 
  {
		bool diff_tstep;	
		Parameter::get("fem.localdg.diffusion_timestep",diff_tstep);
		tstep_eps = (diff_tstep)? epsilon:0;
	}
  inline  void analyticalFlux(typename Traits::EntityType& en,
            double time,  
            const typename Traits::DomainType& x,
            const RangeType& u, 
            FluxRangeType& f) const {
    f = 0;
    f[0] = u*u*0.5;
  }
  inline  void diffusion(typename Traits::EntityType& en,
             double time, 
             const DomainType& x,
             const RangeType& u, 
             DiffusionRangeType& a) const {
    a=0;
    for (int i=0;i<dimDomain;i++)
      a[i][i]=u;
  }
  inline double diffusion(typename Traits::EntityType& en,
        double time, 
        const DomainType& x,
        const RangeType& u, 
        const GradientType& v,
        FluxRangeType& A) const {
    A = 0;
    A[0] = epsilon*v[0];
    return tstep_eps;
  }
  
  // needed for limitation
  inline  void velocity(
             const typename Traits::EntityType& en,
             double time, 
             const DomainType& x,
             const RangeType& u, 
             DomainType& velocity) const 
  {
    velocity = 0;
    velocity[0] = u[0];
    // problem_.velocity(en.geometry().global(x),velocity);
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
// forward declaration of UpwindFlux 
template <class Model>
class UpwindFlux;


template <class GridPartType,class ProblemType>
class AdvectionDiffusionModel {
 public:
  enum { ConstantVelocity = ProblemType :: ConstantVelocity };
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
    velocity_(0), 
    epsilon(problem.epsilon) 
  {
		bool diff_tstep;	
		Parameter::get("fem.localdg.diffusion_timestep",diff_tstep);
		tstep_eps = (diff_tstep)? epsilon:0;

    if(ConstantVelocity) 
    {
      // get constant velocity 
      problem_.velocity(velocity_,velocity_);
    } 
  }
  inline  void analyticalFlux(typename Traits::EntityType& en,
            double time,  
            const typename Traits::DomainType& x,
            const RangeType& u, 
            FluxRangeType& f) const 
  {
    // evaluate velocity 
    problem_.velocity(en.geometry().global(x),f[0]);
    // multiply with u
    f *= u;
  }

  inline  void velocity(
             const typename Traits::EntityType& en,
             double time, 
             const DomainType& x,
             const RangeType& u, 
             DomainType& v) const 
  {
    problem_.velocity(en.geometry().global(x),v);
  }
  
  inline  void diffusion(typename Traits::EntityType& en,
             double time, 
             const DomainType& x,
             const RangeType& u, 
             DiffusionRangeType& a) const 
  {
    a = 0;
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
             const typename Traits::FaceDomainType& x) const 
  {
    return true;
    /*
    const DomainType normal = it.integrationOuterNormal(x);
    int boundaryId = (std::abs(normal[0])>1e-10)? 1:2;
    return (boundaryId==1);
    */
  }

  inline double boundaryFlux(typename Traits::IntersectionIterator& it,
           double time, 
           const typename Traits::FaceDomainType& x,
           const RangeType& uLeft, 
           const GradientType& vLeft, 
           RangeType& gLeft) const  
  {
    gLeft = 0.;
    return 0.;
  }
  inline double boundaryFlux(typename Traits::IntersectionIterator& it,
           double time, 
           const typename Traits::FaceDomainType& x,
           const RangeType& uLeft, 
           RangeType& gLeft) const  
  {
    gLeft = 0.;
    return 0.;
  }
  inline  void boundaryValue(typename Traits::IntersectionIterator& it,
           double time, 
           const typename Traits::FaceDomainType& x,
           const RangeType& uLeft, 
           RangeType& uRight) const 
  {
    // Dirichlet
    DomainType xgl=it->intersectionGlobal().global(x);
    problem_.evaluate(time,xgl,uRight);
  }
  inline double diffusionTimeStep() const {
      return 2.*tstep_eps;
  }
  inline void maxSpeed(const typename Traits::DomainType& normal,
           double time,  
           const typename Traits::DomainType& x,
           const RangeType& u,
           double& advspeed,double& totalspeed) const 
  {
    problem_.velocity(x,velocity_);
    advspeed=std::abs(normal*velocity_);
    totalspeed=advspeed; // +tstep_eps;
  }
 protected:
  const ProblemType& problem_;
 public: 
  mutable DomainType velocity_;
 protected:
  double epsilon;
  double tstep_eps;
  friend class UpwindFlux<AdvectionDiffusionModel<GridPartType,ProblemType> >;
};

// Numerical Upwind-Flux 
template <class GridPartType,class ProblemType>
class UpwindFlux<AdvectionDiffusionModel<GridPartType,ProblemType> > {

 public:
  typedef AdvectionDiffusionModel<GridPartType,ProblemType> Model;
  typedef typename Model::Traits Traits;
  enum { dimRange = Model::dimRange };
  typedef typename Model::DomainType DomainType;
  typedef typename Model::RangeType RangeType;
  typedef typename Model::FluxRangeType FluxRangeType;
  typedef typename Model::DiffusionRangeType DiffusionRangeType;
 protected: 
  template <class Model, bool constVelo>
  struct Velocity
  {
    static inline double upwind(const Model& model,
                                typename Traits::IntersectionIterator& it,
                                double time, 
                                const typename Traits::FaceDomainType& x,
                                const RangeType& uLeft)
    {
      const typename Traits::DomainType normal = it.integrationOuterNormal(x);    
      DomainType velocity;
      model.velocity(*it.inside(),time,
                     it.intersectionSelfLocal().global(x),
                     uLeft,velocity);
      return normal*velocity;
    }
  };
  
  template <class Model>
  struct Velocity<Model,true>
  {
    static inline double upwind(const Model& model,
                                typename Traits::IntersectionIterator& it,
                                double time, 
                                const typename Traits::FaceDomainType& x,
                                const RangeType& uLeft)
    {
      const typename Traits::DomainType normal = it.integrationOuterNormal(x);    
      return normal * model.velocity_;
    }
  };
  
 public:
  UpwindFlux(const Model& mod) : model_(mod) {}
  const Model& model() const {return model_;}
  inline  double numericalFlux(typename Traits::IntersectionIterator& it,
             double time, 
             const typename Traits::FaceDomainType& x,
             const RangeType& uLeft, 
             const RangeType& uRight,
             RangeType& gLeft,
             RangeType& gRight) const 
  {
    const double upwind = Velocity<Model,Model::ConstantVelocity>::upwind(model_,it,time,x,uLeft);
    if (upwind>0)
      gLeft = uLeft;
    else
      gLeft = uRight;
    gLeft *= upwind;
    gRight = gLeft;
    return std::abs(upwind); // +model_.tstep_eps;
  }
 private:
  mutable DomainType velocity_;
  const Model& model_;
};
#endif
