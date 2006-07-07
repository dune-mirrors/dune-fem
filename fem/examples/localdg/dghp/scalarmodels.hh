#include "modeldefault.hh"

template <class GridType,class ProblemType>
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
  BurgersModel(GridType& grid,
	       const ProblemType& problem) :
    problem_(problem),
    epsilon(problem.epsilon), 
    tstep_eps((problem.diff_tstep)?problem.epsilon:0) {}
  inline  void analyticalFlux(const typename Traits::EntityType& en,
			      double time,  
			      const typename Traits::DomainType& x,
			      const RangeType& u, 
			      FluxRangeType& f) const {
    f = 0;
    for (int i=0;i<dimDomain;++i) {
      f[0][i] = u*u*0.5;
    }
  }
  inline  void jacobian(const typename Traits::EntityType& en,
			double time,  
			const typename Traits::DomainType& x,
			const RangeType& u, 
			const FluxRangeType& du,
			RangeType& A) const {
    A = 0.;
    DomainType velocity(u[0]);
    du.umv(velocity,A);
  }
  inline  void diffusion(const typename Traits::EntityType& en,
			 double time, 
			 const DomainType& x,
			 const RangeType& u, 
			 DiffusionRangeType& a) const {
    a*=0;
    for (int i=0;i<dimDomain;i++)
      a[i][i]=u;
  }
  inline double diffusion(const typename Traits::EntityType& en,
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
  inline bool hasBoundaryValue(const typename Traits::IntersectionIterator& it,
			       double time, 
			       const typename Traits::FaceDomainType& x) const {
    return true;
  }
  inline double boundaryFlux(const typename Traits::IntersectionIterator& it,
			     double time, 
			     const typename Traits::FaceDomainType& x,
			     const RangeType& uLeft, 
			     const GradientType& vLeft, 
			     RangeType& gLeft) const  {
    return 0.;
  }
  inline double boundaryFlux(const typename Traits::IntersectionIterator& it,
			     double time, 
			     const typename Traits::FaceDomainType& x,
			     const RangeType& uLeft, 
			     RangeType& gLeft) const  {
    return 0.;
  }
  inline  void boundaryValue(const typename Traits::IntersectionIterator& it,
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
  inline const ProblemType& problem() const {
    return problem_;
  }
 protected:
  const ProblemType& problem_;
  double epsilon;
  double tstep_eps;
};

// ***********************
template <class Model>
class UpwindFlux;
template <class GridType,class ProblemType>
class BuckLevModel {
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
  BuckLevModel(GridType& grid,
			  const ProblemType& problem) :
    problem_(problem),
    epsilon(problem.epsilon), 
    tstep_eps((problem.diff_tstep)?problem.epsilon:0) {}
  inline  void analyticalFlux(const typename Traits::EntityType& en,
			      double time,  
			      const typename Traits::DomainType& x,
			      const RangeType& u, 
			      FluxRangeType& f) const {
    f = 0.;
    f[0][0] = problem_.f(u[0]);
  }
  inline  void jacobian(const typename Traits::EntityType& en,
			double time,  
			const typename Traits::DomainType& x,
			const RangeType& u, 
			const FluxRangeType& du,
			RangeType& A) const {
    A = 0.;
    DomainType velocity(0.0);
    velocity[0] = problem_.f1(u[0]);
    du.umv(velocity,A);
  }
  inline  void diffusion(const typename Traits::EntityType& en,
			 double time, 
			 const DomainType& x,
			 const RangeType& u, 
			 DiffusionRangeType& a) const {
    a*=0;
    for (int i=0;i<dimDomain;i++)
      a[i][i]=u;
  }
  inline double diffusion(const typename Traits::EntityType& en,
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
  inline bool hasBoundaryValue(const typename Traits::IntersectionIterator& it,
			       double time, 
			       const typename Traits::FaceDomainType& x) const {
    return true;
  }
  inline double boundaryFlux(const typename Traits::IntersectionIterator& it,
			     double time, 
			     const typename Traits::FaceDomainType& x,
			     const RangeType& uLeft, 
			     const GradientType& vLeft, 
			     RangeType& gLeft) const  {
    return 0.;
  }
  inline double boundaryFlux(const typename Traits::IntersectionIterator& it,
			     double time, 
			     const typename Traits::FaceDomainType& x,
			     const RangeType& uLeft, 
			     RangeType& gLeft) const  {
    return 0.;
  }
  inline  void boundaryValue(const typename Traits::IntersectionIterator& it,
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
    advspeed=std::abs(normal[0]*problem_.f1(u[0]));
    totalspeed=advspeed; // +tstep_eps;
  }
  inline const ProblemType& problem() const {
    return problem_;
  }
 protected:
  const ProblemType& problem_;
  double epsilon;
  double tstep_eps;
  friend class UpwindFlux<BuckLevModel<GridType,ProblemType> >;
};

// ***********************
template <class GridType,class ProblemType>
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
  AdvectionDiffusionModel(GridType& grid,
			  const ProblemType& problem) :
    problem_(problem),
    velocity(problem.velocity_), epsilon(problem.epsilon), 
    tstep_eps((problem.diff_tstep)?problem.epsilon:0) {}
  inline  void analyticalFlux(const typename Traits::EntityType& en,
			      double time,  
			      const typename Traits::DomainType& x,
			      const RangeType& u, 
			      FluxRangeType& f) const {
    DomainType vel;
    problem_.velocity(time,en.geometry().global(x),vel);
    f[0] = vel;
    f *= u;
  }
  inline  void jacobian(const typename Traits::EntityType& en,
			double time,  
			const typename Traits::DomainType& x,
			const RangeType& u, 
			const FluxRangeType& du,
			RangeType& A) const {
    A = 0.;
    DomainType vel;
    problem_.velocity(time,en.geometry().global(x),vel);
    du.umv(vel,A);
  }
  inline  void diffusion(const typename Traits::EntityType& en,
			 double time, 
			 const DomainType& x,
			 const RangeType& u, 
			 DiffusionRangeType& a) const {
    a*=0;
    for (int i=0;i<dimDomain;i++)
      a[i][i]=u;
  }
  inline double diffusion(const typename Traits::EntityType& en,
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
  inline const ProblemType& problem() const {
    return problem_;
  }
 protected:
  const ProblemType& problem_;
  DomainType velocity;
  double epsilon;
  double tstep_eps;
  friend class UpwindFlux<AdvectionDiffusionModel<GridType,ProblemType> >;
};
// Numerical Upwind-Flux
template <class GridType,class ProblemType>
class UpwindFlux<AdvectionDiffusionModel<GridType,ProblemType> > {
 public:
  typedef AdvectionDiffusionModel<GridType,ProblemType> Model;
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
    double upwind;
    typename Model::DomainType vel;
    model_.problem_.velocity(time,it.intersectionGlobal().global(x),vel);

    upwind = normal*vel;
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

// Numerical Upwind-Flux
template <class GridType,class ProblemType>
class UpwindFlux<BuckLevModel<GridType,ProblemType> > {
 public:
  typedef BuckLevModel<GridType,ProblemType> Model;
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
    double upwind;
    typename Model::DomainType vel;
    model_.problem_.velocity(time,it.intersectionGlobal().global(x),vel);

    upwind = normal*vel;
    if (upwind>0){
      gLeft[0] = model_.problem_.f(uLeft[0]);
      gLeft *= upwind;
      upwind *= model_.problem_.f1(uLeft[0]);
    }
    else{
      gLeft[0] = model_.problem_.f(uRight[0]);
      gLeft *= upwind;
      upwind *= model_.problem_.f1(uRight[0]);
    }
    gRight = gLeft;
    double wave = std::abs(upwind);
    double wavem = std::abs(normal*vel)
      *model_.problem_.f1(0.5*(uLeft[0]+uRight[0]));
    wave = (wave >wavem)?wave:wavem;
    return wave; // +model_.tstep_eps;
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
