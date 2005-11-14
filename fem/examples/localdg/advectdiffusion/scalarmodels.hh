/* TYPES:
   DomainType: vector in world-space FV<dimDomain>
   RangeType: vector in phase-space FV<dimRange>
   FluxRangeType: matrix for analytical flux FM<Range,Domain>
   GradientType: vector in gradient phase-space FV<Range*Domain>
   DiffusionRangeType: matrix for diffusion flux FM<Range*Domain,Domain>
*/
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
  static double epsilon;
  inline static void analyticalFlux(typename Traits::EntityType& en,
				    double time,  
				    const typename Traits::DomainType& x,
				    const RangeType& u, 
				    FluxRangeType& f) {
    f[0] = u*u*0.5;
  }
  inline static void diffusion(typename Traits::EntityType& en,
			       double time, 
			       const DomainType& x,
			       const RangeType& u, 
			       DiffusionRangeType& a) {
    a*=0;
    for (int i=0;i<dimDomain;i++)
      a[i][i]=u;
  }
  inline static void diffusion(typename Traits::EntityType& en,
			       double time, 
			       const DomainType& x,
			       const RangeType& u, 
			       const GradientType& v,
			       FluxRangeType& A) {
    A[0] = v;
    A *= epsilon;
  }
  inline static void boundaryValue(typename Traits::IntersectionIterator& it,
				   double time, 
				   const typename Traits::FaceDomainType& x,
				   const RangeType& uLeft, 
				   RangeType& uRight) {
    
    uRight=uLeft;
  }
  inline static double maxSpeed(const typename Traits::DomainType& normal,
				double time,  
				const typename Traits::DomainType& x,
				const RangeType& u) {
    return abs(normal[0]*u);
  }
};
// ***********************
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
  static DomainType velocity;
  static double epsilon;
  inline static void analyticalFlux(typename Traits::EntityType& en,
				    double time,  
				    const typename Traits::DomainType& x,
				    const RangeType& u, 
				    FluxRangeType& f) {
    f[0] = velocity;
    f *= u;
  }
  inline static void diffusion(typename Traits::EntityType& en,
			       double time, 
			       const DomainType& x,
			       const RangeType& u, 
			       DiffusionRangeType& a) {
    a*=0;
    for (int i=0;i<dimDomain;i++)
      a[i][i]=u;
  }
  inline static void diffusion(typename Traits::EntityType& en,
			       double time, 
			       const DomainType& x,
			       const RangeType& u, 
			       const GradientType& v,
			       FluxRangeType& A) {
    
    A[0] = v;
    A *= epsilon;
  }
  inline static void boundaryValue(typename Traits::IntersectionIterator& it,
				   double time, 
				   const typename Traits::FaceDomainType& x,
				   const RangeType& uLeft, 
				   RangeType& uRight) {
    
    uRight=uLeft;
  }
  inline static double maxSpeed(const typename Traits::DomainType& normal,
				double time,  
				const typename Traits::DomainType& x,
				const RangeType& u) {
    return abs(normal*velocity);
  }
};
// Numerical Upwind-Flux
template <class Model>
class UpwindFlux;
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
  inline static double numericalFlux(typename Traits::IntersectionIterator& it,
				     double time, 
				     const typename Traits::FaceDomainType& x,
				     const RangeType& uLeft, 
				     const RangeType& uRight,
				     RangeType& gLeft,
				     RangeType& gRight) {
    const typename Traits::DomainType normal = it.integrationOuterNormal(x);    
    double upwind = normal*Model::velocity;
    if (upwind>0)
      gLeft = uLeft;
    else
      gLeft = uRight;
    gLeft *= upwind;
    gRight = gLeft;
    return std::abs(upwind);
  }
  inline static double boundaryFlux(typename Traits::IntersectionIterator& it,
				    double time, 
				    const typename Traits::FaceDomainType& x,
				    const RangeType& uLeft, 
				    RangeType& gLeft) {
    RangeType uRight,gRight;
    Model::boundaryValue(it,time,x,uLeft,uRight);
    return numericalFlux(it,time,x,uLeft,uRight,gLeft,gRight);
  }
  inline static void boundaryDiffusionFlux
                    (typename Traits::IntersectionIterator& it,
		     double time, 
		     const typename Traits::FaceDomainType& x,
		     const RangeType& uLeft, 
		     DiffusionRangeType& aLeft) {
    RangeType uRight;
    Model::boundaryValue(it,time,x,uLeft,uRight);
    Model::diffusion(*it.inside(),time,
		     it.intersectionSelfLocal().global(x),
		     uRight,aLeft);
  }
};
