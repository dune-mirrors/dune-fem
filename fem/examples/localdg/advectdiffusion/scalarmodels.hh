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
  BurgersModel(double eps) : epsilon(eps) {}
  inline  void analyticalFlux(typename Traits::EntityType& en,
				    double time,  
				    const typename Traits::DomainType& x,
				    const RangeType& u, 
				    FluxRangeType& f) const {
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
  inline  void diffusion(typename Traits::EntityType& en,
			       double time, 
			       const DomainType& x,
			       const RangeType& u, 
			       const GradientType& v,
			       FluxRangeType& A) const {
    A[0] = v;
    A *= epsilon;
  }
  inline  void boundaryValue(typename Traits::IntersectionIterator& it,
				   double time, 
				   const typename Traits::FaceDomainType& x,
				   const RangeType& uLeft, 
				   RangeType& uRight) const {
    
    uRight=uLeft;
  }
  inline  double maxSpeed(const typename Traits::DomainType& normal,
				double time,  
				const typename Traits::DomainType& x,
				const RangeType& u) const {
    return abs(normal[0]*u);
  }
 protected:
  double epsilon;
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
  AdvectionDiffusionModel(DomainType& velo,double eps) :
    velocity(velo), epsilon(eps) {}
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
  inline  void diffusion(typename Traits::EntityType& en,
			       double time, 
			       const DomainType& x,
			       const RangeType& u, 
			       const GradientType& v,
			       FluxRangeType& A) const {
    
    A[0] = v;
    A *= epsilon;
  }
  inline  void boundaryValue(typename Traits::IntersectionIterator& it,
				   double time, 
				   const typename Traits::FaceDomainType& x,
				   const RangeType& uLeft, 
				   RangeType& uRight) const {
    
    uRight=uLeft;
  }
  inline  double maxSpeed(const typename Traits::DomainType& normal,
				double time,  
				const typename Traits::DomainType& x,
				const RangeType& u) const {
    return abs(normal*velocity);
  }
 protected:
  DomainType velocity;
  double epsilon;
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
    return std::abs(upwind);
  }
  inline  double boundaryFlux(typename Traits::IntersectionIterator& it,
				    double time, 
				    const typename Traits::FaceDomainType& x,
				    const RangeType& uLeft, 
				    RangeType& gLeft) const {
    RangeType uRight,gRight;
    model_.boundaryValue(it,time,x,uLeft,uRight);
    return numericalFlux(it,time,x,uLeft,uRight,gLeft,gRight);
  }
  inline  void boundaryDiffusionFlux
                    (typename Traits::IntersectionIterator& it,
		     double time, 
		     const typename Traits::FaceDomainType& x,
		     const RangeType& uLeft, 
		     DiffusionRangeType& aLeft) const {
    RangeType uRight;
    model_.boundaryValue(it,time,x,uLeft,uRight);
    model_.diffusion(*it.inside(),time,
		     it.intersectionSelfLocal().global(x),
		     uRight,aLeft);
  }
 private:
  const Model& model_;
};
