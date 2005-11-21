/**************************************************************
 **** Structure of model class:
 **************************************************************
 **** Traits class: parameters are dimRange and dimRange1
   DomainType: vector in world-space FV<dimDomain>
   RangeType: vector in phase-space FV<dimRange>
   FluxRangeType: matrix for analytical flux FM<Range,Domain>
   GradientType: vector in gradient phase-space FV<dimRange1>
   DiffusionRangeType: matrix for diffusion flux FM<Range1,Domain>
 **** Description class for    V + div a(U)          = 0
                            dt U + div (F(U)+A(U,V)) = 0
      template <class GridType>
      class Model {
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
        Model(...) {}
        // F
        inline void analyticalFlux(typename Traits::EntityType& en,
				   double time,  
				   const typename Traits::DomainType& x,
				   const RangeType& u, 
				   FluxRangeType& f) const {}
        // a
        inline void diffusion(typename Traits::EntityType& en,
		              double time, 
			      const DomainType& x,
			      const RangeType& u, 
			      DiffusionRangeType& a) const {}
        // A
        inline double diffusion(typename Traits::EntityType& en,
			        double time, 
			        const DomainType& x,
			        const RangeType& u, 
			        const GradientType& v,
			        FluxRangeType& A) const {}
        // for boundary treatment: !!!! has to be implemented !!!!
        // return true:  numerical flux is used with uRight given by boundaryValue
        // return false: boundaryFlux is used for flux computation
        inline bool hasBoundaryValue(int boundaryId) {}
        inline void boundaryValue(typename Traits::IntersectionIterator& it,
				  double time, 
			          const typename Traits::FaceDomainType& x,
			          const RangeType& uLeft, 
			          RangeType& uRight) const {}
        inline void boundaryFlux(typename Traits::IntersectionIterator& it,
				  double time, 
			          const typename Traits::FaceDomainType& x,
			          const RangeType& uLeft, 
                                  RangeType& gLeft) {} 
        // only neaded for LLFFlux for artificial diffusion
        inline double maxSpeed(const typename Traits::DomainType& normal,
       			       double time,  
	       		       const typename Traits::DomainType& x,
		       	       const RangeType& u) const {}
      };

/**************************************************************/
  // Model-Traits
  template <class Grid,int dimRange2,
	    int dimRange1=dimRange2*Grid::dimensionworld>
  class ModelTraits {
  public:
    typedef Grid GridType;
    enum { dimDomain = GridType::dimensionworld };
    enum { dimRange = dimRange2, dimGradRange = dimRange1 };
    typedef FieldVector<double, dimDomain> DomainType;
    typedef FieldVector<double, dimDomain-1> FaceDomainType;
    typedef FieldVector<double,dimRange> RangeType;
    typedef FieldVector<double,dimGradRange> GradientType;
    typedef FieldMatrix<double,dimRange,dimDomain> FluxRangeType;
    typedef FieldMatrix<double,dimGradRange,dimDomain> DiffusionRangeType;
    typedef typename GridType::Traits::IntersectionIterator IntersectionIterator;
    typedef typename GridType::template Codim<0>::Entity EntityType;
   };
  // The Local Lax-Friedrichs Flux Function
  template <class Model>
  class LLFFlux {
  public:
    typedef typename Model::Traits Traits;
    enum { dimRange = Model::dimRange };
    typedef typename Model::RangeType RangeType;
    typedef typename Model::FluxRangeType FluxRangeType;
    typedef typename Model::DiffusionRangeType DiffusionRangeType;
  public:
    LLFFlux(const Model& mod) : model_(mod) {}
    const Model& model() const {return model_;}
    // Return value: maximum wavespeed*length of integrationOuterNormal
    // gLeft,gRight are fluxed * length of integrationOuterNormal
    inline double numericalFlux(typename Traits::IntersectionIterator& it,
				       double time, 
				       const typename Traits::FaceDomainType& x,
				       const RangeType& uLeft, 
				       const RangeType& uRight,
				       RangeType& gLeft,
				       RangeType& gRight) const {
      const typename Traits::DomainType normal = it.integrationOuterNormal(x);  
      typename Traits::RangeType visc;
      typename Traits::FluxRangeType anaflux;
      model_.analyticalFlux(*it.inside(), time,
			    it.intersectionSelfLocal().global(x),
			    uLeft, anaflux);
      gLeft*=0;
      anaflux.umv(normal, gLeft);
      if (it.neighbor())
	model_.analyticalFlux(*it.outside(), time,
			      it.intersectionNeighborLocal().global(x),
			      uRight, anaflux);
      else
	model_.analyticalFlux(*it.inside(), time,
			      it.intersectionSelfLocal().global(x),
			      uRight, anaflux);
      anaflux.umv(normal,gLeft);
      
      double maxspeedl=
	model_.maxSpeed(normal,time,it.intersectionSelfLocal().global(x),uLeft);
      double maxspeedr=
	model_.maxSpeed(normal,time,it.intersectionSelfLocal().global(x),uRight);
      double maxspeed=(maxspeedl>maxspeedr)?maxspeedl:maxspeedr;
      visc=uRight;
      visc-=uLeft;
      visc*=maxspeed;
      gLeft-=visc;
      
      gLeft*=0.5;
      gRight=gLeft;
      return maxspeed;
    }
    // !!!! Should be removed !!!!
    inline double boundaryFlux(typename Traits::IntersectionIterator& it,
				      double time, 
				      const typename Traits::FaceDomainType& x,
				      const RangeType& uLeft, 
				      RangeType& gLeft) const {
      RangeType uRight,gRight;
      model_.boundaryValue(it,time,x,uLeft,uRight);
      return numericalFlux(it,time,x,uLeft,uRight,gLeft,gRight);
    }
    // !!!! Should be removed (what to do if hasBoundaryValue=false?)
    inline void boundaryDiffusionFlux
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
