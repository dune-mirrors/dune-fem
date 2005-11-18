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
    inline double boundaryFlux(typename Traits::IntersectionIterator& it,
				      double time, 
				      const typename Traits::FaceDomainType& x,
				      const RangeType& uLeft, 
				      RangeType& gLeft) const {
      RangeType uRight,gRight;
      model_.boundaryValue(it,time,x,uLeft,uRight);
      return numericalFlux(it,time,x,uLeft,uRight,gLeft,gRight);
    }
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
