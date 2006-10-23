#ifndef DUNE_EXAMPLEDISCRETEMODELS_HH
#define DUNE_EXAMPLEDISCRETEMODELS_HH

#include <dune/fem/pass/dgpass.hh>
#include <dune/fem/pass/discretemodel.hh>
#include <dune/fem/pass/selection.hh>
#include <dune/fem/misc/timeutility.hh>

// Dune includes
#include <dune/common/utility.hh>
#include <dune/grid/common/gridpart.hh>

//- local includes 
#include <dune/fem/space/dgspace.hh>
#include <dune/fem/space/combinedspace.hh>

#include <dune/fem/discretefunction/dfadapt.hh>
#include <dune/fem/discretefunction/adaptivefunction.hh>
#include <dune/fem/quadrature/cachequad.hh>


#include <dune/fem/space/dgspace/dgleafindexset.hh>

//*************************************************************
namespace Dune {  
  template <class Model,class NumFlux,int polOrd >
  class TransportDiffusionDiscreteModel1;
  template <class Model,int polOrd >
  class LimiterDiscreteModel1;
  template <class Model,class NumFlux,int polOrd,
	    bool withDiffusion,bool withAdvection >
  class TransportDiffusionDiscreteModel2;

  // MethodOrderTraits
  template <class Model,int dimRange>
  class PassTraits {
  public:
    typedef typename Model::Traits ModelTraits;
    typedef typename ModelTraits::GridType GridType;
    enum { dimDomain = Model::Traits::dimDomain };
    // choose leaf level for iteration  
    typedef DGAdaptiveLeafGridPart<GridType> GridPartType;
    //typedef DGAdaptiveLeafIndexSet<GridType> DGIndexSetType;
    // typedef typename LeafGridPart<GridType>::IndexSetType DGIndexSetType;
    //typedef DefaultGridPart<GridType,DGIndexSetType> GridPartType;

    typedef FunctionSpace<double, double, dimDomain, dimRange> FunctionSpaceType; 
    typedef ElementQuadrature<GridType,1> FaceQuadratureType;
    typedef CachingQuadrature<GridType,0> VolumeQuadratureType;
    // typedef CachingQuadrature<GridType,1> FaceQuadratureType;
  };
  // DiscreteModelTraits
  template <class Model,class NumFlux,int polOrd >
  struct TransportDiffusionTraits1 
  {
    typedef typename Model::Traits ModelTraits;
    typedef typename ModelTraits::GridType GridType;

    enum { dimRange = ModelTraits::dimGradRange };
    enum { dimDomain = ModelTraits::dimDomain };

    typedef PassTraits<Model,dimRange> Traits;
    typedef typename Traits::FunctionSpaceType FunctionSpaceType;

    typedef typename ModelTraits::DomainType DomainType;
    typedef typename FunctionSpaceType::RangeType RangeType;
    typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;

    typedef typename Traits::VolumeQuadratureType VolumeQuadratureType;
    typedef typename Traits::FaceQuadratureType FaceQuadratureType;
    typedef typename Traits::GridPartType GridPartType;
    //typedef typename Traits::SingleFunctionSpaceType SingleFunctionSpaceType;
    //typedef typename Traits::FunctionSpaceType FunctionSpaceType;

    // typedef DiscontinuousGalerkinSpace<SingleFunctionSpaceType, GridPartType, polOrd> SingleSpaceType;
    // typedef CombinedSpace<SingleSpaceType, dimRange> DiscreteFunctionSpaceType;
    typedef DiscontinuousGalerkinSpace<FunctionSpaceType, GridPartType, polOrd,CachingStorage> DiscreteFunctionSpaceType;
    typedef DiscreteFunctionSpaceType SpaceType;
    //typedef AdaptiveDiscreteFunction<DiscreteFunctionSpaceType> DestinationType;
    typedef DFAdapt<DiscreteFunctionSpaceType> DestinationType;
    typedef DestinationType DiscreteFunctionType;

    typedef TransportDiffusionDiscreteModel1<Model,NumFlux,polOrd> DiscreteModelType;
  };
  template <class Model,class NumFlux,int polOrd,
	    bool withDiffusion,bool withAdvection>
  struct TransportDiffusionTraits2
  {
    typedef typename Model::Traits ModelTraits;
    typedef typename ModelTraits::GridType GridType;

    enum { dimRange = ModelTraits::dimRange };
    enum { dimDomain = ModelTraits::dimDomain };

    typedef PassTraits<Model,dimRange> Traits;
    typedef typename Traits::FunctionSpaceType FunctionSpaceType;
    typedef typename FunctionSpaceType::DomainType DomainType;
    typedef typename FunctionSpaceType::RangeType RangeType;
    typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;

    typedef typename Traits::VolumeQuadratureType VolumeQuadratureType;
    typedef typename Traits::FaceQuadratureType FaceQuadratureType;
    typedef typename Traits::GridPartType GridPartType;
    typedef DiscontinuousGalerkinSpace<FunctionSpaceType, GridPartType, polOrd,CachingStorage> DiscreteFunctionSpaceType;
    typedef DiscreteFunctionSpaceType SpaceType;
    //typedef AdaptiveDiscreteFunction<DiscreteFunctionSpaceType> DestinationType;
    typedef DFAdapt<DiscreteFunctionSpaceType> DestinationType;
    typedef DestinationType DiscreteFunctionType;

    typedef TransportDiffusionDiscreteModel2<Model,NumFlux,polOrd,withDiffusion,withAdvection> DiscreteModelType;
  };
  template <class Model,int polOrd >
  struct LimiterTraits1 
  {
    typedef typename Model::Traits ModelTraits;
    typedef typename ModelTraits::GridType GridType;

    enum { dimRange = ModelTraits::dimRange };
    enum { dimDomain = ModelTraits::dimDomain };

    typedef PassTraits<Model,dimRange> Traits;
    typedef typename Traits::FunctionSpaceType FunctionSpaceType;
    typedef typename ModelTraits::DomainType DomainType;
    typedef typename FunctionSpaceType::RangeType RangeType;
    typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;

    typedef typename Traits::VolumeQuadratureType VolumeQuadratureType;
    typedef typename Traits::FaceQuadratureType FaceQuadratureType;
    // typedef typename Traits::IndexSetType IndexSetType;
    typedef typename Traits::GridPartType GridPartType;

    typedef DiscontinuousGalerkinSpace<FunctionSpaceType, GridPartType, polOrd,CachingStorage> DiscreteFunctionSpaceType;
    typedef DiscreteFunctionSpaceType SpaceType;
    //typedef AdaptiveDiscreteFunction<DiscreteFunctionSpaceType> DestinationType;
    typedef DFAdapt<DiscreteFunctionSpaceType> DestinationType;
    typedef DestinationType DiscreteFunctionType;

    typedef LimiterDiscreteModel1<Model,polOrd> DiscreteModelType;
  };
  // Passes
  template <class Model,class NumFlux,int polOrd>
  class TransportDiffusionDiscreteModel1 :
    public DiscreteModelDefault<TransportDiffusionTraits1<Model,NumFlux,polOrd> > {
  public:
    typedef TransportDiffusionTraits1<Model,NumFlux,polOrd> Traits;
    
    typedef Selector<0> SelectorType;
    typedef FieldVector<double, Traits::dimDomain> DomainType;
    typedef FieldVector<double, Traits::dimDomain-1> FaceDomainType;
    typedef typename Traits::RangeType RangeType;
    typedef typename Traits::GridType GridType;
    typedef typename Traits::JacobianRangeType JacobianRangeType;
    typedef typename GridType::Traits::IntersectionIterator IntersectionIterator;
    typedef typename GridType::template Codim<0>::Entity EntityType;
    
  public:
    TransportDiffusionDiscreteModel1(const DomainType& upwind,
				     const Model& mod,
				     const NumFlux& numf) :
      upwind_(upwind),
      model_(mod),
      numflux_(numf),
      cflDiffinv_(2.*(polOrd+1.))
    {}
    bool hasSource() const { return false; }
    bool hasFlux() const { return true; }
    
    template <class ArgumentTuple>
    double numericalFlux(IntersectionIterator& it,
                         double time, const FaceDomainType& x,
                         const ArgumentTuple& uLeft, 
                         const ArgumentTuple& uRight,
                         RangeType& gLeft,
                         RangeType& gRight)
    { 
      const DomainType normal = it.integrationOuterNormal(x);
      typedef typename ElementType<0, ArgumentTuple>::Type UType;
      const UType& argULeft = Element<0>::get(uLeft);
      const UType& argURight = Element<0>::get(uRight);
      JacobianRangeType diffmatrix;
      RangeType diffflux(0.);
      /* central differences      
      model_.diffusion(*it.inside(),time,it.intersectionSelfLocal().global(x),
		       argULeft,diffmatrix);
      diffmatrix.umv(normal,diffflux);
      model_.diffusion(*it.outside(),time,it.intersectionNeighborLocal().global(x),
		       argURight,diffmatrix);
      diffmatrix.umv(normal,diffflux);
      diffflux*=0.5;
      */
      if (upwind_*normal>0)
	model_.diffusion(*it.inside(),time,
			 it.intersectionSelfLocal().global(x),
			 argULeft,diffmatrix);
      else
	model_.diffusion(*it.outside(),time,
			 it.intersectionNeighborLocal().global(x),
			 argURight,diffmatrix);
      diffmatrix.umv(normal,diffflux);
      gLeft = diffflux;
      gRight = diffflux;
      return model_.diffusionTimeStep()*cflDiffinv_;
    }

    template <class ArgumentTuple>
    double boundaryFlux(IntersectionIterator& it,
                        double time, const FaceDomainType& x,
                        const ArgumentTuple& uLeft,
                        RangeType& gLeft)
    {
      const DomainType normal = it.integrationOuterNormal(x);
      typedef typename ElementType<0, ArgumentTuple>::Type UType;
      const UType& argULeft = Element<0>::get(uLeft);
      JacobianRangeType diffmatrix;
      gLeft*=0.;
      if (upwind_*normal>0)
	model_.diffusion(*it.inside(),time,
			 it.intersectionSelfLocal().global(x),
			 argULeft,diffmatrix);
      else if (model_.hasBoundaryValue(it,time,x)) {
	UType uRight;
	model_.boundaryValue(it,time,x,argULeft,uRight);
	model_.diffusion(*it.inside(),time,
			 it.intersectionSelfLocal().global(x),
			 uRight,diffmatrix);
      } else {
	model_.diffusion(*it.inside(),time,
			 it.intersectionSelfLocal().global(x),
			 argULeft,diffmatrix);
      }
      diffmatrix.umv(normal,gLeft);
      return model_.diffusionTimeStep()*cflDiffinv_; 
    }

    template <class ArgumentTuple>
    void analyticalFlux(EntityType& en,
                        double time, const DomainType& x,
                        const ArgumentTuple& u, JacobianRangeType& f)
    { 
      typedef typename ElementType<0, ArgumentTuple>::Type UType;
      const UType& argU = Element<0>::get(u);
      model_.diffusion(en,time,x,argU,f);
    }
  private:
    const DomainType& upwind_;
    const Model& model_;
    const NumFlux& numflux_;
    const double cflDiffinv_;
  };
  template <class Model,class NumFlux,int polOrd>
  class TransportDiffusionDiscreteModel2<Model,NumFlux,polOrd,true,true> : 
    public DiscreteModelDefault<TransportDiffusionTraits2<Model,NumFlux,polOrd,true,true> > {
  public:
    typedef TransportDiffusionTraits2<Model,NumFlux,polOrd,true,true> Traits;
    typedef Selector<0, 1> SelectorType;
    typedef FieldVector<double, Traits::dimDomain> DomainType;
    typedef FieldVector<double, Traits::dimDomain-1> FaceDomainType;
    
    typedef typename Traits::RangeType RangeType;
    typedef typename Traits::GridType GridType;
    typedef typename Traits::JacobianRangeType JacobianRangeType;
    typedef typename GridType::Traits::IntersectionIterator IntersectionIterator;
    typedef typename GridType::template Codim<0>::Entity EntityType;

  public:
    TransportDiffusionDiscreteModel2(const DomainType& upwind,
				     const Model& mod,const NumFlux& numf) :
      upwind_(upwind),
      model_(mod),
      numflux_(numf)
    {}

    bool hasSource() const { return false; }
    bool hasFlux() const { return true; }

    template <class ArgumentTuple>
    double numericalFlux(IntersectionIterator& it,
                         double time, const FaceDomainType& x,
                         const ArgumentTuple& uLeft, 
                         const ArgumentTuple& uRight,
                         RangeType& gLeft,
                         RangeType& gRight)
    { 
      const DomainType normal = it.integrationOuterNormal(x);
      typedef typename ElementType<0, ArgumentTuple>::Type UType;
      typedef typename ElementType<1, ArgumentTuple>::Type WType;
      const UType& argULeft = Element<0>::get(uLeft);
      const WType& argWLeft = Element<1>::get(uLeft);
      const UType& argURight = Element<0>::get(uRight);
      const WType& argWRight = Element<1>::get(uRight);
      // Advection
      double ldt=numflux_.
	numericalFlux(it,time,x,argULeft,argURight,gLeft,gRight);
      // Diffusion
      JacobianRangeType diffmatrix;
      RangeType diffflux(0.);
      /* central 
      DomainType normal = it.integrationOuterNormal(x);
      model_.
	diffusion(*it.inside(),time,it.intersectionSelfLocal().global(x),
		  argULeft,argWLeft,diffmatrix);
      diffmatrix.umv(normal,diffflux);
      model_.
	diffusion(*it.outside(),time,it.intersectionNeighborLocal().global(x),
		  argURight,argWRight,diffmatrix);
      diffmatrix.umv(normal,diffflux);
      diffflux*=0.5;
      */
      if (upwind_*normal<0)
	model_.diffusion(*it.inside(),time,
			 it.intersectionSelfLocal().global(x),
			 argULeft,argWLeft,diffmatrix);
      else
	model_.diffusion(*it.outside(),time,
			 it.intersectionNeighborLocal().global(x),
			 argURight,argWRight,diffmatrix);
      diffmatrix.umv(normal,diffflux);
      gLeft += diffflux;
      gRight += diffflux;
      return ldt;
    }

    template <class ArgumentTuple>
    double boundaryFlux(IntersectionIterator& it,
                        double time, const FaceDomainType& x,
                        const ArgumentTuple& uLeft,
                        RangeType& gLeft)
    {
      const DomainType normal = it.integrationOuterNormal(x);
      typedef typename ElementType<0, ArgumentTuple>::Type UType;
      typedef typename ElementType<1, ArgumentTuple>::Type WType;
      const UType& argULeft = Element<0>::get(uLeft);
      const WType& argWLeft = Element<1>::get(uLeft);
      // Advection
      double ldt=0.;
      if (model_.hasBoundaryValue(it,time,x)) {
	RangeType uRight,gRight;
	model_.boundaryValue(it,time,x,argULeft,uRight);
	ldt = numflux_.numericalFlux(it,time,x,argULeft,uRight,gLeft,gRight);
	// Diffusion
	JacobianRangeType diffmatrix;
	model_.diffusion(*it.inside(),time,
			 it.intersectionSelfLocal().global(x),
			 argULeft,argWLeft,diffmatrix);
	diffmatrix.umv(normal,gLeft);
      } else {
        ldt = model_.boundaryFlux(it,time,x,argULeft,argWLeft,gLeft);
      }
      return ldt;
    }

    template <class ArgumentTuple>
    void analyticalFlux(EntityType& en,
                        double time, const DomainType& x,
                        const ArgumentTuple& u, JacobianRangeType& f)
    { 
      typedef typename ElementType<0, ArgumentTuple>::Type UType;
      typedef typename ElementType<1, ArgumentTuple>::Type WType;
      const UType& argU = Element<0>::get(u);
      const WType argW = Element<1>::get(u);
      // Advection
      model_.analyticalFlux(en,time,x,argU,f);
      // Diffusion
      JacobianRangeType diffmatrix;
      model_.diffusion(en,time,x,argU,argW,diffmatrix);
      f += diffmatrix;
    }
  private:
    const DomainType& upwind_;
    const Model& model_;
    const NumFlux& numflux_;
  };
  template <class Model,class NumFlux,int polOrd>
  class TransportDiffusionDiscreteModel2<Model,NumFlux,polOrd,false,true> : 
    public DiscreteModelDefault<TransportDiffusionTraits2<Model,NumFlux,polOrd,false,true> > {
  public:
    typedef TransportDiffusionTraits2<Model,NumFlux,polOrd,false,true> Traits;
    typedef Selector<0> SelectorType;
    typedef FieldVector<double, Traits::dimDomain> DomainType;
    typedef FieldVector<double, Traits::dimDomain-1> FaceDomainType;
    
    typedef typename Traits::RangeType RangeType;
    typedef typename Traits::GridType GridType;
    typedef typename Traits::JacobianRangeType JacobianRangeType;
    typedef typename GridType::Traits::IntersectionIterator IntersectionIterator;
    typedef typename GridType::template Codim<0>::Entity EntityType;

  public:
    TransportDiffusionDiscreteModel2(const Model& mod,const NumFlux& numf) :
      model_(mod),
      numflux_(numf) {}

    bool hasSource() const { return false; }
    bool hasFlux() const { return true; }

    template <class ArgumentTuple>
    double numericalFlux(IntersectionIterator& it,
                         double time, const FaceDomainType& x,
                         const ArgumentTuple& uLeft, 
                         const ArgumentTuple& uRight,
                         RangeType& gLeft,
                         RangeType& gRight)
    { 
      const DomainType normal = it.integrationOuterNormal(x);
      typedef typename ElementType<0, ArgumentTuple>::Type UType;
      const UType& argULeft = Element<0>::get(uLeft);
      const UType& argURight = Element<0>::get(uRight);
      // Advection
      double ldt=numflux_.
	numericalFlux(it,time,x,argULeft,argURight,gLeft,gRight);
      return ldt;
    }

    template <class ArgumentTuple>
    double boundaryFlux(IntersectionIterator& it,
                        double time, const FaceDomainType& x,
                        const ArgumentTuple& uLeft,
                        RangeType& gLeft)
    {
      const DomainType normal = it.integrationOuterNormal(x);
      typedef typename ElementType<0, ArgumentTuple>::Type UType;
      const UType& argULeft = Element<0>::get(uLeft);
      // Advection
      double ldt=0.;
      if (model_.hasBoundaryValue(it,time,x)) {
	RangeType uRight,gRight;
	model_.boundaryValue(it,time,x,argULeft,uRight);
	ldt = numflux_.numericalFlux(it,time,x,argULeft,uRight,gLeft,gRight);
      } else {
        ldt = model_.boundaryFlux(it,time,x,argULeft,gLeft);
      }
      return ldt;
    }

    template <class ArgumentTuple>
    void analyticalFlux(EntityType& en,
                        double time, const DomainType& x,
                        const ArgumentTuple& u, JacobianRangeType& f)
    { 
      typedef typename ElementType<0, ArgumentTuple>::Type UType;
      const UType& argU = Element<0>::get(u);
      // Advection
      model_.analyticalFlux(en,time,x,argU,f);
    }
  private:
    const Model& model_;
    const NumFlux& numflux_;
  };
  template <class Model,class NumFlux,int polOrd>
  class TransportDiffusionDiscreteModel2<Model,NumFlux,polOrd,true,false> : 
    public DiscreteModelDefault<TransportDiffusionTraits2<Model,NumFlux,polOrd,true,false> > {
  public:
    typedef TransportDiffusionTraits2<Model,NumFlux,polOrd,true,false> Traits;
    typedef Selector<0, 1> SelectorType;
    typedef FieldVector<double, Traits::dimDomain> DomainType;
    typedef FieldVector<double, Traits::dimDomain-1> FaceDomainType;
    
    typedef typename Traits::RangeType RangeType;
    typedef typename Traits::GridType GridType;
    typedef typename Traits::JacobianRangeType JacobianRangeType;
    typedef typename GridType::Traits::IntersectionIterator IntersectionIterator;
    typedef typename GridType::template Codim<0>::Entity EntityType;

  public:
    TransportDiffusionDiscreteModel2(const DomainType& upwind,
				     const Model& mod,const NumFlux& numf) :
      upwind_(upwind),
      model_(mod),
      numflux_(numf) {}

    bool hasSource() const { return false; }
    bool hasFlux() const { return true; }

    template <class ArgumentTuple>
    double numericalFlux(IntersectionIterator& it,
                         double time, const FaceDomainType& x,
                         const ArgumentTuple& uLeft, 
                         const ArgumentTuple& uRight,
                         RangeType& gLeft,
                         RangeType& gRight)
    { 
      const DomainType normal = it.integrationOuterNormal(x);
      typedef typename ElementType<0, ArgumentTuple>::Type UType;
      typedef typename ElementType<1, ArgumentTuple>::Type WType;
      const UType& argULeft = Element<0>::get(uLeft);
      const WType& argWLeft = Element<1>::get(uLeft);
      const UType& argURight = Element<0>::get(uRight);
      const WType& argWRight = Element<1>::get(uRight);
      double ldt=0;
      // Diffusion
      gLeft*=0.;
      JacobianRangeType diffmatrix;
      if (upwind_*normal<0)
	ldt=model_.diffusion(*it.inside(),time,
			     it.intersectionSelfLocal().global(x),
			     argULeft,argWLeft,diffmatrix);
      else
	ldt=model_.diffusion(*it.outside(),time,
			     it.intersectionNeighborLocal().global(x),
			     argURight,argWRight,diffmatrix);
      diffmatrix.umv(normal,gLeft);
      gRight = gLeft;
      return ldt;
    }

    template <class ArgumentTuple>
    double boundaryFlux(IntersectionIterator& it,
                        double time, const FaceDomainType& x,
                        const ArgumentTuple& uLeft,
                        RangeType& gLeft)
    {
      const DomainType normal = it.integrationOuterNormal(x);
      typedef typename ElementType<0, ArgumentTuple>::Type UType;
      typedef typename ElementType<1, ArgumentTuple>::Type WType;
      const UType& argULeft = Element<0>::get(uLeft);
      const WType& argWLeft = Element<1>::get(uLeft);
      double ldt=0.;
      // Diffusion
      if (model_.hasBoundaryValue(it,time,x)) {
	gLeft*=0.;
	RangeType uRight,gRight;
	model_.boundaryValue(it,time,x,argULeft,uRight);
	JacobianRangeType diffmatrix;
	ldt=model_.diffusion(*it.inside(),time,
			     it.intersectionSelfLocal().global(x),
			     argULeft,argWLeft,diffmatrix);
	diffmatrix.umv(normal,gLeft);
      } else {
        ldt = model_.boundaryFlux(it,time,x,argULeft,argWLeft,gLeft);
      }
      return ldt;
    }

    template <class ArgumentTuple>
    void analyticalFlux(EntityType& en,
                        double time, const DomainType& x,
                        const ArgumentTuple& u, JacobianRangeType& f)
    { 
      typedef typename ElementType<0, ArgumentTuple>::Type UType;
      typedef typename ElementType<1, ArgumentTuple>::Type WType;
      const UType& argU = Element<0>::get(u);
      const WType argW = Element<1>::get(u);
      // Diffusion
      JacobianRangeType diffmatrix;
      model_.diffusion(en,time,x,argU,argW,f);
    }
  private:
    const DomainType& upwind_;
    const Model& model_;
    const NumFlux& numflux_;
  };
  // **********************************************
  // **********************************************
  // **********************************************
  template <class Model,int polOrd>
  class LimiterDiscreteModel1 :
    public DiscreteModelDefault<LimiterTraits1<Model,polOrd> > 
  {
  public:
    typedef LimiterTraits1<Model,polOrd> Traits;
    
    typedef Selector<0> SelectorType;
    typedef FieldVector<double, Traits::dimDomain> DomainType;
    typedef FieldVector<double, Traits::dimDomain-1> FaceDomainType;
    typedef typename Traits::RangeType RangeType;
    typedef typename Traits::GridType GridType;
    typedef typename Traits::JacobianRangeType JacobianRangeType;
    typedef typename GridType::Traits::IntersectionIterator IntersectionIterator;
    typedef typename GridType::template Codim<0>::Entity EntityType;
    
  public:
    LimiterDiscreteModel1(const Model& mod) :
      model_(mod) {}

    bool hasSource() const { return false; }
    bool hasFlux() const { return false; }
    
  private:
    const Model& model_;
  };
}

#endif
