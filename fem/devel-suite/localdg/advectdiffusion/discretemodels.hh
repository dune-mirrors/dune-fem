#ifndef DUNE_EXAMPLEDISCRETEMODELS_HH
#define DUNE_EXAMPLEDISCRETEMODELS_HH

#include <dune/fem/pass/dgpass.hh>
#include <dune/fem/pass/discretemodel.hh>
#include <dune/fem/pass/selection.hh>
#include <dune/fem/solver/timeprovider.hh>
#include <dune/fem/pass/limitpass.hh>

// Dune includes
#include <dune/common/utility.hh>
#include <dune/fem/gridpart/gridpart.hh>
#include <dune/fem/space/dgspace/dgadaptiveleafgridpart.hh>

//- local includes 
#include <dune/fem/space/dgspace.hh>
#include <dune/fem/space/combinedspace.hh>

#include <dune/fem/function/adaptivefunction.hh>
#if HAVE_DUNE_ISTL
#include <dune/fem/function/blockvectorfunction.hh>
#endif
#include <dune/fem/quadrature/cachequad.hh>

//*************************************************************
namespace Dune {  
  

  template < class Model , class NumFlux , int polOrd 
             , int passId1 = -1
             >
  class AdvDiffDModel1;
  

  template < class Model , int polOrd
             , int passId1 = -1
             >
  class LimiterDiscreteModel1;
  

  template < class Model , class NumFlux , int polOrd
             , int passId1 = -1
             , int passId2 = -1
             >
  class AdvDiffDModel2; // true,true
  

  template < class Model , class NumFlux , int polOrd
             , int passId1 = -1
             >
  class AdvDiffDModel3; // false,true

  
  template < class Model , class NumFlux , int polOrd
             , int passId1 = -1
             , int passId2 = -1
             >
  class AdvDiffDModel4; // true,false

  
  // MethodOrderTraits
  template <class Model,int dimRange,int polOrd>
  class PassTraits 
  {
  public:
    typedef typename Model::Traits ModelTraits;
    typedef typename ModelTraits::GridPartType GridPartType;
    typedef typename GridPartType :: GridType GridType;
    enum { dimDomain = Model::Traits::dimDomain };

    typedef CachingQuadrature<GridPartType,0> VolumeQuadratureType;
    typedef CachingQuadrature<GridPartType,1> FaceQuadratureType;
    //typedef ElementQuadrature<GridPartType,0> VolumeQuadratureType;
    //typedef ElementQuadrature<GridPartType,1> FaceQuadratureType;
    #if 0
    typedef FunctionSpace<double, double, dimDomain, dimRange> FunctionSpaceType; 
    typedef DiscontinuousGalerkinSpace<FunctionSpaceType, GridPartType, 
               polOrd,CachingStorage> DiscreteFunctionSpaceType;
    #else
    typedef FunctionSpace<double, double, dimDomain, 1> FunctionSpaceType; 
    typedef DiscontinuousGalerkinSpace<FunctionSpaceType, GridPartType, 
               polOrd,CachingStorage> SingleDiscreteFunctionSpaceType;
    typedef CombinedSpace<SingleDiscreteFunctionSpaceType, dimRange> 
            DiscreteFunctionSpaceType; 
    #endif
    typedef AdaptiveDiscreteFunction<DiscreteFunctionSpaceType> DestinationType;
    //typedef StaticDiscreteFunction<DiscreteFunctionSpaceType> DestinationType;
  };


  // DiscreteModelTraits
  template <class Model,class NumFlux,int polOrd 
             , int passId1 = -1
             >
  struct AdvDiffTraits1 
  {
    typedef typename Model::Traits ModelTraits;
    typedef typename ModelTraits::GridType GridType;

    enum { dimRange = ModelTraits::dimGradRange };
    enum { dimDomain = ModelTraits::dimDomain };

    typedef PassTraits<Model,dimRange,polOrd> Traits;
    typedef typename Traits::FunctionSpaceType FunctionSpaceType;


    typedef typename Traits::VolumeQuadratureType VolumeQuadratureType;
    typedef typename Traits::FaceQuadratureType FaceQuadratureType;
    typedef typename Traits::GridPartType GridPartType;

    typedef typename Traits::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
    typedef typename Traits::DestinationType DestinationType;
    typedef DestinationType DiscreteFunctionType;

    typedef typename ModelTraits::DomainType DomainType;
    typedef typename DestinationType::RangeType RangeType;
    typedef typename DestinationType::JacobianRangeType JacobianRangeType;

    typedef AdvDiffDModel1<Model,NumFlux,polOrd 
             , passId1
             > DiscreteModelType;
  };

  
  template< class Model,class NumFlux,int polOrd
             , int passId1 = -1
             , int passId2 = -1
             >
  struct AdvDiffTraits2
  {
    typedef typename Model::Traits ModelTraits;
    typedef typename ModelTraits::GridType GridType;

    enum { dimRange = ModelTraits::dimRange };
    enum { dimDomain = ModelTraits::dimDomain };

    typedef PassTraits<Model,dimRange,polOrd> Traits;
    typedef typename Traits::FunctionSpaceType FunctionSpaceType;
    // typedef typename Traits::SingleDiscreteFunctionSpaceType SingleDiscreteFunctionSpaceType;
    
    typedef typename Traits::VolumeQuadratureType VolumeQuadratureType;
    typedef typename Traits::FaceQuadratureType FaceQuadratureType;
    typedef typename Traits::GridPartType GridPartType;
    typedef typename Traits::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
    typedef typename Traits::DestinationType DestinationType;
    typedef DestinationType DiscreteFunctionType;

    typedef typename DestinationType::DomainType DomainType;
    typedef typename DestinationType::RangeType RangeType;
    typedef typename DestinationType::JacobianRangeType JacobianRangeType;

    typedef AdvDiffDModel2<Model,NumFlux,polOrd
             , passId1
             , passId2
             > DiscreteModelType;
  };


  template< class Model,class NumFlux,int polOrd
             , int passId1 = -1
             >
  struct AdvDiffTraits3
  {
    typedef typename Model::Traits ModelTraits;
    typedef typename ModelTraits::GridType GridType;

    enum { dimRange = ModelTraits::dimRange };
    enum { dimDomain = ModelTraits::dimDomain };

    typedef PassTraits<Model,dimRange,polOrd> Traits;
    typedef typename Traits::FunctionSpaceType FunctionSpaceType;
    // typedef typename Traits::SingleDiscreteFunctionSpaceType SingleDiscreteFunctionSpaceType;
    
    typedef typename Traits::VolumeQuadratureType VolumeQuadratureType;
    typedef typename Traits::FaceQuadratureType FaceQuadratureType;
    typedef typename Traits::GridPartType GridPartType;
    typedef typename Traits::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
    typedef typename Traits::DestinationType DestinationType;
    typedef DestinationType DiscreteFunctionType;

    typedef typename DestinationType::DomainType DomainType;
    typedef typename DestinationType::RangeType RangeType;
    typedef typename DestinationType::JacobianRangeType JacobianRangeType;

    typedef AdvDiffDModel3<Model,NumFlux,polOrd
             , passId1
             > DiscreteModelType;
  };

  
  template< class Model,class NumFlux,int polOrd
             , int passId1 = -1
             , int passId2 = -1
             >
  struct AdvDiffTraits4
  {
    typedef typename Model::Traits ModelTraits;
    typedef typename ModelTraits::GridType GridType;

    enum { dimRange = ModelTraits::dimRange };
    enum { dimDomain = ModelTraits::dimDomain };

    typedef PassTraits<Model,dimRange,polOrd> Traits;
    typedef typename Traits::FunctionSpaceType FunctionSpaceType;
    // typedef typename Traits::SingleDiscreteFunctionSpaceType SingleDiscreteFunctionSpaceType;
    
    typedef typename Traits::VolumeQuadratureType VolumeQuadratureType;
    typedef typename Traits::FaceQuadratureType FaceQuadratureType;
    typedef typename Traits::GridPartType GridPartType;
    typedef typename Traits::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
    typedef typename Traits::DestinationType DestinationType;
    typedef DestinationType DiscreteFunctionType;

    typedef typename DestinationType::DomainType DomainType;
    typedef typename DestinationType::RangeType RangeType;
    typedef typename DestinationType::JacobianRangeType JacobianRangeType;

    typedef AdvDiffDModel4<Model,NumFlux,polOrd
             , passId1
             , passId2
             > DiscreteModelType;
  };

  
  template < class Model,int polOrd
             , int passId1 = -1
             , int passId2 = -1
             >
  struct LimiterTraits1 
  {
    typedef typename Model::Traits ModelTraits;
    typedef typename ModelTraits::GridType GridType;

    enum { dimRange = ModelTraits::dimRange };
    enum { dimDomain = ModelTraits::dimDomain };

    typedef PassTraits<Model,dimRange,polOrd> Traits;
    typedef typename Traits::FunctionSpaceType FunctionSpaceType;

    //typedef DiscontinuousGalerkinSpace<FunctionSpaceType, GridPartType, polOrd,CachingStorage> DiscreteFunctionSpaceType;
    //typedef AdaptiveDiscreteFunction<DiscreteFunctionSpaceType> DestinationType;
    //typedef DestinationType DiscreteFunctionType;

    typedef typename Traits::VolumeQuadratureType VolumeQuadratureType;
    typedef typename Traits::FaceQuadratureType FaceQuadratureType;
    typedef typename Traits::GridPartType GridPartType;
    typedef typename Traits::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
    typedef typename Traits::DestinationType DestinationType;
    typedef DestinationType DiscreteFunctionType;

    typedef typename DestinationType::DomainType DomainType;
    typedef typename DestinationType::RangeType RangeType;
    typedef typename DestinationType::JacobianRangeType JacobianRangeType;

    typedef LimiterDiscreteModel1< Model , polOrd 
             , passId1
             > DiscreteModelType;
  };

  
  // Passes
  template < class Model , class NumFlux , int polOrd
             , int passId1
             >
  class AdvDiffDModel1 :
    public DiscreteModelDefault
      < AdvDiffTraits1< Model , NumFlux , polOrd , passId1 > , passId1 >
  {
    Int2Type< passId1 > uVar;
  public:
    typedef AdvDiffTraits1< Model , NumFlux , polOrd 
             , passId1
             > Traits;
    
    //typedef Selector< 0 > OldSelectorType;
    //typedef Selector< passId1 > SelectorType;
    typedef FieldVector<double, Traits::dimDomain> DomainType;
    typedef FieldVector<double, Traits::dimDomain-1> FaceDomainType;
    typedef typename Traits::RangeType RangeType;
    typedef typename Traits::GridType GridType;
    typedef typename Traits::JacobianRangeType JacobianRangeType;
    typedef typename Traits::GridPartType::IntersectionIteratorType IntersectionIterator;
    typedef typename GridType::template Codim<0>::Entity EntityType;
    
  public:
    AdvDiffDModel1(const DomainType& upwind,
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
      
      //typedef typename ElementType<0, ArgumentTuple>::Type UType;
      //const UType& argULeft = Element<0>::get(uLeft);
      //const UType& argURight = Element<0>::get(uRight);

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
	model_.diffusion( *it.inside() , time ,
			 it.intersectionSelfLocal().global(x)
       , uLeft[ uVar ] , diffmatrix );
      else
	model_.diffusion( *it.outside() , time ,
			 it.intersectionNeighborLocal().global(x) ,
			 uRight[ uVar ] , diffmatrix );
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
      
      typedef typename ArgumentTuple::template Get< passId1 >::Type UType;
      //const UType& argULeft = Element<0>::get(uLeft);
      
      JacobianRangeType diffmatrix;
      gLeft*=0.;
      if (upwind_*normal>0)
	model_.diffusion(*it.inside(),time,
			 it.intersectionSelfLocal().global(x),
			 uLeft[ uVar ] , diffmatrix );
      else if (model_.hasBoundaryValue(it,time,x)) {
	UType uRight;
	model_.boundaryValue(it,time,x,uLeft[ uVar ] , uRight);
	model_.diffusion(*it.inside(),time,
			 it.intersectionSelfLocal().global(x),
			 uRight,diffmatrix);
      } else {
	model_.diffusion(*it.inside(),time,
			 it.intersectionSelfLocal().global(x),
			 uLeft[ uVar ] , diffmatrix );
      }
      diffmatrix.umv(normal,gLeft);
      return model_.diffusionTimeStep()*cflDiffinv_; 
    }

    template <class ArgumentTuple>
    void analyticalFlux(EntityType& en,
                        double time, const DomainType& x,
                        const ArgumentTuple& u, JacobianRangeType& f)
    { 
      //typedef typename ElementType<0, ArgumentTuple>::Type UType;
      //const UType& argU = Element<0>::get(u);
      
      model_.diffusion(en,time,x,u[ uVar ],f);
    }
  private:
    const DomainType& upwind_;
    const Model& model_;
    const NumFlux& numflux_;
    const double cflDiffinv_;
  };

  
  template <class Model,class NumFlux,int polOrd
             , int passId1
             , int passId2
             >
  class AdvDiffDModel2 :
    public DiscreteModelDefault
      < AdvDiffTraits2< Model , NumFlux , polOrd , passId1 , passId2 > 
        , passId1 , passId2 >
  {
    Int2Type< passId1 > uVar;
    Int2Type< passId2 > vVar;
  public:
    typedef AdvDiffTraits2<Model,NumFlux,polOrd
             , passId1
             , passId2
             > Traits;
    
    //typedef Selector<0, 1> OldSelectorType;
    //typedef Selector< passId1 , passId2 > SelectorType;
    typedef FieldVector<double, Traits::dimDomain> DomainType;
    typedef FieldVector<double, Traits::dimDomain-1> FaceDomainType;
    
    typedef typename Traits::RangeType RangeType;
    typedef typename Traits::GridType GridType;
    typedef typename Traits::JacobianRangeType JacobianRangeType;
    typedef typename Traits::GridPartType::IntersectionIteratorType IntersectionIterator;
    typedef typename GridType::template Codim<0>::Entity EntityType;

  public:
    AdvDiffDModel2(const DomainType& upwind,
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
      
      //typedef typename ElementType<0, ArgumentTuple>::Type UType;
      //typedef typename ElementType<1, ArgumentTuple>::Type WType;
      //const UType& argULeft = Element<0>::get(uLeft);
      //const WType& argWLeft = Element<1>::get(uLeft);
      //const UType& argURight = Element<0>::get(uRight);
      //const WType& argWRight = Element<1>::get(uRight);
      
      // Advection
      double ldt=numflux_.
	numericalFlux( it , time , x , uLeft[ uVar ] , uRight[ uVar ] , gLeft , gRight );
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
			 uLeft[ uVar ] , uLeft[ vVar ] , diffmatrix );
      else
	model_.diffusion(*it.outside(),time,
			 it.intersectionNeighborLocal().global(x),
			 uRight[ uVar ] , uRight[ vVar ] , diffmatrix );
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
      
      //typedef typename ElementType<0, ArgumentTuple>::Type UType;
      //typedef typename ElementType<1, ArgumentTuple>::Type WType;
      //const UType& argULeft = Element<0>::get(uLeft);
      //const WType& argWLeft = Element<1>::get(uLeft);
      
      // Advection
      double ldt=0.;
      if (model_.hasBoundaryValue(it,time,x)) {
	RangeType uRight,gRight;
	model_.boundaryValue( it , time , x , uLeft[ uVar ] , uRight );
	ldt = numflux_.numericalFlux( it , time , x , uLeft[ uVar ], uRight , gLeft , gRight );
	// Diffusion
	JacobianRangeType diffmatrix;
	model_.diffusion(*it.inside(),time,
			 it.intersectionSelfLocal().global(x),
			 uLeft[ uVar ] , uLeft[ vVar ] , diffmatrix );
	diffmatrix.umv(normal,gLeft);
      } else {
        ldt = model_.boundaryFlux( it , time , x , uLeft[ uVar ] , uLeft[ vVar ] , gLeft );
      }
      return ldt;
    }

    template <class ArgumentTuple>
    void analyticalFlux(EntityType& en,
                        double time, const DomainType& x,
                        const ArgumentTuple& u, JacobianRangeType& f)
    { 
      //typedef typename ElementType<0, ArgumentTuple>::Type UType;
      //typedef typename ElementType<1, ArgumentTuple>::Type WType;
      //const UType& argU = Element<0>::get(u);
      //const WType& argW = Element<1>::get(u);
      
      // Advection
      model_.analyticalFlux( en , time , x , u[ uVar ] , f );
      // Diffusion
      JacobianRangeType diffmatrix;
      model_.diffusion( en , time , x , u[ uVar ] , u[ vVar ] , diffmatrix );
      f += diffmatrix;
    }
  private:
    const DomainType& upwind_;
    const Model& model_;
    const NumFlux& numflux_;
  };

  
  template <class Model,class NumFlux,int polOrd
             , int passId1
             >
  class AdvDiffDModel3 :
    public DiscreteModelDefault
      < AdvDiffTraits3< Model , NumFlux , polOrd , passId1 > , passId1 >
  {
    Int2Type< passId1 > uVar;

  public:
    typedef AdvDiffTraits3< Model , NumFlux , polOrd , passId1 > Traits;
    //typedef Selector<0> OldSelectorType;
    //typedef Selector< passId1 > SelectorType;
    typedef FieldVector<double, Traits::dimDomain> DomainType;
    typedef FieldVector<double, Traits::dimDomain-1> FaceDomainType;
    
    typedef typename Traits::RangeType RangeType;
    typedef typename Traits::GridType GridType;
    typedef typename Traits::JacobianRangeType JacobianRangeType;
    typedef typename Traits::GridPartType::IntersectionIteratorType IntersectionIterator;
    typedef typename GridType::template Codim<0>::Entity EntityType;

  public:
    AdvDiffDModel3(const Model& mod,const NumFlux& numf) :
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
      /*
      typedef typename ElementType<0, ArgumentTuple>::Type UType;
      const UType& argULeft = Element<0>::get(uLeft);
      const UType& argURight = Element<0>::get(uRight);
      */

      //typedef typename ArgumentTuple :: template Get< 0 > ::Type UType;
      //const UType& argULeft = uLeft.template get< 0 >();
      //const UType& argURight = uRight[ uVar ];
      
      // Advection
      double ldt=numflux_.
	numericalFlux( it , time , x , uLeft[ uVar ] , uRight[ uVar ] , gLeft , gRight );
      return ldt;
    }

    template <class ArgumentTuple>
    double boundaryFlux(IntersectionIterator& it,
                        double time, const FaceDomainType& x,
                        const ArgumentTuple& uLeft,
                        RangeType& gLeft)
    {
      const DomainType normal = it->integrationOuterNormal(x);
      
      //typedef typename ElementType<0, ArgumentTuple>::Type UType;
      //const UType& argULeft = Element<0>::get(uLeft);
      
      // Advection
      double ldt=0.;
      if (model_.hasBoundaryValue(it,time,x)) {
	RangeType uRight,gRight;
	model_.boundaryValue( it , time , x , uLeft[ uVar ] , uRight );
	ldt = numflux_.numericalFlux( it , time , x , uLeft[ uVar ] , uRight , gLeft , gRight );
      } else {
        ldt = model_.boundaryFlux( it , time , x , uLeft[ uVar ] , gLeft );
      }
      return ldt;
    }

    template <class ArgumentTuple>
    void analyticalFlux(EntityType& en,
                        double time, const DomainType& x,
                        const ArgumentTuple& u, JacobianRangeType& f)
    { 
      //typedef typename ElementType<0, ArgumentTuple>::Type UType;
      //const UType& argU = Element<0>::get(u);
      
      // Advection
      model_.analyticalFlux( en , time , x , u[ uVar ] , f );
    }
  private:
    const Model& model_;
    const NumFlux& numflux_;
  };
  
  template <class Model,class NumFlux,int polOrd
             , int passId1
             , int passId2
             >
  class AdvDiffDModel4 :
    public DiscreteModelDefault< AdvDiffTraits4< Model , NumFlux , polOrd
             , passId1
             , passId2
             >
             , passId1
             , passId2
             >
  {
    Int2Type< passId1 > uVar;
    Int2Type< passId2 > vVar;
  public:
    typedef AdvDiffTraits4< Model , NumFlux , polOrd 
             , passId1
             , passId2
             > Traits;
    //typedef Selector<0, 1> OldSelectorType;
    //typedef Selector< passId1 , passId2 > SelectorType;
    typedef FieldVector<double, Traits::dimDomain> DomainType;
    typedef FieldVector<double, Traits::dimDomain-1> FaceDomainType;
    
    typedef typename Traits::RangeType RangeType;
    typedef typename Traits::GridType GridType;
    typedef typename Traits::JacobianRangeType JacobianRangeType;
    typedef typename Traits::GridPartType::IntersectionIteratorType IntersectionIterator;
    typedef typename GridType::template Codim<0>::Entity EntityType;

  public:
    AdvDiffDModel4(const DomainType& upwind,
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
      
      //typedef typename ElementType<0, ArgumentTuple>::Type UType;
      //typedef typename ElementType<1, ArgumentTuple>::Type WType;
      //const UType& argULeft = Element<0>::get(uLeft);
      //const WType& argWLeft = Element<1>::get(uLeft);
      //const UType& argURight = Element<0>::get(uRight);
      //const WType& argWRight = Element<1>::get(uRight);
      
      double ldt=0;
      // Diffusion
      gLeft*=0.;
      JacobianRangeType diffmatrix;
      if (upwind_*normal<0)
	ldt=model_.diffusion(*it.inside(),time,
			     it.intersectionSelfLocal().global(x),
			     uLeft[ uVar ] , uLeft[ vVar ] , diffmatrix );
      else
	ldt=model_.diffusion(*it.outside(),time,
			     it.intersectionNeighborLocal().global(x),
			     uRight[ uVar ] , uRight[ vVar ] , diffmatrix );
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
      
      //typedef typename ElementType<0, ArgumentTuple>::Type UType;
      //typedef typename ElementType<1, ArgumentTuple>::Type WType;
      //const UType& argULeft = Element<0>::get(uLeft);
      //const WType& argWLeft = Element<1>::get(uLeft);
      
      double ldt=0.;
      // Diffusion
      if (model_.hasBoundaryValue(it,time,x)) {
	gLeft*=0.;
	RangeType uRight,gRight;
	model_.boundaryValue( it , time , x , uLeft[uVar ] , uRight );
	JacobianRangeType diffmatrix;
	ldt=model_.diffusion(*it.inside(),time,
			     it.intersectionSelfLocal().global(x),
			     uLeft[ uVar ] , uLeft[ vVar ] , diffmatrix );
	diffmatrix.umv(normal,gLeft);
      } else {
        ldt = model_.boundaryFlux( it , time , x , uLeft[ uVar ] , uLeft[ vVar ] , gLeft );
      }
      return ldt;
    }

    template <class ArgumentTuple>
    void analyticalFlux(EntityType& en,
                        double time, const DomainType& x,
                        const ArgumentTuple& u, JacobianRangeType& f)
    { 
      //typedef typename ElementType<0, ArgumentTuple>::Type UType;
      //typedef typename ElementType<1, ArgumentTuple>::Type WType;
      //const UType& argU = Element<0>::get(u);
      //const WType argW = Element<1>::get(u);
      
      // Diffusion
      JacobianRangeType diffmatrix;
      model_.diffusion( en , time , x , u[ uVar ] , u[ vVar ] , f );
    }
  private:
    const DomainType& upwind_;
    const Model& model_;
    const NumFlux& numflux_;
  };
  // **********************************************
  // **********************************************
  // **********************************************
  template <class Model,int polOrd
             , int passId1
             >
  class LimiterDiscreteModel1 :
    public DiscreteModelDefault
      < LimiterTraits1< Model , polOrd , passId1 > , passId1 >
  {
    Int2Type< passId1 > uVar;
    public:
    typedef LimiterTraits1< Model , polOrd
             , passId1
             > Traits;
    
    //typedef Selector<0> OldSelectorType;
    //typedef Selector< passId1 > SelectorType;
    typedef FieldVector<double, Traits::dimDomain> DomainType;
    typedef FieldVector<double, Traits::dimDomain-1> FaceDomainType;
    typedef typename Traits::RangeType RangeType;
    typedef typename Traits::GridType GridType;
    typedef typename Traits::JacobianRangeType JacobianRangeType;
    typedef typename Traits::GridPartType::IntersectionIteratorType IntersectionIterator;
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
