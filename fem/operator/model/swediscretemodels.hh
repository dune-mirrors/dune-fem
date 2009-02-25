#ifndef DUNE_FEM_LDG_DISCRETEMODELS_HH
#define DUNE_FEM_LDG_DISCRETEMODELS_HH

//#include <dune/fem/pass/dgpass.hh>
#include <dune/swe/pass/wdswedgpass.hh>
#include <dune/fem/pass/discretemodel.hh>
#include <dune/fem/pass/selection.hh>
#include <dune/fem/solver/timeprovider.hh>
//#include <dune/fem/pass/limitpass.hh>
#include <dune/swe/pass/wdswelimitpass.hh>

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
  
  //declarations

  template < class Model , class NumFlux , int polOrd
             , int passId1 = -1
             >
  class ShallowWaterDiscreteModel; // shallow water discrete model


  template < class GridPartType , class ProblemType >
  class ThinLayerShallowWaterModel;  //thin layer shallow water model


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

  // shallow water discrete model traits, 
  // also used by shallow water discrete model with wetting-drying treatment using thin water layer approach
  template< class Model,class NumFlux,int polOrd
             , int passId1 = -1
             >
  struct ShallowWaterDiscreteModelTraits
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

    typedef ShallowWaterDiscreteModel< Model, NumFlux, polOrd, passId1 > DiscreteModelType;
  }; //end of ShallowWaterDiscreteModelTraits


  // shallow water discrete model
  template < class Model,class NumFlux,int polOrd, int passId1 >
  class ShallowWaterDiscreteModel :
    public DiscreteModelDefault
      < ShallowWaterDiscreteModelTraits< Model, NumFlux, polOrd, passId1 > , passId1 >
  {
    Int2Type< passId1 > uVar;

  public:
    typedef ShallowWaterDiscreteModelTraits< Model, NumFlux, polOrd, passId1 > Traits;
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
    ShallowWaterDiscreteModel(const Model& mod,const NumFlux& numf) :
      model_(mod),
      numflux_(numf) {}

    bool hasSource() const { return true; }
    template <class ArgumentTuple, class JacobianTuple>
    void source(EntityType& en,
                double time,
                const DomainType& x,
                const ArgumentTuple& u,
                const JacobianTuple& jac,
                RangeType& s)
    {
//      typedef typename ElementType<0, ArgumentTuple>::Type UType;
//      const UType& argU = Element<0>::get(u);
      model_.source( en , time , x , u[ uVar ] , s );
    }

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
  }; //end of ShallowWaterDiscreteModel


  // shallow water discrete model with wetting-drying treatment using thin water layer approach
  // only used for analytical Model = ThinLayerShallowWaterModel
  template <class GridPartType, class ProblemType,class NumFlux,int polOrd, int passId1 >
  class ShallowWaterDiscreteModel< ThinLayerShallowWaterModel< GridPartType, ProblemType >, NumFlux, polOrd, passId1 > :
    public DiscreteModelDefault
      < ShallowWaterDiscreteModelTraits< ThinLayerShallowWaterModel< GridPartType, ProblemType >, NumFlux , polOrd , passId1 > , passId1 >
  {
    Int2Type< passId1 > uVar;
    typedef ThinLayerShallowWaterModel< GridPartType, ProblemType > Model;

  public:
    typedef ShallowWaterDiscreteModelTraits< Model, NumFlux, polOrd, passId1 > Traits;
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
    ShallowWaterDiscreteModel(const Model& mod,const NumFlux& numf) :
      model_(mod),
      numflux_(numf) {}

    bool hasSource() const { return true; }
    template <class ArgumentTuple, class JacobianTuple>
    void source(EntityType& en,
                double time,
                const DomainType& x,
                const ArgumentTuple& u,
                const JacobianTuple& jac,
                RangeType& s)
    {
//      typedef typename ElementType<0, ArgumentTuple>::Type UType;
//      const UType& argU = Element<0>::get(u);
      model_.source( en , time , x , u[ uVar ] , s );
    }

    bool hasFlux() const { return true; }

    template <class ArgumentTuple>
    double numericalFlux(IntersectionIterator& it,
                         double time, const FaceDomainType& x,
                         const ArgumentTuple& uLeft, 
                         const ArgumentTuple& uRight,
                         bool reflection,
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
      double ldt=numflux_.numericalFlux( it , time , x , uLeft[ uVar ] , uRight[ uVar ] , reflection, gLeft , gRight );
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
	ldt = numflux_.numericalFlux( it , time , x , uLeft[ uVar ] , uRight , false, gLeft , gRight );
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
  // public access is needed to access WDflag in applyLocalNeighbor in swe/pass/wdswedgpass.hh
  public:
    const Model& model_;
  private:
//    const Model& model_;
    const NumFlux& numflux_;
  }; //end of ShallowWaterDiscreteModel with wetting-drying treatment

}
#endif
