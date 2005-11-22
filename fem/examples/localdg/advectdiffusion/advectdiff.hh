// Gitterdefinition?
#ifndef DUNE_DGOPERATORS_HH
#define DUNE_DGOPERATORS_HH

#include "../../../pass/dgpass.hh"
#include "../../../pass/discretemodel.hh"
#include "../../../pass/selection.hh"
#include "../../../misc/timeutility.hh"
#include "../../../space/dgspace.hh"

// Dune includes
#include <dune/common/utility.hh>
#include <dune/fem/space/combinedspace.hh>
#include <dune/fem/dfadapt.hh>
#include <dune/fem/discretefunction/adaptivefunction.hh>
#include <dune/grid/common/gridpart.hh>
#include <dune/quadrature/fixedorder.hh>

#include "discretemodels.hh"

using namespace Dune;

//*************************************************************
namespace Dune {  
  template <class Model,template<class M> class NumFlux,int polOrd >
  class DGAdvectionDiffusionOperator : 
    public Operator<double,double,
	   typename TransportDiffusionDiscreteModel2<Model,NumFlux<Model>,polOrd,true,true>::Traits::DiscreteFunctionType,
	   typename TransportDiffusionDiscreteModel2<Model,NumFlux<Model>,polOrd,true,true>::Traits::DiscreteFunctionType> {
  public:
    enum { dimRange = Model::dimRange };
    enum { dimDomain = Model::Traits::dimDomain };

    typedef NumFlux<Model> NumFluxType;
    typedef typename Model::Traits::GridType GridType;

    typedef TransportDiffusionDiscreteModel1<Model,NumFluxType,polOrd> DiscreteModel1Type;
    typedef TransportDiffusionDiscreteModel2<Model,NumFluxType,polOrd,true,true> DiscreteModel2Type;
    typedef typename DiscreteModel1Type::Traits Traits1;
    typedef typename DiscreteModel2Type::Traits Traits2;

    typedef typename Traits2::DomainType DomainType;
    typedef typename Traits2::DiscreteFunctionType DiscreteFunction2Type;
    typedef StartPass<DiscreteFunction2Type> Pass0Type;
    typedef LocalDGPass<DiscreteModel1Type, Pass0Type> Pass1Type;
    typedef LocalDGPass<DiscreteModel2Type, Pass1Type> Pass2Type;

    //typedef typename Traits1::SingleSpaceType SingleSpace1Type;
    //typedef typename Traits2::SingleSpaceType SingleSpace2Type;
    typedef typename Traits1::DiscreteFunctionSpaceType 
      Space1Type;
    typedef typename Traits2::DiscreteFunctionSpaceType  
      Space2Type;
    typedef typename Traits1::DestinationType Destination1Type;
    typedef typename Traits2::DestinationType Destination2Type;
    typedef Destination2Type DestinationType;
    typedef Space2Type SpaceType;

    typedef typename Traits1::IndexSetType IndexSet1Type;
    typedef typename Traits2::IndexSetType IndexSet2Type;
    typedef typename Traits1::GridPartType GridPart1Type;
    typedef typename Traits2::GridPartType GridPart2Type;

  public:
    DGAdvectionDiffusionOperator(GridType& grid,
				 const NumFluxType& numf,
				 DomainType& upwind) :
      upwind_(upwind),
      grid_(grid),
      model_(numf.model()),
      numflux_(numf),
      iset1_(grid_,grid_.maxLevel()),
      iset2_(grid_,grid_.maxLevel()),
      gridPart1_(grid_, iset1_),
      gridPart2_(grid_, iset2_),
      //singleSpace1_(gridPart1_),
      //singleSpace2_(gridPart2_),
      //space1_(singleSpace1_),
      //space2_(singleSpace2_),
      space1_(gridPart1_),
      space2_(gridPart2_),
      problem1_(upwind_,model_,numflux_),
      problem2_(upwind_,model_,numflux_),
      pass1_(problem1_, pass0_, space1_),
      pass2_(problem2_, pass1_, space2_) 
    {}
    void timeProvider(TimeProvider* timeprovider)
    {
      pass2_.timeProvider(timeprovider);
    }
    void operator()(const DestinationType& arg, DestinationType& dest) const {
      pass2_(arg,dest);
      // upwind_*=(-1.);
    }
    void switchupwind() {upwind_*=(-1.);}
    const SpaceType& space() const {
      return space2_;
    }
  private:
    mutable DomainType upwind_;
    GridType& grid_;
    const Model& model_;
    const NumFluxType& numflux_;
    IndexSet1Type iset1_;
    IndexSet2Type iset2_;
    GridPart1Type gridPart1_;
    GridPart2Type gridPart2_;
    //SingleSpace1Type singleSpace1_;
    //SingleSpace2Type singleSpace2_;
    Space1Type space1_;
    Space2Type space2_;
    DiscreteModel1Type problem1_;
    DiscreteModel2Type problem2_;
    Pass0Type pass0_;
    Pass1Type pass1_;
    Pass2Type pass2_;
  };
/**************************************************************/
  template <class Model,template<class M> class NumFlux,int polOrd >
  class DGAdvectionOperator : 
    public Operator<double,double,
		    typename TransportDiffusionDiscreteModel2<Model,NumFlux<Model>,polOrd,false,true>::Traits::DiscreteFunctionType,
		    typename TransportDiffusionDiscreteModel2<Model,NumFlux<Model>,polOrd,false,true>::Traits::DiscreteFunctionType> {
  public:
    enum { dimRange = Model::dimRange };
    enum { dimDomain = Model::Traits::dimDomain };

    typedef NumFlux<Model> NumFluxType;
    typedef typename Model::Traits::GridType GridType;

    typedef TransportDiffusionDiscreteModel2<Model,NumFluxType,polOrd,false,true> DiscreteModel2Type;
     typedef typename DiscreteModel2Type::Traits Traits2;

    typedef typename Traits2::DomainType DomainType;
    typedef typename Traits2::DiscreteFunctionType DiscreteFunction2Type;
    typedef StartPass<DiscreteFunction2Type> Pass0Type;
    typedef LocalDGPass<DiscreteModel2Type, Pass0Type> Pass2Type;

    // typedef typename Traits2::SingleSpaceType SingleSpace2Type;
    typedef typename Traits2::DiscreteFunctionSpaceType Space2Type;
    typedef typename Traits2::DestinationType Destination2Type;
    typedef Destination2Type DestinationType;
    typedef Space2Type SpaceType;

    typedef typename Traits2::IndexSetType IndexSet2Type;
    typedef typename Traits2::GridPartType GridPart2Type;

  public:
    DGAdvectionOperator(GridType& grid,
				 const NumFluxType& numf) :
      grid_(grid),
      model_(numf.model()),
      numflux_(numf),
      iset2_(grid_,grid_.maxLevel()),
      gridPart2_(grid_, iset2_),
      //singleSpace2_(gridPart2_),
      //space2_(singleSpace2_),
      space2_(gridPart2_),
      problem2_(model_,numflux_),
      pass2_(problem2_, pass0_, space2_) 
    {}
    void timeProvider(TimeProvider* timeprovider)
    {
      pass2_.timeProvider(timeprovider);
    }
    void operator()(const DestinationType& arg, DestinationType& dest) const {
      pass2_(arg,dest);
    }
    const SpaceType& space() const {
      return space2_;
    }
  private:
    GridType& grid_;
    const Model& model_;
    const NumFluxType& numflux_;
    IndexSet2Type iset2_;
    GridPart2Type gridPart2_;
    // SingleSpace2Type singleSpace2_;
    Space2Type space2_;
    DiscreteModel2Type problem2_;
    Pass0Type pass0_;
    Pass2Type pass2_;
  };
/**************************************************************/
  template <class Model,template<class M> class NumFlux,int polOrd >
  class DGDiffusionOperator : 
    public Operator<double,double,
		    typename TransportDiffusionDiscreteModel2<Model,NumFlux<Model>,polOrd,true,false>::Traits::DiscreteFunctionType,
		    typename TransportDiffusionDiscreteModel2<Model,NumFlux<Model>,polOrd,true,false>::Traits::DiscreteFunctionType> {
  public:
    enum { dimRange = Model::dimRange };
    enum { dimDomain = Model::Traits::dimDomain };

    typedef NumFlux<Model> NumFluxType;
    typedef typename Model::Traits::GridType GridType;

    typedef TransportDiffusionDiscreteModel1<Model,NumFluxType,polOrd> DiscreteModel1Type;
    typedef TransportDiffusionDiscreteModel2<Model,NumFluxType,polOrd,true,false> DiscreteModel2Type;
    typedef typename DiscreteModel1Type::Traits Traits1;
    typedef typename DiscreteModel2Type::Traits Traits2;

    typedef typename Traits2::DomainType DomainType;
    typedef typename Traits2::DiscreteFunctionType DiscreteFunction2Type;
    typedef StartPass<DiscreteFunction2Type> Pass0Type;
    typedef LocalDGPass<DiscreteModel1Type, Pass0Type> Pass1Type;
    typedef LocalDGPass<DiscreteModel2Type, Pass1Type> Pass2Type;

    //typedef typename Traits1::SingleSpaceType SingleSpace1Type;
    //typedef typename Traits2::SingleSpaceType SingleSpace2Type;
    typedef typename Traits1::DiscreteFunctionSpaceType 
      Space1Type;
    typedef typename Traits2::DiscreteFunctionSpaceType  
      Space2Type;
    typedef typename Traits1::DestinationType Destination1Type;
    typedef typename Traits2::DestinationType Destination2Type;
    typedef Destination2Type DestinationType;
    typedef Space2Type SpaceType;

    typedef typename Traits1::IndexSetType IndexSet1Type;
    typedef typename Traits2::IndexSetType IndexSet2Type;
    typedef typename Traits1::GridPartType GridPart1Type;
    typedef typename Traits2::GridPartType GridPart2Type;

  public:
    DGDiffusionOperator(GridType& grid,
			const NumFluxType& numf,
			DomainType& upwind) :
      upwind_(upwind),
      grid_(grid),
      model_(numf.model()),
      numflux_(numf),
      iset1_(grid_,grid_.maxLevel()),
      iset2_(grid_,grid_.maxLevel()),
      gridPart1_(grid_, iset1_),
      gridPart2_(grid_, iset2_),
      //singleSpace1_(gridPart1_),
      //singleSpace2_(gridPart2_),
      //space1_(singleSpace1_),
      //space2_(singleSpace2_),
      space1_(gridPart1_),
      space2_(gridPart2_),
      problem1_(upwind_,model_,numflux_),
      problem2_(upwind_,model_,numflux_),
      pass1_(problem1_, pass0_, space1_),
      pass2_(problem2_, pass1_, space2_) 
    {}
    void timeProvider(TimeProvider* timeprovider)
    {
      pass2_.timeProvider(timeprovider);
    }
    void operator()(const DestinationType& arg, DestinationType& dest) const {
      pass2_(arg,dest);
      // upwind_*=(-1.);
    }
    void switchupwind() {upwind_*=(-1.);}
    const SpaceType& space() const {
      return space2_;
    }
  private:
    mutable DomainType upwind_;
    GridType& grid_;
    const Model& model_;
    const NumFluxType& numflux_;
    IndexSet1Type iset1_;
    IndexSet2Type iset2_;
    GridPart1Type gridPart1_;
    GridPart2Type gridPart2_;
    //SingleSpace1Type singleSpace1_;
    //SingleSpace2Type singleSpace2_;
    Space1Type space1_;
    Space2Type space2_;
    DiscreteModel1Type problem1_;
    DiscreteModel2Type problem2_;
    Pass0Type pass0_;
    Pass1Type pass1_;
    Pass2Type pass2_;
  };
}
#endif
