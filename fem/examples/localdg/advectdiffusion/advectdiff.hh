// Gitterdefinition?
#ifndef DUNE_EXAMPLEPROBLEM_HH
#define DUNE_EXAMPLEPROBLEM_HH

#include "../../../pass/dgpass.hh"
#include "../../../pass/discretemodel.hh"
#include "../../../pass/selection.hh"
#include "../../../misc/timeutility.hh"

// Dune includes
#include <dune/common/utility.hh>
#include "../../../space/dgspace.hh"
#include <dune/fem/space/combinedspace.hh>
#include <dune/fem/dfadapt.hh>
#include <dune/fem/discretefunction/adaptivefunction.hh>
#include <dune/grid/common/gridpart.hh>
#include <dune/quadrature/fixedorder.hh>

using namespace Dune;

//*************************************************************
namespace Dune {  
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

  template <class Model,class NumFlux,int polOrd >
  class TransportDiffusionDiscreteModel1;
  template <class Model,class NumFlux,int polOrd >
  class TransportDiffusionDiscreteModel2;

  // MethodOrderTraits
  template <class Model,int polOrd>
  class PassTraits {
  public:
    // * Need to do: adapt quadrature order
    typedef FixedOrderQuad<
      double, FieldVector<double, 2>, 5> VolumeQuadratureType;
    typedef FixedOrderQuad<
      double, FieldVector<double, 1>, 5> FaceQuadratureType;

    typedef typename Model::Traits::GridType GridType;
    enum { dimDomain = Model::Traits::dimDomain };
    typedef DefaultGridIndexSet<GridType, LevelIndex> IndexSetType;
    typedef DefaultGridPart<GridType, IndexSetType> GridPartType;
    typedef FunctionSpace<double, double, dimDomain, 1> SingleFunctionSpaceType; 
  };
  // DiscreteModelTraits
  template <class Model,class NumFlux,int polOrd >
  struct TransportDiffusionTraits1 
  {
    typedef typename Model::Traits ModelTraits;
    typedef typename ModelTraits::GridType GridType;

    enum { dimRange = ModelTraits::dimGradRange };
    enum { dimDomain = ModelTraits::dimDomain };

    typedef FunctionSpace<
      double, double, dimDomain, dimRange> FunctionSpaceType;
    typedef typename ModelTraits::DomainType DomainType;
    typedef typename FunctionSpaceType::RangeType RangeType;
    typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;

    typedef PassTraits<Model,polOrd> Traits;
    typedef typename Traits::VolumeQuadratureType VolumeQuadratureType;
    typedef typename Traits::FaceQuadratureType FaceQuadratureType;
    typedef typename Traits::IndexSetType IndexSetType;
    typedef typename Traits::GridPartType GridPartType;
    typedef typename Traits::SingleFunctionSpaceType SingleFunctionSpaceType;

    typedef DiscontinuousGalerkinSpace<SingleFunctionSpaceType, GridPartType, polOrd> SingleSpaceType;
    typedef CombinedSpace<SingleSpaceType, dimRange> DiscreteFunctionSpaceType;
    typedef DiscreteFunctionSpaceType SpaceType;
    typedef AdaptiveDiscreteFunction<DiscreteFunctionSpaceType> DestinationType;
    typedef DestinationType DiscreteFunctionType;

    typedef TransportDiffusionDiscreteModel1<Model,NumFlux,polOrd> DiscreteModelType;
  };
  template <class Model,class NumFlux,int polOrd >
  struct TransportDiffusionTraits2
  {
    typedef typename Model::Traits ModelTraits;
    typedef typename ModelTraits::GridType GridType;

    enum { dimRange = ModelTraits::dimRange };
    enum { dimDomain = ModelTraits::dimDomain };

    typedef FunctionSpace<
      double, double, dimDomain, dimRange> FunctionSpaceType;
    typedef typename FunctionSpaceType::DomainType DomainType;
    typedef typename FunctionSpaceType::RangeType RangeType;
    typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;

    typedef PassTraits<Model,polOrd> Traits;
    typedef typename Traits::VolumeQuadratureType VolumeQuadratureType;
    typedef typename Traits::FaceQuadratureType FaceQuadratureType;
    typedef typename Traits::IndexSetType IndexSetType;
    typedef typename Traits::GridPartType GridPartType;
    typedef typename Traits::SingleFunctionSpaceType SingleFunctionSpaceType;

    typedef DiscontinuousGalerkinSpace<SingleFunctionSpaceType, GridPartType, polOrd> SingleSpaceType;
    typedef CombinedSpace<SingleSpaceType, dimRange> DiscreteFunctionSpaceType;
    typedef DiscreteFunctionSpaceType SpaceType;
    typedef AdaptiveDiscreteFunction<DiscreteFunctionSpaceType> DestinationType;
    typedef DestinationType DiscreteFunctionType;

    typedef TransportDiffusionDiscreteModel2<Model,NumFlux,polOrd> DiscreteModelType;
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
      return 1e-10;
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
      else
	numflux_.boundaryDiffusionFlux(it,time,x,
				       argULeft,diffmatrix);
      diffmatrix.umv(normal,gLeft);
      return 1e-10; 
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
  };
  template <class Model,class NumFlux,int polOrd >
  class TransportDiffusionDiscreteModel2 : 
    public DiscreteModelDefault<TransportDiffusionTraits2<Model,NumFlux,polOrd> > {
  public:
    typedef TransportDiffusionTraits2<Model,NumFlux,polOrd> Traits;
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
      double ldt=numflux_.
	boundaryFlux(it,time,x,argULeft,gLeft);
      JacobianRangeType diffmatrix;
      model_.diffusion(*it.inside(),time,
		       it.intersectionSelfLocal().global(x),
		       argULeft,argWLeft,diffmatrix);
      diffmatrix.umv(normal,gLeft);
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
  // The actual operator
  template <class Model,template<class M> class NumFlux,int polOrd >
  class DGAdvectionDiffusionOperator : 
    public Operator<double,double,
		    typename TransportDiffusionDiscreteModel2<Model,NumFlux<Model>,polOrd>::Traits::DiscreteFunctionType,
		    typename TransportDiffusionDiscreteModel2<Model,NumFlux<Model>,polOrd>::Traits::DiscreteFunctionType> {
  public:
    enum { dimRange = Model::dimRange };
    enum { dimDomain = Model::Traits::dimDomain };

    typedef NumFlux<Model> NumFluxType;
    typedef typename Model::Traits::GridType GridType;

    typedef TransportDiffusionDiscreteModel1<Model,NumFluxType,polOrd> DiscreteModel1Type;
    typedef TransportDiffusionDiscreteModel2<Model,NumFluxType,polOrd> DiscreteModel2Type;
    typedef typename DiscreteModel1Type::Traits Traits1;
    typedef typename DiscreteModel2Type::Traits Traits2;

    typedef typename Traits2::DomainType DomainType;
    typedef typename Traits2::DiscreteFunctionType DiscreteFunction2Type;
    typedef StartPass<DiscreteFunction2Type> Pass0Type;
    typedef LocalDGPass<DiscreteModel1Type, Pass0Type> Pass1Type;
    typedef LocalDGPass<DiscreteModel2Type, Pass1Type> Pass2Type;

    typedef typename Traits1::SingleSpaceType
      SingleSpace1Type;
    typedef typename Traits2::SingleSpaceType  
      SingleSpace2Type;
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
    DGAdvectionDiffusionOperator(GridType& grid,TimeProvider& timeprovider,
				 const NumFluxType& numf,
				 DomainType& upwind) :
      upwind_(upwind),
      grid_(grid),
      model_(numflux_.model()),
      numflux_(numf),
      iset1_(grid_,grid_.maxLevel()),
      iset2_(grid_,grid_.maxLevel()),
      gridPart1_(grid_, iset1_),
      gridPart2_(grid_, iset2_),
      singleSpace1_(gridPart1_),
      singleSpace2_(gridPart2_),
      space1_(singleSpace1_),
      space2_(singleSpace2_),
      problem1_(upwind_,model_,numflux_),
      problem2_(upwind_,model_,numflux_),
      pass1_(problem1_, pass0_, space1_),
      pass2_(problem2_, pass1_, space2_)      
    {
      pass2_.timeProvider(&timeprovider);
    }
    void operator()(const DestinationType& arg, DestinationType& dest) const {
      pass2_(arg,dest);
      upwind_*=(-1.);
    }
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
    GridPart1Type gridPart2_;
    SingleSpace1Type singleSpace1_;
    SingleSpace2Type singleSpace2_;
    Space1Type space1_;
    Space2Type space2_;
    DiscreteModel1Type problem1_;
    DiscreteModel2Type problem2_;
    Pass0Type pass0_;
    Pass1Type pass1_;
    Pass2Type pass2_;
  };
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
      model_.analyticalFlux(*it.outside(), time,
			    it.intersectionNeighborLocal().global(x),
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
}

#endif
