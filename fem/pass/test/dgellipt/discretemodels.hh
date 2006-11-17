#ifndef DUNE_EXAMPLEDISCRETEMODELS_HH
#define DUNE_EXAMPLEDISCRETEMODELS_HH

// FR passes 
#include <dune/fem/pass/dgpass.hh>
#include <dune/fem/pass/discretemodel.hh>
#include <dune/fem/pass/selection.hh>
#include <dune/fem/misc/timeutility.hh>
#include <dune/fem/space/dgspace.hh>

// Dune includes
#include <dune/common/utility.hh>
#include <dune/fem/space/combinedspace.hh>
#include <dune/fem/discretefunction/adaptivefunction.hh>
#include <dune/grid/common/gridpart.hh>
//#include "dune/fem/quadrature/cachingquad.hh"
#include <dune/grid/common/referenceelements.cc>
#include <dune/fem/quadrature/cachequad.hh>

#include <dune/fem/misc/boundaryidentifier.hh>

#include <dune/fem/operator/matrix/matrixhandler.hh>
#if HAVE_BLAS
#include <dune/fem/solver/oemsolver/oemsolver.hh>
#endif
#include <dune/fem/operator/inverseoperators.hh>

//*************************************************************
namespace LDGExample {  

  using namespace Dune;

  // MethodOrderTraits
  template <class Model,int polOrd,int dimRange>
  class ElliptPassTraits {
  public:
    typedef typename Model::Traits ModelTraits;
    typedef typename ModelTraits::GridType GridType;
    enum { dimDomain = Model::Traits::dimDomain };
    
    //typedef LeafGridPart<GridType, All_Partition > GridPartType;
    typedef LeafGridPart<GridType> GridPartType;
    //typedef HierarchicGridPart<GridType> GridPartType;

    typedef CachingQuadrature<GridType,0> VolumeQuadratureType;
    typedef CachingQuadrature<GridType,1> FaceQuadratureType;
    
    // typical tpye of space 
    typedef FunctionSpace<double, double, dimDomain, dimRange > SingleFunctionSpaceType; 
    typedef SingleFunctionSpaceType FunctionSpaceType;
    typedef DiscontinuousGalerkinSpace<SingleFunctionSpaceType, GridPartType, polOrd> ContainedSpaceType;
    typedef ContainedSpaceType DiscreteFunctionSpaceType;
    
    //typedef DFAdapt<DiscreteFunctionSpaceType> DiscreteFunctionType;
    typedef AdaptiveDiscreteFunction<DiscreteFunctionSpaceType> DiscreteFunctionType;
    
    template <class RowSpaceType, class ColSpaceType> 
    struct MatrixHandler
    { 
      typedef MatrixHandlerSPMat<RowSpaceType,ColSpaceType> MatrixHandlerType; 
    };
    
    template <class DiscreteFunctionImp, class OperatorImp> 
    struct InverseOperator
    { 
#if HAVE_BLAS
      //typedef GMRESOp    <DiscreteFunctionImp, OperatorImp> InverseOperatorType;
      typedef OEMBICGSTABOp <DiscreteFunctionImp, OperatorImp> InverseOperatorType;
      //typedef OEMGMRESOp    <DiscreteFunctionImp, OperatorImp> InverseOperatorType;
      //typedef OEMCGOp       <DiscreteFunctionImp, OperatorImp> InverseOperatorType;
      //typedef CGInverseOp   <DiscreteFunctionImp, OperatorImp> InverseOperatorType;
#else 
      typedef CGInverseOp   <DiscreteFunctionImp, OperatorImp> InverseOperatorType;
#endif
    };
  };

  ///////////////////////////////////////////////////////
  //
  //  --Gradient Traits 
  //
  ///////////////////////////////////////////////////////
  template <class Model,class NumFlux,int polOrd>
  class GradientDiscreteModel ; 
    
  template <class Model,class NumFlux,int polOrd>
  struct GradientTraits
  {
    typedef typename Model::Traits ModelTraits;
    typedef typename ModelTraits::GridType GridType;

    enum { dimDomain = ModelTraits::dimDomain };
    enum { dimRange  = dimDomain };

    typedef ElliptPassTraits<Model,polOrd,dimDomain> Traits;
    typedef typename Traits::FunctionSpaceType FunctionSpaceType;
    typedef typename FunctionSpaceType::DomainType DomainType;
    typedef typename FunctionSpaceType::RangeType RangeType;
    typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;

    typedef typename Traits::VolumeQuadratureType VolumeQuadratureType;
    typedef typename Traits::FaceQuadratureType FaceQuadratureType;
    typedef typename Traits::GridPartType GridPartType;

    typedef typename Traits::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
    typedef typename Traits::ContainedSpaceType ContainedSpaceType;
    typedef typename Traits::DiscreteFunctionType DiscreteFunctionType;
    typedef DiscreteFunctionType DestinationType;

    typedef GradientDiscreteModel<Model,NumFlux,polOrd> DiscreteModelType;
  };

  template <class Model,class NumFlux,int polOrd>
  class GradientDiscreteModel : 
    public DiscreteModelDefault<GradientTraits<Model,NumFlux,polOrd> >
  { 
  public:
    enum { polynomialOrder = polOrd };

    typedef GradientTraits<Model,NumFlux,polOrd> Traits;
    typedef Dune::Selector<0> SelectorType;
    typedef FieldVector<double, Traits::dimDomain> DomainType;
    typedef FieldVector<double, Traits::dimDomain-1> FaceDomainType;

    typedef typename Traits::RangeType RangeType;
    typedef typename Traits::GridType GridType;
    typedef typename Traits::JacobianRangeType JacobianRangeType;
    typedef typename GridType::Traits::IntersectionIterator IntersectionIterator;
    typedef typename GridType::template Codim<0>::Entity EntityType;

    typedef BoundaryIdentifier BoundaryIdentifierType ;

    enum { dimRange = Traits::dimRange };

  public:
    GradientDiscreteModel(const Model& mod,const NumFlux& numf, bool preCon = false) :
      model_(mod),
      numflux_(numf),
      preCon_(preCon)
    {}

    const Model & data () const { return model_; }

    bool preconditioning () const { return preCon_; }
    bool hasSource() const { return true; }
    bool hasFlux() const { return false; }

    template <class ArgumentTuple, class ReturnType>
    double numericalFlux(IntersectionIterator& it,
                         double time, const FaceDomainType& x,
                         const ArgumentTuple& uLeft,
                         const ArgumentTuple& uRight,
                         ReturnType& gLeft,
                         ReturnType& gRight) const 
    {
      return 0.0;
    }

    template <class ArgumentTuple, class ReturnType >
    double boundaryFlux(IntersectionIterator& it,
                        double time, const FaceDomainType& x,
                        const ArgumentTuple& uLeft,
                        ReturnType& gLeft) const
    {
      return 0.0;
    }

    // returns true, when Dirichlet boundary, false when other boundary
    // (Neumann) 
    template <class IntersectionIteratorType, class FaceDomType, 
              class ReturnType >
    BoundaryIdentifierType 
    boundaryValue(const IntersectionIteratorType& it,
                  double time, const FaceDomType& x,
                  ReturnType& bndVal) const
    {
      return BoundaryIdentifierType::undefined; 
    }
    
    template <class ArgumentTuple>
    void analyticalFlux(EntityType& en,
                        double time, const DomainType& x,
                        const ArgumentTuple& u, JacobianRangeType& f)
    {
      model_.diffusion(en,time,x,f);
      double tmp = f[0][0];
      tmp = sqrt(tmp);

      for(int i=0;i<dimRange; ++i) f[i][i] = tmp;
      
      //std::cout << f << " f \n";
    }

    template <class EntityType , class DomainType,
              class ArgumentTuple, class JacobianTuple>
    void source(EntityType &en,
                const double time,
                const DomainType& local,
                const ArgumentTuple& u,
                const JacobianTuple& jac,
                RangeType & s)
    {
      model_.diffusion(en,time,local,s);
    }

    template <class EntityType , class DomainType,
              class ArgumentTuple, class JacobianTuple, class RanType>
    void rightHandSide(EntityType &en,
                       const double time,
                       const DomainType& local,
                       const ArgumentTuple& u,
                       const JacobianTuple& jac,
                       RanType& s)
    {
      abort();
    }

  private:
    const Model& model_;
    const NumFlux& numflux_;
    const bool preCon_;
  };

  ///////////////////////////////////////////////////////////
  //
  //  --Laplace Traits 
  //
  ///////////////////////////////////////////////////////////
  template <class Model,class NumFlux,int polOrd>
  class LaplaceDiscreteModel;

  template <class Model,class NumFlux,int polOrd>
  struct LaplaceTraits
  {
    typedef typename Model::Traits ModelTraits;
    typedef typename ModelTraits::GridType GridType;

    enum { dimDomain = ModelTraits::dimDomain };
    enum { dimRange = ModelTraits::dimRange };

    typedef ElliptPassTraits<Model,polOrd,dimRange> Traits;
    typedef typename Traits::FunctionSpaceType FunctionSpaceType;
    typedef typename FunctionSpaceType::DomainType DomainType;
    typedef typename FunctionSpaceType::RangeType RangeType;
    typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;

    typedef typename Traits::VolumeQuadratureType VolumeQuadratureType;
    typedef typename Traits::FaceQuadratureType FaceQuadratureType;
    typedef typename Traits::GridPartType GridPartType;

    typedef typename Traits::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
    typedef typename Traits::ContainedSpaceType ContainedSpaceType;
    typedef typename Traits::DiscreteFunctionType DiscreteFunctionType;
    typedef DiscreteFunctionType DestinationType;

    typedef LaplaceDiscreteModel<Model,NumFlux,polOrd> DiscreteModelType;
  };

  template <class Model,class NumFlux,int polOrd>
  class LaplaceDiscreteModel : 
    public DiscreteModelDefault<LaplaceTraits<Model,NumFlux,polOrd> >
  { 
  public:
    enum { polynomialOrder = polOrd };

    typedef LaplaceTraits<Model,NumFlux,polOrd> Traits;
    typedef Dune::Selector<0> SelectorType;
    typedef FieldVector<double, Traits::dimDomain> DomainType;
    typedef FieldVector<double, Traits::dimDomain-1> FaceDomainType;

    typedef typename Traits::RangeType RangeType;
    typedef typename Traits::GridType GridType;
    typedef typename Traits::JacobianRangeType JacobianRangeType;
    typedef typename GridType::Traits::IntersectionIterator IntersectionIterator;
    typedef typename GridType::template Codim<0>::Entity EntityType;

    typedef BoundaryIdentifier BoundaryIdentifierType ;

  public:
    LaplaceDiscreteModel(const Model& mod,const NumFlux& numf, bool preCon = false) :
      model_(mod),
      numflux_(numf),
      preCon_(preCon)
    {}

    const Model & data () const { return model_; }

    bool preconditioning () const { return preCon_; }
    bool hasSource() const { return false; }
    bool hasFlux() const { return false; }

    template <class ArgumentTuple, class ReturnType>
    double numericalFlux(IntersectionIterator& it,
                         double time, const FaceDomainType& x,
                         const ArgumentTuple& uLeft,
                         const ArgumentTuple& uRight,
                         ReturnType& gLeft,
                         ReturnType& gRight) const 
    {
      return 0.0;
    }

    template <class ArgumentTuple, class ReturnType >
    double boundaryFlux(IntersectionIterator& it,
                        double time, const FaceDomainType& x,
                        const ArgumentTuple& uLeft,
                        ReturnType& gLeft) const
    {
      return 0.0;
    }

    // returns true, when Dirichlet boundary, false when other boundary
    // (Neumann) 
    template <class IntersectionIteratorType, class FaceDomType, 
              class ReturnType >
    BoundaryIdentifierType 
    boundaryValue(const IntersectionIteratorType& it,
                  double time, const FaceDomType& x,
                  ReturnType& bndVal) const
    {
      DomainType p = it.intersectionGlobal().global(x);
      bool dirich = boundaryDataFunction(&p[0],bndVal[0]);

      if(!dirich)
      {
        const DomainType normal = it.integrationOuterNormal(x);
        DomainType grad;
        rhsNeumann(&p[0],&grad[0]);
        bndVal = grad * normal;
      }
      return (dirich) ? BoundaryIdentifierType::DirichletNonZero : BoundaryIdentifierType::NeumannNonZero;
    }
    
    template <class ArgumentTuple>
    void analyticalFlux(EntityType& en,
                        double time, const DomainType& x,
                        const ArgumentTuple& u, JacobianRangeType& f)
    {
      model_.diffusion(en,time,x,f);
      double tmp = f[0][0];
      tmp = sqrt(tmp);
      f = tmp;
    }

    template <class EntityType , class DomainType,
              class ArgumentTuple, class JacobianTuple, class RanType>
    void source(EntityType &en,
                const double time,
                const DomainType& local,
                const ArgumentTuple& u,
                const JacobianTuple& jac,
                RanType& s)
    {
      abort();
    }

    template <class EntityType , class DomainType,
              class ArgumentTuple, class JacobianTuple, class RanType>
    void rightHandSide(EntityType &en,
                       const double time,
                       const DomainType& local,
                       const ArgumentTuple& u,
                       const JacobianTuple& jac,
                       RanType& s)
    {
      model_.rhsData().evaluate(en.geometry().global(local),s);
    }

  private:
    const Model& model_;
    const NumFlux& numflux_;
    const bool preCon_;
  };

}
#endif
