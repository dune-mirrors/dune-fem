#ifndef DUNE_EXAMPLEDISCRETEMODELS_HH
#define DUNE_EXAMPLEDISCRETEMODELS_HH

// FR passes 
#include <dune/fem/pass/dgpass.hh>
#include <dune/fem/pass/selection.hh>
#include <dune/fem/pass/dgdiscretemodel.hh>
#include <dune/fem/space/dgspace.hh>
#include <dune/fem/space/combinedspace.hh>

// Dune includes
#include <dune/fem/space/combinedspace.hh>
#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/gridpart/gridpart.hh>
#include <dune/fem/quadrature/cachingquadrature.hh>

#include <dune/fem/misc/boundaryidentifier.hh>

#include <dune/fem/operator/matrix/spmatrix.hh>
#if HAVE_BLAS
#include <dune/fem/solver/oemsolver/oemsolver.hh>
#endif


#include <dune/fem/operator/matrix/istlmatrix.hh>
#include <dune/fem/solver/istlsolver.hh>

#include <dune/fem/solver/inverseoperators.hh>

#include <dune/fem/gridpart/adaptiveleafgridpart.hh>

#include <dune/fem/operator/2order/dgprimaloperator.hh>

#include <dune/fem/pass/dgelliptpass.hh>

#if HAVE_DUNE_ISTL
#include <dune/istl/bvector.hh>
#include <dune/fem/function/blockvectorfunction.hh>
#endif


#define USE_DUNE_ISTL HAVE_DUNE_ISTL
//#define USE_DUNE_ISTL 0

//*************************************************************
namespace LDGExample
{

  using namespace Dune;

  // MethodOrderTraits
  template< class Model, int polOrd, int dimRange >
  struct ElliptPassTraits
  {
    typedef typename Model::Traits ModelTraits;
    typedef typename ModelTraits::GridType GridType;
    enum { dimDomain = Model::Traits::dimDomain };
    
    //typedef LeafGridPart<GridType> GridPartType;
    typedef DGAdaptiveLeafGridPart<GridType> GridPartType;

    typedef CachingQuadrature<GridPartType,0> VolumeQuadratureType;
    typedef CachingQuadrature<GridPartType,1> FaceQuadratureType;
    
    //typedef CachingQuadrature<GridPartType,0> VolumeQuadratureType;
    //typedef ElementQuadrature<GridPartType,1> FaceQuadratureType;
    
    // typical tpye of space 
    typedef FunctionSpace< typename Model::DomainFieldType, typename Model::RangeFieldType , 
                           dimDomain, dimRange > FunctionSpaceType; 
    typedef DiscontinuousGalerkinSpace<FunctionSpaceType, GridPartType, polOrd, CachingStorage > ContainedSpaceType;
    typedef ContainedSpaceType DiscreteFunctionSpaceType;
    
    typedef AdaptiveDiscreteFunction<DiscreteFunctionSpaceType> DiscreteFunctionType;
  };

  ///////////////////////////////////////////////////////////
  //
  //  --Laplace Traits 
  //
  ///////////////////////////////////////////////////////////
  template <class Model,int polOrd, int passId = -1 >
  class LaplaceDiscreteModel;

  template <class Model,int polOrd, int passId = -1 >
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
    typedef typename FunctionSpaceType::RangeFieldType RangeFieldType;

    typedef typename Traits::VolumeQuadratureType VolumeQuadratureType;
    typedef typename Traits::FaceQuadratureType FaceQuadratureType;
    typedef typename Traits::GridPartType GridPartType;

    typedef typename Traits::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
    typedef typename Traits::ContainedSpaceType ContainedSpaceType;

#if USE_DUNE_ISTL
    //! only working for DGSpace     
    typedef BlockVectorDiscreteFunction<DiscreteFunctionSpaceType> DiscreteFunctionType;
#else 
    typedef typename Traits::DiscreteFunctionType DiscreteFunctionType;
#endif
    typedef DiscreteFunctionType DestinationType;

    typedef LaplaceDiscreteModel< Model, polOrd, passId > DGDiscreteModelType;

    template <class PreviousPassType>
    struct LocalOperatorSelector
    {
#if USE_DUNE_ISTL                                
      typedef ISTLMatrixTraits<DiscreteFunctionSpaceType> MatrixObjectTraits;
#else 
      typedef SparseRowMatrixTraits<DiscreteFunctionSpaceType,
                                    DiscreteFunctionSpaceType> MatrixObjectTraits; 
      //typedef BlockMatrixTraits<DiscreteFunctionSpaceType,
      //                          DiscreteFunctionSpaceType> MatrixObjectTraits; 
#endif
      //! The pass id for DGPrimalOperator-pass should be something different than, 
      //! template-given passId, which is given for this template class.
      //! Since this is the last pass, doesn't play role, 
      //! but it's pass id can not be -1 or omitted
      typedef DGPrimalOperator< DGDiscreteModelType, PreviousPassType, MatrixObjectTraits, passId >
        LocalOperatorType;

#if USE_DUNE_ISTL
      typedef ISTLBICGSTABOp <DiscreteFunctionType, LocalOperatorType> InverseOperatorType;
#elif HAVE_BLAS
      //typedef GMRESOp    <DiscreteFunctionType, LocalOperatorType> InverseOperatorType;
      typedef OEMBICGSTABOp <DiscreteFunctionType, LocalOperatorType> InverseOperatorType;
      //typedef OEMGMRESOp    <DiscreteFunctionType, LocalOperatorType> InverseOperatorType;
      //typedef OEMCGOp       <DiscreteFunctionType, LocalOperatorType> InverseOperatorType;
      //typedef CGInverseOp   <DiscreteFunctionType, LocalOperatorType> InverseOperatorType;
      //typedef BICGSTABOp    <DiscreteFunctionType, LocalOperatorType> InverseOperatorType;
#else 
      typedef CGInverseOp   <DiscreteFunctionType, LocalOperatorType> InverseOperatorType;
#endif
    };
  };

  template< class Model, int polOrd, int passId >
  struct LaplaceDiscreteModel
  : public DGDiscreteModelDefaultWithInsideOutside
    < LaplaceTraits< Model, polOrd, passId >, passId >
  { 
    enum { polynomialOrder = polOrd };

    typedef LaplaceTraits< Model, polOrd, passId > Traits;
    typedef FieldVector<double, Traits::dimDomain> DomainType;
    typedef FieldVector<double, Traits::dimDomain-1> FaceDomainType;

    typedef typename Traits::RangeType RangeType;
    typedef typename Traits::GridType GridType;
    typedef typename Traits::JacobianRangeType JacobianRangeType;
    typedef typename Traits::GridPartType GridPartType; 
    typedef typename GridPartType :: IntersectionIteratorType IntersectionIterator; 
    typedef typename IntersectionIterator :: Intersection  Intersection; 
    typedef typename GridType::template Codim<0>::EntityPointer  EntityPointerType;
    typedef typename GridType::template Codim<0>::Entity EntityType;

    typedef BoundaryIdentifier BoundaryIdentifierType ;
    enum { dimDomain = Traits::dimDomain }; 

  public:
    LaplaceDiscreteModel(const Model& mod) :
      model_(mod)
    {}

    const Model & data () const { return model_; }

    bool constantCoefficient () const { return false; }

    bool hasSource() const { return false; }
    bool hasFlux() const { return false; }
    bool hasCoefficient() const { return true; }
    bool hasRHS() const { return true; }

    // returns true, when Dirichlet boundary, false when other boundary
    // (Neumann) 
    template <class IntersectionIteratorType, class FaceDomType, 
              class ArgumentTuple>
    BoundaryIdentifierType 
    boundaryValue(const IntersectionIteratorType& it,
                  double time, const FaceDomType& x,
                  const ArgumentTuple& u,               
                  RangeType& bndVal) const
    {
      DomainType p = it.geometry().global(x);
      bool dirich = model_.boundaryValue(p,bndVal);

      if(!dirich)
      {
        JacobianRangeType grad;
        model_.neumann(p,grad);
        bndVal = 0.0;
        DomainType normal (it.integrationOuterNormal(x));
        for(int i=0; i<dimDomain; ++i)
        {
          bndVal += grad[0][i] * normal[i];
        }
      }
      return (dirich) ? BoundaryIdentifierType::DirichletNonZero : BoundaryIdentifierType::NeumannNonZero;
    }
    
    template <class ArgumentTuple>
    void analyticalFlux(const EntityType& en,
                        double time, const DomainType& x,
                        const ArgumentTuple& u, JacobianRangeType& f) const
    {
      model_.diffusion(en,time,x,f);
    }

    template <class ArgumentTuple,class CoefficientType>
    void coefficient(const EntityType& en,
                     double time, const DomainType& x,
                     const ArgumentTuple& u, CoefficientType& coeff) const
    {
      model_.diffusion(en,time,x,coeff);
    }

    template <class ArgumentTuple, class CoefficientType> 
    void coefficientFace(const Intersection& it,
                         const double time, 
                         const FaceDomainType& local,
                         const ArgumentTuple& uLeft,
                         const ArgumentTuple& uRight,
                         CoefficientType & coeffLeft, 
                         CoefficientType & coeffRight) const 
    {
      assert( it.neighbor() );
      model_.diffusion( this->inside(), time, it.geometryInInside().global( local ), coeffLeft );
      model_.diffusion( this->outside(), time, it.geometryInOutside().global( local ), coeffRight );
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
  };

  template <class ModelImp, int polOrd, int passId = -1 >
  class VelocityDiscreteModel;

  // DiscreteModelTraits
  template <class ModelImp, int polOrd, int passId = -1 >
  struct VelocityTraits
  {
    enum { myPolOrd = polOrd-1 };

    typedef typename ModelImp::Traits ModelTraits;
    typedef typename ModelTraits::GridType GridType;

    enum { dimRange = ModelTraits::dimGradRange };
    enum { dimDomain = ModelTraits::dimDomain };

    typedef ElliptPassTraits<ModelImp,myPolOrd,dimRange> Traits;
    typedef typename Traits::FunctionSpaceType FunctionSpaceType;
    typedef typename FunctionSpaceType::RangeType RangeType;
    typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;
    typedef typename  FunctionSpaceType:: RangeFieldType RangeFieldType; 
    typedef typename  FunctionSpaceType:: DomainFieldType DomainFieldType; 
    typedef typename Traits::VolumeQuadratureType VolumeQuadratureType;
    typedef typename Traits::FaceQuadratureType FaceQuadratureType;
    typedef typename Traits::GridPartType GridPartType;

    typedef FunctionSpace<DomainFieldType, RangeFieldType, 
             dimDomain, 1 > SingleFunctionSpaceType;
#if defined YASPGRID || defined ALUGRID_CUBE || defined SGRID
    typedef LegendreDiscontinuousGalerkinSpace< SingleFunctionSpaceType,
            GridPartType, polOrd > ContainedFunctionSpaceType;
#else 
    typedef DiscontinuousGalerkinSpace< SingleFunctionSpaceType,
            GridPartType, polOrd > ContainedFunctionSpaceType;
#endif
    typedef typename ContainedFunctionSpaceType :: template ToNewDimRange< dimRange > ::
      Type DiscreteFunctionSpaceType;

    typedef typename ModelTraits::DomainType DomainType;

    //typedef typename Traits::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
    //typedef typename Traits::DiscreteFunctionType DiscreteFunctionType;
    typedef AdaptiveDiscreteFunction< DiscreteFunctionSpaceType > DiscreteFunctionType;
    typedef DiscreteFunctionType DestinationType;


    typedef VelocityDiscreteModel< ModelImp, polOrd, passId > DGDiscreteModelType;
  };


  template <class ModelImp,int polOrd, int passId >
  class VelocityDiscreteModel
  : public DGDiscreteModelDefaultWithInsideOutside
    < VelocityTraits< ModelImp, polOrd, passId >, passId >
  {
    // do not copy this class 
    VelocityDiscreteModel(const VelocityDiscreteModel&);

  public:
    typedef VelocityTraits< ModelImp, polOrd, passId > Traits;
    
    // select Pressure, which comes from pass before 
    typedef FieldVector<double, Traits::dimDomain> DomainType;
    typedef FieldVector<double, Traits::dimDomain-1> FaceDomainType;
    typedef typename Traits::RangeType RangeType;
    typedef typename Traits::GridType GridType;
    typedef typename Traits::JacobianRangeType JacobianRangeType;
    typedef typename Traits::GridPartType::IntersectionIteratorType IntersectionIteratorType;
    typedef typename IntersectionIteratorType :: Intersection Intersection;
    typedef typename GridType::template Codim<0>::Entity EntityType;

    enum { polynomialOrder = polOrd };

  public:
    VelocityDiscreteModel ( const ModelImp &model, bool methodIP )
    : model_( model ),
      methodIP_( methodIP )
    {}

    bool hasSource() const { return false; }
    bool hasFlux() const   { return true; }

    template <class ArgumentTuple>
    double numericalFlux(const Intersection& it,
                         double time, const FaceDomainType& x,
                         const ArgumentTuple& uLeft,
                         const ArgumentTuple& uRight,
                         RangeType& gLeft,
                         RangeType& gRight)
    {
      const DomainType normal = it.integrationOuterNormal(x);

      // get saturation 
      typedef typename ElementType<0, ArgumentTuple>::Type UType;
      const UType& argULeft  = Element<0>::get(uLeft);
      const UType& argURight = Element<0>::get(uRight);

      JacobianRangeType diffmatrix;

      RangeType average ( argULeft );
      average += argURight ;
      average *= 0.5;

      // + n * [ u ] = uL - uR (NIPG, BO)
      if ( ! methodIP_ )
      {
        average += argULeft ;
        average -= argURight ;
      }

      // set diffmatrix 
      model_.gradient( this->inside(),time,
            it.geometryInInside().global(x),
            average,diffmatrix);

      diffmatrix.mv( normal, gLeft );
      gRight = gLeft ;

      return 0.;
    }

    template< class ArgumentTuple >
    double boundaryFlux ( const Intersection &it,
                          double time, const FaceDomainType &x,
                          const ArgumentTuple &uLeft,
                          RangeType &gLeft )
    {
      const DomainType normal = it.integrationOuterNormal(x);
      // get saturation 
      typedef typename ElementType<0, ArgumentTuple>::Type SType;
      const SType& argSLeft  = Element<0>::get(uLeft);

      JacobianRangeType diffmatrix;

      {
        model_.gradient( this->inside(),time,
            it.geometryInInside().global(x),
            argSLeft,diffmatrix);
      }

      diffmatrix.mv(normal,gLeft);
      return 0.0;
    }

    template< class ArgumentTuple >
    void analyticalFlux ( const EntityType &entity,
                          double time, const DomainType &x,
                          const ArgumentTuple &u, JacobianRangeType &f )
    {
      typedef typename ElementType< 0, ArgumentTuple >::Type UType;
      const UType &argU = Element< 0 >::get( u );
      // get saturation 
      model_.gradient( entity, time, x, argU, f );
    }

  private:
    const ModelImp &model_;
    const bool methodIP_;
  };


}
#endif
