#ifndef DUNE_EXAMPLEDISCRETEMODELS_HH
#define DUNE_EXAMPLEDISCRETEMODELS_HH

// FR passes 
#include <dune/fem/pass/dgpass.hh>
#include <dune/fem/pass/selection.hh>
#include <dune/fem/pass/discretemodel.hh>
#include <dune/fem/space/dgspace.hh>
#include <dune/fem/space/combinedspace.hh>

// Dune includes
#include <dune/common/utility.hh>
#include <dune/fem/space/combinedspace.hh>
#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/gridpart/gridpart.hh>
#include <dune/fem/quadrature/cachequad.hh>

#include <dune/fem/misc/boundaryidentifier.hh>

#include <dune/fem/operator/matrix/matrixhandler.hh>
#if HAVE_BLAS
#include <dune/fem/solver/oemsolver/oemsolver.hh>
#endif


#include <dune/fem/operator/matrix/istlmatrix.hh>
#include <dune/fem/solver/istlsolver.hh>

#include <dune/fem/solver/inverseoperators.hh>

#include <dune/fem/space/dgspace/dgadaptiveleafgridpart.hh>
#include <dune/fem/space/common/adaptiveleafgridpart.hh>

#include <dune/fem/operator/2order/dgprimaloperator.hh>

#include <dune/fem/pass/dgelliptpass.hh>

#include <dune/fem/misc/double.hh>

#if HAVE_DUNE_ISTL
#include <dune/istl/bvector.hh>
#include <dune/fem/function/blockvectorfunction.hh>
#endif


#define USE_DUNE_ISTL HAVE_DUNE_ISTL
//#define USE_DUNE_ISTL 0

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
    
    //typedef LeafGridPart<GridType> GridPartType;
    typedef DGAdaptiveLeafGridPart<GridType> GridPartType;

    typedef CachingQuadrature<GridPartType,0> VolumeQuadratureType;
    typedef CachingQuadrature<GridPartType,1> FaceQuadratureType;
    
    //typedef CachingQuadrature<GridPartType,0> VolumeQuadratureType;
    //typedef ElementQuadrature<GridPartType,1> FaceQuadratureType;
    
    // typical tpye of space 
    typedef FunctionSpace< typename Model::DomainFieldType, typename Model::RangeFieldType , 
                           dimDomain, dimRange > FunctionSpaceType; 
#if 0
    typedef FunctionSpace< typename Model::DomainFieldType, typename Model::RangeFieldType , 
                           dimDomain, 1 > SingleFunctionSpaceType; 
    typedef DiscontinuousGalerkinSpace<SingleFunctionSpaceType, GridPartType, polOrd, CachingStorage > ContainedSpaceType;
    typedef CombinedSpace< ContainedSpaceType, dimRange, PointBased > DiscreteFunctionSpaceType;
#else 
    typedef DiscontinuousGalerkinSpace<FunctionSpaceType, GridPartType, polOrd, CachingStorage > ContainedSpaceType;
    typedef ContainedSpaceType DiscreteFunctionSpaceType;
#endif
    
    typedef AdaptiveDiscreteFunction<DiscreteFunctionSpaceType> DiscreteFunctionType;
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
    typedef typename Traits::GridPartType GridPartType; 
    typedef typename GridPartType :: IntersectionIteratorType IntersectionIterator; 
    typedef typename GridType::template Codim<0>::EntityPointer  EntityPointerType;
    typedef typename GridType::template Codim<0>::Entity EntityType;

    typedef BoundaryIdentifier BoundaryIdentifierType ;

    enum { dimRange = Traits::dimRange };

  public:
    GradientDiscreteModel(const Model& mod,const NumFlux& numf) :
      model_(mod),
      numflux_(numf),
      one_(0.0)
    {
      for(int i=0; i<dimRange; ++i) one_[i][i] = 1.0;
    }

    const Model & data () const { return model_; }

    bool hasSource() const { return true; }
    bool hasFlux() const { return false; }

    template <class ArgumentTuple> 
    double numericalFlux(const IntersectionIterator& it,
                         const double time, 
                         const FaceDomainType& local,
                         const ArgumentTuple& uLeft,
                         const ArgumentTuple& uRight,
                         RangeType & sigmaLeft, 
                         RangeType & sigmaRight,
                         RangeType & gLeft,
                         RangeType & gRight) const 
    {
      // calculate unit normal and face volume  
      DomainType unitNormal = it.integrationOuterNormal(local);
      const double faceVol = unitNormal.two_norm();
      unitNormal *= 1.0/faceVol;

      if( hasFlux() )
      {
        DomainType dom = it.intersectionGlobal().global(local);
        // default is id matrix 
        JacobianRangeType anaFluxLeft(one_);
        JacobianRangeType anaFluxRight(one_);

        RangeType tmpLeft(0.0);
        RangeType tmpRight(0.0);
        
        {
          EntityPointerType ep = it.inside();
          const EntityType & en = *ep;
          analyticalFlux(en,time,dom,
                         uLeft,anaFluxLeft);
        }
                      
        {
          assert( it.neighbor() );
          EntityPointerType ep = it.outside();
          const EntityType & en = *ep;
          analyticalFlux(en,time,dom,
                         uRight,anaFluxRight);
        }

        anaFluxLeft.umv(unitNormal,tmpLeft); 
        anaFluxRight.umv(unitNormal,tmpRight); 
        
        // evaluate u part of ldg flux 
        numflux_.uFlux(unitNormal,faceVol,
                       tmpLeft,tmpRight, 
                       sigmaLeft,sigmaRight,
                       gLeft,gRight);
      }
      else 
      {
        // evaluate u part of ldg flux 
        numflux_.uFlux(unitNormal,faceVol,
                       unitNormal,unitNormal, 
                       sigmaLeft,sigmaRight,
                       gLeft,gRight);
      }
      return 0.0;
    }

    template <class ArgumentTuple>
    void analyticalFlux(const EntityType& en,
                        const double time, 
                        const DomainType& local,
                        const ArgumentTuple& u, 
                        JacobianRangeType& f) const
    {
      model_.diffusion(en,time,local,f);
      double tmp = f[0][0];
      tmp = sqrt(tmp);

      for(int i=0;i<dimRange; ++i) f[i][i] = tmp;
      
      //std::cout << f << " f \n";
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
              class ArgumentTuple>
    BoundaryIdentifierType 
    boundaryValue(const IntersectionIteratorType& it,
                  double time, const FaceDomType& x,
                  const ArgumentTuple& u,               
                  RangeType& bndVal) const
    {
      return BoundaryIdentifierType::undefined; 
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
    JacobianRangeType one_;
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

    typedef LaplaceDiscreteModel<Model,NumFlux,polOrd> DiscreteModelType;
    typedef DiscreteModelType ThisType;

    template <class PreviousPassType>
    struct LocalOperatorSelector
    {
      //! select one pass before fake gradient pass 
      typedef typename PreviousPassType :: PreviousPassType ElliptPrevPassType;
        
      // get type of gradient pass 
      typedef typename PreviousPassType :: DiscreteFunctionSpaceType
        PrevSpaceType;

      // old type 
      typedef MatrixHandlerSPMat<DiscreteFunctionSpaceType,
                                 PrevSpaceType> MatrixHandlerType; 
      
      // new type 
#if USE_DUNE_ISTL                                
      typedef ISTLMatrixTraits<DiscreteFunctionSpaceType> MatrixObjectTraits;
#else 
      //typedef SparseRowMatrixTraits<DiscreteFunctionSpaceType,
      //                              DiscreteFunctionSpaceType> MatrixObjectTraits; 
      typedef BlockMatrixTraits<DiscreteFunctionSpaceType,
                                DiscreteFunctionSpaceType> MatrixObjectTraits; 
#endif
      typedef DGPrimalOperator<ThisType,PreviousPassType,
                ElliptPrevPassType, MatrixObjectTraits> LocalOperatorType;

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

  template <class Model,class NumFlux,int polOrd>
  class LaplaceDiscreteModel : 
    public DiscreteModelDefaultWithInsideOutSide<LaplaceTraits<Model,NumFlux,polOrd> >
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
    typedef typename Traits::GridPartType GridPartType; 
    typedef typename GridPartType :: IntersectionIteratorType IntersectionIterator; 
    typedef typename GridType::template Codim<0>::EntityPointer  EntityPointerType;
    typedef typename GridType::template Codim<0>::Entity EntityType;

    typedef BoundaryIdentifier BoundaryIdentifierType ;
    enum { dimDomain = Traits::dimDomain }; 

  public:
    LaplaceDiscreteModel(const Model& mod,const NumFlux& numf) :
      model_(mod),
      numflux_(numf)
    {}

    const Model & data () const { return model_; }

    bool hasSource() const { return false; }
    bool hasFlux() const { return false; }
    bool hasCoefficient() const { return true; }
    bool hasRHS() const { return true; }

    template <class ArgumentTuple> 
    double numericalFlux(const IntersectionIterator& it,
                         const double time, 
                         const FaceDomainType& local,
                         const ArgumentTuple& uLeft,
                         const ArgumentTuple& uRight,
                         RangeType & sigmaLeft, 
                         RangeType & sigmaRight,
                         RangeType & gLeft,
                         RangeType & gRight) const 
    {
      // calculate unit normal and face volume  
      DomainType unitNormal = it.integrationOuterNormal(local);
      const double faceVol = unitNormal.two_norm();
      unitNormal *= 1.0/faceVol;

      RangeType tmpLeft(1.0);
      RangeType tmpRight(1.0);

      if( hasFlux() )
      {
        DomainType dom = it.intersectionGlobal().global(local);
        // default is id matrix 
        JacobianRangeType anaFluxLeft;
        anaFluxLeft[0] = unitNormal;
        JacobianRangeType anaFluxRight;
        anaFluxRight[0] = unitNormal;

        {
          EntityPointerType ep = it.inside();
          const EntityType & en = *ep;
          analyticalFlux(en,time,dom,
                         uLeft,anaFluxLeft);
          anaFluxLeft.umv(unitNormal,tmpLeft); 
        }
                      
        if( it.neighbor() )
        {
          EntityPointerType ep = it.outside();
          const EntityType & en = *ep;
          analyticalFlux(en,time,dom,
                         uRight,anaFluxRight);
          anaFluxRight.umv(unitNormal,tmpRight); 
        }
        
        tmpLeft  = unitNormal * anaFluxLeft[0];
        tmpRight = unitNormal * anaFluxRight[0];
      }
      
      numflux_.sigmaFlux(unitNormal,faceVol,
                         tmpLeft,tmpRight, 
                         sigmaLeft,sigmaRight,
                         gLeft,gRight);
      return 0.0;
    }


    template <class ArgumentTuple> 
    double boundaryFlux(const IntersectionIterator& it,
                        const double time, 
                        const FaceDomainType& local,
                        const ArgumentTuple& uLeft,
                        RangeType & sigmaLeft ) const 
                      //, RangeType & gLeft) const 
    {
      // calculate unit normal and face volume  
      DomainType unitNormal = it.integrationOuterNormal(local);
      const double faceVol = unitNormal.two_norm();
      unitNormal *= 1.0/faceVol;

      RangeType sigmaRight;
      RangeType gRight; 

      RangeType tmpLeft(1.0);
      
      if( hasFlux() )
      {
        DomainType dom = it.intersectionGlobal().global(local);
        // default is id matrix 
        JacobianRangeType anaFluxLeft;
        anaFluxLeft[0] = unitNormal;

        {
          EntityPointerType ep = it.inside();
          const EntityType & en = *ep;
          analyticalFlux(en,time,dom,
                         uLeft,anaFluxLeft);
          anaFluxLeft.umv(unitNormal,tmpLeft); 
        }
                      
        tmpLeft  = unitNormal * anaFluxLeft[0];
      }
          
      // don't apply beta stabilization at boundary 
      numflux_.sigmaFluxBetaZero(unitNormal,faceVol,
                                 tmpLeft,tmpLeft, 
                                 sigmaLeft,sigmaRight);
      
      return 0.0;
    }

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
      DomainType p = it.intersectionGlobal().global(x);
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
    void coefficientFace(const IntersectionIterator& it,
                         const double time, 
                         const FaceDomainType& local,
                         const ArgumentTuple& uLeft,
                         const ArgumentTuple& uRight,
                         CoefficientType & coeffLeft, 
                         CoefficientType & coeffRight) const 
    {
      assert( it.neighbor() );
      model_.diffusion(this->inside(),time,it.intersectionSelfLocal().global(local),coeffLeft);
      model_.diffusion(this->outside(),time,it.intersectionNeighborLocal().global(local),coeffRight);
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
  };

  template <class ModelImp, class NumFluxImp, int polOrd>
  class VelocityDiscreteModel;

  // DiscreteModelTraits
  template <class ModelImp,class NumFluxImp, int polOrd >
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
    typedef CombinedSpace< ContainedFunctionSpaceType, dimRange >
      DiscreteFunctionSpaceType;

    typedef typename ModelTraits::DomainType DomainType;

    //typedef typename Traits::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
    //typedef typename Traits::DiscreteFunctionType DiscreteFunctionType;
    typedef AdaptiveDiscreteFunction< DiscreteFunctionSpaceType > DiscreteFunctionType;
    typedef DiscreteFunctionType DestinationType;


    typedef VelocityDiscreteModel<ModelImp,NumFluxImp,polOrd> DiscreteModelType;
  };
  template <class ModelImp,class NumFluxImp,int polOrd>
  class VelocityDiscreteModel :
    public DiscreteModelDefaultWithInsideOutSide<VelocityTraits<ModelImp,NumFluxImp,polOrd> >
  {
    // do not copy this class 
    VelocityDiscreteModel(const VelocityDiscreteModel&);
  public:
    typedef VelocityTraits<ModelImp,NumFluxImp,polOrd> Traits;

    // select Pressure, which comes from pass before 
    typedef Dune::Selector<1> SelectorType;
    typedef FieldVector<double, Traits::dimDomain> DomainType;
    typedef FieldVector<double, Traits::dimDomain-1> FaceDomainType;
    typedef typename Traits::RangeType RangeType;
    typedef typename Traits::GridType GridType;
    typedef typename Traits::JacobianRangeType JacobianRangeType;
    typedef typename Traits::GridPartType::IntersectionIteratorType IntersectionIteratorType;
    typedef typename GridType::template Codim<0>::Entity EntityType;
    typedef NumFluxImp NumFluxType;

    enum { polynomialOrder = polOrd };

  public:
    VelocityDiscreteModel(const ModelImp& mod, const NumFluxType& flux) 
      : model_(mod)
      , numFlux_(flux) 
    {}

    bool hasSource() const { return false; }
    bool hasFlux() const   { return true; }

    template <class ArgumentTuple>
    double numericalFlux(IntersectionIteratorType& it,
                         double time, const FaceDomainType& x,
                         const ArgumentTuple& uLeft,
                         const ArgumentTuple& uRight,
                         RangeType& gLeft,
                         RangeType& gRight)
    {
      const DomainType normal = it.integrationOuterNormal(x);
      const double faceVol = normal.two_norm();

      // get saturation 
      typedef typename ElementType<0, ArgumentTuple>::Type UType;
      const UType& argULeft  = Element<0>::get(uLeft);
      const UType& argURight = Element<0>::get(uRight);

      JacobianRangeType diffmatrix;

      RangeType diffflux(0.);
      gLeft  = 0.0;
      gRight = 0.0;

      UType result;
        
      // left value 
      {
        // eval num flux 
        numFlux_.uFlux(faceVol,argULeft,argURight,result);

        // set diffmatrix 
        model_.gradient( this->inside(),time,
              it.intersectionSelfLocal().global(x),
              result,diffmatrix);

        diffmatrix.umv(normal,gLeft);
      }
      
      {
        // eval num flux 
        numFlux_.uFlux(faceVol,argURight,argULeft,result);
        
        // set diffmatrix 
        model_.gradient( this->inside(),time,
              it.intersectionSelfLocal().global(x),
              result,diffmatrix);

        diffmatrix.umv(normal,gRight);
      }
      return 0.;
    }

    template <class ArgumentTuple>
    double boundaryFlux(IntersectionIteratorType& it,
                        double time, const FaceDomainType& x,
                        const ArgumentTuple& uLeft,
                        RangeType& gLeft)
    {
      const DomainType normal = it.integrationOuterNormal(x);
      // get saturation 
      typedef typename ElementType<0, ArgumentTuple>::Type SType;
      const SType& argSLeft  = Element<0>::get(uLeft);

      JacobianRangeType diffmatrix;
      gLeft = 0.0;

      {
        model_.gradient( this->inside(),time,
            it.intersectionSelfLocal().global(x),
            argSLeft,diffmatrix);
      }

      diffmatrix.umv(normal,gLeft);
      return 0.0;
    }

    template <class ArgumentTuple>
    void analyticalFlux(EntityType& en,
                        double time, const DomainType& x,
                        const ArgumentTuple& u, JacobianRangeType& f)
    {
      typedef typename ElementType<0, ArgumentTuple>::Type UType;
      const UType& argU = Element<0>::get(u);
      // get saturation 
      model_.gradient(en,time,x,argU,f);
    }
  private:
    const ModelImp & model_;
    const NumFluxType& numFlux_;
  };


}
#endif
