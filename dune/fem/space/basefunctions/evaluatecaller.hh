#ifndef DUNE_EVALUATECALLER_HH
#define DUNE_EVALUATECALLER_HH

#include <cstdlib>
#include <vector>

////////////////////////////////////////////
//
// pre-define these values for faster compilation 
//
////////////////////////////////////////////
#ifndef MAX_NUMBER_OF_QUAD_POINTS
#define MAX_NUMBER_OF_QUAD_POINTS 20 
#endif

#ifndef MAX_NUMBER_OF_BASE_FCT 
#define MAX_NUMBER_OF_BASE_FCT 10 
#endif

#ifndef MIN_NUMBER_OF_QUAD_POINTS
#define MIN_NUMBER_OF_QUAD_POINTS 1 
#endif

#ifndef MIN_NUMBER_OF_BASE_FCT 
#define MIN_NUMBER_OF_BASE_FCT 10
#endif

namespace Dune
{

namespace Fem 
{
  // empty class for specialization of evaluation classes in basefunctionsets.hh 
  class EmptyGeometry {};

  // forward declaration 
  template <class Traits,
            int quadNop,
            int numBaseFct >
  class EvaluateCaller;

  template< class QuadratureImp, 
            class FactorImp, 
            class LocalDofVectorImp,
            class GeometryImp = EmptyGeometry > 
  struct EvaluateCallerInterfaceTraits 
  {
    typedef QuadratureImp     QuadratureType;
    typedef FactorImp         FactorType;
    typedef LocalDofVectorImp LocalDofVectorType;
    typedef GeometryImp       Geometry;        
  };

  template <class Traits,
            class BaseFunctionSet,
            class RangeVectorImp> 
  struct EvaluateCallerTraits 
  {
    typedef Traits  BaseTraits;
    typedef typename Traits :: QuadratureType      QuadratureType ;
    typedef typename Traits :: FactorType          FactorType ;
    typedef typename Traits :: LocalDofVectorType  LocalDofVectorType ;
    typedef typename Traits :: Geometry            Geometry ;

    typedef BaseFunctionSet   BaseFunctionSetType;
    typedef RangeVectorImp    RangeVectorType;
  };

  //- base function evaluation interface 
  template <class Traits>  
  class EvaluateCallerInterface
  {
    typedef EvaluateCallerInterface< Traits >  ThisType;
  protected:
    enum { maxNumBaseFunctions = MAX_NUMBER_OF_BASE_FCT };
    enum { maxQuadratures = 50 };
    enum { maxQuadNop = MAX_NUMBER_OF_QUAD_POINTS };

    class EvaluatorStorage
    {
    protected:  
      typedef ThisType* ValueType ;
      std::vector< ValueType > storage_;
    public:
      EvaluatorStorage() : 
        storage_( maxQuadratures, (ThisType* ) 0 ) 
      {}

      ~EvaluatorStorage() 
      {
        for( size_t i = 0; i<storage_.size(); ++i ) 
          delete storage_[ i ];
      }

      ValueType& operator [] ( const int i ) { return storage_[ i ]; }
      const ValueType& operator [] ( const int i ) const { return storage_[ i ]; }
    };


    EvaluateCallerInterface() {}

  public:
    typedef typename Traits :: QuadratureType      QuadratureType ;
    typedef typename Traits :: FactorType          FactorType ;
    typedef typename Traits :: LocalDofVectorType  LocalDofVectorType ;
    typedef typename Traits :: Geometry            Geometry ;

    virtual ~EvaluateCallerInterface() {}

    virtual void axpyRanges( const QuadratureType&,
                             const FactorType& ,
                             LocalDofVectorType & ) const = 0;

    virtual void evaluateRanges( const QuadratureType& quad,
                                 const LocalDofVectorType & dofs,
                                 FactorType& factors) const = 0;

    virtual void axpyJacobians( const QuadratureType&,
                                const Geometry&,
                                const FactorType& ,
                                LocalDofVectorType & ) const = 0;

    virtual void evaluateJacobians( const QuadratureType&,
                                    const Geometry&,
                                    const LocalDofVectorType&,
                                    FactorType&) const = 0;

    template < class BaseFunctionSet, class Storage >
    static const ThisType& storage(const BaseFunctionSet& baseSet,
                                   const Storage& dataCache,
                                   const QuadratureType& quad ) 
    {
      // assert that max numbers are big enough 
      assert( baseSet.numDifferentBaseFunctions() <= maxNumBaseFunctions );
      assert( quad.nop()  <= maxQuadNop );
      assert( quad.id()   < maxQuadratures );

      // static vector holding all evaluator instances 
      static EvaluatorStorage evaluators; 

      // check if object already created 
      const size_t quadId = quad.id();
      if( evaluators[ quadId ] == 0 )
      {
        typedef EvaluateCallerTraits< Traits, BaseFunctionSet, Storage> NewTraits;
        // create appropriate evaluator 
        evaluators[ quadId ] = 
          EvaluateCaller< NewTraits, maxQuadNop, maxNumBaseFunctions > 
            :: create( dataCache , quad.nop(), baseSet.numDifferentBaseFunctions() );
      }

      // return reference to evaluator 
      return *(evaluators[ quadId ]);
    }
  };

  template <class Traits,
            int quadNop,
            int numBaseFct >
  class EvaluateCallerDefault 
    : public EvaluateCallerInterface< typename Traits :: BaseTraits >
  {
  protected:  
    typedef typename Traits :: BaseFunctionSetType BaseFunctionSetType;
    typedef typename Traits :: QuadratureType      QuadratureType ;
    typedef typename Traits :: FactorType          FactorType ;
    typedef typename Traits :: LocalDofVectorType  LocalDofVectorType ;
    typedef typename Traits :: Geometry            Geometry ;
    typedef typename Traits :: RangeVectorType     RangeVectorType ;

    enum { dimRange = BaseFunctionSetType :: dimRange };
    typedef EvaluateCallerInterface< typename Traits :: BaseTraits >   BaseType;
    // type of interface class
    typedef BaseType InterfaceType;

    typedef EvaluateCallerDefault< Traits, quadNop, numBaseFct > ThisType;
    const RangeVectorType& rangeStorage_;
  public:  
    EvaluateCallerDefault( const RangeVectorType& rangeStorage ) 
      : rangeStorage_( rangeStorage )
    {}

    virtual void axpyRanges( const QuadratureType& quad,
                             const FactorType& rangeFactors,
                             LocalDofVectorType & dofs ) const
    {
      BaseFunctionSetType :: template AxpyRanges 
        < BaseFunctionSetType, Geometry, dimRange, quadNop, numBaseFct > :: axpy 
        ( quad, rangeStorage_, rangeFactors, dofs );
    }

    virtual void axpyJacobians( const QuadratureType& quad,
                                const Geometry& geometry,
                                const FactorType& jacFactors,
                                LocalDofVectorType& dofs) const 
    {
      BaseFunctionSetType :: template AxpyJacobians 
        < BaseFunctionSetType, Geometry, dimRange, quadNop, numBaseFct > :: axpy 
        ( quad, geometry, rangeStorage_, jacFactors, dofs );
    }

    virtual void evaluateRanges( const QuadratureType& quad,
                                 const LocalDofVectorType & dofs,
                                 FactorType& rangeFactors) const
    {
      BaseFunctionSetType :: template EvaluateRanges
        < BaseFunctionSetType, Geometry, dimRange, quadNop, numBaseFct >
        :: eval ( quad, rangeStorage_, dofs, rangeFactors );
    }

    virtual void evaluateJacobians( const QuadratureType& quad,
                                    const Geometry& geometry,
                                    const LocalDofVectorType& dofs,
                                    FactorType& jacFactors) const 
    {
      BaseFunctionSetType :: template EvaluateJacobians 
        < BaseFunctionSetType, Geometry, dimRange, quadNop, numBaseFct > :: eval 
        ( quad, geometry, rangeStorage_, dofs, jacFactors );
    } 

  };
   
  template <class Traits,
            int quadNop,
            int numBaseFct >
  class EvaluateCaller 
    : public EvaluateCallerDefault< Traits, quadNop, numBaseFct >
  {
  protected:  
    typedef typename Traits :: RangeVectorType     RangeVectorType ;
    typedef EvaluateCallerDefault< Traits, quadNop, numBaseFct >   BaseType;
    typedef EvaluateCaller< Traits, quadNop, numBaseFct > ThisType;
    typedef typename BaseType :: InterfaceType InterfaceType;
  public:  
    explicit EvaluateCaller( const RangeVectorType& rangeStorage ) 
      : BaseType( rangeStorage )
    {}

    static InterfaceType* createObj( const RangeVectorType& rangeStorage, 
                                const size_t numbase ) 
    {
      if( numBaseFct == numbase ) 
        return new ThisType( rangeStorage );
      else 
        return EvaluateCaller< Traits, quadNop, numBaseFct - 1 > :: createObj( rangeStorage, numbase );
    }

    static InterfaceType* create( const RangeVectorType& rangeStorage, 
                             const size_t quadnop, const size_t numbase ) 
    {
      if( quadNop == quadnop ) 
        return EvaluateCaller< Traits, quadNop, numBaseFct > :: createObj( rangeStorage, numbase );
      else  
        return EvaluateCaller< Traits, quadNop - 1, numBaseFct > :: create( rangeStorage, quadnop, numbase );
    }
  };
   
  template <class Traits,
            int numBaseFct >
  class EvaluateCaller< Traits, MIN_NUMBER_OF_QUAD_POINTS, numBaseFct >
    : public EvaluateCallerDefault< Traits, MIN_NUMBER_OF_QUAD_POINTS, numBaseFct >
  {
  protected:  
    enum { quadNop = MIN_NUMBER_OF_QUAD_POINTS };
    typedef typename Traits :: RangeVectorType     RangeVectorType ;
    typedef EvaluateCallerDefault< Traits, quadNop, numBaseFct >   BaseType;
    typedef typename BaseType :: InterfaceType InterfaceType;
    typedef EvaluateCaller< Traits, quadNop, numBaseFct > ThisType;
  public:  
    explicit EvaluateCaller( const RangeVectorType& rangeStorage ) 
      : BaseType( rangeStorage )
    {}

    static InterfaceType* createObj( const RangeVectorType& rangeStorage, 
                                const size_t numbase ) 
    {
      if( numBaseFct == numbase ) 
        return new ThisType( rangeStorage );
      else 
        return EvaluateCaller< Traits, quadNop, numBaseFct - 1 > :: createObj( rangeStorage, numbase );
    }

    static InterfaceType* create( const RangeVectorType& rangeStorage, 
                             const size_t quadnop, const size_t numbase ) 
    {
      if( quadNop == quadnop ) 
        return EvaluateCaller< Traits, quadNop, numBaseFct > :: createObj( rangeStorage, numbase );
      else  
      {
        abort();
      }
    }
  };
   
  template <class Traits,
            int quadNop>
  class EvaluateCaller< Traits, quadNop, MIN_NUMBER_OF_BASE_FCT >
    : public EvaluateCallerDefault< Traits, quadNop, MIN_NUMBER_OF_BASE_FCT >
  {
  protected:  
    enum { numBaseFct = MIN_NUMBER_OF_BASE_FCT };
    typedef typename Traits :: RangeVectorType     RangeVectorType ;
    typedef EvaluateCallerDefault< Traits, quadNop, numBaseFct >   BaseType;
    typedef EvaluateCaller< Traits, quadNop, numBaseFct > ThisType;
    typedef typename BaseType :: InterfaceType InterfaceType;
  public:  
    explicit EvaluateCaller( const RangeVectorType& rangeStorage ) 
      : BaseType( rangeStorage )
    {}

    static InterfaceType* createObj( const RangeVectorType& rangeStorage, 
                                const size_t numbase ) 
    {
      if( numBaseFct == numbase ) 
        return new ThisType( rangeStorage );
      else 
        abort();
    }

    static InterfaceType* create( const RangeVectorType& rangeStorage, 
                             const size_t quadnop, const size_t numbase ) 
    {
      if( quadNop == quadnop ) 
        return EvaluateCaller< Traits, quadNop, numBaseFct > :: createObj( rangeStorage, numbase );
      else  
      {
        return EvaluateCaller< Traits, quadNop - 1, numBaseFct > :: create( rangeStorage, quadnop, numbase );
      }
    }
  };
   
  template <class Traits>
  class EvaluateCaller< Traits, MIN_NUMBER_OF_QUAD_POINTS, MIN_NUMBER_OF_BASE_FCT> 
    : public EvaluateCallerDefault< Traits, MIN_NUMBER_OF_QUAD_POINTS, MIN_NUMBER_OF_BASE_FCT >
  {
  protected:  
    enum { quadNop = MIN_NUMBER_OF_QUAD_POINTS };
    enum { numBaseFct = MIN_NUMBER_OF_BASE_FCT };
    typedef typename Traits :: RangeVectorType     RangeVectorType ;
    typedef EvaluateCallerDefault< Traits, quadNop, numBaseFct >   BaseType;
    typedef EvaluateCaller< Traits, quadNop, numBaseFct > ThisType;
    typedef typename BaseType :: InterfaceType InterfaceType;
  public:  
    explicit EvaluateCaller( const RangeVectorType& rangeStorage ) 
      : BaseType( rangeStorage )
    {}

    static InterfaceType* createObj( const RangeVectorType& rangeStorage, 
                                const size_t numbase ) 
    {
      if( numBaseFct == numbase ) 
        return new ThisType( rangeStorage );
      else 
        abort();
    }

    static InterfaceType* create( const RangeVectorType& rangeStorage, 
                             const size_t quadnop, const size_t numbase ) 
    {
      if( quadNop == quadnop ) 
        return EvaluateCaller< Traits, quadNop, numBaseFct > :: createObj( rangeStorage, numbase );
      else  
      {
        abort();
      }
    }
  };

} // end namespace Fem 

} // end namespace Dune
#endif
