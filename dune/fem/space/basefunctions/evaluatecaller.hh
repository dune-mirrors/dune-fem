#ifndef DUNE_EVALUATECALLER_HH
#define DUNE_EVALUATECALLER_HH

#include <vector>

namespace Dune
{

namespace Fem 
{

  // forward declaration 
  template <class BaseFunctionSet,
            class QuadratureImp, 
            class RangeVectorType,
            class RangeFactorType, 
            class LocalDofVectorType,
            int quadNop,
            int numBaseFct >
  class EvaluateCaller;

  //- base function evaluation interface 
  template <class QuadratureImp, 
            class RangeFactorType, 
            class LocalDofVectorType>
  class EvaluateCallerInterface
  {
    typedef EvaluateCallerInterface
      < QuadratureImp, RangeFactorType, LocalDofVectorType > ThisType;
  protected:
    enum { maxNumBaseFunctions = 10 };
    enum { maxQuadratures = 100 };
    enum { maxQuadNop = 20 };

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
    typedef QuadratureImp QuadratureType ;

    virtual ~EvaluateCallerInterface() {}

    virtual void axpyRanges( const QuadratureType&,
                             const RangeFactorType& ,
                             LocalDofVectorType & ) const = 0;

    virtual void evaluateRanges( const QuadratureType& quad,
                                 const LocalDofVectorType & dofs,
                                 RangeFactorType& rangeFactors) const = 0;

    template < class BaseFunctionSet >
    static const ThisType& storage(const BaseFunctionSet& baseSet,
                                   const QuadratureType& quad, 
                                   const size_t numDiffBase ) 
    {
      assert( numDiffBase < maxNummBaseFunctions );
      assert( quad.nop()  < maxQuadNop );
      assert( quad.id()   < maxQuadratures );

      static EvaluatorStorage evaluators; 

      typedef typename BaseFunctionSet :: StorageType :: RangeVectorType RangeVectorType;

      const size_t quadId = quad.id();
      if( evaluators[ quadId ] == 0 )
      {
        evaluators[ quadId ] = 
          EvaluateCaller< BaseFunctionSet, 
                                      QuadratureType, 
                                      RangeVectorType,
                                      RangeFactorType, 
                                      LocalDofVectorType, maxQuadNop, maxNumBaseFunctions > 
            :: create( baseSet.storage().getRangeStorage( quad ), 
                       quad.nop(),
                       baseSet.numDifferentBaseFunctions() );
      }

      return *(evaluators[ quadId ]);
    }
  };

  template <class BaseFunctionSet,
            class QuadratureImp, 
            class RangeVectorType,
            class RangeFactorType, 
            class LocalDofVectorType,
            int quadNop,
            int numBaseFct >
  class EvaluateCaller : public EvaluateCallerInterface
         < QuadratureImp, RangeFactorType, LocalDofVectorType >
  {
    typedef EvaluateCallerInterface
               < QuadratureImp, RangeFactorType, LocalDofVectorType >  BaseType;

    typedef EvaluateCaller< BaseFunctionSet, QuadratureImp, RangeVectorType, RangeFactorType,
            LocalDofVectorType , quadNop, numBaseFct > ThisType;
    const RangeVectorType& rangeStorage_;
  public:  
    EvaluateCaller( const RangeVectorType& rangeStorage ) 
      : rangeStorage_( rangeStorage )
    {}

    virtual void axpyRanges( const QuadratureImp& quad,
                             const RangeFactorType& rangeFactors,
                             LocalDofVectorType & dofs ) const
    {
      BaseFunctionSet ::template 
        doAxpyRanges < quadNop, numBaseFct > 
        ( quad, rangeStorage_, rangeFactors, dofs );
    }

    virtual void evaluateRanges( const QuadratureImp& quad,
                                 const LocalDofVectorType & dofs,
                                 RangeFactorType& rangeFactors) const
    {
      BaseFunctionSet ::template 
        doEvaluateRanges < quadNop, numBaseFct > 
        ( quad, rangeStorage_, dofs, rangeFactors );
    }

    static BaseType* createObj( const RangeVectorType& rangeStorage, 
                                const size_t numbase ) 
    {
      if( numBaseFct == numbase ) 
        return new ThisType( rangeStorage );
      else 
        return EvaluateCaller< BaseFunctionSet, QuadratureImp, RangeVectorType, RangeFactorType,
                    LocalDofVectorType , quadNop, numBaseFct - 1 > :: createObj( rangeStorage, numbase );
    }

    static BaseType* create( const RangeVectorType& rangeStorage, 
                             const size_t quadnop, const size_t numbase ) 
    {
      if( quadNop == quadnop ) 
        return EvaluateCaller< BaseFunctionSet, QuadratureImp, RangeVectorType, RangeFactorType,
                    LocalDofVectorType , quadNop, numBaseFct > :: createObj( rangeStorage, numbase );
      else  
        return EvaluateCaller< BaseFunctionSet, QuadratureImp, RangeVectorType, RangeFactorType,
                    LocalDofVectorType , quadNop - 1, numBaseFct > :: create( rangeStorage, quadnop, numbase );
    }
  };
   
  template <class BaseFunctionSet,
            class QuadratureImp, 
            class RangeVectorType,
            class RangeFactorType, 
            class LocalDofVectorType,
            int numBaseFct >
  class EvaluateCaller
    < BaseFunctionSet, QuadratureImp, RangeVectorType, RangeFactorType,
  LocalDofVectorType, 1, numBaseFct> : public EvaluateCallerInterface
         < QuadratureImp, RangeFactorType, LocalDofVectorType >
  {
    enum { quadNop = 1 };
    typedef EvaluateCallerInterface
               < QuadratureImp, RangeFactorType, LocalDofVectorType > BaseType;
    typedef EvaluateCaller< BaseFunctionSet, QuadratureImp, RangeVectorType, RangeFactorType,
            LocalDofVectorType , quadNop, numBaseFct > ThisType;
    const RangeVectorType& rangeStorage_;
  public:  
    EvaluateCaller( const RangeVectorType& rangeStorage ) 
      : rangeStorage_( rangeStorage )
    {}

    virtual void axpyRanges( const QuadratureImp& quad,
                             const RangeFactorType& rangeFactors,
                             LocalDofVectorType & dofs ) const
    {
      BaseFunctionSet ::template 
        doAxpyRanges < quadNop, numBaseFct > 
        ( quad, rangeStorage_, rangeFactors, dofs );
    }

    virtual void evaluateRanges( const QuadratureImp& quad,
                                 const LocalDofVectorType & dofs,
                                 RangeFactorType& rangeFactors) const
    {
      BaseFunctionSet ::template 
        doEvaluateRanges < quadNop, numBaseFct > 
        ( quad, rangeStorage_, dofs, rangeFactors );
    }

    static BaseType* createObj( const RangeVectorType& rangeStorage, 
                                const size_t numbase ) 
    {
      if( numBaseFct == numbase ) 
        return new ThisType( rangeStorage );
      else 
        return EvaluateCaller< BaseFunctionSet, QuadratureImp, RangeVectorType, RangeFactorType,
                    LocalDofVectorType , quadNop, numBaseFct - 1 > :: createObj( rangeStorage, numbase );
    }

    static BaseType* create( const RangeVectorType& rangeStorage, 
                             const size_t quadnop, const size_t numbase ) 
    {
      if( quadNop == quadnop ) 
        return EvaluateCaller< BaseFunctionSet, QuadratureImp, RangeVectorType, RangeFactorType,
                    LocalDofVectorType , quadNop, numBaseFct > :: createObj( rangeStorage, numbase );
      else  
      {
        abort();
      }
    }
  };
   
  template <class BaseFunctionSet,
            class QuadratureImp, 
            class RangeVectorType,
            class RangeFactorType, 
            class LocalDofVectorType,
            int quadNop >
  class EvaluateCaller
    < BaseFunctionSet, QuadratureImp, RangeVectorType, RangeFactorType,
  LocalDofVectorType, quadNop, 1> : public EvaluateCallerInterface
         < QuadratureImp, RangeFactorType, LocalDofVectorType >
  {
    enum { numBaseFct = 1 };
    typedef EvaluateCallerInterface
               < QuadratureImp, RangeFactorType, LocalDofVectorType > BaseType;

    typedef EvaluateCaller< BaseFunctionSet, QuadratureImp, RangeVectorType, RangeFactorType,
            LocalDofVectorType , quadNop, numBaseFct > ThisType;
    const RangeVectorType& rangeStorage_;
  public:  
    EvaluateCaller( const RangeVectorType& rangeStorage ) 
      : rangeStorage_( rangeStorage )
    {}

    virtual void axpyRanges( const QuadratureImp& quad,
                             const RangeFactorType& rangeFactors,
                             LocalDofVectorType & dofs ) const
    {
      BaseFunctionSet ::template 
        doAxpyRanges < quadNop, numBaseFct > 
        ( quad, rangeStorage_, rangeFactors, dofs );
    }

    virtual void evaluateRanges( const QuadratureImp& quad,
                                 const LocalDofVectorType & dofs,
                                 RangeFactorType& rangeFactors) const
    {
      BaseFunctionSet ::template 
        doEvaluateRanges < quadNop, numBaseFct > 
        ( quad, rangeStorage_, dofs, rangeFactors );
    }

    static BaseType* createObj( const RangeVectorType& rangeStorage, 
                                const size_t numbase ) 
    {
      if( numBaseFct == numbase ) 
        return new ThisType( rangeStorage );
      else 
        abort();
    }

    static BaseType* create( const RangeVectorType& rangeStorage, 
                             const size_t quadnop, const size_t numbase ) 
    {
      if( quadNop == quadnop ) 
        return EvaluateCaller< BaseFunctionSet, QuadratureImp, RangeVectorType, RangeFactorType,
                    LocalDofVectorType , quadNop, numBaseFct > :: createObj( rangeStorage, numbase );
      else  
      {
        return EvaluateCaller< BaseFunctionSet, QuadratureImp, RangeVectorType, RangeFactorType,
                    LocalDofVectorType , quadNop - 1, numBaseFct > :: create( rangeStorage, quadnop, numbase );
      }
    }
  };
   
  template <class BaseFunctionSet,
            class QuadratureImp, 
            class RangeVectorType,
            class RangeFactorType, 
            class LocalDofVectorType>
  class EvaluateCaller
    < BaseFunctionSet, QuadratureImp, RangeVectorType, RangeFactorType,
  LocalDofVectorType, 1, 1> : public EvaluateCallerInterface
         < QuadratureImp, RangeFactorType, LocalDofVectorType >
  {
    enum { quadNop = 1 };
    enum { numBaseFct = 1 };
    typedef EvaluateCallerInterface
               < QuadratureImp, RangeFactorType, LocalDofVectorType > BaseType;

    typedef EvaluateCaller< BaseFunctionSet, QuadratureImp, RangeVectorType, RangeFactorType,
            LocalDofVectorType , quadNop, numBaseFct > ThisType;
    const RangeVectorType& rangeStorage_;
  public:  
    EvaluateCaller( const RangeVectorType& rangeStorage ) 
      : rangeStorage_( rangeStorage )
    {}

    virtual void axpyRanges( const QuadratureImp& quad,
                             const RangeFactorType& rangeFactors,
                             LocalDofVectorType & dofs ) const
    {
      BaseFunctionSet ::template 
        doAxpyRanges < quadNop, numBaseFct > 
        ( quad, rangeStorage_, rangeFactors, dofs );
    }

    virtual void evaluateRanges( const QuadratureImp& quad,
                                 const LocalDofVectorType & dofs,
                                 RangeFactorType& rangeFactors) const
    {
      BaseFunctionSet ::template 
        doEvaluateRanges < quadNop, numBaseFct > 
        ( quad, rangeStorage_, dofs, rangeFactors );
    }

    static BaseType* createObj( const RangeVectorType& rangeStorage, 
                                const size_t numbase ) 
    {
      if( numBaseFct == numbase ) 
        return new ThisType( rangeStorage );
      else 
        abort();
    }

    static BaseType* create( const RangeVectorType& rangeStorage, 
                             const size_t quadnop, const size_t numbase ) 
    {
      if( quadNop == quadnop ) 
        return EvaluateCaller< BaseFunctionSet, QuadratureImp, RangeVectorType, RangeFactorType,
                    LocalDofVectorType , quadNop, numBaseFct > :: createObj( rangeStorage, numbase );
      else  
      {
        abort();
      }
    }
  };

} // end namespace Fem 

} // end namespace Dune

#endif
