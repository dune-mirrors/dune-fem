#ifndef DUNE_BASEFUNCTIONSTORAGE_HH
#define DUNE_BASEFUNCTIONSTORAGE_HH

#define DUNE_FEM_BASEFUNC_USE_SSE

#include <map>
#include <vector>
#include <list>

#include <dune/common/typetraits.hh>
#include <dune/grid/common/grid.hh>

#include <dune/fem/space/common/arrays.hh>
#include <dune/fem/space/basefunctions/basefunctionfactory.hh>
#include <dune/fem/space/basefunctions/storageinterface.hh>
#include <dune/fem/quadrature/quadrature.hh>

#include <dune/fem/storage/ssematrixdynamic.hh>

namespace Dune
{

  //! forward deklaration for CachingInterface 
  class CachingInterface;
   
  template <class FunctionSpaceImp>
  class StorageBase : public StorageInterface<FunctionSpaceImp::dimDomain> 

  {
  public:
    typedef BaseFunctionFactory<FunctionSpaceImp> FactoryType;
    typedef typename FunctionSpaceImp::DomainType DomainType;
    typedef typename FunctionSpaceImp::RangeType RangeType;
    typedef typename FunctionSpaceImp::JacobianRangeType JacobianRangeType;

  public:
    //! Constructor
    explicit StorageBase( const FactoryType& factory );

    //! Destructor (must be virtual)
    virtual ~StorageBase();

    //! Number of base functions
    inline
    int numBaseFunctions() const;

    //! Evaluation of base functions and their derivatives
    template <int diffOrd>
    inline
    void evaluate(int baseFunct,
                  const FieldVector<int, diffOrd>& diffVar,
                  const DomainType& xLocal, 
                  RangeType& result) const;

    //! Evaluation of jacobians
    inline
    void jacobian(int baseFunct, 
                  const DomainType& xLocal, 
                  JacobianRangeType& result) const;

    inline void cacheQuadrature(const size_t id, 
                                const size_t codim, 
                                const size_t quadSize) const 
    {}
    
    //! return geometry type 
    GeometryType geometryType () const { return elementGeometry_; }

    // default implementation to fill cache with base function evaluation 
    template <class StorageType,
              class QuadratureType,
              class RangeMatrixType>
    static inline 
    void fillRangeCache(const StorageType& storage,
                        const QuadratureType& quad, 
                        RangeMatrixType& ranges ) 
    {
      // for all quad points evaluate all basis functions 
      const size_t quadNop = quad.nop();
      const int numBase = storage.numBaseFunctions();

      FieldVector<int, 0> diffVar;
#ifndef DUNE_FEM_BASEFUNC_USE_SSE
      ranges.resize( quadNop );
#endif
      for( size_t qp = 0; qp < quadNop; ++qp ) 
      {
#ifndef DUNE_FEM_BASEFUNC_USE_SSE
        ranges[ qp ].resize( numBase );
#endif
        for( int i = 0 ; i< numBase; ++i ) 
        {
#ifndef DUNE_FEM_BASEFUNC_USE_SSE
          storage.evaluate( i, diffVar, quad[ qp ], ranges[ qp ][ i ] );
#else 
          RangeType val;
          storage.evaluate( i, diffVar, quad[ qp ], val );
          ranges[ qp ][ i ]  = val[ 0 ];
#endif
        }
      }
    }

    // default implementation to fill cache with base function evaluation 
    template <class StorageType,
              class QuadratureType,
              class JacobianRangeVectorType>
    static inline 
    void fillJacobianCache(const StorageType& storage,
                           const QuadratureType& quad, 
                           JacobianRangeVectorType& jacobians)
    {
      // for all quad points evaluate all basis functions 
      const size_t quadNop = quad.nop();
      const int numBase = storage.numBaseFunctions();

      FieldVector<int, 0> diffVar;
      jacobians.resize( quadNop );
      for( size_t qp = 0; qp < quadNop; ++qp ) 
      {
        jacobians[ qp ].resize( numBase );
        for( int i = 0 ; i< numBase; ++i ) 
        {
          storage.jacobian( i, quad[ qp ], jacobians[ qp ][ i ] );
        }
      }
    }

    // default method simply returns quadrature point 
    template <class QuadratureType> 
    int applyCaching (const QuadratureType&, const int quadraturePoint ) const 
    {
      return quadraturePoint;
    } 

  private:
    typedef typename FactoryType::BaseFunctionType BaseFunctionType;

  private:
    const int storageSize_;
    BaseFunctionType** storage_;
    mutable FieldVector<int, 1> diffVar1_;

  protected:
    const GeometryType elementGeometry_;
  };

  /** \brief simple base function storage
   *
   *  This base function storage just forwards any evaluations
   *  (quadrature point or not) to the underlying base functions.
   */
  template< class FunctionSpaceImp >
  class SimpleStorage
  : public StorageBase< FunctionSpaceImp >
  {
  public:
    typedef FunctionSpaceImp FunctionSpaceType;
    
  private:
    typedef SimpleStorage< FunctionSpaceType > ThisType;
    typedef StorageBase< FunctionSpaceType > BaseType;

  public:
    typedef BaseFunctionFactory< FunctionSpaceType > FactoryType;
    typedef typename FunctionSpaceType :: DomainType DomainType;
    typedef typename FunctionSpaceType :: RangeType RangeType;
    typedef typename FunctionSpaceType :: JacobianRangeType JacobianRangeType;

    typedef MutableArray< MutableArray<RangeType> >         RangeVectorType;
    typedef MutableArray< MutableArray<JacobianRangeType> > JacobianRangeVectorType;

  public:
    using BaseType :: evaluate;
    using BaseType :: jacobian;

  public:
    //! Constructor
    inline explicit SimpleStorage ( const FactoryType &factory )
    : BaseType( factory )
    {
    }

    template< int diffOrder, class QuadratureType >
    inline void evaluate ( const int baseFunction,
                           const FieldVector< int, diffOrder > &diffVariable,
                           const QuadraturePointWrapper< QuadratureType > &x,
                           RangeType &ret ) const
    {
      evaluate( baseFunction, diffVariable, coordinate( x ), ret );
    }
    
    template< class QuadratureType >
    inline void jacobian( const int baseFunction,
                          const QuadraturePointWrapper< QuadratureType > &x,
                          JacobianRangeType &ret ) const
    {
      jacobian( baseFunction, coordinate( x ), ret );
    }

    template <class QuadratureType> 
    const RangeVectorType& getRangeStorage( const QuadratureType& quad ) const 
    {
      this->fillRangeCache( *this, quad, ranges_ );
      return ranges_;
    }

    template <class QuadratureType> 
    const JacobianRangeVectorType& getJacobianStorage( const QuadratureType& quad ) const 
    {
      this->fillJacobianCache( *this, quad, jacobianTmp_ );
      return jacobianTmp_;
    }

  protected:  
    mutable RangeVectorType ranges_;
    mutable JacobianRangeVectorType jacobianTmp_;
  };


  //! \brief Storage scheme which caches evaluations of base function values
  //! and derivatives.
  //! This storage scheme works in conjunction with the CacheQuadrature.
  //! \warning This works only for conforming grids so far!!!!!
  //! \todo Implement switch for non-conforming grids!!!!
  template< class FunctionSpaceImp >
  class CachingStorage
  : public StorageBase< FunctionSpaceImp >
  {
  public:
    typedef FunctionSpaceImp FunctionSpaceType;

  private:
    typedef CachingStorage< FunctionSpaceType > ThisType;
    typedef StorageBase< FunctionSpaceType > BaseType;
    
    friend class StorageInterface< FunctionSpaceType :: dimDomain >;
    
  public:
    typedef BaseFunctionFactory< FunctionSpaceType > FactoryType;
    typedef typename FunctionSpaceType :: DomainType DomainType;
    typedef typename FunctionSpaceType :: RangeType RangeType;
    typedef typename FunctionSpaceType :: JacobianRangeType JacobianRangeType;

    typedef typename FunctionSpaceImp :: RangeFieldType  RealType;

    typedef MutableArray< MutableArray<RangeType> >         RangeVectorType;
    typedef MutableArray< MutableArray<JacobianRangeType> > JacobianRangeVectorType;

    typedef Fem :: SSEMatrix< RealType > RangeCacheMatrixType ;
    typedef std::vector< std::pair< RangeCacheMatrixType*, RangeCacheMatrixType* > > RangeCacheMatrixContainerType ;

  private:
    typedef std::map<size_t, bool> RangeStoredType;
    typedef std::map<size_t, bool> JacobianRangeStoredType;
    typedef typename RangeStoredType::iterator RangeIteratorType;
    typedef typename JacobianRangeStoredType::iterator JacobianRangeIteratorType;
    typedef std::pair<
      RangeIteratorType, JacobianRangeIteratorType> ReturnPairType;


    typedef MutableArray< RangeVectorType >         RangeContainerType;
    typedef MutableArray< JacobianRangeVectorType > JacobianRangeContainerType;

    //! evaluation, calls generic evaluate method
    //! used quad.point(quadPoint) for evaluation 
    template <class QuadratureType, bool isCachingQuadrature > 
    struct Evaluate
    {
      static void evaluate(const ThisType& storage,
                           const int baseFunct, 
                           const FieldVector<int, 0>& diffVar,
                           const QuadratureType& quad, 
                           const int quadPoint,
                           const RangeContainerType& ranges,
                           RangeType& result)
      {
        storage.evaluate(baseFunct,diffVar,quad.point(quadPoint),result);
      }

      static void evaluate(const ThisType& storage,
                           const int baseFunct, 
                           const FieldVector<int, 1>& diffVar,
                           const QuadratureType& quad, 
                           const int quadPoint,
                           const JacobianRangeContainerType& jacobians,
                           RangeType& result)
      {
        storage.evaluate(baseFunct,diffVar,quad.point(quadPoint),result);
      }

      static void jacobian(const ThisType& storage,
                           const int baseFunct, 
                           const QuadratureType& quad, 
                           const int quadPoint, 
                           const JacobianRangeContainerType& jacobians,
                           JacobianRangeType& result)
      {
        storage.jacobian(baseFunct,quad.point(quadPoint),result); 
      }

      static const RangeCacheMatrixType& 
      evaluate(const ThisType& storage,
               const QuadratureType& quad, 
               RangeCacheMatrixContainerType& rangeMatrices)
      {
        assert( rangeMatrices[ quad.id() ].second );
        RangeCacheMatrixType& matrix = * (rangeMatrices[ quad.id() ].second);
        storage.fillRangeCache( storage, quad, matrix );
        return matrix; 
      }

      static const RangeVectorType& 
      evaluate(const ThisType& storage,
               const QuadratureType& quad, 
               const RangeContainerType&,
               RangeVectorType& rangeTmp )
      {
        storage.fillRangeCache( storage, quad, rangeTmp );
        return rangeTmp;
      }

      static const JacobianRangeVectorType& 
      jacobian(const ThisType& storage,
               const QuadratureType& quad, 
               const JacobianRangeContainerType&,
               const JacobianRangeVectorType& jacobianTmp )
      {
        storage.fillJacobianCache( storage, quad, jacobianTmp );
        return jacobianTmp;
      }
      
    };
    
    //! evaluation for CachingQuadratures 
    template <class QuadratureType> 
    struct Evaluate<QuadratureType,true> 
    {
      static void evaluate(const ThisType& storage,
                           const int baseFunct, 
                           const FieldVector<int, 0>& diffVar,
                           const QuadratureType& quad, 
                           const int quadPoint,
                           const RangeContainerType& ranges,
                           RangeType& result)
      {
        result = ranges[quad.id()][quad.cachingPoint(quadPoint)][baseFunct];
      }

      static void evaluate(const ThisType& storage,
                           const int baseFunct, 
                           const FieldVector<int, 1>& diffVar,
                           const QuadratureType& quad, 
                           const int quadPoint,
                           const JacobianRangeContainerType& jacobians,
                           RangeType& result)
      {
        const JacobianRangeType& jResult =
          jacobians[quad.id()][quad.cachingPoint(quadPoint)][baseFunct];

        for (size_t i = 0; i < RangeType::dimension; ++i)
        {
          result[i] = jResult[i][diffVar[0]];
        }
      }
      
      static void jacobian(const ThisType& storage,
                           const int baseFunct, 
                           const QuadratureType& quad, 
                           const int quadPoint, 
                           const JacobianRangeContainerType& jacobians,
                           JacobianRangeType& result)
      {
        result = jacobians[quad.id()][quad.cachingPoint(quadPoint)][baseFunct];
      }

      static const RangeCacheMatrixType& 
      evaluate(const ThisType& storage,
               const QuadratureType& quad, 
               const RangeCacheMatrixContainerType& rangeMatrices)
      {
        assert( rangeMatrices[ quad.id() ].first );
        return * (rangeMatrices[ quad.id() ].first);
      }

      static const RangeVectorType& 
      evaluate(const ThisType& storage,
               const QuadratureType& quad, 
               const RangeContainerType& ranges,
               RangeVectorType& )
      {
        return ranges[ quad.id() ];
      }

      static const JacobianRangeVectorType& 
      jacobian(const ThisType& storage,
               const QuadratureType& quad, 
               const JacobianRangeContainerType& jacobians,
               const JacobianRangeVectorType& )
      {
        return jacobians[ quad.id() ];
      }
      
    };

  public:
    using BaseType :: cacheExistingQuadratures;
    using BaseType :: evaluate;
    using BaseType :: jacobian;
   
  public:
    //! Constructor
    inline explicit CachingStorage ( const FactoryType &factory )
    : BaseType( factory )
    {
      cacheExistingQuadratures( *this );
    }
    
    //! Destructor
    ~CachingStorage ()
    {
    }

    //! evaulate base function 
    template< class QuadratureType >
    inline void evaluate ( const int baseFunction,
                           const FieldVector< int, 0 > &diffVariable,
                           const QuadraturePointWrapper< QuadratureType > &x, 
                           RangeType &ret ) const;

    //! evaluate derivative of base function 
    template< class QuadratureType >
    inline void evaluate ( const int baseFunction,
                           const FieldVector< int, 1 > &diffVariable,
                           const QuadraturePointWrapper< QuadratureType > &x, 
                           RangeType &ret ) const;

    //! get derivative of base function 
    template< class QuadratureType >
    inline void jacobian ( const int baseFunction,
                           const QuadraturePointWrapper< QuadratureType > &x,
                           JacobianRangeType &ret ) const;

    template <class QuadratureType> 
    const RangeCacheMatrixType& getRangeMatrix( const QuadratureType& quad ) const 
    {
      enum { cachable = Conversion< QuadratureType, CachingInterface > :: exists };
      return Evaluate< QuadratureType, cachable > :: evaluate( *this, quad, rangeMatrices_ );
    }

    template <class QuadratureType> 
    const RangeVectorType& getRangeStorage( const QuadratureType& quad ) const 
    {
      enum { cachable = Conversion< QuadratureType, CachingInterface > :: exists };
      return Evaluate< QuadratureType, cachable > :: evaluate( *this, quad, ranges_, rangeTmp_ );
    }

    template <class QuadratureType> 
    const JacobianRangeVectorType& getJacobianStorage( const QuadratureType& quad ) const 
    {
      enum { cachable = Conversion< QuadratureType, CachingInterface > :: exists };
      return Evaluate< QuadratureType, cachable > :: 
        jacobian( *this, quad, jacobians_, jacobianTmp_ );
    }

    template < class QuadratureType, bool caching > 
    struct ApplyCaching
    {
      static int apply(const QuadratureType& quad, const int i ) 
      {
        return i;
      }
    };

    template < class QuadratureType > 
    struct ApplyCaching< QuadratureType, true >
    {
      static int apply(const QuadratureType& quad, const int i ) 
      {
        return quad.cachingPoint( i );
      }
    };

    template <class QuadratureType> 
    int applyCaching (const QuadratureType& quad, const int i ) const 
    {
      enum { applyCache = Conversion< QuadratureType, CachingInterface > :: exists };
      return ApplyCaching< QuadratureType, applyCache > :: apply( quad, i );
    } 
  private:
    // caches the quadrature, see also addEntry.. 
    inline void cacheQuadrature(const size_t id, 
                                const size_t codim,
                                const size_t quadSize) const;
 
    // here a switch-case for codim is done and then addEntry called
    inline ReturnPairType addEntryInterface(const size_t id, 
                                            const size_t codim,
                                            const size_t quadSize) const;
    
    // here the real caching is done 
    template <int codim>
    inline ReturnPairType addEntry(const size_t id, const size_t quadSize) const;

  protected:
    mutable RangeContainerType ranges_;
    mutable JacobianRangeContainerType jacobians_;
    mutable RangeStoredType rangestored_;
    mutable JacobianRangeStoredType jacobianstored_;

    mutable RangeCacheMatrixContainerType rangeMatrices_;

    mutable RangeVectorType rangeTmp_;
    mutable JacobianRangeVectorType jacobianTmp_;
  };

} // end namespace Dune

#include "basefunctionstorage.cc"
#endif
