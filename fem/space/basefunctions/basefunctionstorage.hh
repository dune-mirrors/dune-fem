#ifndef DUNE_BASEFUNCTIONSTORAGE_HH
#define DUNE_BASEFUNCTIONSTORAGE_HH

//- system includes 
#include <map>
#include <vector>
#include <list>

//- Dune includes 
#include <dune/common/typetraits.hh>
#include <dune/grid/common/grid.hh>

//- local includes 
#include <dune/fem/space/common/arrays.hh>
#include <dune/fem/space/common/basefunctionfactory.hh>
#include "storageinterface.hh"

namespace Dune {

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

    inline void cacheQuadrature(size_t id, int codim) const {}
    
    //! return geometry type 
    GeometryType geometryType () const { return elementGeometry_; }

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

  private:
    typedef typename FunctionSpaceImp::DomainFieldType RealType;
    typedef MutableArray< MutableArray<RangeType> >         RangeVectorType;
    typedef MutableArray< MutableArray<JacobianRangeType> > JacobianRangeVectorType;
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

  private:
    // caches the quadrature, see also addEntry.. 
    inline void cacheQuadrature(size_t id, int codim) const;
 
    // here a switch-case for codim is done and then addEntry called
    inline ReturnPairType addEntryInterface(size_t id , int codim) const;
    
    // here the real caching is done 
    template <int codim>
    inline ReturnPairType addEntry(size_t id) const;

  private:
    mutable RangeContainerType ranges_;
    mutable JacobianRangeContainerType jacobians_;
    mutable RangeStoredType rangestored_;
    mutable JacobianRangeStoredType jacobianstored_;
  };

} // end namespace Dune

#include "basefunctionstorage.cc"
#endif
