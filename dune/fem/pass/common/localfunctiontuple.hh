#ifndef DUNE_FEM_PASS_COMMON_LOCALFUNCTIONTUPLE_HH
#define DUNE_FEM_PASS_COMMON_LOCALFUNCTIONTUPLE_HH

#include <cassert>
#include <cstddef>

#include <dune/common/exceptions.hh>
#include <dune/common/nullptr.hh>
#include <dune/common/tuples.hh>
#include <dune/common/tupleutility.hh>

#include <dune/fem/quadrature/quadrature.hh>

#include "localfunctionselector.hh"

namespace
{

  // LocalFunctionEvaluator
  // ----------------------

  /*
   * \brief Choose local function for given discrete function.
   *        Use LocalFunctionSelector to define the type of 
   *        the local function.
   */
  template< class DiscreteFunction >
  struct LocalFunctionEvaluator
  {
    typedef typename Dune::Fem::LocalFunctionSelector< 
        typename Dune::TypeTraits< DiscreteFunction >::ReferredType
      >::Type Type; 

    static Type apply ( const DiscreteFunction &discreteFunction )
    {
      return Type( discreteFunction );
    }
  };



  // RangeTypeEvaluator
  // ------------------

  /*
   * \brief Get range type from local function type.
   */
  template< class LocalFunction >
  struct RangeTypeEvaluator
  {
    typedef typename LocalFunction::RangeType Type;
  };



  // JacobianRangeTypeEvaluator
  // --------------------------

  /*
   * \brief Get jacobian range type from local function type.
   */
  template< class LocalFunction >
  struct JacobianRangeTypeEvaluator
  {
    typedef typename LocalFunction::JacobianRangeType Type;
  };



  // HessianRangeTypeEvaluator
  // --------------------------

  /*
   * \brief Get hessian range type from local function type.
   */
  template< class LocalFunction >
  struct HessianRangeTypeEvaluator
  {
    typedef typename LocalFunction::HessianRangeType Type;
  };

} // namespace 



namespace Dune
{

  namespace Fem
  {

    /** \brief wrapper class to convert a vector of tuples of RangeTypes into something 
               that behaves like a vector< RangeType >
    */           
    template <class VectorTupleType, int passId >
    class TupleToVectorConverter
    {
      // no copying
      TupleToVectorConverter(const TupleToVectorConverter&);

      //! standard case 
      template <int pos, class Tuple> 
      struct TupleElement
      {
        typedef typename tuple_element< pos, Tuple> :: type type ;
      };

      //! specialization to extract correct tuple component
      template <int pos, class Tuple, class Types> 
      struct TupleElement< pos, TypeIndexedTuple< Tuple, Types > > 
      {
        typedef typename tuple_element< pos, Tuple> :: type type ;
      };

    public:
      typedef typename VectorTupleType :: value_type TupleType;
      typedef typename TupleElement< passId, TupleType > :: type  ValueType;
      typedef ValueType value_type ;

      //! constructor
      TupleToVectorConverter(VectorTupleType& vec)
        : vector_( vec )
      {}

      //! return reference to i-th entry of vector and passId's tuple component
      ValueType& operator [] (const size_t i)
      {
        assert( i < vector_.size() );
        return get< passId >( vector_[ i ] );
      }

      //! return reference to i-th entry of vector and passId's tuple component
      const ValueType& operator [] (const size_t i) const
      {
        assert( i < vector_.size() );
        return get< passId >( vector_[ i ] );
      }

      //! return size of vector 
      size_t size() const
      {
        return vector_.size();
      }

    protected:
      VectorTupleType& vector_;
    };

    /** \brief ForEach class that is used with the evaluateQuadrature call below */
    template <class TupleType, class VectorType>
    class ForEachValueVector {
    public:
      //! Constructor
      //! \param t1 First tuple.
      //! \param t2 Second tuple.
      ForEachValueVector(TupleType& tuple, VectorType& vec ) :
        tuple_( tuple ),
        vector_( vec )
      {}

      //! Applies the function object f to the pair of tuples.
      //! \param f The function object to apply on the pair of tuples.
      template <class Functor>
      void apply(Functor& f)
      {
        // iterate over tuple elements for 0 to tuple_size-1 
        Apply<0, tuple_size< TupleType >::value >::apply(f, tuple_, vector_ );
      }

    private:
      template <int passId, int size>  
      struct Apply
      {
        template <class Functor, class Tuple, class VectorOfTuples> 
        static void apply( Functor& f, Tuple& tuple, VectorOfTuples& vectorOfTuples )
        {
          TupleToVectorConverter< VectorOfTuples, passId > vector ( vectorOfTuples );
          // call functor (here evaluateQuadrature on local function)
          f.visit( get< passId >(tuple), vector );
          // got to next tuple element
          Apply< passId+1, size> :: apply( f, tuple, vectorOfTuples );
        }
      };

      //! termination of tuple iteration when size is reached 
      template <int size>  
      struct Apply< size, size >
      {
        template <class Functor, class Tuple, class VectorOfTuples> 
        static void apply( Functor& f, Tuple& tuple, VectorOfTuples& vector )
        {
          // do nothing here, since this is the terminating call
        }
      };

    private:
      TupleType&   tuple_;
      VectorType&  vector_;
    };



    // LocalFunctionTuple
    // ------------------

    /*
     * \brief A wrapper for a tuple of localfunctions. It's interface
     *        mimicks the LocalFunction interface
     *        (see dune/fem/function/localfunction/localfunction.hh).
     */
    template< class DiscreteFunctionTuple, class Entity >
    class LocalFunctionTuple
    {
      typedef LocalFunctionTuple< DiscreteFunctionTuple, Entity > ThisType;

      struct SetEntity;
      struct Evaluate;
      template <class Quadrature> 
      struct EvaluateQuadrature ;  
      struct Jacobian;
      struct Hessian;
   
    public:
      //! \brief discrete function tuple
      typedef DiscreteFunctionTuple DiscreteFunctionTupleType;
      //! \brief entity type
      typedef Entity EntityType;

      //! \brief type of local function tuple
      typedef typename Dune::ForEachType< LocalFunctionEvaluator, DiscreteFunctionTupleType >::Type LocalFunctionTupleType;

      // ! \brief type of range type tuple
      typedef typename Dune::ForEachType< RangeTypeEvaluator, LocalFunctionTupleType >::Type RangeTupleType;
      // ! \brief type of jacobian range type tuple
      typedef typename Dune::ForEachType< JacobianRangeTypeEvaluator, LocalFunctionTupleType >::Type JacobianRangeTupleType;
      // ! \brief type of hessian range type tuple
      typedef typename Dune::ForEachType< HessianRangeTypeEvaluator, LocalFunctionTupleType >::Type HessianRangeTupleType;

    protected:
      typedef typename EntityType::Geometry GeometryType;
      typedef typename GeometryType::LocalCoordinate LocalCoordinateType;

    public:
      template< class Factory >
      LocalFunctionTuple ( Factory factory )
      : localFunctionTuple_( Dune::transformTuple< LocalFunctionEvaluator, Factory >( factory ) )
      {}

      /** \brief set local functions to given entity
       *
       *  \param[in]  entity  grid part entity
       */
      void init ( const EntityType &entity )
      {
        entity_ = &entity;
        init( entity, localFunctionTuple_ );
      }

      /** \brief return entity */
      const EntityType &entity () const
      {
        assert( entity_ );
        return *entity_;
      }

      /** \brief evaluate local functions
       *
       *  \param[in]  x       quadrature point or local coordinate
       *  \param[in]  values  values of local functions
       */
      template< class PointType >
      void evaluate ( const PointType &x, RangeTupleType &values ) const
      {
        Dune::ForEachValuePair< LocalFunctionTupleType, RangeTupleType > forEach( localFunctionTuple_, values );
        Evaluate functor( coordinate( x ) );
        forEach.apply( functor );
      }

      /** \brief evaluate jacobians of local functions
       *
       *  \param[in]  x          quadrature point or local coordinate
       *  \param[in]  jacobians  jacobians of local functions
       */
      template< class PointType >
      void jacobian ( const PointType &x, JacobianRangeTupleType &jacobians ) const
      {
        Dune::ForEachValuePair< LocalFunctionTupleType, JacobianRangeTupleType > forEach( localFunctionTuple_, jacobians );
        Jacobian functor( coordinate( x ) );
        forEach.apply( functor );
      }

      /** \brief evaluate hessians of local functions
       *
       *  \param[in]  x         quadrature point or local coordinate
       *  \param[in]  hessians  hessians of local functions
       */
      template< class PointType >
      void hessian ( const PointType &x, HessianRangeTupleType &hessians ) const
      {
        Dune::ForEachValuePair< LocalFunctionTupleType, HessianRangeTupleType > forEach( localFunctionTuple_, hessians );
        Hessian functor( coordinate( x ) );
        forEach.apply( functor );
      }

      /** \brief evaluate local functions for quadrature
       *
       *  \param[in]  quadrature  quadrature
       *  \param[in]  vector      a (dynamic) vector of range tuples
       */
      template< class QuadratureType, class TupleVectorType >
      void evaluateQuadrature ( const QuadratureType &quadrature, TupleVectorType &vector ) const
      {
        assert( vector.size() >= quadrature.nop() );
        ForEachValueVector< LocalFunctionTupleType, TupleVectorType > forEach( localFunctionTuple_, vector );
        EvaluateQuadrature< QuadratureType > functor( quadrature );
        forEach.apply( functor );
      }

    protected:
      void init ( const EntityType &entity, LocalFunctionTupleType &localFunctions )
      {
        ForEachValue< LocalFunctionTupleType > forEach( localFunctionTuple_ );
        SetEntity functor( entity );
        forEach.apply( functor );
      }

      LocalFunctionTupleType &localFunctions () { return localFunctionTuple_; }

      const LocalFunctionTupleType &localFunctions () const { return localFunctionTuple_; }

    private:
      LocalFunctionTuple ( const ThisType & );
      ThisType &operator= ( const ThisType & );

      mutable LocalFunctionTupleType localFunctionTuple_;
      const EntityType *entity_;
    };



    // Implementation of LocalFunctionTuple::SetEntity
    // -----------------------------------------------

    template< class DiscreteModel, class Entity >
    struct LocalFunctionTuple< DiscreteModel, Entity>::SetEntity
    {
      SetEntity ( const EntityType &entity ) : entity_( entity ) {}

      template< class LocalFunction >
      void visit ( LocalFunction &localFunction ) const
      {
        localFunction.init( entity_ );
      }

    private:
      const EntityType &entity_;
    };



    // Implementation of LocalFunctionTuple::Evaluate
    // ----------------------------------------------

    template< class DiscreteModel, class Entity >
    struct LocalFunctionTuple< DiscreteModel, Entity >::Evaluate
    {
      Evaluate ( const LocalCoordinateType &x ) : x_( x ) {}

      template< class LocalFunction, class Range >
      void visit ( LocalFunction &localFunction, Range &value ) const
      {
        localFunction.evaluate( x_, value );
      }

    private:
      const LocalCoordinateType &x_;
    };




    // Implementation of LocalFunctionTuple::EvaluateQuadrature
    // --------------------------------------------------------

    template< class DiscreteModel, class Entity >
    template< class Quadrature >
    struct LocalFunctionTuple< DiscreteModel, Entity >::EvaluateQuadrature
    {
      EvaluateQuadrature ( const Quadrature& quadrature ) : quadrature_( quadrature ) {}

      template< class LocalFunction, class ValueVector >
      void visit ( LocalFunction &localFunction, ValueVector &values ) const
      {
        localFunction.evaluateQuadrature( quadrature_, values );
      }

    private:
      const Quadrature& quadrature_;
    };



    // Implementation of LocalFunctionTuple::Jacobian
    // ----------------------------------------------

    template< class DiscreteModel, class Entity >
    struct LocalFunctionTuple< DiscreteModel, Entity >::Jacobian
    {
      Jacobian ( const LocalCoordinateType &x ) : x_( x ) {}

      template< class LocalFunction, class JacobianRangeType >
      void visit ( LocalFunction &localFunction, JacobianRangeType &jacobian ) const
      {
        localFunction.jacobian( x_, jacobian );
      }

    private:
      const LocalCoordinateType &x_;
    };



    // Implementation of LocalFunctionTuple::Hessian
    // ---------------------------------------------

    template< class DiscreteModel, class Entity >
    struct LocalFunctionTuple< DiscreteModel, Entity>::Hessian
    {
      Hessian ( const LocalCoordinateType &x ) : x_( x ) {}

      template< class LocalFunction, class HessianRangeType >
      void visit ( LocalFunction &localFunction, HessianRangeType &hessian ) const
      {
        localFunction.hessian( x_, hessian );
      }

    private:
      const LocalCoordinateType &x_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_PASS_COMMON_LOCALFUNCTIONTUPLE_HH
