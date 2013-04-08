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
    public:
      typedef typename VectorTupleType :: value_type TupleType;
      typedef typename tuple_element< passId, TupleType > :: type  ValueType;
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
      template < int passId > 
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
        // loop over local function tuple and call EvaluateQuadrature::apply
        ForLoop< EvaluateQuadrature, 0, tuple_size< LocalFunctionTupleType >::value-1 > 
          :: apply( quadrature, localFunctionTuple_, vector );
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
    template< int passId >
    struct LocalFunctionTuple< DiscreteModel, Entity >::EvaluateQuadrature
    {
      template< class Quadrature, class LocalFunctionTuple, class VectorOfTuples >
      static void apply( const Quadrature& quadrature, 
                         LocalFunctionTuple &localFunctionTuple, 
                         VectorOfTuples &vectorOfTuples )
      {
        TupleToVectorConverter< VectorOfTuples, passId > vector ( vectorOfTuples );
        // get passId-th local function and call evaluateQuadrature  
        get< passId > ( localFunctionTuple ).evaluateQuadrature( quadrature, vector );
      }
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
