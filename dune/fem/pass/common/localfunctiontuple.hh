#ifndef DUNE_FEM_PASS_COMMON_LOCALFUNCTIONTUPLE_HH
#define DUNE_FEM_PASS_COMMON_LOCALFUNCTIONTUPLE_HH

#include <cassert>
#include <cstddef>

#include <dune/common/exceptions.hh>
#include <dune/common/forloop.hh>
#include <dune/common/nullptr.hh>
#include <dune/common/tuples.hh>
#include <dune/common/tupleutility.hh>

#include <dune/fem/common/tupleutility.hh>

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

    // LocalFunctionTuple
    // ------------------

    /*
     * \brief A wrapper for a tuple of localfunctions. It's interface
     *        mimicks the LocalFunction interface
     *        (see dune/fem/function/localfunction/localfunction.hh).
     */
    template< class DiscreteFunctionTuple, class Entity,
              size_t TupleSize = tuple_size< DiscreteFunctionTuple >::value >
    class LocalFunctionTuple
    {
      typedef LocalFunctionTuple< DiscreteFunctionTuple, Entity, TupleSize > ThisType;

      template< int pos > struct SetEntity;
      template< int pos > struct Evaluate;
      template< int pos > struct EvaluateQuadrature;
      template< int pos > struct Jacobian;
      template< int pos > struct Hessian;
   
    public:
      //! \brief discrete function tuple
      typedef DiscreteFunctionTuple DiscreteFunctionTupleType;
      //! \brief entity type
      typedef Entity EntityType;

    protected:
      typedef typename Dune::ForEachType< LocalFunctionEvaluator, DiscreteFunctionTupleType >::Type LocalFunctionTupleType;

      typedef typename EntityType::Geometry GeometryType;
      typedef typename GeometryType::LocalCoordinate LocalCoordinateType;

    public:
      // ! \brief type of range type tuple
      typedef typename Dune::ForEachType< RangeTypeEvaluator, LocalFunctionTupleType >::Type RangeTupleType;
      // ! \brief type of jacobian range type tuple
      typedef typename Dune::ForEachType< JacobianRangeTypeEvaluator, LocalFunctionTupleType >::Type JacobianRangeTupleType;
      // ! \brief type of hessian range type tuple
      typedef typename Dune::ForEachType< HessianRangeTypeEvaluator, LocalFunctionTupleType >::Type HessianRangeTupleType;

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
        ForLoop< SetEntity, 0, TupleSize-1 >::apply( localFunctions(), entity );
      }

      /** \brief return entity */
      const EntityType &entity () const
      {
        return Dune::get< 0 >( localFunctions() ).entity();
      }

      /** \brief evaluate local functions
       *
       *  \param[in]  x       quadrature point or local coordinate
       *  \param[in]  values  values of local functions
       */
      template< class Point >
      void evaluate ( const Point &x, RangeTupleType &values ) const
      {
        ForLoop< Evaluate, 0, TupleSize-1 >::apply( localFunctions(), x, values );
      }

      /** \brief evaluate jacobians of local functions
       *
       *  \param[in]  x          quadrature point or local coordinate
       *  \param[in]  jacobians  jacobians of local functions
       */
      template< class Point >
      void jacobian ( const Point &x, JacobianRangeTupleType &jacobians ) const
      {
        ForLoop< Jacobian, 0, TupleSize-1 >::apply( localFunctions(), x, jacobians );
      }

      /** \brief evaluate hessians of local functions
       *
       *  \param[in]  x         quadrature point or local coordinate
       *  \param[in]  hessians  hessians of local functions
       */
      template< class Point >
      void hessian ( const Point &x, HessianRangeTupleType &hessians ) const
      {
        ForLoop< Hessian, 0, TupleSize-1 >::apply( localFunctions(), x, hessians );
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
        ForLoop< EvaluateQuadrature, 0, TupleSize-1 >::apply( quadrature, localFunctions(), vector );
      }

    protected:
      LocalFunctionTupleType &localFunctions () { return localFunctionTuple_; }
      const LocalFunctionTupleType &localFunctions () const { return localFunctionTuple_; }

    private:
      mutable LocalFunctionTupleType localFunctionTuple_;
    };



    // LocalFunctionTuple for Empty Tuples
    // -----------------------------------

    template< class DiscreteFunctionTuple, class Entity >
    class LocalFunctionTuple< DiscreteFunctionTuple, Entity, 0 >
    {
      typedef LocalFunctionTuple< DiscreteFunctionTuple, Entity, 0 > ThisType;

    public:
      typedef DiscreteFunctionTuple DiscreteFunctionTupleType;
      typedef Entity EntityType;

      typedef Dune::tuple<> RangeTupleType;
      typedef Dune::tuple<> JacobianRangeTupleType;
      typedef Dune::tuple<> HessianRangeTupleType;

      template< class Factory >
      LocalFunctionTuple ( Factory factory )
      : entity_( nullptr )
      {}

      void init ( const EntityType &entity ) { entity_ = &entity; }

      const EntityType &entity () const
      {
        assert( entity_ );
        return *entity_;
      }

      template< class Point >
      void evaluate ( const Point &x, RangeTupleType &values ) const
      {}

      template< class Point >
      void jacobian ( const Point &x, JacobianRangeTupleType &jacobians ) const
      {}

      template< class Point >
      void hessian ( const Point &x, HessianRangeTupleType &hessians ) const
      {}

      template< class QuadratureType, class TupleVectorType >
      void evaluateQuadrature ( const QuadratureType &quadrature, TupleVectorType &vector ) const
      {
        assert( vector.size() >= quadrature.nop() );
      }

    private:
      EntityType *entity_;
    };



    // Implementation of LocalFunctionTuple::SetEntity
    // -----------------------------------------------

    template< class DiscreteFunctionTuple, class Entity, size_t TupleSize >
    template< int pos >
    struct LocalFunctionTuple< DiscreteFunctionTuple, Entity, TupleSize >::SetEntity
    {
      static void apply ( LocalFunctionTupleType &localFunctions,
                          const EntityType &entity )
      {
        Dune::get< pos >( localFunctions ).init( entity );
      }
    };



    // Implementation of LocalFunctionTuple::Evaluate
    // ----------------------------------------------

    template< class DiscreteFunctionTuple, class Entity, size_t TupleSize >
    template< int pos >
    struct LocalFunctionTuple< DiscreteFunctionTuple, Entity, TupleSize >::Evaluate
    {
      template< class Point >
      static void apply ( const LocalFunctionTupleType &localFunctions,
                          const Point &x,
                          RangeTupleType &values )
      {
        Dune::get< pos >( localFunctions ).evaluate( x, Dune::get< pos >( values ) );
      }
    };



    // Implementation of LocalFunctionTuple::EvaluateQuadrature
    // --------------------------------------------------------

    template< class DiscreteFunctionTuple, class Entity, size_t TupleSize >
    template< int pos >
    struct LocalFunctionTuple< DiscreteFunctionTuple, Entity, TupleSize >::EvaluateQuadrature
    {
      template< class Quadrature, class VectorOfTuples >
      static void apply ( const Quadrature &quadrature,
                          const LocalFunctionTupleType &localFunctions,
                          VectorOfTuples &vectorOfTuples )
      {
        TupleToVectorConverter< VectorOfTuples, pos > vector( vectorOfTuples );
        Dune::get< pos >( localFunctions ).evaluateQuadrature( quadrature, vector );
      }
    };



    // Implementation of LocalFunctionTuple::Jacobian
    // ----------------------------------------------

    template< class DiscreteFunctionTuple, class Entity, size_t TupleSize >
    template< int pos >
    struct LocalFunctionTuple< DiscreteFunctionTuple, Entity, TupleSize >::Jacobian
    {
      template< class Point >
      static void apply ( const LocalFunctionTupleType &localFunctions,
                          const Point &x,
                          JacobianRangeTupleType &jacobians )
      {
        Dune::get< pos >( localFunctions ).jacobian( x, Dune::get< pos >( jacobians ) );
      }
    };



    // Implementation of LocalFunctionTuple::Hessian
    // ---------------------------------------------

    template< class DiscreteFunctionTuple, class Entity, size_t TupleSize >
    template< int pos >
    struct LocalFunctionTuple< DiscreteFunctionTuple, Entity, TupleSize >::Hessian
    {
      template< class Point >
      static void apply ( const LocalFunctionTupleType &localFunctions,
                          const Point &x,
                          HessianRangeTupleType &hessians )
      {
        Dune::get< pos >( localFunctions ).hessian( x, Dune::get< pos >( hessians ) );
      }
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_PASS_COMMON_LOCALFUNCTIONTUPLE_HH
