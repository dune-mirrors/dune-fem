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

      template < int passId > 
      struct SetEntity;
      template < int passId > 
      struct Evaluate;
      template < int passId > 
      struct EvaluateQuadrature ;  
      template< int passId > 
      struct Jacobian;
      template< int passId > 
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
        ForLoop< Evaluate, 0, tuple_size< LocalFunctionTupleType >::value-1 > 
          :: apply( localFunctionTuple_, coordinate( x ), values );
      }

      /** \brief evaluate jacobians of local functions
       *
       *  \param[in]  x          quadrature point or local coordinate
       *  \param[in]  jacobians  jacobians of local functions
       */
      template< class PointType >
      void jacobian ( const PointType &x, JacobianRangeTupleType &jacobians ) const
      {
        ForLoop< Jacobian, 0, tuple_size< LocalFunctionTupleType >::value-1 > 
          :: apply( localFunctionTuple_, coordinate( x ), jacobians );
      }

      /** \brief evaluate hessians of local functions
       *
       *  \param[in]  x         quadrature point or local coordinate
       *  \param[in]  hessians  hessians of local functions
       */
      template< class PointType >
      void hessian ( const PointType &x, HessianRangeTupleType &hessians ) const
      {
        ForLoop< Hessian, 0, tuple_size< LocalFunctionTupleType >::value-1 >
          :: apply( localFunctionTuple_, coordinate( x ), hessians );
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
        ForLoop< SetEntity, 0, tuple_size< LocalFunctionTupleType >::value-1 >
          :: apply( localFunctions, entity );
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
    template< int passId >
    struct LocalFunctionTuple< DiscreteModel, Entity>::SetEntity
    {
      template< class LocalFunctionTuple >
      static void apply ( LocalFunctionTuple &localFunctionTuple, const EntityType& entity )
      {
        get< passId >( localFunctionTuple ).init( entity );
      }
    };



    // Implementation of LocalFunctionTuple::Evaluate
    // ----------------------------------------------

    template< class DiscreteModel, class Entity >
    template< int passId >
    struct LocalFunctionTuple< DiscreteModel, Entity >::Evaluate
    {
      template< class LocalFunctionTuple, class Range >
      static void apply ( LocalFunctionTuple &localFunctionTuple, 
                          const LocalCoordinateType& x, 
                          Range &value )
      {
        get< passId >( localFunctionTuple ).evaluate( x, value );
      }
    };



    // Implementation of LocalFunctionTuple::EvaluateQuadrature
    // --------------------------------------------------------

    template< class DiscreteModel, class Entity >
    template< int pos >
    struct LocalFunctionTuple< DiscreteModel, Entity >::EvaluateQuadrature
    {
      template< class Quadrature, class LocalFunctionTuple, class VectorOfTuples >
      static void apply ( const Quadrature &quadrature, 
                          LocalFunctionTuple &localFunctionTuple, 
                          VectorOfTuples &vectorOfTuples )
      {
        TupleToVectorConverter< VectorOfTuples, pos > vector( vectorOfTuples );
        Dune::get< pos >( localFunctionTuple ).evaluateQuadrature( quadrature, vector );
      }
    };



    // Implementation of LocalFunctionTuple::Jacobian
    // ----------------------------------------------

    template< class DiscreteModel, class Entity >
    template< int passId >
    struct LocalFunctionTuple< DiscreteModel, Entity >::Jacobian
    {
      template< class LocalFunctionTuple, class JacobianRangeType >
      static void apply ( LocalFunctionTuple &localFunctionTuple, 
                          const LocalCoordinateType& x, 
                          JacobianRangeType &jacobian )
      {
        get< passId >( localFunctionTuple ).jacobian( x, jacobian );
      }
    };



    // Implementation of LocalFunctionTuple::Hessian
    // ---------------------------------------------

    template< class DiscreteModel, class Entity >
    template< int passId >
    struct LocalFunctionTuple< DiscreteModel, Entity>::Hessian
    {
      template< class LocalFunctionTuple, class HessianRangeType >
      static void visit ( LocalFunctionTuple &localFunctionTuple, 
                          const LocalCoordinateType& x,
                          HessianRangeType &hessian )
      {
        get< passId >( localFunctionTuple ).hessian( x, hessian );
      }

    private:
      const LocalCoordinateType &x_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_PASS_COMMON_LOCALFUNCTIONTUPLE_HH
