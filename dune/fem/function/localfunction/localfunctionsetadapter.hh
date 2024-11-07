#ifndef DUNE_FEM_FUNCTION_LOCALFUNCTION_LOCALFUNCTIONSETADAPTER_HH
#define DUNE_FEM_FUNCTION_LOCALFUNCTION_LOCALFUNCTIONSETADAPTER_HH

#include <cassert>
#include <cstddef>

#include <dune/fem/common/coordinate.hh>
#include <dune/fem/storage/entitygeometry.hh>

namespace Dune
{

  namespace Fem
  {

    // LocalFunctionSetAdapter
    // -----------------------

    /** \class LocalFunctionSetAdapter
     *
     *  \brief convert (global) function set to local function set
     *
     *  \tparam  Entity       entity type
     *  \tparam  FunctionSet  implementation of FunctionSet
     */
    template< class Entity, class FunctionSet >
    struct LocalFunctionSetAdapter
      : public EntityGeometryStorage< Entity >
    {
    protected:
      typedef EntityGeometryStorage< Entity > BaseType;

    public:
      //! \brief entity type
      typedef typename BaseType::EntityType EntityType;

      //! \brief geometry
      typedef typename BaseType::Geometry  Geometry;

      //! \brief function set type
      typedef FunctionSet FunctionSetType;

      //! \brief function space type
      typedef typename FunctionSet::FunctionSpaceType FunctionSpaceType;

      //! \brief domain type
      typedef typename FunctionSpaceType::DomainType DomainType;
      //! \brief range type
      typedef typename FunctionSpaceType::RangeType RangeType;
      //! \brief jacobian range type
      typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;
      //! \brief hessian range type
      typedef typename FunctionSpaceType::HessianRangeType HessianRangeType;

      explicit LocalFunctionSetAdapter ( const FunctionSet &functionSet = FunctionSet() )
        : BaseType(),
          functionSet_( functionSet )
      {}

      explicit LocalFunctionSetAdapter ( const Entity &entity, const FunctionSet &functionSet = FunctionSet() )
        : BaseType( entity ), functionSet_( functionSet )
      {}

      /** \copydoc Dune::Fem::LocalFunctionSet::order() const */
      int order () const { return functionSet_.order(); }

      using BaseType :: entity;
      using BaseType :: geometry;
      using BaseType :: bind;
      using BaseType :: unbind;

      /** \copydoc Dune::Fem::LocalFunctionSet::size() const */
      std::size_t size () const { return functionSet().size(); }

      /** \copydoc Dune::Fem::LocalFunctionSet::evaluateEach() const */
      template< class Point, class Functor >
      void evaluateEach ( const Point &x, Functor functor ) const
      {
        DomainType y = geometry().global( coordinate( x ) );
        functionSet().evaluateEach( y, functor );
      }

      /** \copydoc Dune::Fem::LocalFunctionSet::jacobianEach() const */
      template< class Point, class Functor >
      void jacobianEach ( const Point &x, Functor functor ) const
      {
        DomainType y = geometry().global( coordinate( x ) );
        functionSet().jacobianEach( y, functor );
      }

      /** \copydoc Dune::Fem::LocalFunctionSet::hessianEach() const */
      template< class Point, class Functor >
      void hessianEach ( const Point &x, Functor functor ) const
      {
        DomainType y = geometry().global( coordinate( x ) );
        functionSet().hessianEach( y, functor );
      }

      // Non-interface methods
      // ---------------------

      //! \brief return function set
      const FunctionSetType functionSet () const { return functionSet_; }

    private:
      FunctionSetType functionSet_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_FUNCTION_LOCALFUNCTION_LOCALFUNCTIONSETADAPTER_HH
