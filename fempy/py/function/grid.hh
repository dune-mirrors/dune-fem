#ifndef DUNE_FEMPY_PY_FUNCTION_GRID_HH
#define DUNE_FEMPY_PY_FUNCTION_GRID_HH

#include <string>
#include <utility>

#include <dune/fem/function/common/discretefunction.hh>
#include <dune/fem/space/common/functionspace.hh>

namespace Dune
{

  namespace FemPy
  {

    // PyGridFunction
    // --------------

    template< class GridFunction >
    class PyGridFunction
      : public Fem::Function< typename GridFunction::FunctionSpaceType, PyGridFunction< GridFunction > >,
        public Fem::HasLocalFunction
    {
      typedef Fem::Function< typename GridFunction::FunctionSpaceType, PyGridFunction< GridFunction > > Base;

    public:
      typedef typename GridFunction::GridPartType GridPartType;

      typedef typename GridFunction::EntityType EntityType;

      typedef typename Base::DomainType DomainType;
      typedef typename Base::RangeType RangeType;
      typedef typename Base::JacobianRangeType JacobianRangeType;

      class LocalFunctionType
      {
        typedef typename GridFunction::LocalFunctionType Impl;

      public:
        typedef typename Impl::EntityType EntityType;

        typedef typename Impl::FunctionSpaceType FunctionSpaceType;

        static const int dimDomain = FunctionSpaceType::dimDomain;
        static const int dimRange = FunctionSpaceType::dimRange;

        typedef typename FunctionSpaceType::DomainType DomainType;
        typedef typename FunctionSpaceType::RangeType RangeType;
        typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;

        explicit LocalFunctionType ( const PyGridFunction &gf ) : impl_( *gf.impl_ ), pyObj_( gf.pyObj_ ) {}

        LocalFunctionType ( Impl impl, pybind11::object pyObj ) : impl_( std::move( impl ) ), pyObj_( std::move( pyObj ) ) {}

        void init ( const EntityType &entity ) { impl_.init( entity ); }

        template< class Point >
        void evaluate ( const Point &x, RangeType &value ) const
        {
          impl_.evaluate( x, value );
        }

        template< class Quadrature, class Values >
        void evaluateQuadrature ( const Quadrature &quadrature, Values &values ) const
        {
          impl_.evaluateQuadrature( quadrature, values );
        }

        template< class Point >
        void jacobian ( const Point &x, JacobianRangeType &jacobian ) const
        {
          impl_.jacobian( x, jacobian );
        }

        template< class Quadrature, class Jacobians >
        void jacobianQuadrature ( const Quadrature &quadrature, Jacobians &jacobians ) const
        {
          impl_.jacobianQuadrature( quadrature, jacobians );
        }

        int order () const { return impl_.order(); }

        const EntityType &entity () const { return impl_.entity(); }

      private:
        Impl impl_;
        pybind11::object pyObj_;
      };

    public:
      PyGridFunction ( const GridFunction &impl, pybind11::object pyObj ) : impl_( &impl ), pyObj_( std::move( pyObj ) ) {}
      PyGridFunction ( const GridFunction &impl ) : impl_( &impl ), pyObj_( pybind11::detail::get_object_handle( impl_ ), true ) {}

      LocalFunctionType localFunction ( const EntityType &entity ) const { return LocalFunctionType( impl_->localFunction( entity ), pyObj_ ); }

      std::string name () const { return impl_->name(); }

      const GridPartType &gridPart () const { return impl_->gridPart(); }

      void evaluate ( const DomainType &x, RangeType &value ) const { return impl_->evaluate( x, value ); }
      void jacobian ( const DomainType &x, JacobianRangeType &jacobian ) const { return impl_->jacobian( x, jacobian ); }

    protected:
      const GridFunction *impl_;
      pybind11::object pyObj_;
    };



    // pyGridFunction
    // --------------

    template< class GridFunction >
    inline static PyGridFunction< GridFunction > pyGridFunction ( const GridFunction &gridFunction ) noexcept
    {
      return PyGridFunction< GridFunction >( gridFunction );
    }

    template< class GridFunction >
    inline static PyGridFunction< GridFunction > pyGridFunction ( const GridFunction &gridFunction, pybind11::object pyObj ) noexcept
    {
      return PyGridFunction< GridFunction >( gridFunction, std::move( pyObj ) );
    }

  } // namespace FemPy

} // namespace Dune

#endif // #ifndef DUNE_FEMPY_PY_FUNCTION_GRID_HH
