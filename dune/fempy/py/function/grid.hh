#ifndef DUNE_FEMPY_PY_FUNCTION_GRID_HH
#define DUNE_FEMPY_PY_FUNCTION_GRID_HH

// this class is only used in dune/fempy/py/discretefunction.hh - should be moved

#include <string>
#include <utility>

#include <dune/common/visibility.hh>

#include <dune/fem/function/common/discretefunction.hh>
#include <dune/fem/space/common/functionspace.hh>

namespace Dune
{

  namespace FemPy
  {

    // PyGridFunction
    // --------------

    template< class GridFunction >
    class DUNE_PRIVATE PyGridFunction
      : public Fem::BindableGridFunctionWithSpace<typename GridFunction::GridPartType,
                         Dim<GridFunction::FunctionSpaceType::RangeType::dimension>>
    {
      // typedef Fem::Function< typename GridFunction::FunctionSpaceType, PyGridFunction< GridFunction > > Base;
      typedef Fem::BindableGridFunctionWithSpace<typename GridFunction::GridPartType,
                         Dim<GridFunction::FunctionSpaceType::RangeType::dimension>> Base;

    public:
      typedef typename Base::DiscreteFunctionSpaceType::FunctionSpaceType FunctionSpaceType;
      typedef typename GridFunction::GridPartType GridPartType;

      typedef typename GridFunction::EntityType EntityType;

      typedef typename FunctionSpaceType::DomainType DomainType;
      typedef typename FunctionSpaceType::RangeType RangeType;
      typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;
      typedef typename FunctionSpaceType::HessianRangeType HessianRangeType;

#if 0
      class DUNE_PRIVATE LocalFunctionType
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

        template< class Point >
        void hessian ( const Point &x, HessianRangeType &hessian ) const
        {
          impl_.hessian( x, hessian );
        }

        int order () const { return impl_.order(); }

        const EntityType &entity () const { return impl_.entity(); }

      private:
        Impl impl_;
        pybind11::object pyObj_;
      };
#endif

    public:
      PyGridFunction ( const GridFunction &impl, pybind11::object pyObj )
        : impl_( &impl ), pyObj_( std::move( pyObj ) ),
          lf_(impl)
      {}

      PyGridFunction ( const GridFunction &impl )
        : Base(impl.gridPart(), impl.name(), impl.order()),
          impl_( &impl ),
          // pyObj_( pybind11::detail::get_object_handle( impl_, pybind11::detail::get_type_info( typeid( GridFunction ) ) ), true ),
          pyObj_( pybind11::reinterpret_borrow<pybind11::object>(
                  pybind11::detail::get_object_handle( impl_, pybind11::detail::get_type_info( typeid( GridFunction ) ) )
                ) ),
          lf_(impl)
      {}

      template< class Point >
      void evaluate ( const Point &x, RangeType &value ) const
      {
        lf_.evaluate( x, value );
      }
#if 0
      template< class Quadrature, class Values >
      void evaluateQuadrature ( const Quadrature &quadrature, Values &values ) const
      {
        lf_.evaluateQuadrature( quadrature, values );
      }
#endif
      template< class Point >
      void jacobian ( const Point &x, JacobianRangeType &value ) const
      {
        lf_.jacobian( x, value );
      }
#if 0
      template< class Quadrature, class Jacobians >
      void jacobianQuadrature ( const Quadrature &quadrature, Jacobians &jacobians ) const
      {
        lf_.jacobianQuadrature( quadrature, jacobians );
      }
#endif
      template< class Point >
      void hessian ( const Point &x, HessianRangeType &value ) const
      {
        lf_.hessian( x, value );
      }

#if 0
      LocalFunctionType localFunction ( const EntityType &entity ) const { return LocalFunctionType( impl_->localFunction( entity ), pyObj_ ); }

      std::string name () const { return impl_->name(); }

      const GridPartType &gridPart () const { return impl_->gridPart(); }

      // !!!! void evaluate ( const DomainType &x, RangeType &value ) const { return impl_->evaluate( x, value ); }
      // !!!! void jacobian ( const DomainType &x, JacobianRangeType &jacobian ) const { return impl_->jacobian( x, jacobian ); }
#endif
      void bind(const EntityType &entity) { lf_.bind(entity); }
    protected:
      const GridFunction *impl_;
      pybind11::object pyObj_;
      Fem::ConstLocalFunction<GridFunction> lf_;
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
