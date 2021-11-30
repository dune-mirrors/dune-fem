#warning "This header should not be needed anymore. Remove it from the include list!"

#ifndef DUNE_FEMPY_PY_FUNCTION_GRID_HH
#define DUNE_FEMPY_PY_FUNCTION_GRID_HH

// this class is only used in dune/fempy/py/discretefunction.hh - should be moved

#include <string>
#include <utility>

#include <dune/python/pybind11/pybind11.h>
#include <dune/common/visibility.hh>

#include <dune/fem/function/localfunction/bindable.hh>
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

    public:
      PyGridFunction ( const GridFunction &impl, pybind11::object pyObj )
        : impl_( &impl ), pyObj_( std::move( pyObj ) ),
          lf_(impl)
      {}

      PyGridFunction ( const GridFunction &impl )
        : Base(impl.gridPart(), impl.name(), impl.order()),
          impl_( &impl ),
          pyObj_( pybind11::reinterpret_borrow<pybind11::object>(
                  pybind11::detail::get_object_handle( impl_, pybind11::detail::get_type_info( typeid( GridFunction ) ) )
                ) ),
          lf_(impl)
      {}

      template< class Point >
      auto operator() ( const Point &x ) const
      {
        RangeType value;
        lf_.evaluate( x, value );
        return value;
      }
      template< class Point >
      void evaluate ( const Point &x, RangeType &value ) const
      {
        lf_.evaluate( x, value );
      }
      template< class Point >
      void jacobian ( const Point &x, JacobianRangeType &value ) const
      {
        lf_.jacobian( x, value );
      }
      template< class Point >
      void hessian ( const Point &x, HessianRangeType &value ) const
      {
        lf_.hessian( x, value );
      }

      template <class Entity>
      void bind(const Entity &entity) { Base::bind(entity); lf_.bind(entity); }
      void unbind() { Base::unbind(); lf_.unbind(); }
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
