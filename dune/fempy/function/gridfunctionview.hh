#ifndef DUNE_FEMPY_FUNCTION_GRIDFUNCTIONVIEW_HH
#define DUNE_FEMPY_FUNCTION_GRIDFUNCTIONVIEW_HH

#include <dune/python/function/simplegridfunction.hh>
#include <dune/fem/function/common/discretefunction.hh>
#include <dune/fem/function/localfunction/const.hh>

namespace Dune
{

  namespace FemPy
  {

    // GridFunctionView
    // ----------------

    template< class GF, bool isDiscreteFunction = std::is_base_of< Fem::IsDiscreteFunction, GF >::value >
    struct GridFunctionView;

    template< class GF >
    struct GridFunctionView< GF, false >
    {
      typedef typename GF::EntityType Entity;
      typedef typename GF::RangeType Value;

      typedef typename Entity::Geometry::LocalCoordinate LocalCoordinate;

      GridFunctionView ( const GF &gf ) : localFunction_( gf ) {}

      Value operator() ( const LocalCoordinate &x ) const
      {
        Value value;
        localFunction_.evaluate( x, value );
        return value;
      }

      void bind ( const Entity &entity ) { localFunction_.bind( entity ); }
      void unbind () {}

    private:
      Dune::Fem::ConstLocalFunction<GF> localFunction_;
    };

    template< class GF >
    struct GridFunctionView< GF, true >
    {
      typedef typename GF::EntityType Entity;
      typedef typename GF::RangeType Value;

      typedef typename Entity::Geometry::LocalCoordinate LocalCoordinate;

    private:
      typedef typename GF::DiscreteFunctionSpaceType DiscreteFunctionSpace;
      typedef typename DiscreteFunctionSpace::BasisFunctionSetType BasisFunctionSet;

    public:
      GridFunctionView ( const GF &gf )
        : gf_( gf )
      {
        localDofVector_.reserve( DiscreteFunctionSpace::localBlockSize * space().blockMapper().maxNumDofs() );
      }

      Value operator() ( const LocalCoordinate &x ) const
      {
        Value value;
        basisFunctionSet_.evaluateAll( x, localDofVector_, value );
        return value;
      }

      void bind ( const Entity &entity )
      {
        basisFunctionSet_ = space().basisFunctionSet( entity );
        localDofVector_.resize( basisFunctionSet_.size() );
        gf_.getLocalDofs( entity, localDofVector_ );
      }

      void unbind ()
      {
        basisFunctionSet_ = BasisFunctionSet();
        localDofVector_.resize( 0u );
      }

    private:
      const DiscreteFunctionSpace &space () const { return gf_.space(); }

      const GF &gf_;
      BasisFunctionSet basisFunctionSet_;
      DynamicVector< typename GF::DofType > localDofVector_;
    };



    // localFunction
    // -------------

    template< class GF, std::enable_if_t< std::is_base_of< Fem::HasLocalFunction, GF >::value, int > = 0 >
    inline static GridFunctionView< GF > localFunction ( const GF &gf )
    {
      return GridFunctionView< GF >( gf );
    }

  } // namespace FemPy


  namespace Fem
  {

    using FemPy::localFunction;

    namespace Impl
    {
      template< class GV, class Evaluate >
      struct ConstLocalFunction< Dune::Python::SimpleGridFunction<GV,Evaluate> >
      {
        struct Type
        {
          typedef Dune::Python::SimpleGridFunction<GV,Evaluate> GridFunctionType;
          typedef decltype(localFunction(std::declval<const GridFunctionType&>())) LocalFunctionType;
          typedef typename LocalFunctionType::Element EntityType;

          typedef typename LocalFunctionType::LocalCoordinate DomainType;
          typedef typename LocalFunctionType::Value RangeType;
          // typedef typename LocalFunctionType::JacobianRangeType JacobianRangeType;
          // typedef typename LocalFunctionType::HessianRangeType HessianRangeType;

          explicit Type ( const GridFunctionType &gridFunction )
              : localFunction_( localFunction(gridFunction) ),
                gridFunction_( gridFunction )
            {}
            explicit Type ( const EntityType &entity, const GridFunctionType &gridFunction )
              : localFunction_( localFunction(gridFunction) ),
                gridFunction_( gridFunction )
            { bind(entity); }

            //! evaluate local function
            template< class Point >
            RangeType operator() ( const Point &p ) const
            {
              return localFunction_( Fem::coordinate(p) );
            }
            template< class Point >
            RangeType evaluate ( const Point &p ) const
            {
              return localFunction_( Fem::coordinate(p) );
            }
            //! evaluate local function
            template< class Point >
            void evaluate ( const Point &p, RangeType &val ) const
            {
              val = localFunction_( Fem::coordinate(p) );
            }
            template< class Point, class JacobianRangeType >
            void jacobian ( const Point &p, JacobianRangeType &val ) const
            {
              DUNE_THROW(NotImplemented,"SimpleLocalFunction does not provide a jocobian!");
            }
            template< class Point, class HessianRangeType >
            void hessian ( const Point &p, HessianRangeType &val ) const
            {
              DUNE_THROW(NotImplemented,"SimpleLocalFunction does not provide a hessian!");
            }

            void bind ( const EntityType &entity ) { localFunction_.bind( entity ); }
            void unbind () {}

            const GridFunctionType &gridFunction () const { return gridFunction_; }

          private:
            const GridFunctionType &gridFunction_;
            LocalFunctionType localFunction_;
        };
      };
    } // namespace Impl

  } // namespace Fem

  namespace Python
  {

    using FemPy::localFunction;

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEMPY_FUNCTION_GRIDFUNCTIONVIEW_HH
