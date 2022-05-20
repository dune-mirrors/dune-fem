#ifndef DUNE_FEM_FUNCTION_LOCALFUNCTION_BINDABLE_HH
#define DUNE_FEM_FUNCTION_LOCALFUNCTION_BINDABLE_HH

#include <dune/fem/space/common/discretefunctionspace.hh>
#include <dune/fem/function/common/discretefunction.hh>
#include <dune/fem/common/coordinate.hh>
#include <dune/fem/quadrature/quadrature.hh> // shouldn't be here (but otherwise the coordinate method doesn't work)
#include <dune/fem/common/intersectionside.hh>

namespace Dune
{
  namespace Fem
  {
    struct BindableFunction : public HasLocalFunction {};

    template <class GridPart, class Range>
    struct BindableGridFunction : public BindableFunction
    {
      typedef GridPart GridPartType;
      typedef typename GridPart::template Codim<0>::EntityType   EntityType;
      typedef typename GridPart::IntersectionType                IntersectionType;
      typedef typename EntityType::Geometry                      Geometry;
      typedef typename Geometry::GlobalCoordinate                DomainType;
      typedef Dune::Fem::GridFunctionSpace<GridPartType, Range>  FunctionSpaceType;
      typedef typename FunctionSpaceType::RangeFieldType         RangeFieldType;
      typedef typename FunctionSpaceType::RangeType              RangeType;
      typedef typename FunctionSpaceType::JacobianRangeType      JacobianRangeType;
      typedef typename FunctionSpaceType::HessianRangeType       HessianRangeType;
      BindableGridFunction(const GridPart &gridPart)
      : gridPart_(gridPart)
#ifdef TESTTHREADING
      , thread_(-1)
#endif
      {}

      void bind(const EntityType &entity)
      {
#ifdef TESTTHREADING
        if (thread_==-1) thread_ = MPIManager::thread();
        if (thread_ != MPIManager::thread())
        {
          std::cout << "wrong thread number\n";
          assert(0);
          std::abort();
        }
        if (entity_ || geometry_)
        {
          std::cout << "BindableGF: bind called on object before unbind was called\n";
          std::abort();
        }
        assert(!entity_ && !geometry_); // this will fail with dune-fem-dg
#endif
        unbind();

        entity_.emplace( entity );
        geometry_.emplace( this->entity().geometry() );
      }

      void unbind()
      {
#ifdef TESTTHREADING
        if (thread_ != MPIManager::thread())
        {
          std::cout << "wrong thread number\n";
          assert(0);
          std::abort();
        }
#endif
        geometry_.reset();
        entity_.reset();
      }

      void bind(const IntersectionType &intersection, IntersectionSide side)
      {
        // store local copy to avoid problems with casting to temporary types
        const EntityType entity = side==IntersectionSide::in? intersection.inside(): intersection.outside();
        bind( entity );
      }

      bool continuous() const { return true; }
      template <class Point>
      DomainType global(const Point &x) const
      {
        return geometry_.value().global( Dune::Fem::coordinate(x) );
      }

      // this method needs to be overloaded in the derived class
      template <class Point>
      void evaluate( const Point& x, RangeType& ret ) const;

      template <class Quadrature, class RangeArray>
      void evaluate( const Quadrature& quadrature, RangeArray& values ) const
      {
        const unsigned int nop = quadrature.nop();
        values.resize( nop );
        for( unsigned int qp=0; qp<nop; ++qp)
        {
          evaluate( quadrature[ qp ], values[ qp ]);
        }
      }

      const GridPart& gridPart() const { return gridPart_; }
      const EntityType &entity() const { return entity_.value(); }
      const Geometry& geometry() const { return geometry_.value(); }

    protected:
      std::optional< EntityType > entity_;
      std::optional< Geometry > geometry_;
      const GridPart &gridPart_;
#ifdef TESTTHREADING
      int thread_;
#endif
    };

    template <class GridPart, class Range>
    struct BindableGridFunctionWithSpace : public BindableGridFunction<GridPart,Range>
    {
      typedef BindableGridFunction<GridPart,Range> Base;
      typedef GridPart GridPartType;
      typedef typename GridPart::template Codim<0>::EntityType EntityType;
      typedef typename EntityType::Geometry::GlobalCoordinate DomainType;
      typedef Dune::Fem::GridFunctionSpace<GridPartType, Range> FunctionSpaceType;
      typedef typename FunctionSpaceType::RangeFieldType RangeFieldType;
      typedef typename FunctionSpaceType::RangeType RangeType;
      typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;
      typedef typename FunctionSpaceType::HessianRangeType HessianRangeType;
      typedef DiscreteFunctionSpaceAdapter< FunctionSpaceType, GridPartType > DiscreteFunctionSpaceType;
      BindableGridFunctionWithSpace(const GridPart &gridPart, const std::string &name, int order)
      : Base(gridPart),
        space_( gridPart, order ),
        name_(name)
      {}
      //! return the order of the space
      unsigned int order() const
      {
        return space().order();
      }
      const std::string &name() const
      {
        return name_;
      }
      const DiscreteFunctionSpaceType &space () const
      {
        return space_;
      }
      private:
      DiscreteFunctionSpaceType space_;
      const std::string name_;
    };

    namespace detail
    {
      template <class,class,class>
      struct canBind
         : std::false_type {};
      template <class GP,class LF>
      struct canBind<GP,LF,
             std::void_t< decltype( std::declval<LF>().
                            bind(std::declval<const typename GP::template Codim<0>::EntityType&>())) >>
        : std::true_type {};
      template <class GP,class LF>
      using canBind_t = canBind<GP,LF,void>;
    }

    template <class GP,class LF>
    constexpr detail::canBind_t<GP,LF> checkGridPartValid() { return {}; }

  } // namespace Fem
} // namespace Dune
#endif // DUNE_FEM_FUNCTION_LOCALFUNCTION_BINDABLE_HH
