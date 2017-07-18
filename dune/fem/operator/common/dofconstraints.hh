#ifndef DUNE_FEM_OPERATOR_COMMON_DOFCONSTRAINTS
#define DUNE_FEM_OPERATOR_COMMON_DOFCONSTRAINTS

#include <dune/fem/function/common/rangegenerators.hh>
#include <dune/fem/common/bindguard.hh>
#include <dune/fem/function/common/localcontribution.hh>
#include <dune/fem/space/common/interpolate.hh>
#include <dune/fem/function/localfunction/const.hh>

namespace Dune
{
  namespace Fem
  {
    /***********************************************************
     * Constrain class interface
     *   bind(Entity);
     *   unbind();
     *   bool operator()(unsinged int localDof);
     ***********************************************************/

    template <class DiscreteFunctionSpace, class Constrain>
    struct DofConstraints
    {
      typedef Constrain ConstrainType;
      typedef DiscreteFunctionSpace DiscreteFunctionSpaceType;
      typedef Dune::Fem::AdaptiveDiscreteFunction<DiscreteFunctionSpace> InternalStorageType;
      DofConstraints(const DiscreteFunctionSpaceType &space, Constrain &constrain)
      : space_( space ),
        constrain_( constrain ),
        values_( "costraint_values", space_ ),
        mask_( "dof_mask", space_ ),
        sequence_( -1 )
      { update(); }

      template <class GridFunction>
      void set(const GridFunction &gf)
      {
        checkUpdate();
        interpolate( gf, values_, constrain_, mask_ );
      }
      template <class DiscreteFunction>
      void operator()(DiscreteFunction &df)
      {
        (*this)(values_,df);
      }
      template <class LinearOperator>
      void applyToOperator(LinearOperator &lo) // could use operator() with some SFINAE
      {
        checkUpdate();
        std::vector<unsigned int> rows;
        rows.reserve(mask_.size()); // need something better
        unsigned int r=0;
        for ( const auto &m : dofs(mask_) )
        {
          if (m>0)
            rows.push_back(r);
          ++r;
        }
        lo.setUnitRows( rows );
      }
      template <class From, class To>
      void operator()(const From &from, To &to)
      {
        checkUpdate();

        auto fromIt = from.dbegin();
        auto maskIt = mask_.dbegin();
        int idx = 0;
        for ( auto&& dof : dofs(to) )
        {
          if ( (*maskIt) > 0)
            dof = (*fromIt);
          assert( maskIt != mask_.dend() && fromIt != from.dend() );
          ++maskIt;
          ++fromIt;
          ++idx;
        }
      }
      void clear()
      {
        checkUpdate();
      }

      void update()
      {
        mask_.clear();
        // should use a ZeroDF here
        interpolate( values_, values_, constrain_, mask_ );
        values_.clear();
        sequence_ = space_.sequence();
      }
      private:
      void checkUpdate()
      {
        if (space_.sequence() != sequence_)
          update();
      }
      const DiscreteFunctionSpaceType &space_;
      Constrain &constrain_;
      InternalStorageType values_;
      InternalStorageType mask_;
      int sequence_;
    };

    ///////////////////////////////////////////////////////////

    template <unsigned int dimR>
    struct EmptyModel
    {
      template <class I>
      std::pair< unsigned int, std::bitset< dimR > >
      isDirichlet ( const I& )
      {
        std::bitset<dimR> components;
        return std::make_pair(1, components.set());
      }
      template <class E>
      void bind ( const E& ) {}
      void unbind ( ) {}
    };
    template <class Space, class Model>
    struct ConstrainOnBoundary
    {
      typedef typename Space::EntityType EntityType;
      ConstrainOnBoundary(const Space& space, Model &model)
      : space_(space), model_(model) {}

      template <class DofVector>
      void operator()(const EntityType &entity, DofVector &localMask)
      {
        auto modelGuard = bindGuard( model_, entity );

        const unsigned int numDofs = space_.blockMapper().numDofs( entity );
        std::fill( localMask.begin(), localMask.end(), 0. );

        // vector for subentity filter
        std::vector< bool > onSubEntityFilter( numDofs );

        for( const auto &intersection : intersections( space_.gridPart(), entity ) )
        {
          // only need to do something if intersection is with boundary
          if( !intersection.boundary() )
            continue;
          // and model says that it is part of the dirichlet boundary
          auto bnd = model_.isDirichlet( intersection );
          if ( bnd.second.none() )
            continue;

          // get face number of boundary intersection
          const int face = intersection.indexInInside();
          space_.blockMapper().onSubEntity( entity, face, 1, onSubEntityFilter );
          for( unsigned int i = 0; i < numDofs; ++i )
            for( unsigned int r = 0; r < Space::localBlockSize; ++r )
              localMask[ i*Space::localBlockSize+r ] =
                ( (localMask[ i*Space::localBlockSize+r ]==1.) |
                  (onSubEntityFilter[ i ] & bnd.second[r]) ) ? 1.:0.;
        }
      }

      private:
      const Space &space_;
      Model &model_;
    };
    template <class Space>
    struct ConstrainOnFullBoundary
    : public ConstrainOnBoundary< Space, EmptyModel<Space::dimRange> >
    {
      ConstrainOnFullBoundary(const Space& space)
      : ConstrainOnBoundary< Space, EmptyModel<Space::dimRange> >(space,model_) {}
      EmptyModel<Space::dimRange> model_;
    };

    /////////////////////////////////////////////////////////////////

#if HAVE_DUNE_LOCALFUNCTIONS
    template <class Model>
    struct PiecewiseGridFunction
    {
      typedef typename Model::GridPartType GridPartType;
      typedef typename Model::EntityType EntityType;
      typedef typename Model::RangeType RangeType;
      typedef typename Model::FunctionSpaceType FunctionSpaceType;
      typedef Dune::Fem::LagrangeSpace< FunctionSpaceType, GridPartType > DiscreteFunctionSpaceType;
      typedef std::vector<double> TemporaryStorage;
      typedef typename DiscreteFunctionSpaceType::BasisFunctionSetType BasisFunctionSetType;

      PiecewiseGridFunction(Model &model, GridPartType &gridPart, unsigned int order)
        : model_(model),
          space_(gridPart, order),
          localFunction_( space_.maxNumDofs() ),
          ltmp_( space_.maxNumDofs() )
      {}

      void init(const EntityType &entity)
      {
        bind(entity);
      }
      void bind(const EntityType &entity)
      {
        basisFunctionSet_ = space_.basisFunctionSet( entity );
        model_.bind(entity);

        const int localBlocks = space_.blockMapper().numDofs( entity );
        std::vector< bool > onSubEntity( localBlocks );

        const auto interpolation = space_.interpolation( entity );
        for( const auto &intersection : intersections( space_.gridPart(), entity ) )
        {
          // if intersection is not with boundary, skip
          if( !intersection.boundary() )
            continue;

          // get dirichlet information from model
          // Remark: assuming for now that block[i] is zero or always has the same value
          auto bnd = model_.isDirichlet( intersection );
          if( bnd.second.none() )
            continue;

          // get face number of boundary intersection
          const int face = intersection.indexInInside();
          space_.blockMapper().onSubEntity( entity, face, 1, onSubEntity );
          interpolation( BoundaryWrapper( model_, bnd.first ), ltmp_ );

          int localDof = 0;
          for( int localBlock = 0; localBlock < localBlocks; ++localBlock )
          {
            for( int l = 0; l < DiscreteFunctionSpaceType::localBlockSize; ++l, ++localDof )
            {
              if (bnd.second[l]==0) continue;
              localFunction_[ localDof ] = ltmp_[ localDof ];
            }
          }
        }
      }
      void unbind()
      {
        model_.unbind();
      }

      template <class Point>
      void evaluate(const Point &x, RangeType &value) const
      {
        basisFunctionSet_.evaluateAll( x, localFunction_, value );
      }

      private:
      struct BoundaryWrapper
      {
        typedef typename DiscreteFunctionSpaceType::DomainType DomainType;
        typedef typename DiscreteFunctionSpaceType::RangeType RangeType;
        typedef typename DiscreteFunctionSpaceType::JacobianRangeType JacobianRangeType;
        BoundaryWrapper( const Model& impl, int bndId )
        : impl_( impl ), bndId_(bndId) {}

        //! evaluate function
        template <class Point>
        void evaluate( const Point& x, RangeType& ret ) const
        {
          ret = impl_.dirichlet(bndId_,Dune::Fem::coordinate(x));
        }
        private:
        const Model& impl_;
        int bndId_;
      };

      Model &model_;
      const DiscreteFunctionSpaceType space_;
      TemporaryStorage localFunction_;
      BasisFunctionSetType basisFunctionSet_;
      TemporaryStorage ltmp_;
    };
#endif // HAVE_DUNE_LOCALFUNCTIONS

    template <class Model, class Space>
    using DirichletConstraints =
      DofConstraints< Space, ConstrainOnBoundary< Space, Model > >;

  }
}

#endif // DUNE_FEM_OPERATOR_COMMON_DOFCONSTRAINTS
