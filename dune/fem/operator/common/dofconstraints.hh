#ifndef DUNE_FEM_OPERATOR_COMMON_DOFCONSTRAINTS
#define DUNE_FEM_OPERATOR_COMMON_DOFCONSTRAINTS

#include <dune/fem/function/common/rangegenerators.hh>
#include <dune/fem/common/bindguard.hh>
#include <dune/fem/function/common/localcontribution.hh>
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
      : space_(space),
        constrain_(constrain),
        values_( "costraint_values", space_),
        mask_( "dof_mask", space_),
        sequence_( -1 )
      { update(); }

      template <class GridFunction>
      void set(const GridFunction &gf)
      {
        checkUpdate();
        values_.clear();
        {
          Dune::Fem::ConstLocalFunction< GridFunction > lgf( gf );
          Dune::Fem::AddLocalContribution< InternalStorageType > lvalues( values_ );

          for ( const auto& entity : space_ )
          {
            auto gfGuard = bindGuard( lgf, entity );
            auto valuesGuard = bindGuard( lvalues, entity );
            const auto interpolation = space_.interpolation( entity );
            interpolation( lgf, lvalues );
          }
        }

        typedef typename InternalStorageType::DofType DofType;
        // divide DoFs by the mask
        std::transform( values_.dbegin(), values_.dend(), mask_.dbegin(), values_.dbegin(),
            [] ( DofType u, DofType w ) {
            using std::abs;
            typename Dune::FieldTraits< DofType >::field_type weight = abs( w );
            return (weight > 1e-12 ? u / weight : 0.);
            } );
      }
      template <class DiscreteFunction>
      void operator()(DiscreteFunction &df)
      {
        (*this)(values_,df);
      }
      template <class LinearOperator>
      void applyToOperator(LinearOperator &df)
      {
        checkUpdate();
#if 0
        Dune::Fem::LocalContribution<LinearOperator> localMatrix;
        // if Dirichlet Dofs have been found, treat them
        for( const auto &entity : space_ )
        {
          auto localMatrixGuard = bindGuard( localMatrix, entity, entity );
          // get slave dof structure (for parallel runs)
          const auto &slaveDofs = linearOperator.rangeSpace().slaveDofs();
          // get number of basis functions
          const int localBlocks = space_.blockMapper().numDofs( entity );
          // map local to global dofs
          std::vector< std::size_t > globalBlockDofs( localBlocks );
          // obtain all DofBlocks for this element
          space_.blockMapper().map( entity, globalBlockDofs );
          // counter for all local dofs (i.e. localBlockDof * localBlockSize + ... )
          int localDof = 0;
          // iterate over face dofs and set unit row
          for( int localBlockDof = 0; localBlockDof < localBlocks; ++localBlockDof )
          {
            int global = globalBlockDofs[ localBlockDof ];
            for( int l = 0; l < localBlockSize; ++l, ++localDof )
            {
              if( !dirichletBlocks_[ global ][ l ] )
                continue;
              // clear all other columns
              localMatrix.clearRow( localDof );
              // set diagonal to 1
              double value = slaveDofs.isSlave( global ) ? 0.0 : 1.0;
              localMatrix.set( localDof, localDof, value );
            }
          }
        }
#endif
      }
      template <class From, class To>
      void operator()(const From &from, To &to)
      {
        checkUpdate();

        auto valueIt = from.dbegin();
        auto maskIt = mask_.dbegin();
        int idx = 0;
        for ( auto& dof : dofs(to) )
        {
          if ( (*maskIt) > 0)
            dof = (*valueIt);
          assert( maskIt != mask_.dend() && valueIt != values_.dend() );
          ++maskIt;
          ++valueIt;
          ++idx;
        }
      }
      void clear()
      {
        checkUpdate();
        values_.clear();
      }

      void update()
      {
        mask_.clear();
        values_.clear();

        Dune::Fem::AddLocalContribution< InternalStorageType > lmask( mask_ );

        for ( const auto& entity : space_ )
        {
          // initialize local function and local weight
          auto maskGuard = bindGuard( lmask, entity );
          auto constraintGuard = bindGuard( constrain_, entity );
          for (unsigned int dof = 0; dof < lmask.size(); ++dof)
            if ( constrain_(dof) )
              lmask[dof] = 1.;
        }
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

    struct EmptyModel
    {
      static const unsigned int dimRange = 1;
      template <class I, class V>
      bool isDirichletIntersection (const I&, const V&) { return true; }
      template <class E>
      void bind ( const E& ) {}
      void unbind ( ) {}
    };
    template <class Space, class Model>
    struct ConstrainOnBoundary
    {
      typedef typename Space::EntityType EntityType;
      ConstrainOnBoundary(const Space& space, Model &model)
      : space_(space), model_(model), localMask_(0) {}

      bool operator()(unsigned int localDof) const
      {
        assert( localDof < localMask_.size() );
        return localMask_[localDof];
      }

      void bind(const EntityType &entity)
      {
        auto modelGuard = bindGuard( model_, entity );

        const unsigned int numDofs = space_.blockMapper().numDofs( entity );
        localMask_.resize( numDofs*Space::localBlockSize );
        std::fill( localMask_.begin(), localMask_.end(), false );

        // vector for subentity filter
        std::vector< bool > onSubEntityFilter( numDofs );

        FieldVector< int, Model::dimRange > block( 0 );
        for( const auto &intersection : intersections( space_.gridPart(), entity ) )
        {
          // only need to do something if intersection is with boundary
          if( !intersection.boundary() )
            continue;
          // and model says that it is part of the dirichlet boundary
          if( !model_.isDirichletIntersection( intersection, block ) )
            continue;

          // get face number of boundary intersection
          const int face = intersection.indexInInside();
          space_.blockMapper().onSubEntity( entity, face, 1, onSubEntityFilter );
          // Q: is this general enough?
          for( unsigned int i = 0; i < numDofs; ++i )
            for( unsigned int r = 0; r < Space::localBlockSize; ++r )
              localMask_[ i*Space::localBlockSize+r ] =
                localMask_[ i*Space::localBlockSize+r ] |
                (onSubEntityFilter[ i ] & block[i]>0);
        }
      }
      void unbind() {}

      private:
      const Space &space_;
      Model &model_;
      std::vector<bool> localMask_;
    };
    template <class Space>
    struct ConstrainOnFullBoundary : public ConstrainOnBoundary<Space,EmptyModel>
    {
      ConstrainOnFullBoundary(const Space& space)
      : ConstrainOnBoundary<Space,EmptyModel>(space,model_) {}
      EmptyModel model_;
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
          FieldVector< int, RangeType::dimension > block( 0 );
          if( !model_.isDirichletIntersection( intersection, block ) )
            continue;

          // get face number of boundary intersection
          const int face = intersection.indexInInside();
          space_.blockMapper().onSubEntity( entity, face, 1, onSubEntity );
          interpolation( BoundaryWrapper( model_, block[0] ), ltmp_ );

          int localDof = 0;
          for( int localBlock = 0; localBlock < localBlocks; ++localBlock )
          {
            for( int l = 0; l < DiscreteFunctionSpaceType::localBlockSize; ++l, ++localDof )
            {
              if (block[l]==0) continue;
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
