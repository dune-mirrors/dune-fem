#ifndef DUNE_FEM_OPERATOR_COMMON_DOFCONSTRAINTS
#define DUNE_FEM_OPERATOR_COMMON_DOFCONSTRAINTS

#include <dune/fem/function/common/rangegenerators.hh>
#include <dune/fem/common/bindguard.hh>
#include <dune/fem/function/common/localcontribution.hh>

/***********************************************************
 * Constrain class interface
 *   bind(Entity);
 *   unbind();
 *   bool operator()(unsinged int localDof);
 ***********************************************************/

template <class DiscreteFunctionSpace, class Constrain>
struct DofConstraints
{
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
    checkUpdate();

    typedef typename InternalStorageType::ConstDofIteratorType ConstDofIteratorType;
    ConstDofIteratorType valueIt = values_.dbegin();
    ConstDofIteratorType maskIt = mask_.dbegin();
    int idx = 0;
    for ( auto& dof : dofs(df) )
    {
      if ( (*maskIt) > 0)
        dof = (*valueIt);
      assert( maskIt != mask_.dend() && valueIt != values_.dend() );
      ++maskIt;
      ++valueIt;
      ++idx;
    }
  }
  template <class LinearOperator>
  void operator()(const LinearOperator &df)
  {
    checkUpdate();
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

template <class Space>
struct ConstrainOnBoundary
{
  typedef typename Space::EntityType EntityType;
  ConstrainOnBoundary(const Space& space)
  : space_(space), localMask_(0) {}

  bool operator()(unsigned int localDof) const
  {
    assert( localDof < localMask_.size() );
    return localMask_[localDof];
  }

  void bind(const EntityType &entity)
  {
    const unsigned int numDofs = space_.blockMapper().numDofs( entity );
    localMask_.resize( numDofs*Space::localBlockSize );
    std::fill( localMask_.begin(), localMask_.end(), false );

    // vector for subentity filter
    std::vector< bool > globalBlockDofsFilter( numDofs );

    for( const auto &intersection : intersections( space_.gridPart(), entity ) )
    {
      // only need to do something if intersection is with boundary
      if( !intersection.boundary() )
        continue;

      // get face number of boundary intersection
      const int face = intersection.indexInInside();
      space_.blockMapper().onSubEntity( entity, face, 1, globalBlockDofsFilter );
      // Q: is this general enough?
      for( unsigned int i = 0; i < numDofs; ++i )
        for( unsigned int r = 0; r < Space::localBlockSize; ++r )
          localMask_[ i*Space::localBlockSize+r ] =
            localMask_[ i*Space::localBlockSize+r ] |
            globalBlockDofsFilter[ i ];
    }
  }
  void unbind() {}

  private:
  const Space &space_;
  std::vector<bool> localMask_;
};

#endif // DUNE_FEM_OPERATOR_COMMON_DOFCONSTRAINTS
