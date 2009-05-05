#ifndef DUNE_ADAPTCALLBACKHANDLE_HH
#define DUNE_ADAPTCALLBACKHANDLE_HH

#include <dune/common/version.hh>
#include <dune/grid/common/adaptcallback.hh>

/** \file
 *  \author Martin Nolte
 *  \brief  interfaces and wrappers needed for the callback adaptation provided
 *          by AlbertaGrid and ALUGrid
 */

namespace Dune
{

  // RestrictProlongWrapper
  // ----------------------

#if DUNE_VERSION_NEWER(DUNE_GRID,1,3,0)
  template< class Grid, class DofManager, class RestrictProlongOperator >
  class RestrictProlongWrapper
  : public AdaptDataHandle
    < Grid, RestrictProlongWrapper< Grid, DofManager, RestrictProlongOperator > >
  {
    typedef RestrictProlongWrapper< Grid, DofManager, RestrictProlongOperator > This;
    typedef AdaptDataHandle< Grid, This > Base;

  protected:  
    DofManager &dofManager_;
    RestrictProlongOperator &rpOp_;

  public:
    typedef typename Base::Entity Entity;

    RestrictProlongWrapper ( DofManager &dofManager, RestrictProlongOperator &rpOp )
    : dofManager_( dofManager ),
      rpOp_( rpOp )
    {}

    RestrictProlongWrapper ( const RestrictProlongWrapper& org ) 
    : dofManager_( org.dofManager_ ), 
      rpOp_( org.rpOp_ )
    {}

    void preAdapt ( const size_t estimatedAdditionalElements )
    {
      dofManager_.reserveMemory( estimatedAdditionalElements );
    }

    void postAdapt ()
    {
      dofManager_.compress();
    }

    void preCoarsening ( const Entity &father ) const
    {
      typedef typename Entity::HierarchicIterator HIterator;

      bool initialize = true;
      const int childLevel = father.level() + 1;
      const HIterator end = father.hend( childLevel );
      for( HIterator it = father.hbegin( childLevel ); it != end; ++it )
      {
        restrictLocal( father, *it, initialize );
        initialize = false;
      }
    }
    
    void restrictLocal ( const Entity &father, const Entity &son, bool initialize ) const
    {
      dofManager_.indexSetRestrictProlong().restrictLocal( const_cast< Entity & >( father ), const_cast< Entity & >( son ), initialize );
      rpOp_.restrictLocal( const_cast< Entity & >( father ), const_cast< Entity & >( son ), initialize );
    }

    void postRefinement ( const Entity &father ) const
    {
      typedef typename Entity::HierarchicIterator HIterator;

      bool initialize = true;
      const int childLevel = father.level() + 1;
      const HIterator end = father.hend( childLevel );
      for( HIterator it = father.hbegin( childLevel ); it != end; ++it )
      {
        prolongLocal( father, *it, initialize );
        initialize = false;
      }
    }

    void prolongLocal ( const Entity &father, const Entity &son, bool initialize ) const
    {
      dofManager_.indexSetRestrictProlong().prolongLocal( const_cast< Entity & >( father ), const_cast< Entity & >( son ), initialize );
      rpOp_.prolongLocal( const_cast< Entity & >( father ), const_cast< Entity & >( son ), initialize );
    }
  };
#endif

} // end namespace Dune 

#endif // #ifndef DUNE_ADAPTCALLBACKHANDLE_HH
