#ifndef DUNE_ADAPTCALLBACKHANDLE_HH
#define DUNE_ADAPTCALLBACKHANDLE_HH

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

    // flag that is set to true when at least one entity was coarsend or refined 
    mutable bool wasChanged_ ;

  public:
    typedef typename Base::Entity Entity;

    RestrictProlongWrapper ( DofManager &dofManager, RestrictProlongOperator &rpOp )
    : dofManager_( dofManager ),
      rpOp_( rpOp ),
      wasChanged_( false )
    {}

    RestrictProlongWrapper ( const RestrictProlongWrapper& org ) 
    : dofManager_( org.dofManager_ ), 
      rpOp_( org.rpOp_ ),
      wasChanged_( org.wasChanged_ )
    {}

    bool isValidEntity( const Entity& entity ) const
    {
      // grid was changed, if this method is called
      wasChanged_ = true ;

      // ghosts are not valid for restriction/prolongation
      return entity.partitionType() != GhostEntity ;
    }

    void preAdapt ( const unsigned int estimatedAdditionalElements )
    {
      // unset was changed 
      wasChanged_ = false;
      // reserve memory 
      dofManager_.reserveMemory( estimatedAdditionalElements );
    }

    void postAdapt ()
    {
      // is something was changed we need to call compress 
      if( wasChanged_ )
      {
        // make sure that no communication calls 
        // are done during DofManager::compress
        dofManager_.compress();

        // unset was changed flag
        wasChanged_ = false;
      }

      // make sequence counter globally equal 
      dofManager_.notifySequence();
    }

    void preCoarsening ( const Entity &father ) const
    {
      if( isValidEntity( father ) )
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
    }
    
    void restrictLocal ( const Entity &father, const Entity &son, bool initialize ) const
    {
      if( isValidEntity( father ) )
      {
        dofManager_.indexSetRestrictProlong().restrictLocal( const_cast< Entity & >( father ), const_cast< Entity & >( son ), initialize );
        rpOp_.restrictLocal( const_cast< Entity & >( father ), const_cast< Entity & >( son ), initialize );
      }
    }

    void postRefinement ( const Entity &father ) const
    {
      if( isValidEntity( father ) )
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
    }

    void prolongLocal ( const Entity &father, const Entity &son, bool initialize ) const
    {
      if( isValidEntity( father ) ) 
      {
        dofManager_.indexSetRestrictProlong().prolongLocal( const_cast< Entity & >( father ), const_cast< Entity & >( son ), initialize );
        rpOp_.prolongLocal( const_cast< Entity & >( father ), const_cast< Entity & >( son ), initialize );
      }
    }
  };

} // end namespace Dune 

#endif // #ifndef DUNE_ADAPTCALLBACKHANDLE_HH
