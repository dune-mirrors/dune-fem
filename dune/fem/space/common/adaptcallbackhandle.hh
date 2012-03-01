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
      wasChanged_( false )
    {}

    void preAdapt ( const unsigned int estimatedAdditionalElements )
    {
      // unset was changed 
      wasChanged_ = false;
      dofManager_.reserveMemory( estimatedAdditionalElements );
    }

    void postAdapt ()
    {
      // is something was changed we need to call compress 
      // don't call any communication in DofManager::compress
      if( wasChanged_ )
      {
        dofManager_.compress();
        // unset was changed 
        wasChanged_ = false;
      }

      // call dofmanger finalize to make flags know globally 
      dofManager_.globalFinalize();
    }

    void preCoarsening ( const Entity &father ) const
    {
      wasChanged_ = true ;

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
      wasChanged_ = true ;

      dofManager_.indexSetRestrictProlong().restrictLocal( const_cast< Entity & >( father ), const_cast< Entity & >( son ), initialize );
      rpOp_.restrictLocal( const_cast< Entity & >( father ), const_cast< Entity & >( son ), initialize );
    }

    void postRefinement ( const Entity &father ) const
    {
      wasChanged_ = true ;

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
      wasChanged_ = true ;

      dofManager_.indexSetRestrictProlong().prolongLocal( const_cast< Entity & >( father ), const_cast< Entity & >( son ), initialize );
      rpOp_.prolongLocal( const_cast< Entity & >( father ), const_cast< Entity & >( son ), initialize );
    }
  };

} // end namespace Dune 

#endif // #ifndef DUNE_ADAPTCALLBACKHANDLE_HH
