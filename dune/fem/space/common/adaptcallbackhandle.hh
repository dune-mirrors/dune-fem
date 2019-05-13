#ifndef DUNE_FEM_ADAPTCALLBACKHANDLE_HH
#define DUNE_FEM_ADAPTCALLBACKHANDLE_HH

#include <dune/grid/common/gridenums.hh>
#include <dune/grid/common/adaptcallback.hh>

/** \file
 *  \author Martin Nolte
 *  \brief  interfaces and wrappers needed for the callback adaptation provided
 *          by AlbertaGrid and ALUGrid
 */

namespace Dune
{

  namespace Fem
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
      bool initializeCalled_;
      bool finalizeCalled_;

    public:
      typedef typename Base::Entity Entity;

      RestrictProlongWrapper ( DofManager &dofManager, RestrictProlongOperator &rpOp )
      : dofManager_( dofManager ),
        rpOp_( rpOp ),
        wasChanged_( false ),
        initializeCalled_( false ),
        finalizeCalled_( false )
      {
      }

      RestrictProlongWrapper ( const RestrictProlongWrapper& org )
      : dofManager_( org.dofManager_ ),
        rpOp_( org.rpOp_ ),
        wasChanged_( org.wasChanged_ ),
        initializeCalled_( org.initializeCalled_ ),
        finalizeCalled_( org.finalizeCalled_ )
      {}

      bool isValidEntity( const Entity& entity ) const
      {
        // grid was changed, if this method is called
        wasChanged_ = true ;

        // ghosts are not valid for restriction/prolongation
        assert( entity.partitionType() != GhostEntity );
        return true ;
      }

      // old interface methods
      void preAdapt ( const unsigned int estimatedAdditionalElements ) { initialize (); }
      void postAdapt () { finalize(); }

      /** \brief initialize basically reserves some memory on the DofManager */
      void initialize ( unsigned int estimatedAdditionalElements = 0 )
      {
        // if preAdapt was already called just return
        if( initializeCalled_ ) return ;

        // initialize restrict-prolong object
        rpOp_.initialize();

        // unset was changed
        wasChanged_ = false;
        // reserve memory
        dofManager_.reserveMemory( estimatedAdditionalElements );

        // set initializeCalled_ flag in case method is called again (only dune-grid version)
        initializeCalled_ = true;
        // reset postAdaptCalled flag
        finalizeCalled_ = false ;

      }

      /** \brief finalize calls the compress on the DofManager */
      void finalize ()
      {
        // if method has been called already do nothing
        if( finalizeCalled_ ) return ;

        // notifyGlobalChange make wasChanged equal on all cores
        if( dofManager_.notifyGlobalChange( wasChanged_ ) )
        {
          // make sure that no communication calls
          // are done during DofManager::compress
          dofManager_.compress();

          // unset was changed flag
          wasChanged_ = false;
        }

        // initialize restrict-prolong object
        rpOp_.finalize();

        // set postAdaptCalled flag
        finalizeCalled_ = true ;

        // reset initializeCalled_ flag
        initializeCalled_ = false ;
      }

      /** \copydoc Dune::AdaptDataHandleInterface::preCoarsening */
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
          rpOp_.restrictFinalize( father );
        }
      }

      void restrictLocal ( const Entity &father, const Entity &son, bool initialize ) const
      {
        if( isValidEntity( father ) )
        {
          dofManager_.indexSetRestrictProlong().restrictLocal( father, son, initialize );
          rpOp_.restrictLocal( father, son, initialize );
        }
      }

      /** \copydoc Dune::AdaptDataHandleInterface::postRefinement */
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
          dofManager_.indexSetRestrictProlong().prolongLocal( father, son, initialize );
          rpOp_.prolongLocal( father, son, initialize );
        }
      }
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_ADAPTCALLBACKHANDLE_HH
