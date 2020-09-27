#ifndef DUNE_FEMPY_GRID_ADAPTATION_HH
#define DUNE_FEMPY_GRID_ADAPTATION_HH

#include <utility>
#include <vector>

#include <dune/grid/common/rangegenerators.hh>

#include <dune/fem/space/common/adaptationmanager.hh>
#include <dune/fem/space/common/communicationmanager.hh>
#include <dune/fem/space/common/restrictprolonginterface.hh>

#include <dune/python/grid/enums.hh>

#include <dune/fempy/grid/discretefunctionmanager.hh>
#include <dune/fempy/grid/virtualizedrestrictprolong.hh>
#include <dune/fempy/parameter.hh>

namespace Dune
{

  namespace FemPy
  {

    // RestrictProlong
    // ---------------

    template< class Grid >
    class RestrictProlong
      : public Fem::RestrictProlongInterface< Fem::RestrictProlongTraits< RestrictProlong< Grid >, typename Grid::ctype > >
    {
      typedef Fem::RestrictProlongInterface< Fem::RestrictProlongTraits< RestrictProlong< Grid >, typename Grid::ctype > > BaseType;

    public:
      typedef typename BaseType::DomainFieldType DomainFieldType;

      typedef typename Grid::template Codim< 0 >::Entity ElementType;
      typedef typename Grid::template Codim< 0 >::LocalGeometry LocalGeometryType;

      explicit RestrictProlong ( Grid &grid )
        : discreteFunctionManager_( grid )
      {}

      void initialize ()
      {
        for( auto &rp : restrictProlongs_ )
          rp.initialize();
      }

      void finalize ()
      {
        for( auto &rp : restrictProlongs_ )
          rp.finalize();
      }

      void setFatherChildWeight ( const DomainFieldType &weight ) const
      {
        for( const auto &rp : restrictProlongs_ )
          rp.setFatherChildWeight( weight );
      }

      void restrictLocal ( const ElementType &father, const ElementType &child, bool initialize ) const
      {
        for( const auto &rp : restrictProlongs_ )
          rp.restrictLocal( father, child, initialize );
      }
      void restrictFinalize ( const ElementType &father ) const
      {
        for( const auto &rp : restrictProlongs_ )
          rp.restrictFinalize( father );
      }

      void prolongLocal ( const ElementType &father, const ElementType &child, bool initialize ) const
      {
        for( const auto &rp : restrictProlongs_ )
          rp.prolongLocal( father, child, initialize );
      }

      void addToList ( Fem::CommunicationManagerList &commList )
      {
        if( commList_ )
          DUNE_THROW( InvalidStateException, "Only one communication list supported." );
        commList_ = &commList;
        addToCommList();
      }

      void removeFromList ( Fem::CommunicationManagerList &commList )
      {
        if( commList_ != &commList )
          DUNE_THROW( InvalidStateException, "Only one communication list supported." );
        removeFromCommList();
        commList_ = nullptr;
      }

      void addToLoadBalancer ( Fem::LoadBalancer< Grid > &loadBalancer ) { discreteFunctionManager_.registerToLoadBalancer( loadBalancer ); }

      void clear ()
      {
        discreteFunctionManager_.clear();
        removeFromCommList();
        restrictProlongs_.clear();
      }

      template< class Iterator >
      void assign ( Iterator begin, Iterator end )
      {
        discreteFunctionManager_.clear();
        removeFromCommList();
        restrictProlongs_.assign( std::move( begin ), std::move( end ) );
        addToCommList();
        for( auto &rp : restrictProlongs_ )
          rp.addToLoadBalancer( discreteFunctionManager_ );
      }

      const Grid &grid () const { return discreteFunctionManager_.grid(); }
      Grid &grid () { return discreteFunctionManager_.grid(); }

    private:
      void addToCommList ()
      {
        if( commList_ )
          for( auto &rp : restrictProlongs_ )
            rp.addToList( *commList_ );
      }

      void removeFromCommList ()
      {
        if( commList_ )
          for( auto &rp : restrictProlongs_ )
            rp.removeFromList( *commList_ );
      }

      std::vector< VirtualizedRestrictProlong< Grid > > restrictProlongs_;
      DiscreteFunctionManager< Grid > discreteFunctionManager_;
      Fem::CommunicationManagerList *commList_ = nullptr;
    };



    // GridAdaptation
    // --------------

    template< class G >
    struct GridAdaptation
    {
      typedef G Grid;

      typedef Fem::AdaptationManager< Grid, RestrictProlong< Grid > > AdaptationManager;

      typedef typename Grid::template Codim< 0 >::Entity Element;

      explicit GridAdaptation ( Grid &grid )
        : restrictProlong_( grid ),
          adaptationManager_( grid, restrictProlong_ )
      {}

      template< class Iterator >
      void adapt ( Iterator begin, Iterator end )
      {
        restrictProlong_.assign( begin, end );
        adaptationManager_.adapt();
        restrictProlong_.clear();
      }

      void markAll ( )
      {
        for ( const auto &e : elements(grid().leafGridView()) )
          grid().mark(1,e);
      }

      void globalRefine ( int level )
      {
        Fem::GlobalRefine::apply( grid(), level );
      }

      template< class Iterator >
      void loadBalance ( Iterator begin, Iterator end )
      {
        restrictProlong_.assign( begin, end );
        adaptationManager_.loadBalance();
        restrictProlong_.clear();
      }

      const Grid &grid () const { return restrictProlong_.grid(); }
      Grid &grid () { return restrictProlong_.grid(); }

    private:
      RestrictProlong< Grid > restrictProlong_;
      AdaptationManager adaptationManager_;
    };

  } // namespace FemPy

} // namespace Dune

#endif // #ifndef DUNE_FEMPY_GRID_ADAPTATION_HH
