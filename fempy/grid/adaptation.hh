#ifndef DUNE_FEMPY_GRID_ADAPTATION_HH
#define DUNE_FEMPY_GRID_ADAPTATION_HH

#include <cstddef>

#include <functional>
#include <memory>
#include <string>
#include <typeindex>
#include <utility>
#include <unordered_map>
#include <vector>

#include <dune/common/dynvector.hh>

#include <dune/grid/common/rangegenerators.hh>

#include <dune/fem/space/common/adaptmanager.hh>
#include <dune/fem/space/common/communicationmanager.hh>
#include <dune/fem/space/common/restrictprolonginterface.hh>

#include <dune/fempy/grid/virtualizedrestrictprolong.hh>
#include <dune/fempy/parameter.hh>

namespace Dune
{

  namespace FemPy
  {

    // LoadBalanceContainsCheck
    // ------------------------

    template< class Grid, class Predicate >
    struct LoadBalanceContainsCheck
    {
      explicit LoadBalanceContainsCheck ( Predicate predicate )
        : predicate_( std::move( predicate ) )
      {}

      bool contains ( const typename Grid::template Codim< 0 >::Entity &element ) const
      {
        return predicate_( element );
      }

    private:
      Predicate predicate_;
    };



    // loadBalanceContainsCheck
    // ------------------------

    template< class Grid, class Predicate >
    inline static LoadBalanceContainsCheck< Grid, Predicate > loadBalanceContainsCheck ( Predicate predicate ) noexcept
    {
      return LoadBalanceContainsCheck< Grid, Predicate >( std::move( predicate ) );
    }



    // FakeGridPart
    // ------------

    template< class Grid >
    class FakeGridPart
    {
      typedef typename Grid::template Codim< 0 >::Entity Element;

    public:
      explicit FakeGridPart ( Grid &grid ) : grid_( grid ) {}

      struct IndexSetType
      {
        bool contains ( const Element & ) const { return true; }
      };

      typedef Grid GridType;

      Grid &grid () const { return grid_; }

      const Element &convert ( const Element &entity ) const { return entity; }

    private:
      Grid &grid_;
    };



    // FakeDiscreteFunctionSpace
    // -------------------------

    template< class Grid, class NumLocalDofs = std::function< std::size_t( typename Grid::template Codim< 0 >::Entity ) > >
    struct FakeDiscreteFunctionSpace
    {
      typedef FakeGridPart< Grid > GridPartType;

      typedef typename GridPartType::IndexSetType IndexSetType;
      typedef typename GridPartType::GridType GridType;

      typedef typename GridType::template Codim< 0 >::Entity EntityType;

      struct BasisFunctionSetType
      {
        BasisFunctionSetType ( NumLocalDofs numLocalDofs, const EntityType &entity ) : numLocalDofs_( std::move( numLocalDofs ) ), entity_( entity ) {}

        std::size_t size () const { return numLocalDofs_( entity_ ); }

      private:
        NumLocalDofs numLocalDofs_;
        const EntityType &entity_;
      };

      explicit FakeDiscreteFunctionSpace ( GridType &grid, NumLocalDofs numLocalDofs ) : gridPart_( grid ), numLocalDofs_( std::move( numLocalDofs ) ) {}
      explicit FakeDiscreteFunctionSpace ( GridPartType gridPart, NumLocalDofs numLocalDofs ) : gridPart_( std::move( gridPart ) ), numLocalDofs_( std::move( numLocalDofs ) ) {}

      const GridPartType &gridPart () const { return gridPart_; }
      IndexSetType indexSet () const { return IndexSetType(); }

      BasisFunctionSetType basisFunctionSet ( const EntityType &entity ) const { return BasisFunctionSetType( numLocalDofs_, entity ); }

    private:
      GridPartType gridPart_;
      NumLocalDofs numLocalDofs_;
    };



    // fakeDiscreteFunctionSpace
    // -------------------------

    template< class Grid, class NumLocalDofs >
    inline static FakeDiscreteFunctionSpace< Grid, NumLocalDofs > fakeDiscreteFunctionSpace ( const Grid &grid, NumLocalDofs numLocalDofs )
    {
      return FakeDiscreteFunctionSpace< Grid, NumLocalDofs >( grid, std::move( numLocalDofs ) );
    }



    // DiscreteFunctionList
    // --------------------

    template< class Grid, class D = double >
    struct DiscreteFunctionList
      : public Fem::IsDiscreteFunction
    {
      typedef D DofType;

      typedef FakeDiscreteFunctionSpace< Grid > DiscreteFunctionSpaceType;

      typedef typename DiscreteFunctionSpaceType::ElementType ElementType;
      typedef typename DiscreteFunctionSpaceType::GridPartType GridPartType;

      typedef std::allocator< DofType > LocalDofVectorAllocatorType;

    private:
      struct DofVector
      {
        virtual ~DofVector () = default;

        virtual void enableDofCompression () = 0;

        virtual std::size_t numLocalDofs ( const ElementType &element ) const = 0;

        virtual DofType *getLocalDofs ( const ElementType &element, DofType *localDofs ) const = 0;
        virtual const DofType *setLocalDofs ( const ElementType &element, const DofType *localDofs ) = 0;
      };

      template< class Impl >
      struct DiscreteFunction final
        : public DofVector
      {
        DiscreteFunction ( Impl &impl ) : impl_( impl ) {}

        virtual void enableDofCompression () override { impl_.enableDofCompression(); }

        virtual std::size_t numLocalDofs ( const ElementType &element ) const override
        {
          return blockMapper().numDofs( element ) * Impl::DiscreteFunctionSpaceType::localBlockSize;
        }

        virtual DofType *getLocalDofs ( const ElementType &element, DofType *localDofs ) const override
        {
          //using Fem::dofBlockFunctor;

          Fem::AssignFunctor< DofType * > assignFunctor( localDofs );
          blockMapper().mapEach( element, dofBlockFunctor( impl_.dofVector(), assignFunctor ) );

          return localDofs + numLocalDofs( element );
        }

        virtual const DofType *setLocalDofs ( const ElementType &element, const DofType *localDofs ) override
        {
          //using Fem::dofBlockFunctor;

          Fem::LeftAssign< const DofType * > assignFunctor( localDofs );
          blockMapper().mapEach( element, dofBlockFunctor( impl_.dofVector(), assignFunctor ) );

          return localDofs + numLocalDofs( element );
        }

      private:
        const typename Impl::DiscreteFunctionSpaceType::BlockMapperType &blockMapper () const { return impl_.space().blockMapper(); }

        Impl &impl_;
      };

    public:
      explicit DiscreteFunctionList ( Grid &grid )
        : space_( grid, [ this ] ( const ElementType &element ) {
              std::size_t size = 0;
              for( const auto &dv : dofVectors_ )
                size += dv->numLocalDofs( element );
              return size;
            } )
      {}

      const GridPartType &gridPart () const { return space_.gridPart(); }

      const DiscreteFunctionSpaceType &space () const { return space_; }

      void enableDofCompression ()
      {
        for( const auto &dv : dofVectors_ )
          dv->enableDofCompression();
      }

      LocalDofVectorAllocatorType localDofVectorAllocator () const { return LocalDofVectorAllocatorType(); }

      template< class A >
      void getLocalDofs ( const ElementType &entity, DynamicVector< DofType, A > &localDofs ) const
      {
        DofType *it = &localDofs[ 0 ];
        for( const auto &dv : dofVectors_ )
          it = dv->getLocalDofs( entity, it );
      }

      template< class A >
      void setLocalDofs ( const ElementType &entity, const DynamicVector< DofType, A > &localDofs )
      {
        const DofType *it = &localDofs[ 0 ];
        for( const auto &dv : dofVectors_ )
          it = dv->setLocalDofs( entity, it );
      }

      void clear () { dofVectors_.clear(); }

      template< class DF >
      void add ( DF &df )
      {
        dofVectors_.push_back( new DiscreteFunction< DF >( df ) );
      }

    private:
      DiscreteFunctionSpaceType space_;
      std::vector< std::unique_ptr< DofVector > > dofVectors_;
    };

  } // namespace FemPy



  namespace Fem
  {

    // DiscreteFunctionTraits for DiscreteFunctionList
    // -----------------------------------------------

    template< class Grid, class D >
    struct DiscreteFunctionTraits< FemPy::DiscreteFunctionList< Grid, D > >
    {
      typedef D DofType;
      typedef std::allocator< DofType > LocalDofVectorAllocatorType;
    };

  } // namespace Fem



  namespace FemPy
  {

    // DiscreteFunctionManager
    // -----------------------

    template< class Grid >
    class DiscreteFunctionManager
    {
      struct DFList
      {
        virtual ~DFList () = default;
        virtual void clear () = 0;
        virtual void registerToLoadBalancer ( Fem::LoadBalancer< Grid > &loadBalancer ) = 0;
      };

      template< class D >
      struct DFListImpl final
        : public DFList
      {
        DFListImpl ( Grid &grid ) : list( grid ) {}

        virtual void clear () override { list.clear(); }

        virtual void registerToLoadBalancer ( Fem::LoadBalancer< Grid > &loadBalancer ) override
        {
          // todo: implement this check
          auto contains = [] ( const typename Grid::template Codim< 0 >::Entity &e ) { return false; };
          loadBalancer.addDiscreteFunction( list, loadBalanceContainsCheck< Grid >( contains ) );
        }

        DiscreteFunctionList< Grid, D > list;
      };

    public:
      DiscreteFunctionManager ( Grid &grid ) : grid_( grid ) {}

      void clear ()
      {
        for( auto &dfList : dfLists_ )
          dfList.second->clear();
      }

      template< class DF >
      void addDiscreteFunction ( DF &df )
      {
        dfList< typename DF::DofType >.add( df );
      }

      template< class DF, class ContainsCheck >
      void addDiscreteFunction ( DF &df, const ContainsCheck &containsCheck )
      {
        addDiscreteFunction( df );
      }

      template< class DF >
      void addToLoadBalancer ( DF &df )
      {
        addDiscreteFunction( df );
      }

      const Grid &grid () const { return grid_; }
      Grid &grid () { return grid_; }

      void registerToLoadBalancer ( Fem::LoadBalancer< Grid > &loadBalancer )
      {
        loadBalancer_ = &loadBalancer;
        for( auto &dfList : dfLists_ )
          dfList.second->registerToLoadBalancer( loadBalancer );
      }

    private:
      template< class D >
      DiscreteFunctionList< Grid, D > &dfList ()
      {
        auto result = dfLists_.insert( std::make_pair( std::type_index( typeid( D ) ), nullptr ) );
        if( result.second )
        {
          result.first->second.reset( new DFListImpl< D >( grid() ) );
          if( loadBalancer_ )
            result.first->second->registerToLoadBalancer( *loadBalancer_ );
        }
        return static_cast< DFListImpl< D > * >( result.first->second )->list;
      }

      Grid &grid_;
      Fem::LoadBalancer< Grid > *loadBalancer_ = nullptr;
      std::unordered_map< std::type_index, std::unique_ptr< DFList > > dfLists_;
    };



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

      void restrictLocal ( const ElementType &father, const ElementType &child, const LocalGeometryType &geometryInFather, bool initialize ) const
      {
        for( const auto &rp : restrictProlongs_ )
          rp.restrictLocal( father, child, geometryInFather, initialize );
      }

      void prolongLocal ( const ElementType &father, const ElementType &child, bool initialize ) const
      {
        for( const auto &rp : restrictProlongs_ )
          rp.prolongLocal( father, child, initialize );
      }

      void prolongLocal ( const ElementType &father, const ElementType &child, const LocalGeometryType &geometryInFather, bool initialize ) const
      {
        for( const auto &rp : restrictProlongs_ )
          rp.prolongLocal( father, child, geometryInFather, initialize );
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

      enum class Marker { Coarsen = -1, Keep = 0, Refine = 1 };

      typedef typename Grid::template Codim< 0 >::Entity Element;

      explicit GridAdaptation ( Grid &grid )
        : restrictProlong_( grid ),
          adaptationManager_( grid, restrictProlong_, noParameter() )
      {}

      template< class Marking >
      void mark ( Marking marking )
      {
        for( const Element &element : elements( grid().leafGridView() ) )
        {
          Marker marker = marking( element );
          grid().mark( static_cast< int >( marker ), element );
        }
      }

      template< class Iterator >
      void adapt ( Iterator begin, Iterator end )
      {
        restrictProlong_.assign( begin, end );
        adaptationManager_.adapt();
        restrictProlong_.clear();
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
