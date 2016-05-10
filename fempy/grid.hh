#ifndef PYTHON_DUNEGRID_HH
#define PYTHON_DUNEGRID_HH

#include <cassert>

#include <memory>
#include <string>
#include <utility>
#include <vector>

#include <dune/common/dynvector.hh>

#include <dune/grid/common/rangegenerators.hh>

#include <dune/fem/gridpart/adaptiveleafgridpart.hh>
#include <dune/fem/gridpart/leafgridpart.hh>
#include <dune/fem/io/file/vtkio.hh>
#include <dune/fem/space/common/adaptmanager.hh>
#include <dune/fem/space/common/communicationmanager.hh>
#include <dune/fem/space/common/restrictprolonginterface.hh>
#include <dune/fem/storage/envelope.hh>

#include <dune/fempy/gridfunction.hh>

namespace Dune
{

  namespace FemPy
  {

    /*!
        @file
        @brief Contains C++ template classes for grids.
        \ingroup Grids
    */

    // CornerIterator
    // --------------

    template< class Geometry >
    class CornerIterator
      : public std::iterator< std::random_access_iterator_tag, typename Geometry::GlobalCoordinate, int, Fem::Envelope< typename Geometry::GlobalCoordinate >, typename Geometry::GlobalCoordinate >
    {
      typedef std::iterator< std::random_access_iterator_tag, typename Geometry::GlobalCoordinate, int, Fem::Envelope< typename Geometry::GlobalCoordinate >, typename Geometry::GlobalCoordinate > Base;

    public:
      typedef typename Base::reference reference;
      typedef typename Base::pointer pointer;

      CornerIterator () = default;

      CornerIterator ( const Geometry &geometry, int index )
        : geometry_( &geometry ), index_( index )
      {}

      CornerIterator ( const Geometry &geometry )
        : geometry_( &geometry ), index_( geometry.corners() )
      {}

      reference operator* () const { return geometry().corner( index_ ); }
      pointer operator-> () const { return pointer( geometry().corner( index_ ) ); }

      reference operator[] ( int n ) { return geometry().corner( index_ + n ); }

      CornerIterator &operator++ () { ++index_; return *this; }
      CornerIterator &operator-- () { --index_; return *this; }

      CornerIterator &operator+= ( int n ) { index_ += n; return *this; }
      CornerIterator &operator-= ( int n ) { index_ -= n; return *this; }

      CornerIterator operator++ ( int ) { CornerIterator copy( *this ); ++(*this); return copy; }
      CornerIterator operator-- ( int ) { CornerIterator copy( *this ); --(*this); return copy; }

      CornerIterator operator+ ( int n ) const { CornerIterator copy( *this ); copy += n; return copy; }
      CornerIterator operator- ( int n ) const { CornerIterator copy( *this ); copy -= n; return copy; }

      int operator- ( const CornerIterator &other ) const { return (index_ - other.index_); }

      friend CornerIterator operator+ ( int n, const CornerIterator &it ) { CornerIterator copy( it ); copy += n; return copy; }

      bool operator== ( const CornerIterator &other ) const { return (index_ == other.index_); }
      bool operator!= ( const CornerIterator &other ) const { return (index_ != other.index_); }

      bool operator< ( const CornerIterator &other ) const { return (index_ < other.index_); }
      bool operator<= ( const CornerIterator &other ) const { return (index_ <= other.index_); }
      bool operator> ( const CornerIterator &other ) const { return (index_ > other.index_); }
      bool operator>= ( const CornerIterator &other ) const { return (index_ >= other.index_); }

      const Geometry &geometry () const { assert( geometry_ ); return *geometry_; }

    private:
      const Geometry *geometry_ = nullptr;
      int index_ = 0;
    };



    // noParameter
    // -----------

    inline static Fem::ParameterReader noParameter ()
    {
      return Fem::ParameterReader( [] ( const std::string &, const std::string *def ) { return def; } );
    }



    // readDGF
    // -------

    template< class Grid >
    static std::shared_ptr< Grid > readDGF ( const std::string &dgf )
    {
      GridPtr< Grid > gridPtr( dgf );
      gridPtr->loadBalance();
      return std::shared_ptr< Grid >( gridPtr.release() );
    }



    // DiscreteFunctionList
    // --------------------

    template< class Grid >
    struct DiscreteFunctionList
      : public Fem::IsDiscreteFunction
    {
      typedef AdaptiveGridFunction< Grid > DiscreteFunction;

      typedef typename DiscreteFunction::DofType DofType;

      typedef typename Grid::template Codim< 0 >::Entity ElementType;

      struct GridPartType
      {
        explicit GridPartType ( std::shared_ptr< Grid > grid ) : grid_( std::move( grid ) ) {}

        struct IndexSetType
        {
          bool contains ( const ElementType & ) const { return true; }
        };

        typedef Grid GridType;

        const Grid &grid () const { return *grid_; }

        const ElementType &convert ( const ElementType &entity ) const { return entity; }

      private:
        std::shared_ptr< Grid > grid_;
      };

      struct DiscreteFunctionSpaceType
      {
        typedef std::vector< std::shared_ptr< DiscreteFunction > > DiscreteFunctions;

        typedef DiscreteFunctionList::GridPartType GridPartType;

        typedef typename GridPartType::IndexSetType IndexSetType;
        typedef typename GridPartType::GridType GridType;
        typedef ElementType EntityType;

        struct BasisFunctionSetType
        {
          BasisFunctionSetType ( const DiscreteFunctions &discreteFunctions, const ElementType &entity )
            : discreteFunctions_( discreteFunctions ), entity_( entity )
          {}

          std::size_t size () const
          {
            std::size_t size = 0;
            for( const auto &df : discreteFunctions_ )
              size += df->numLocalDofs( entity_ );
            return size;
          }

        private:
          const DiscreteFunctions &discreteFunctions_;
          const ElementType &entity_;
        };

        explicit DiscreteFunctionSpaceType ( GridPartType gridPart ) : gridPart_( std::move( gridPart ) ) {}

        const GridPartType &gridPart () const { return gridPart_; }
        IndexSetType indexSet () const { return IndexSetType(); }

        BasisFunctionSetType basisFunctionSet ( const ElementType &entity ) const { return BasisFunctionSetType( discreteFunctions_, entity ); }

        const DiscreteFunctions &discreteFunctions () const { return discreteFunctions_; }
        DiscreteFunctions &discreteFunctions () { return discreteFunctions_; }

      private:
        GridPartType gridPart_;
        DiscreteFunctions discreteFunctions_;
      };

      typedef typename DiscreteFunctionSpaceType::DiscreteFunctions::const_iterator ConstIterator;
      typedef typename DiscreteFunctionSpaceType::DiscreteFunctions::iterator Iterator;

      typedef std::allocator< DofType > LocalDofVectorAllocatorType;

      explicit DiscreteFunctionList ( std::shared_ptr< Grid > grid )
        : space_( GridPartType( std::move( grid ) ) )
      {}

      const GridPartType &gridPart () const { return space_.gridPart(); }

      const DiscreteFunctionSpaceType &space () const { return space_; }

      void enableDofCompression ()
      {
        for( const auto &df : space_.discreteFunctions() )
          df->enableDofCompression();
      }

      void assign () { space_.discreteFunctions().clear(); }

      template< class Iterator >
      void assign ( Iterator begin, Iterator end )
      {
        space_.discreteFunctions().assign( begin, end );
      }

      ConstIterator begin () const { return space_.discreteFunctions().begin(); }
      Iterator begin () { return space_.discreteFunctions().begin(); }
      ConstIterator end () const { return space_.discreteFunctions().end(); }
      Iterator end () { return space_.discreteFunctions().end(); }

      LocalDofVectorAllocatorType localDofVectorAllocator () const { return LocalDofVectorAllocatorType(); }

      template< class A >
      void getLocalDofs ( const ElementType &entity, DynamicVector< DofType, A > &localDofs ) const
      {
        DofType *it = &localDofs[ 0 ];
        for( const auto &df : space_.discreteFunctions() )
          it = df->getLocalDofs( entity, it );
      }

      template< class A >
      void setLocalDofs ( const ElementType &entity, const DynamicVector< DofType, A > &localDofs )
      {
        const DofType *it = &localDofs[ 0 ];
        for( auto &df : space_.discreteFunctions() )
          it = df->setLocalDofs( entity, it );
      }

    private:
      DiscreteFunctionSpaceType space_;
    };

  } // namespace FemPy



  namespace Fem
  {

    // DiscreteFunctionTraits for DiscreteFunctionList
    // -----------------------------------------------

    template< class Grid >
    struct DiscreteFunctionTraits< FemPy::DiscreteFunctionList< Grid > >
    {
      typedef typename FemPy::AdaptiveGridFunction< Grid >::DofType DofType;
      typedef std::allocator< DofType > LocalDofVectorAllocatorType;
    };

  } // namespace Fem



  namespace FemPy
  {

    // RestrictProlong
    // ---------------

    template< class Grid >
    class RestrictProlong
      : public Fem::RestrictProlongInterface< Fem::RestrictProlongTraits< RestrictProlong< Grid >, typename Grid::ctype > >
    {
      typedef Fem::RestrictProlongInterface< Fem::RestrictProlongTraits< RestrictProlong< Grid >, typename Grid::ctype > > BaseType;

      struct LoadBalanceContainsCheck;

    public:
      typedef typename BaseType::DomainFieldType DomainFieldType;

      typedef typename Grid::template Codim< 0 >::Entity ElementType;
      typedef typename Grid::template Codim< 0 >::LocalGeometry LocalGeometryType;

      explicit RestrictProlong ( std::shared_ptr< Grid > grid )
        : discreteFunctions_( std::move( grid ) )
      {}

      void setFatherChildWeight ( const DomainFieldType &weight ) const
      {
        for( const auto &df : discreteFunctions_ )
          df->restrictProlong().setFatherChildWeight( weight );
      }

      void restrictLocal ( const ElementType &father, const ElementType &child, bool initialize ) const
      {
        for( const auto &df : discreteFunctions_ )
          df->restrictProlong().restrictLocal( father, child, initialize );
      }

      void restrictLocal ( const ElementType &father, const ElementType &child, const LocalGeometryType &geometryInFather, bool initialize ) const
      {
        for( const auto &df : discreteFunctions_ )
          df->restrictProlong().restrictLocal( father, child, geometryInFather, initialize );
      }

      void prolongLocal ( const ElementType &father, const ElementType &child, bool initialize ) const
      {
        for( const auto &df : discreteFunctions_ )
          df->restrictProlong().prolongLocal( father, child, initialize );
      }

      void prolongLocal ( const ElementType &father, const ElementType &child, const LocalGeometryType &geometryInFather, bool initialize ) const
      {
        for( const auto &df : discreteFunctions_ )
          df->restrictProlong().prolongLocal( father, child, geometryInFather, initialize );
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

      template< class LoadBalancer >
      void addToLoadBalancer ( LoadBalancer &loadBalancer );

      void assign ()
      {
        removeFromCommList();
        discreteFunctions_.assign();
        // no need to add to comm list; discrete functions is empty
      }

      template< class Iterator >
      void assign ( Iterator begin, Iterator end )
      {
        removeFromCommList();
        discreteFunctions_.assign( std::move( begin ), std::move( end ) );
        addToCommList();
      }

    private:
      void addToCommList ()
      {
        if( commList_ )
          for( const auto &df : discreteFunctions_ )
            df->restrictProlong().addToList( *commList_ );
      }

      void removeFromCommList ()
      {
        if( commList_ )
          for( const auto &df : discreteFunctions_ )
            df->restrictProlong().removeFromList( *commList_ );
      }

      DiscreteFunctionList< Grid > discreteFunctions_;
      Fem::CommunicationManagerList *commList_ = nullptr;
    };



    // RestrictProlong::LoadBalanceContainsCheck
    // -----------------------------------------

    template< class Grid >
    struct RestrictProlong< Grid >::LoadBalanceContainsCheck
    {
      explicit LoadBalanceContainsCheck ( const DiscreteFunctionList< Grid > &discreteFunctions )
        : discreteFunctions_( discreteFunctions )
      {}

      bool contains ( const ElementType &element ) const
      {
        // todo: implement this predicate
        return false;
      }

    private:
      const DiscreteFunctionList< Grid > &discreteFunctions_;
    };



    // Implementation of RestrictProlong
    // ---------------------------------

    template< class Grid >
    template< class LoadBalancer >
    inline void RestrictProlong< Grid >::addToLoadBalancer ( LoadBalancer &loadBalancer )
    {
      loadBalancer.addDiscreteFunction( discreteFunctions_, LoadBalanceContainsCheck( discreteFunctions_ ) );
    }



    // HierarchicalGrid
    // ----------------

    template< class G >
    struct HierarchicalGrid
    {
      typedef G Grid;

      typedef Fem::AdaptationManager< Grid, RestrictProlong< Grid > > AdaptationManager;

      enum class Marker { Coarsen = -1, Keep = 0, Refine = 1 };

      typedef typename Grid::template Codim< 0 >::Entity Element;

      explicit HierarchicalGrid ( const std::string &dgf )
        : grid_( readDGF< Grid >( dgf ) ),
          restrictProlong_( new RestrictProlong< Grid >( grid_ ) ),
          adaptationManager_( new AdaptationManager( *grid_, *restrictProlong_, noParameter() ) )
      {}

      template< class Marking >
      void mark ( Marking marking )
      {
        for( const Element &element : elements( grid_->leafGridView() ) )
        {
          Marker marker = marking( element );
          grid_->mark( static_cast< int >( marker ), element );
        }
      }

      template< class Iterator >
      void adapt ( Iterator begin, Iterator end )
      {
        restrictProlong_->assign( begin, end );
        adaptationManager_->adapt();
        restrictProlong_->assign();
      }

      void globalRefine ( int level )
      {
        Fem::GlobalRefine::apply( *grid_, level );
      }

      template< class Iterator >
      void loadBalance ( Iterator begin, Iterator end )
      {
        restrictProlong_->assign( begin, end );
        adaptationManager_->loadBalance();
        restrictProlong_->assign();
      }

      std::shared_ptr< Grid > grid () const { return grid_; }

    private:
      std::shared_ptr< Grid > grid_;
      std::shared_ptr< RestrictProlong< Grid > > restrictProlong_;
      std::shared_ptr< AdaptationManager > adaptationManager_;
    };



    // LeafGrid
    // --------

    template< class GP >
    struct LeafGrid
    {
      typedef GP GridPart;
      typedef typename GridPart::GridType Grid;

      explicit LeafGrid ( const std::string &dgf )
        : grid_( dgf ),
          gridPart_( new GridPart( *grid_.grid() ) )
      {}

      template< class Mark >
      void adapt ( Mark mark )
      {
        return grid_.adapt( mark );
      }

      void globalRefine ( int level ) { grid_.globalRefine( level ); }

      std::size_t size ( int codim ) const { return gridPart_->indexSet().size( codim ); }

      const std::shared_ptr< Grid > &grid () const { return grid_.grid(); }
      std::shared_ptr< Grid > &grid () { return grid_.grid(); }

      const std::shared_ptr< GridPart > &gridPart () const { return gridPart_; }
      std::shared_ptr< GridPart > &gridPart () { return gridPart_; }

      HierarchicalGrid< Grid > hierarchicalGrid () const { return grid_; }

    private:
      HierarchicalGrid< Grid > grid_;
      std::shared_ptr< GridPart > gridPart_;
    };



    // VTKOutput
    // ---------

    template< class GridPart >
    struct VTKOutput
    {
      explicit VTKOutput ( const LeafGrid< GridPart > &grid )
        : gridPart_( grid.gridPart() ), vtk_( *gridPart_ )
      {}

      template< int R >
      void add ( std::shared_ptr< GridFunction< GridPart, R > > gf )
      {
        gfVector_.push_back( gf );
        vtk_.addCellData( *gf );
      }

      void write ( const char *name ) { vtk_.write( name ); }

    private:
      std::shared_ptr< GridPart > gridPart_;
      Fem::VTKIO< GridPart > vtk_;
      std::vector< std::shared_ptr< GridFunctionBase<GridPart> > > gfVector_;
    };

  } // namespace FemPy

} // namespace Dune

#endif // #ifndef PYTHON_DUNEGRID_HH
