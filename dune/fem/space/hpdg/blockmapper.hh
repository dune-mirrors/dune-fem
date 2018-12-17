#ifndef DUNE_FEM_HPDG_SPACE_DISCONTINUOUSGALERKIN_BLOCKMAPPER_HH
#define DUNE_FEM_HPDG_SPACE_DISCONTINUOUSGALERKIN_BLOCKMAPPER_HH

#include <cassert>
#include <cstddef>

#include <algorithm>
#include <functional>
#include <limits>
#include <type_traits>
#include <utility>
#include <vector>

#include <dune/common/exceptions.hh>

#include <dune/geometry/dimension.hh>

#include <dune/grid/common/gridenums.hh>
#include <dune/grid/utility/persistentcontainer.hh>

#include <dune/fem/gridpart/common/gridpart.hh>
#include <dune/fem/io/streams/streams.hh>
#include <dune/fem/space/common/dofmanager.hh>
#include <dune/fem/space/mapper/dofmapper.hh>

#include <dune/fem/space/common/capabilities.hh>

#include "datahandle.hh"
#include "localdofstorage.hh"

namespace Dune
{

  namespace Fem
  {

    namespace hpDG
    {

      // Internal forward declaration
      // ----------------------------

      template< class GridPart, class LocalKeys >
      struct DiscontinuousGalerkinBlockMapper;



#ifndef DOXYGEN

      // DiscontinuousGalerkinBlockMapperTraits
      // --------------------------------------

      template< class GridPart, class LocalKeys >
      struct DiscontinuousGalerkinBlockMapperTraits
      {
        using DofMapperType = DiscontinuousGalerkinBlockMapper< GridPart, LocalKeys >;
        using ElementType = typename GridPart::template Codim< 0 >::EntityType;
        using SizeType = std::size_t;
      };

#endif // #ifndef DOXYGEN



      // DiscontinuousGalerkinBlockMapper
      // --------------------------------

      /** \brief An \f$hp\f$-adaptive Dune::Fem::DofMapper
       *
       *  \tparam  GridPart  a Dune::Fem::GridPart type
       *  \tparam  LocalKeys  see documentation below
       *
       *  \note The second template parameter is required to provide the
       *        following member methods:
       *  \code
       *  struct LocalKeys
       *  {
       *    // type of key
       *    using KeyType = ImplementationDefined;
       *    // type of data
       *    using DataType = ImplementationDefined;
       *
       *    // return maximum number of dofs associated with an entity
       *    std::size_t maxBlocks () const;
       *
       *    // return maximum number of dofs for given type and key
       *    std::size_t blocks ( GeometryType type, const KeyType &key ) const;
       *
       *    // map key to data type
       *    DataType encode ( const KeyType &key ) const;
       *
       *    // map data to key type
       *    KeyType decode ( const DataType &data ) const;
       *  };
       *  \endcode
       *
       *  \ingroup DiscreteFunctionSpace_Mapper
       */
      template< class GridPart, class LocalKeys >
      struct DiscontinuousGalerkinBlockMapper final
      : public AdaptiveDofMapper< DiscontinuousGalerkinBlockMapperTraits< GridPart, LocalKeys > >
      {
        using ThisType = DiscontinuousGalerkinBlockMapper< GridPart, LocalKeys >;
        using BaseType = AdaptiveDofMapper< DiscontinuousGalerkinBlockMapperTraits< GridPart, LocalKeys > >;

      public:
        /** \brief size type */
        using SizeType = typename BaseType::SizeType;
        /** \brief global key type */
        using GlobalKeyType = typename BaseType::GlobalKeyType;

        /** \brief grid part type */
        using GridPartType = GridPart;
        /** \brief element type */
        using ElementType = typename BaseType::ElementType;

        /** \brief basis function sets type */
        using LocalKeysType = LocalKeys;
        /** \brief key type */
        using KeyType = typename LocalKeysType::KeyType;

      private:
        using DataHandleType = DataHandle< ThisType >;
        friend class DataHandle< DiscontinuousGalerkinBlockMapper< GridPart, LocalKeys > >;

        struct Reserve;
        struct Resize;

        using LocalDofStorageType = LocalDofStorage< GlobalKeyType >;

        using GridType = typename GridPartType::GridType;
        using GridElementType = typename GridType::template Codim< 0 >::Entity;

      public:
        /** \brief dof manager type */
        using DofManagerType = DofManager< GridType >;

        /** \name Construction
         *  \{
         */

        template< class Function >
        DiscontinuousGalerkinBlockMapper ( const GridPartType &gridPart,
                                           const LocalKeysType &localKeys,
                                           const KeyType &value, Function function )
          : gridPart_( gridPart ),
            localKeys_( localKeys ),
            key_( value ),
            keys_( gridPart_.get().grid(), 0, std::make_pair( key_, key_ ) ),
            size_( 0u ),
            dofs_( gridPart_.get().grid(), 0 ),
            dofManager_( DofManagerType::instance( gridPart.grid() ) )
        {
          auto first = gridPart.template begin< 0, All_Partition >();
          auto last = gridPart.template end< 0, All_Partition >();
          for( ; first != last; ++first )
          {
            const ElementType &element = *first;
            const KeyType key = function( element );

            const GridElementType &gridElement = gridEntity( element );
            keys_[ gridElement ] = std::make_pair( key, key );

            auto &dofs = dofs_[ gridElement ];
            resize( dofs, localKeys.blocks( gridElement.type(), key ) );
          }

          dofManager_.get().addIndexSet( *this );
        }

        DiscontinuousGalerkinBlockMapper ( const GridPartType &gridPart,
                                           const LocalKeysType &localKeys,
                                           const KeyType &value )
          : DiscontinuousGalerkinBlockMapper( gridPart, localKeys, value, [&value]( const ElementType & ){ return value; } )
        {}

        /** \} */

        /** \name Copying and assignment
         *  \{
         */

        /** \brief copy constructor */
        DiscontinuousGalerkinBlockMapper ( const ThisType & ) = delete;

        /** \brief move constructor */
        DiscontinuousGalerkinBlockMapper ( ThisType && ) = default;

        /** \brief assignment operator */
        ThisType &operator= ( const ThisType & ) = delete;

        /** \brief move assignment operator */
        ThisType &operator= ( ThisType && ) = default;

        /** \} */

#ifndef DOXYGEN

        ~DiscontinuousGalerkinBlockMapper ()
        {
          dofManager_.get().removeIndexSet( *this );
        }

#endif // #ifndef DOXYGEN

        /** \name Interface methods
         *  \{
         */

        /** \brief return number of dofs */
        SizeType size () const { return size_; }

        /** \brief return upper bound for number of dofs */
        int maxNumDofs () const { return localKeys().maxBlocks(); }

        /** \brief return number of dofs for given element */
        SizeType numDofs ( const ElementType &element ) const
        {
          return numDofs( element, Codim< ElementType::codimension >() );
        }

        /** \brief return number of dofs for given element */
        template< class Entity >
        SizeType numEntityDofs ( const Entity &entity ) const
        {
          return numDofs( entity, Codim< Entity::codimension >() );
        }

        void onSubEntity ( const ElementType &element, int i, int c, std::vector< bool > &indices ) const
        {
          indices.resize( numDofs(element) );
          if (c == 0)
            std::fill(indices.begin(),indices.end(),true);
          else
            std::fill(indices.begin(),indices.end(),false);
        }

        /** \brief return \b true if dofs are associated to codimension */
        static constexpr bool contains ( const int codim ) { return (codim == 0); }

        /** \brief return \b true if number of dofs is fixed for given codimension */
        static constexpr bool fixedDataSize ( int codim ) { return (codim != 0); }

        /** \brief map local dof to global key */
        template< class Function >
        void mapEach ( const ElementType &element, Function function ) const
        {
          mapEach( element, function, Codim< ElementType::codimension >() );
        }

        /** \brief map local dof to global key */
        template< class Entity, class Function >
        void mapEachEntityDof ( const Entity &entity, Function function ) const
        {
          mapEach( entity, function, Codim< Entity::codimension >() );
        }

        /** \brief return number of holes during compression */
        SizeType numberOfHoles ( const int block ) const
        {
          assert( block == 0 );
          return indices_.size();
        }

        /** \brief return old index of given hole during compression */
        GlobalKeyType oldIndex ( const int hole, const int block ) const
        {
          assert( block == 0 );
          GlobalKeyType dof = indices_[ hole ].first;
          assert( dof != std::numeric_limits< GlobalKeyType >::max() );
          return std::move( dof );
        }

        /** \brief return new index of given hole during compression */
        GlobalKeyType newIndex ( const int hole, const int block ) const
        {
          assert( block == 0 );
          return indices_[ hole ].second;
        }

        /** \brief return \b true (this mapper yields a consecutive DOF numbering) */
        static constexpr bool consecutive () { return true; }

        /** \brief return 0 (this mapper has no offset) */
        static constexpr SizeType oldOffSet ( const int block ) { return 0u; }

        /** \brief return 0 (this mapper has no offset) */
        static constexpr SizeType offSet ( const int block ) { return 0u; }

        /** \brief return 1 (this mapper has one block) */
        static constexpr SizeType numBlocks () { return 1u; }

        /** \} */

        /** \name Adaptation interface
         *  \{
         */

        /** \brief get key currently assigned to an entity */
        template< class Element >
        typename std::enable_if< (std::is_same< Element, ElementType >::value || std::is_same< Element, GridElementType >::value), const KeyType & >::type
        key ( const Element &element ) const
        {
          return keys_[ gridEntity( element ) ].first;
        }

        /** \brief set key to be assigned to an entity after next call to adapt() */
        void mark ( const KeyType &key, const ElementType &element )
        {
          assert( gridPart().indexSet().contains( element ) );
          keys_[ gridEntity( element ) ].second = key;
        }

        /** \brief get key to be assigned to an entity after next call to adapt() */
        KeyType getMark ( const ElementType &element ) const
        {
          assert( gridPart().indexSet().contains( element ) );
          return keys_[ gridEntity( element ) ].second;
        }

        /** \brief please doc me */
        template< class Function >
        bool adapt ( Function function );

        /** \brief please doc me */
        bool adapt ()
        {
          auto function = []( const ElementType &,
                              const KeyType &, const KeyType &,
                              const std::vector< std::size_t > &,
                              const std::vector< std::size_t > & )
          {};
          return adapt( function );
        }

        /** \} */

        /** \name Public interface for adaptation managers
         *  \{
         */

        /** \brief return DOF manager */
        const DofManagerType &dofManager () const { return dofManager_.get(); }

        /** \brief add DOFs for element */
        void insertEntity ( const GridElementType &gridElement )
        {
          resize( dofs_ );
          resize( keys_, std::make_pair( key_, key_ ) );

          resize( dofs_[ gridElement ], localKeys().blocks( gridElement.type(), key( gridElement ) ) );
        }

        /** \brief mark DOFs for removal */
        void removeEntity ( const GridElementType &gridElement )
        {
          auto &dofs = dofs_[ gridElement ];

          size_ = dofs.reserve( 0u, Reserve( size_ ) );
          Resize function( holes_ );
          for( auto dof : dofs )
            function( dof );
        }

        /** \brief add DOFs for new element */
        void insertNewEntity ( const GridElementType &gridElement )
        {
          if( gridElement.isNew() )
          {
            resize( dofs_[ gridElement ], localKeys().blocks( gridElement.type(), key( gridElement ) ) );
            assert( gridElement.hasFather() );
            const GridElementType &father = gridElement.father();
            insertNewEntity( father );
            removeEntity( father );
          }
          else
          {
            auto &dofs = dofs_[ gridElement ];
            resize( dofs, localKeys().blocks( gridElement.type(), key( gridElement ) ) );
            for( auto dof : dofs )
            {
              auto iterator = std::lower_bound( holes_.begin(), holes_.end(), dof );
              if( iterator != holes_.end() && *iterator == dof )
                holes_.erase( iterator );
            }
          }
        }

        void resize ()
        {
          resize( dofs_ );
          resize( keys_, std::make_pair( key_, key_ ) );

          auto first = gridPart().template begin< 0, All_Partition >();
          auto last = gridPart().template end< 0, All_Partition >();
          for( ; first != last; ++first )
          {
            const ElementType &element = *first;
            insertNewEntity( gridEntity( element ) );
          }
#if 0
            auto &dofs = dofs_[ gridEntity( element ) ];
            resize( dofs, localKeys().blocks( element.type(), key( element ) ) );
#endif
        }

        /** \brief compress DOF mapping */
        bool compress ();

        /** \brief this mapper has no I/O capabilities */
        template< class Traits >
        void write ( OutStreamInterface< Traits > & )
        {
          DUNE_THROW( NotImplemented, "Method write() not implemented yet" );
        }

        /** \brief this mapper has no I/O capabilities */
        template< class Traits >
        void read ( InStreamInterface< Traits > & )
        {
          DUNE_THROW( NotImplemented, "Method read() not implemented yet" );
        }

        /** \brief this mapper has no I/O capabilities */
        void backup () const
        {
          DUNE_THROW( NotImplemented, "Method backup() not implemented" );
        }

        /** \brief this mapper has no I/O capabilities */
        void restore ()
        {
          DUNE_THROW( NotImplemented, "Method restore() not implemented" );
        }

      private:
        bool compressed () const { return holes_.empty(); }

        template< class Element, class Function >
        void mapEach ( const Element &element, Function function, Codim< 0 > ) const
        {
          const auto &dofs = dofs_[ gridEntity( element ) ];
          std::size_t size = dofs.size();
          for( std::size_t i = 0; i < size; ++i )
            function( i, dofs[ i ] );
        }

        template< class Entity, class Function, int codim >
        void mapEach ( const Entity &entity, Function function, Codim< codim > ) const
        {}

        template< class Element >
        SizeType numDofs ( const Element &element, Codim< 0 > ) const
        {
          return dofs_[ gridEntity( element ) ].size();
        }

        template< class Entity, int codim >
        SizeType numDofs ( const Entity &entity, Codim< codim > ) const
        {
          return 0u;
        }

        void resize ( LocalDofStorageType &dofs, SizeType n )
        {
          size_ = dofs.reserve( n, Reserve( size_ ) );
          dofs.resize( Resize( holes_ ) );
        }

        template< class T >
        static void resize ( PersistentContainer< GridType, T > &container, const T &value = T() )
        {
          container.resize( value );
          container.shrinkToFit();
        }

        const GridPartType &gridPart () const { return gridPart_.get(); }

        const LocalKeysType &localKeys () const { return localKeys_.get(); }

        std::reference_wrapper< const GridPartType > gridPart_;
        std::reference_wrapper< const LocalKeysType > localKeys_;

        const KeyType key_;
        PersistentContainer< GridType, std::pair< KeyType, KeyType > > keys_;

        SizeType size_;
        PersistentContainer< GridType, LocalDofStorageType > dofs_;

        std::vector< GlobalKeyType > holes_;
        std::vector< std::pair< GlobalKeyType, GlobalKeyType > > indices_;

        std::reference_wrapper< DofManagerType > dofManager_;
      };



      // DiscontinuousGalerkinBlockMapper::Reserve
      // -----------------------------------------

      template< class GridPart, class LocalKeys >
      struct DiscontinuousGalerkinBlockMapper< GridPart, LocalKeys >::Reserve
      {
        explicit Reserve ( SizeType &size ) : size_( size ) {}

        GlobalKeyType operator() () { return size_++; }

        operator SizeType () const { return size_; }

      private:
        SizeType &size_;
      };



      // DiscontinuousGalerkinBlockMapper::Resize
      // ----------------------------------------

      template< class GridPart, class LocalKeys >
      struct DiscontinuousGalerkinBlockMapper< GridPart, LocalKeys >::Resize
      {
        explicit Resize ( std::vector< GlobalKeyType > &holes )
          : holes_( holes )
        {
          assert( std::is_sorted( holes_.begin(), holes_.end() ) );
        }

        ~Resize ()
        {
          assert( std::is_sorted( holes_.begin(), holes_.end() ) );
        }

        void operator() ( const GlobalKeyType &dof )
        {
          auto iterator = std::lower_bound( holes_.begin(), holes_.end(), dof );
          if( iterator == holes_.end() || *iterator != dof )
            holes_.insert( iterator, dof );
        }

      private:
        std::vector< GlobalKeyType > &holes_;
      };



      // DiscontinuousGalerkinBlockMapper::adapt
      // ---------------------------------------

      template< class GridPart, class LocalKeys >
      template< class Function >
      bool DiscontinuousGalerkinBlockMapper< GridPart, LocalKeys >::adapt ( Function function )
      {
        DataHandleType dataHandle( *this );
        gridPart().communicate( dataHandle, InteriorBorder_All_Interface, ForwardCommunication );

        std::vector< std::size_t > origin, destination;
        origin.reserve( maxNumDofs() );
        destination.reserve( maxNumDofs() );

        std::size_t count = 0;

        auto first = gridPart().template begin< 0, All_Partition >();
        auto last = gridPart().template end< 0, All_Partition >();
        for( ; first != last; ++first )
        {
          const ElementType &element = *first;
          const GridElementType &gridElement = gridEntity( element );
          auto &keys = keys_[ gridElement ];
          if( keys.first == keys.second )
            continue;

          // remember old state
          const KeyType prior = keys.first;
          auto &dofs = dofs_[ gridElement ];
          origin.resize( dofs.size() );
          std::copy( dofs.begin(), dofs.end(), origin.begin() );

          // reset to new state
          keys.first = keys.second;
          resize( dofs, localKeys().blocks( gridElement.type(), keys.first ) );

          // call function
          destination.resize( dofs.size() );
          std::copy( dofs.begin(), dofs.end(), destination.begin() );
          function( element, prior, keys.first, origin, destination );

          ++count;
        }
        return (count > 0u);
      }



      // DiscontinuousGalerkinBlockMapper::compress
      // ------------------------------------------

      template< class GridPart, class LocalKeys >
      bool DiscontinuousGalerkinBlockMapper< GridPart, LocalKeys >::compress ()
      {
        // resize persistent containers
        resize( dofs_ );
        resize( keys_, std::make_pair( key_, key_ ) );

        for( auto &dofs : dofs_ )
          dofs.resize( Resize( holes_ ) );

        // check if compress is necessary
        if( holes_.empty() )
        {
          indices_.clear();
          return false;
        }

        // get number of dofs after compress
        assert( size_ >= holes_.size() );
        size_ -= holes_.size();

        // remove trailing indices
        auto iterator = std::lower_bound( holes_.begin(), holes_.end(), size_ );
        holes_.erase( iterator, holes_.end() );

        // initialize vector of old and new indices
        auto assign = []( GlobalKeyType newIndex ) -> std::pair< GlobalKeyType, GlobalKeyType >
        {
          GlobalKeyType oldIndex = std::numeric_limits< GlobalKeyType >::max();
          return std::make_pair( oldIndex, newIndex );
        };
        indices_.resize( holes_.size() );
        std::transform( holes_.begin(), holes_.end(), indices_.begin(), assign );

        // update local dof storages
        std::size_t hole = 0u;
        for( auto &dofs : dofs_ )
        {
          for( auto &dof : dofs )
          {
            if( dof >= size_ )
            {
              auto &indices = indices_.at( hole++ );
              indices.first = dof;
              dof = indices.second;
            }
          }
        }

        assert( hole == holes_.size() );
        holes_.clear();

        // replace all but markers currently used by default value
        PersistentContainer< GridType, std::pair< KeyType, KeyType > > tmp( gridPart().grid(), 0, std::make_pair( key_, key_ ) );
        auto first = gridPart().template begin< 0, All_Partition >();
        auto last = gridPart().template end< 0, All_Partition >();
        for( ; first != last; ++first )
        {
          const ElementType &element = *first;
          const GridElementType &gridElement = gridEntity( element );
          tmp[ gridElement ] = keys_[ gridElement ];
        }
        keys_.swap( tmp );

        return true;
      }

    } // namespace hpDG

    namespace Capabilities
    {
      // isConsecutiveIndexSet
      // ---------------------

      template< class GridPart, class LocalKeys >
      struct isConsecutiveIndexSet< hpDG::DiscontinuousGalerkinBlockMapper< GridPart, LocalKeys > >
      {
        static const bool v = true;
      };

      // isAdaptiveDofMapper
      // -------------------

      template< class GridPart, class LocalKeys >
      struct isAdaptiveDofMapper< hpDG::DiscontinuousGalerkinBlockMapper< GridPart, LocalKeys > >
      {
        static const bool v = true;
      };

    } // namespace Capabilities

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_HPDG_SPACE_DISCONTINUOUSGALERKIN_BLOCKMAPPER_HH
