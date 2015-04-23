#ifndef DUNE_FEM_SPACE_COMBINEDSPACE_TUPLEMAPPER_HH
#define DUNE_FEM_SPACE_COMBINEDSPACE_TUPLEMAPPER_HH

#include <array>
#include <tuple>
#include <utility>

#include <dune/common/forloop.hh>
#include <dune/common/tuples.hh>

#include <dune/common/std/utility.hh>
#include <dune/fem/common/utility.hh>

#include <dune/fem/space/common/dofmanager.hh>
#include <dune/fem/space/mapper/dofmapper.hh>
#include <dune/fem/space/mapper/nonblockmapper.hh>


namespace Dune
{

  namespace Fem
  {

    //! forward declaration
    template< class ... Mapper >
    class TupleMapper;


    // TupleMapperTraits
    // -----------------

    template< class ... Mapper >
    struct TupleMapperTraits
    {
      typedef Dune::tuple< Mapper ... > MapperTupleType;

      // we still need entities of codimension 0 here
      static_assert( Std::are_all_same< typename Mapper::ElementType ... >::value,
                     "TupleMapper needs common ElementType" );
      typedef typename Dune::tuple_element< 0, MapperTupleType >::type::ElementType ElementType;

      // type of Size's
      typedef typename Dune::tuple_element< 0, MapperTupleType >::type::SizeType SizeType;

      // type of dofmanager
      typedef TupleMapper< Mapper ... > DofMapperType;
    };


    // TupleMapper
    // -----------

    template< class ... Mapper >
    class TupleMapper
      : public AdaptiveDofMapper< TupleMapperTraits< Mapper ... > >
    {
      typedef TupleMapper< Mapper ... > ThisType;
      typedef AdaptiveDofMapper< TupleMapperTraits< Mapper ... > >  BaseType;

      // helper classes
      template< int >
      struct ComputeOffSet;
      template< int >
      struct MapEach;
      template< int >
      struct MapEachEntityDof;
      template< int >
      struct Offset;
      template< int >
      struct OldOffset;
      template< int >
      struct NewIndex;
      template< int >
      struct OldIndex;
      template< int >
      struct NumberOfHoles;


      // FunctorWrapper
      // --------------

      template< class Functor >
      struct FunctorWrapper
      {
        FunctorWrapper ( Functor functor, int localOffset, int globalOffset )
          : functor_( functor ),
            localOffset_( localOffset ),
            globalOffset_( globalOffset )
        {}

        template< class GlobalKey >
        void operator() ( int localDof, const GlobalKey &globalKey )
        {
          functor_( localDof + localOffset_, globalKey + globalOffset_ );
        }

        template< class GlobalKey >
        void operator() ( const GlobalKey &globalKey )
        {
          functor_( globalKey + globalOffset_ );
        }

      private:
        Functor functor_;
        const int localOffset_;
        const int globalOffset_;
      };

      // size of the Mapper Tuple
      static const int mapperTupleSize = sizeof ... ( Mapper );

      typedef std::array< typename BaseType::Traits::SizeType, mapperTupleSize + 1 > OffsetType;

    public:
      //! Traits class
      typedef typename BaseType::Traits Traits;

      //! Type of Element (e.g. Codim 0 Entiy)
      typedef typename BaseType::ElementType ElementType;

      //! Type of indices
      typedef typename Traits::SizeType SizeType;

      typedef typename Traits::MapperTupleType MapperTupleType;

      //! constructor taking a sequence of mapper
      TupleMapper ( const Mapper & ... mapper )
        : mapperTuple_( mapper ... )
      {
        update();
      }

      //! constructor taking a sequence of mapper
      TupleMapper ( Mapper && ... mapper )
        : mapperTuple_( std::move( mapper ) ... )
      {
        update();
      }

      //! constructor taking a sequence of mapper
      TupleMapper ( const MapperTupleType &mapperTuple )
        : mapperTuple_( mapperTuple )
      {
        update();
      }

      //! constructor taking a sequence of mapper
      TupleMapper ( MapperTupleType &&mapperTuple )
        : mapperTuple_( std::move( mapperTuple ) )
      {
        update();
      }

      //! copy constructor
      TupleMapper ( const ThisType & ) = default;
      TupleMapper ( ThisType && ) = default;

      /***** DofManager Interface Methods < begin > ******/

      /** \copydoc Dune::DofMapper::size const */
      SizeType size () const
      {
        return size( Std::index_sequence_for< Mapper ... >() );
      }

      /** \copydoc Dune::DofMapper::contains(const int codim) const */
      bool contains ( const int codim ) const
      {
        return contains( codim, Std::index_sequence_for< Mapper ... >() );
      }

      /** \copydoc Dune::DofMapper::fixedDataSize(int codim) const */
      bool fixedDataSize ( int codim ) const
      {
        return fixedDataSize( codim, Std::index_sequence_for< Mapper ... >() );
      }

      /** \copydoc Dune::DofMapper::mapEach(const ElementType &element, Functor f) const */
      template< class Functor >
      void mapEach ( const ElementType &element, Functor f ) const
      {
        OffsetType localOffset;
        localOffset[ 0 ] = 0;
        ForLoop< MapEach, 0, mapperTupleSize - 1 >::apply( localOffset, globalOffset_, element, f, mapperTuple_ );
      }

      /** \copydoc Dune::DofMapper::mapEachEntityDof(const Entity &entity,Functor f) const */
      template< class Entity, class Functor >
      void mapEachEntityDof ( const Entity &entity, Functor f ) const
      {
        OffsetType localOffset;
        localOffset[ 0 ] = 0;
        ForLoop< MapEachEntityDof, 0, mapperTupleSize - 1 >::apply( localOffset, globalOffset_, entity, f, mapperTuple_ );
      }

      /** \copydoc Dune::DofMapper::maxNumDofs const */
      int maxNumDofs () const
      {
        return maxNumDofs( Std::index_sequence_for< Mapper ... >() );
      }

      /** \copydoc Dune::DofMapper::numDofs(const ElementType &element) const */
      int numDofs ( const ElementType &element ) const
      {
        return numDofs( element, Std::index_sequence_for< Mapper ... >() );
      }

      /** \copydoc Dune::DofMapper::numEntityDofs(const Entity &entity) const */
      template< class Entity >
      int numEntityDofs ( const Entity &entity ) const
      {
        return numEntityDofs( entity, Std::index_sequence_for< Mapper ... >() );
      }

      /***** DofManager Interface Methods < end > ******/

      /***** AdaptiveDofManager Interface Methods < begin > *****/

      /** \copydoc Dune::DofMapper::numberOfHoles(const int block) const */
      SizeType numberOfHoles ( const int block ) const
      {
        SizeType nHoles = 0;
        int comp = -1;
        ForLoop< NumberOfHoles, 0, mapperTupleSize - 1 >::apply( nHoles, comp, block, mapperTuple_ );
        return nHoles;
      }

      /** \copydoc Dune::DofMapper::oldIndex(const int hole, const int block) const */
      SizeType oldIndex ( const int hole, const int block ) const
      {
        int comp = -1;
        SizeType oIndex = 0;
        ForLoop< OldIndex, 0, mapperTupleSize - 1 >::apply( oIndex, comp, hole, block, mapperTuple_ );
        assert( comp >= 0 );
        return oIndex + globalOffset_[ comp ];
      }

      /** \copydoc Dune::DofMapper::newIndex(const int hole, const int block) const */
      SizeType newIndex ( const int hole, const int block ) const
      {
        int comp = -1;
        SizeType nIndex = 0;
        ForLoop< NewIndex, 0, mapperTupleSize - 1 >::apply( nIndex, comp, hole, block, mapperTuple_ );
        assert( comp > 0 );
        return nIndex + globalOffset_[ comp ];
      }

      /** \copydoc Dune::DofMapper::consecutive const */
      bool consecutive () const
      {
        return true;
      }

      /** \copydoc Dune::DofMapper::oldOffSet(const int block) const */
      SizeType oldOffSet ( const int block ) const
      {
        int comp = -1;
        SizeType oOffset = 0;
        ForLoop< OldOffset, 0, mapperTupleSize - 1 >::apply( oOffset, comp, block, mapperTuple_ );
        assert( comp >= 0 );
        return oOffset + oldGlobalOffset_[ comp ];
      }

      /** \copydoc Dune::DofMapper::offSet(const int block) const */
      SizeType offSet ( const int block ) const
      {
        int comp = -1;
        SizeType offset = 0;
        ForLoop< Offset, 0, mapperTupleSize - 1 >::apply( offset, comp, block, mapperTuple_ );
        assert( comp >= 0 );
        return offset + globalOffset_[ comp ];
      }

      SizeType numBlocks () const
      {
        return numBlocks( Std::index_sequence_for< Mapper ... >() );
      }

      /***** AdaptiveDofManager Interface Methods < end > *****/

      /***** Methods need from the DofManager < begin > *****/

      void resize () { update(); }

      bool compress ()
      {
        update();
        return true;
      }

      void backup () const {}

      void restore () { update(); }

      template< class IOStream >
      void read ( IOStream &in )
      {
        update();
      }

      template< class IOStream >
      void write ( IOStream &out ) const {}

      template< class Entity >
      void insertEntity ( const Entity & )
      {
        update();
      }

      template< class Entity >
      void removeEntity ( const Entity & )
      {
        update();
      }
      /***** Methods need from the DofManager < end > *****/

    protected:

      void update ()
      {
        oldGlobalOffset_ = globalOffset_;

        globalOffset_[ 0 ] = 0;
        // compute new offsets
        ForLoop< ComputeOffSet, 0, mapperTupleSize - 1 >::apply( globalOffset_, mapperTuple_ );
      }


      template< std::size_t ... i >
      SizeType size ( Std::index_sequence< i ... > ) const
      {
        return Std::sum( std::get< i >( mapperTuple_ ).size() ... );
      }

      // hmmm, lets see how this will work
      template< std::size_t ... i >
      bool contains ( const int codim, Std::index_sequence< i ... > ) const
      {
        return Std::Or( std::get< i >( mapperTuple_ ).contains( codim ) ... );
      }

      template< std::size_t ... i >
      int maxNumDofs ( Std::index_sequence< i ... > ) const
      {
        return Std::sum( std::get< i >( mapperTuple_ ).maxNumDofs() ... );
      }

      template< std::size_t ... i >
      SizeType numDofs ( const ElementType &element, Std::index_sequence< i ... > ) const
      {
        return Std::sum( std::get< i >( mapperTuple_ ).numDofs( element ) ... );
      }

      template< class Entity, std::size_t ... i >
      SizeType numEntityDofs ( const Entity &entity, Std::index_sequence< i ... > ) const
      {
        return Std::sum( std::get< i >( mapperTuple_ ).numEntityDofs( entity ) ... );
      }

      template< std::size_t ... i >
      bool fixedDataSize ( const int codim, Std::index_sequence< i ... > ) const
      {
        return Std::And( std::get< i >( mapperTuple_ ).fixedDataSize( codim ) ... );
      }

      template< std::size_t ... i >
      SizeType numBlocks ( Std::index_sequence< i ... > ) const
      {
        return Std::sum( std::get< i >( mapperTuple_ ).numBlocks() ... );
      }

    private:
      typename Traits::MapperTupleType mapperTuple_;
      OffsetType globalOffset_, oldGlobalOffset_;
    };


    // ComputeOffSet
    // -------------

    template< class ... Mapper >
    template< int i >
    struct TupleMapper< Mapper ... >::
    ComputeOffSet
    {
      template< class Tuple >
      static void apply ( OffsetType &offset, const Tuple &tuple )
      {
        offset[ i + 1 ] = offset[ i ] + std::get< i >( tuple ).size();
      }
    };


    // MapEachEntityDof
    // ----------------

    template< class ... Mapper >
    template< int i >
    struct TupleMapper< Mapper ... >::
    MapEachEntityDof
    {
      template< class Entity, class Functor, class Tuple >
      static void apply ( OffsetType &localOffset, const OffsetType &globalOffset, const Entity &entity, Functor f, const Tuple &tuple )
      {
        FunctorWrapper< Functor > wrappedFunctor( f, localOffset[ i ], globalOffset[ i ] );
        std::get< i >( tuple ).mapEachEntityDof( entity, wrappedFunctor );
        localOffset[ i + 1 ] = localOffset[ i ] + std::get< i >( tuple ).numEntityDofs( entity );
      }
    };


    // MapEach
    // -------

    template< class ... Mapper >
    template< int i >
    struct TupleMapper< Mapper ... >::
    MapEach
    {
      template< class Functor, class Tuple >
      static void apply ( OffsetType &localOffset, const OffsetType &globalOffset, const ElementType &element, Functor f, const Tuple &tuple )
      {
        FunctorWrapper< Functor > wrappedFunctor( f, localOffset[ i ], globalOffset[ i ] );
        std::get< i >( tuple ).mapEach( element, wrappedFunctor );
        localOffset[ i + 1 ] = localOffset[ i ] + std::get< i >( tuple ).numDofs( element );
      }
    };


    // NumberOfHoles
    // -------------

    template< class ... Mapper >
    template< int i >
    struct TupleMapper< Mapper ... >::
    NumberOfHoles
    {
      template< class Tuple >
      static void apply ( SizeType &nHoles, int &comp, const int block, const Tuple &tuple )
      {
        if( comp >= 0 )
          return;
        const int localBlock = block - std::get< i >( tuple ).numBlocks();
        if( localBlock >= 0 )
        {
          comp = i;
          nHoles = std::get< i >( tuple ).numberOfHoles( localBlock );
        }
      }
    };



    // OldIndex
    // --------

    template< class ... Mapper >
    template< int i >
    struct TupleMapper< Mapper ... >::
    OldIndex
    {
      template< class Tuple >
      static void apply ( SizeType &oldIndex, int &comp, const int hole, const int block, const Tuple &tuple )
      {
        if( comp >= 0 )
          return;
        const int localBlock = block - std::get< i >( tuple ).numBlocks();
        if( localBlock >= 0 )
        {
          comp = i;
          oldIndex = std::get< i >( tuple ).oldIndex( hole, localBlock );
        }
      }
    };


    // NewIndex
    // --------

    template< class ... Mapper >
    template< int i >
    struct TupleMapper< Mapper ... >::
    NewIndex
    {
      template< class Tuple >
      static void apply ( SizeType &nIndex, int &comp, const int hole, const int block, const Tuple &tuple )
      {
        if( comp >= 0 )
          return;
        const int localBlock = block - std::get< i >( tuple ).numBlocks();
        if( localBlock >= 0 )
        {
          comp = i;
          nIndex = std::get< i >( tuple ).newIndex( hole, localBlock );
        }
      }
    };


    // OldOffset
    // ---------

    template< class ... Mapper >
    template< int i >
    struct TupleMapper< Mapper ... >::
    OldOffset
    {
      template< class Tuple >
      static void apply ( SizeType &offset, int &comp, const int block, const Tuple &tuple )
      {
        if( comp >= 0 )
          return;
        const int localBlock = block - std::get< i >( tuple ).numBlocks();
        if( localBlock >= 0 )
        {
          comp = i;
          offset = std::get< i >( tuple ).oldOffSet( localBlock );
        }
      }
    };


    // Offset
    // ------

    template< class ... Mapper >
    template< int i >
    struct TupleMapper< Mapper ... >::
    Offset
    {
      template< class Tuple >
      static void apply ( SizeType &offset, int &comp, const int block, const Tuple &tuple )
      {
        if( comp >= 0 )
          return;
        const int localBlock = block - std::get< i >( tuple ).numBlocks();
        if( localBlock >= 0 )
        {
          comp = i;
          offset = std::get< i >( tuple ).offSet( localBlock );
        }
      }
    };

  }   // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SPACE_COMBINEDSPACE_TUPLEMAPPER_HH
