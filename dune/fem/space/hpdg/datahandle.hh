#ifndef DUNE_FEM_HPDG_SPACE_DISCONTINUOUSGALERKIN_DATAHANDLE_HH
#define DUNE_FEM_HPDG_SPACE_DISCONTINUOUSGALERKIN_DATAHANDLE_HH

#include <cassert>
#include <cstddef>

#include <functional>

#include <dune/grid/common/datahandleif.hh>

namespace Dune
{

  namespace Fem
  {

    namespace hpDG
    {

      // External forward declaration
      // ----------------------------

      template< class GridPart, class LocalKeys >
      struct DiscontinuousGalerkinBlockMapper;



      // Internal forward declaration
      // ----------------------------

      template< class BlockMapper >
      class DataHandle;



      // DataHandle
      // ----------

      template< class GridPart, class LocalKeys >
      class DataHandle< DiscontinuousGalerkinBlockMapper< GridPart, LocalKeys > >
      : public Dune::CommDataHandleIF< DataHandle< DiscontinuousGalerkinBlockMapper< GridPart, LocalKeys > >, typename LocalKeys::DataType >
      {
        using ThisType = DataHandle< DiscontinuousGalerkinBlockMapper< GridPart, LocalKeys > >;
        using BaseType = Dune::CommDataHandleIF< DataHandle< DiscontinuousGalerkinBlockMapper< GridPart, LocalKeys > >, typename LocalKeys::DataType >;

        using BlockMapperType = DiscontinuousGalerkinBlockMapper< GridPart, LocalKeys >;

      public:
        /** \brief key type */
        using KeyType = typename LocalKeys::KeyType;
        /** \brief data type */
        using DataType = typename BaseType::DataType;

        /** \name Construction
         *  \{
         */

        explicit DataHandle ( BlockMapperType &blockMapper )
          : blockMapper_( blockMapper )
        {}

        /** \} */

        /** \name Copying and assignment
         *  \{
         */

        DataHandle ( const ThisType & ) = default;

        ThisType &operator= ( const ThisType & ) = default;

        DataHandle ( ThisType && ) = default;

        ThisType &operator= ( ThisType && ) = default;

        /** \} */

        /** \name Public member methods
         *  \{
         */

        bool contains ( int dim, int codim ) const { return (codim == 0); }

        bool fixedSize ( int dim, int codim ) const { return true; }

        template< class Entity >
        std::size_t size ( const Entity &entity ) const
        {
          return (Entity::codimension == 0 ? 1u : 0u);
        }

        template< class Buffer, class Entity >
        void gather ( Buffer &buffer, const Entity &entity ) const
        {
          const auto &keys = blockMapper().keys_[ entity ];
          buffer.write( encode( keys.second ) );
        }

        template< class Buffer, class Entity >
        void scatter ( Buffer &buffer, const Entity &entity, std::size_t n )
        {
          assert( n == 1u );
          auto &keys = blockMapper().keys_[ entity ];
          DataType data;
          buffer.read( data );
          keys.second = decode( data );
        }

        /** \} */

      private:
        DataType encode ( const KeyType &key ) const { return localKeys().encode( key ); }

        KeyType decode ( const DataType &data ) const { return localKeys().decode( data ); }

        KeyType reduce ( const KeyType &a, const KeyType &b ) const { return localKeys().reduce( a, b ); }

        const LocalKeys &localKeys () const { return blockMapper().localKeys(); }

        BlockMapperType &blockMapper () { return blockMapper_.get(); }

        const BlockMapperType &blockMapper () const { return blockMapper_.get(); }

        std::reference_wrapper< BlockMapperType > blockMapper_;
      };

    } // namespace hpDG

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_HPDG_SPACE_DISCONTINUOUSGALERKIN_DATAHANDLE_HH
