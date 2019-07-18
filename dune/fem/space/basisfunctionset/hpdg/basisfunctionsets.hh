#ifndef DUNE_FEM_HPDG_SPACE_BASISFUNCTIONSETS_BASISFUNCTIONSETS_HH
#define DUNE_FEM_HPDG_SPACE_BASISFUNCTIONSETS_BASISFUNCTIONSETS_HH

#include <cstddef>

#include <dune/common/bartonnackmanifcheck.hh>

#include <dune/geometry/type.hh>

namespace Dune
{

  namespace Fem
  {

    namespace hpDG
    {

      // BasisFunctionSets
      // -----------------

      /** \brief abstract interface class for a family of local basis function sets
       *
       *  \tparam  T  traits class
       *
       *  \ingroup DiscreteFunctionSpace_API
       */
      template< class T >
      class BasisFunctionSets
      {
      protected:
        using Traits = T;

      public:
        /** \brief grid part type */
        using GridPartType = typename Traits::GridPartType;

        /** \brief key type */
        using KeyType = typename Traits::KeyType;

        /** \brief basis function set type */
        using BasisFunctionSetType = typename Traits::BasisFunctionSetType;
        /** \brief entity type */
        using EntityType = typename BasisFunctionSetType::EntityType;

        /** \brief a range of geometry types */
        using Types = typename Traits::Types;

        /** \brief block size */
        static const int localBlockSize = Traits::localBlockSize;

        /** \brief data type */
        using DataType = typename Traits::DataType;

      protected:
        BasisFunctionSets () = default;

      public:
        /** \name Supported geometry types
         *  \{
         */

        /** \brief return range of supported geometry types */
        Types types () const
        {
          CHECK_INTERFACE_IMPLEMENTATION( impl().types() );
          return impl().types();
        }

        /** \} */

        /** \name Block usage
         *  \{
         */

        /** \brief return maximum number of blocks used */
        std::size_t maxBlocks () const
        {
          CHECK_INTERFACE_IMPLEMENTATION( impl().maxBlocks() );
          return impl().maxBlocks();
        }

        /** \brief return maximum number of blocks used per geometry type */
        std::size_t maxBlocks ( GeometryType type ) const
        {
          CHECK_INTERFACE_IMPLEMENTATION( impl().maxBlocks( type ) );
          return impl().maxBlocks( type );
        }

        /** \brief return number of blocks used */
        std::size_t blocks ( GeometryType type, const KeyType &key ) const
        {
          CHECK_INTERFACE_IMPLEMENTATION( impl().blocks( type, key ) );
          return impl().blocks( type, key );
        }

        /** \} */

        /** \name Communication
         *  \{
         */

        /** \brief map key to data type */
        DataType encode ( const KeyType &key ) const
        {
          CHECK_INTERFACE_IMPLEMENTATION( impl().encode( key ) );
          return impl().encode( key );
        }

        /** \brief map data to key type */
        KeyType decode ( const DataType &data ) const
        {
          CHECK_INTERFACE_IMPLEMENTATION( impl().decode( data ) );
          return impl().decode( data );
        }

        /** \} */

        /** \name Basis function sets
         *  \{
         */

        /** \brief return \b true if basis function sets are always orthogonal, \b false otherwise */
        static constexpr bool orthogonal () noexcept
        {
          return Traits::ImplementationType::orthogonal();
        }

        /** \brief return maximum polynomial order */
        int order () const
        {
          CHECK_INTERFACE_IMPLEMENTATION( impl().order() );
          return impl().order();
        }

        /** \brief return maximum polynomial order per geometry type */
        int order ( GeometryType type ) const
        {
          CHECK_INTERFACE_IMPLEMENTATION( impl().order( type ) );
          return impl().order( type );
        }

        /** \brief return polynomial order */
        int order ( GeometryType type, const KeyType &key ) const
        {
          CHECK_INTERFACE_IMPLEMENTATION( impl().order( type, key ) );
          return impl().order( type, key );
        }

        /** \brief return basis function set */
        BasisFunctionSetType basisFunctionSet ( const EntityType &entity, const KeyType &key ) const
        {
          CHECK_INTERFACE_IMPLEMENTATION( impl().basisFunctionSet( entity, key ) );
          return impl().basisFunctionSet( entity, key );
        }

        /** \} */

      protected:
        const typename Traits::ImplementationType &impl () const
        {
          return static_cast< const typename Traits::ImplementationType & >( *this );
        }
      };

    } // namespace hpDG

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_HPDG_SPACE_BASISFUNCTIONSETS_BASISFUNCTIONSETS_HH
