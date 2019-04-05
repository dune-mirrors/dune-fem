#ifndef DUNE_FEM_HPDG_SPACE_BASISFUNCTIONSETS_VECTORIAL_HH
#define DUNE_FEM_HPDG_SPACE_BASISFUNCTIONSETS_VECTORIAL_HH

#include <cstddef>

#include <dune/common/exceptions.hh>

#include <dune/geometry/type.hh>

#include <dune/fem/space/basisfunctionset/vectorial.hh>

#include "basisfunctionsets.hh"

namespace Dune
{

  namespace Fem
  {

    namespace hpDG
    {

      // Internal forward declaration
      // ----------------------------

      template< class BasisFunctionSets, class Range >
      class VectorialBasisFunctionSets;



#ifndef DOXYGEN

      // VectorialBasisFunctionSetTraits
      // -------------------------------

      template< class BasisFunctionSets, class Range >
      class VectorialBasisFunctionSetsTraits
      {
      public:
        using ImplementationType = VectorialBasisFunctionSets< BasisFunctionSets, Range >;

        using GridPartType = typename BasisFunctionSets::GridPartType;
        using Types = typename BasisFunctionSets::Types;

        using KeyType = typename BasisFunctionSets::KeyType;
        using BasisFunctionSetType = Dune::Fem::VectorialBasisFunctionSet< typename BasisFunctionSets::BasisFunctionSetType, Range >;
        static const int localBlockSize = Range::dimension * BasisFunctionSets::localBlockSize;

        using DataType = void *;
      };

#endif // #ifndef DOXYGEN



      // VectorialBasisFunctionSet
      // -------------------------

      /** \brief A meta implemenation of a family of local basis function sets
       *
       *  \tparam BasisFunctionSets  a scalar family of local basis function sets
       *  \tparam Range  the new range type
       *
       *  \ingroup DiscreteFunctionSpace_API
       */
      template< class BasisFunctionSets, class Range >
      class VectorialBasisFunctionSets
        : public Dune::Fem::hpDG::BasisFunctionSets< VectorialBasisFunctionSetsTraits< BasisFunctionSets, Range > >
      {
        using BaseType = Dune::Fem::hpDG::BasisFunctionSets< VectorialBasisFunctionSetsTraits< BasisFunctionSets, Range > >;

      public:
        /** \copydoc Dune::Fem::BasisFunctionSets::KeyType */
        using KeyType = typename BaseType::KeyType;

        /** \copydoc Dune::Fem::BasisFunctionSets::BasisFunctionSetType */
        using BasisFunctionSetType = typename BaseType::BasisFunctionSetType;
        /** \copydoc Dune::Fem::BasisFunctionSets::EntityType */
        using EntityType = typename BaseType::EntityType;

        /** \copydoc Dune::Fem::BasisFunctionSets::DataType */
        using DataType = typename BaseType::DataType;

        /** \name Construction
         *  \{
         */

        explicit VectorialBasisFunctionSets ( const BasisFunctionSets &basisFunctionSets )
          : basisFunctionSets_( basisFunctionSets )
        {}

        /** \} */

        /** \name Public member methods
         *  \{
         */

        /** \copydoc Dune::Fem::BasisFunctionSets::types */
        typename BaseType::Types types () const { return impl().types(); }

        /** \copydoc Dune::Fem::BasisFunctionSets::maxBlocks */
        std::size_t maxBlocks () const { return impl().maxBlocks(); }

        /** \copydoc Dune::Fem::BasisFunctionSets::maxBlocks */
        std::size_t maxBlocks ( GeometryType type ) const { return impl().maxBlocks(); }

        /** \copydoc Dune::Fem::BasisFunctionSets::blocks */
        std::size_t blocks ( GeometryType type, const KeyType &key ) const { return impl().maxBlocks(); }

        /** \brief map key to data type */
        DataType encode ( const KeyType &key ) const
        {
          DUNE_THROW( NotImplemented, "Method encode() not implemented yet" );
        }

        /** \brief map data to key type */
        KeyType decode ( const DataType &data ) const
        {
          DUNE_THROW( NotImplemented, "Method decode() not implemented yet" );
        }

        /** \copydoc Dune::Fem::BasisFunctionSets::orthogonal */
        static constexpr bool orthogonal () noexcept { return BasisFunctionSets::orthogonal(); }

        /** \copydoc Dune::Fem::BasisFunctionSets::order */
        int order () const { return impl().order(); }

        /** \copydoc Dune::Fem::BasisFunctionSets::order */
        int order ( GeometryType type ) const { return impl().order( type ); }

        /** \copydoc Dune::Fem::BasisFunctionSets::order */
        int order ( GeometryType type, const KeyType &key ) const { return impl().order( type, key ); }

        /** \copydoc Dune::Fem::BasisFunctionSets::size */
        std::size_t size ( GeometryType type, const KeyType &key ) const
        {
          return static_cast< std::size_t >( Range::dimension ) * impl().order( type, key );
        }

        /** \copydoc Dune::Fem::BasisFunctionSets::basisFunctionSet */
        BasisFunctionSetType basisFunctionSet ( const EntityType &entity, const KeyType &key ) const
        {
          return BasisFunctionSetType( impl().basisFunctionSet( entity, key ) );
        }

        /** \} */

        /** \name Non-interface methods
         *  \{
         */

        /** \brief return scalar basis function sets */
        const BasisFunctionSets &impl() const { return basisFunctionSets_; }

        /** \} */

      private:
        BasisFunctionSets basisFunctionSets_;
      };

    } // namespace hpDG

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_HPDG_SPACE_BASISFUNCTIONSETS_VECTORIAL_HH
