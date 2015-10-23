#ifndef DUNE_FEM_SPACE_FINITEVOLUME_BASISFUNCTIONSETS_HH
#define DUNE_FEM_SPACE_FINITEVOLUME_BASISFUNCTIONSETS_HH

#include "basisfunctionset.hh"

namespace Dune
{

  namespace Fem
  {

    // FiniteVolumeBasisFunctionSets
    // -----------------------------

    template< class Entity, class Range >
    class FiniteVolumeBasisFunctionSets
    {
      typedef FiniteVolumeBasisFunctionSets< Entity, Range > ThisType;

    public:
      /** \copydoc Dune::Fem::BasisFunctionSets::BasisFunctionSetType */
      typedef FiniteVolumeBasisFunctionSet< Entity, Range > BasisFunctionSetType;
      /** \copydoc Dune::Fem::BasisFunctionSets::EntityType */
      typedef typename BasisFunctionSetType::EntityType EntityType;

      /** \name Construction
       *  \{
       */

      FiniteVolumeBasisFunctionSets () {}

      /** \} */

      /** \name Copying and assignment
       *  \{
       */

      FiniteVolumeBasisFunctionSets ( const ThisType & ) = default;

      FiniteVolumeBasisFunctionSets &operator= ( const ThisType & ) = default;

      /** \} */

      /** \name Public member methods
       *  \{
       */

      /** \copydoc Dune::Fem::BasisFunctionSets::order */
      static constexpr int order () { return 0; }

      /** \copydoc Dune::Fem::BasisFunctionSets::order */
      static constexpr int order ( const EntityType & ) { return 0; }

      /** \copydoc Dune::Fem::BasisFunctionSets::basisFunctionSet */
      static BasisFunctionSetType basisFunctionSet ( const EntityType &entity )
      {
        return BasisFunctionSetType( entity );
      }

      /** \} */
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SPACE_FINITEVOLUME_BASISFUNCTIONSETS_HH
