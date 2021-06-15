#ifndef DUNE_FEM_SPACE_FINITEVOLUME_INTERPOLATION_HH
#define DUNE_FEM_SPACE_FINITEVOLUME_INTERPOLATION_HH

#include <functional>

#include <dune/fem/function/localfunction/average.hh>

#include "basisfunctionset.hh"

namespace Dune
{

  namespace Fem
  {

    // FiniteVolumeLocalInterpolation
    // ------------------------------

    template< class GridPart, class Range >
    class FiniteVolumeLocalInterpolation
    {
      typedef FiniteVolumeLocalInterpolation< GridPart, Range > ThisType;

    public:
      /** \brief grid part type */
      typedef GridPart GridPartType;
      /** \brief entity type */
      typedef typename GridPartType::template Codim< 0 >::EntityType EntityType;

      /** \brief basis function set type */
      typedef FiniteVolumeBasisFunctionSet< EntityType, Range > BasisFunctionSetType;

      /** \name Construction
       * \{
       */

      FiniteVolumeLocalInterpolation () {}

      void bind( const EntityType &entity ) {}
      void unbind() {}

      /*
      explicit FiniteVolumeLocalInterpolation ( const EntityType &entity )
        : entity_( entity )
      {}
      */

      /** \} */

      /** \name Copying and assignment
       * \{
       */

      FiniteVolumeLocalInterpolation ( const ThisType & ) = default;

      FiniteVolumeLocalInterpolation &operator= ( const ThisType & ) = default;

      /** \} */

      /** \name Public member methods
       *  \{
       */

      /** \brief return basis function set */
      /*
      BasisFunctionSetType basisFunctionSet () const
      {
        return BasisFunctionSetType( entity() );
      }
      */

      /** \brief interpolate local function */
      template< class LocalFunction, class LocalDofVector >
      void operator() ( const LocalFunction &localFunction, LocalDofVector &localDofVector ) const
      {
        apply( localFunction, localDofVector );
      }

      /** \brief interpolate local function */
      template< class LocalFunction, class LocalDofVector >
      void apply ( const LocalFunction &localFunction, LocalDofVector &localDofVector ) const
      {
        Range value;
        LocalAverage< LocalFunction, GridPartType >::apply( localFunction, value );
        for( int i = 0; i < Range::dimension; ++i )
          localDofVector[ i ] = value[ i ];
      }

      /** \} */

    private:
      //const EntityType &entity () const { return entity_.get(); }

      //std::reference_wrapper< const EntityType > entity_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SPACE_FINITEVOLUME_INTERPOLATION_HH
