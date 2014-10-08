#ifndef DUNE_FEM_OPERATOR_PROJECTION_LOCAL_RIESZ_LOCALRIESZPROJECTION_HH
#define DUNE_FEM_OPERATOR_PROJECTION_LOCAL_RIESZ_LOCALRIESZPROJECTION_HH

#include <dune/common/bartonnackmanifcheck.hh>

namespace Dune
{

  namespace Fem
  {

    // LocalRieszProjection
    // --------------------

    /** \class LocalRieszProjection
     *
     *  \brief interface documentation of a local Riesz projection
     *
     *  \tparam  BasisFunctionSet  basis function set type
     */
    template< class BasisFunctionSet, class Implementation >
    class LocalRieszProjection
    {
    public:
      /** \brief basis function set */
      typedef BasisFunctionSet BasisFunctionSetType;

    protected:
#ifndef DOXYGEN

      LocalRieszProjection () = default;

#endif // #ifndef DOXYGEN

      /** \name Copying and assignment
       *  \{
       */

      /** \brief copy constructor */
      LocalRieszProjection ( const LocalRieszProjection & ) = default;

      /** \brief move constructor */
      LocalRieszProjection ( LocalRieszProjection && ) {}

      /** \brief assignment operator */
      LocalRieszProjection &operator= ( const LocalRieszProjection & ) = default;

      /** \brief move assignment operator */
      LocalRieszProjection &operator= ( LocalRieszProjection && ) {}

      /** \} */

      /** \name Public member methods
       *  \{
       */

      /** \brief return basis function set */
      BasisFunctionSet basisFunctionSet () const
      {
        CHECK_INTERFACE_IMPLEMENTATION( impl().basisFunctionSet() );
        return impl().basisFunctionSet();
      }

      /** \brief compute Riesz representation
       *
       *  \tparam  LocalDofVector  local dof vector type
       *
       *  \param[in]  f  please doc me
       *  \param[out]  localDofVector  please doc me
       */
      template< class F, class LocalDofVector >
      void operator() ( const F &f, LocalDofVector &localDofVector ) const
      {
        impl()( f, localDofVector );
      }

      /** \brief compute Riesz representation
       *
       *  \tparam  LocalDofVector  local dof vector type
       *
       *  \param[in]  f  please doc me
       *  \param[out]  localDofVector  please doc me
       */
      template< class F, class LocalDofVector >
      void apply ( const F &f, LocalDofVector &localDofVector ) const
      {
        CHECK_INTERFACE_IMPLEMENTATION( impl().apply( f, localDofVector ) );
        impl().apply( f, localDofVector );
      }

      /** \} */

    protected:
      const Implementation &impl () const
      {
        return static_cast< const Implementation & >( *this );
      }
    };

  } // namespace Fem

} // namepsace Dune

#endif // #ifndef DUNE_FEM_OPERATOR_PROJECTION_LOCAL_RIESZ_LOCALRIESZPROJECTION_HH
