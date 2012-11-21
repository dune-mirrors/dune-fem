#ifndef DUNE_FEM_SPACE_RANNACHERTUREK_LOCALFINITEELEMENT_HH
#define DUNE_FEM_SPACE_RANNACHERTUREK_LOCALFINITEELEMENT_HH

// dune-localfunctions includes
#include <dune/localfunctions/rannacherturek.hh>

// dune-fem includes
#include <dune/fem/space/common/functionspace.hh>


namespace Dune
{

  namespace Fem
  {

    // LocalFiniteElement
    // ------------------

    /**
     * \brief wrapper for dune-localfunctions local finite elements
     *
     * \tparam  LocalFiniteElementImp  local finite element
     */
    template< class LocalFiniteElementImp >
    class LocalFiniteElement
    {
      typedef typename LocalFiniteElementImp::Traits Traits;

    public:
      //! \brief local basis type
      typedef typename Traits::LocalBasisType LocalBasisType;
      //! \brief local coefficients type
      typedef typename Traits::LocalCoefficientsType LocalCoefficientsType;
      //! \brief local interpolation type
      typedef typename Traits::LocalInterpolationType LocalInterpolationType;

      LocalFiniteElement ( const LocalFiniteElementImp &localFiniteElement = LocalFiniteElementImp() )
      : localFiniteElement_( localFiniteElement )
      {}
      
      //! \brief return local basis
      const LocalBasisType &localBasis () const
      {
        return localFiniteElement().localBasis();
      }

      //! \brief return local coefficients
      const LocalCoefficientsType &localCoefficients () const
      {
        return localFiniteElement().localCoefficients();
      }

      //! \brief reuturn local interpolation
      const LocalInterpolationType &localInterpolation () const
      {
        return localFiniteElement().localInterpolation ();
      }

      //! \brief return geometry type
      GeometryType type () const
      {
        return localFiniteElement().type();
      }

    protected:
      const LocalFiniteElementImp &localFiniteElement () const
      {
        return localFiniteElement_;
      }

    private:
      LocalFiniteElementImp localFiniteElement_;
    };



    // RannacherTurekLocalFiniteElement
    // ------------------

    template< class FunctionSpace >
    class RannacherTurekLocalFiniteElement;

    template< class DomainFieldType, class RangeFieldType >
    struct RannacherTurekLocalFiniteElement< Dune::Fem::FunctionSpace< DomainFieldType, RangeFieldType, 2, 1 > >
    : public LocalFiniteElement< Dune::RannacherTurek2DLocalFiniteElement< DomainFieldType, RangeFieldType > >
    {};

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SPACE_RANNACHERTUREK_LOCALFINITEELEMENT_HH
