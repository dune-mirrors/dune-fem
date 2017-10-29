#ifndef DUNE_FEM_FUNCTIONSPACEINTERFACE_HH
#define DUNE_FEM_FUNCTIONSPACEINTERFACE_HH

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>
#include <dune/fem/common/explicitfieldvector.hh>

namespace Dune
{

  namespace Fem
  {

    /** \addtogroup FunctionSpace
     *
     *  This provides the interfaces for analytical function spaces.
     *  A function space is characterized by it's domain and range field type
     *  and the two finite dimensional vector spaces over these fields.
     *
     *  All corresponding types are given in a Traits class. In addition to the
     *  domain and range types also types for the Jacobian and the Hessian must
     *  be provieded.
     *
     *  \remarks The interface for analytical functions spaces is defined by the
     *           class FunctionSpaceInterface.
     */



    /** \class FunctionSpaceInterface
     *  \brief interface for an arbitrary function space
     *  \ingroup FunctionSpace
     *
     *  Base class for specific function spaces.
     *
     *  \interfaceclass
     */
    template< typename FunctionSpaceTraits >
    class FunctionSpaceInterface
    {
    public:
      /** \brief Dimensions of domain and range */
      enum
      {
        //! dimension of domain vector space
        dimDomain = FunctionSpaceTraits :: dimDomain,
        //! dimension of range vector space
        dimRange  = FunctionSpaceTraits :: dimRange
      };

      // for compatibility with GrapDataDisplay (see dune-grid), we also export
      // the following two values:
      enum
      {
        DimDomain = dimDomain,
        DimRange = dimRange
      };

      /** \brief Intrinsic type used for values in the domain field (usually a double) */
      typedef typename FunctionSpaceTraits::DomainFieldType DomainFieldType;

      /** \brief Intrinsic type used for values in the range field (usually a double) */
      typedef typename FunctionSpaceTraits::RangeFieldType RangeFieldType;

      /** \brief Type of domain vector (using type of domain field)
          has a Dune::FieldVector type interface */
      typedef typename FunctionSpaceTraits::DomainType DomainType;

      /** \brief Type of range vector (using type of range field)
          has a Dune::FieldVector type interface */
      typedef typename FunctionSpaceTraits::RangeType RangeType;

      /** \brief Intrinsic type used for the jacobian values
          has a Dune::FieldMatrix type interface */
      typedef typename FunctionSpaceTraits::LinearMappingType JacobianRangeType;

      /** \brief Intrinsic type used for the hessian values
          has a Dune::FieldMatrix type interface */
      typedef ExplicitFieldVector< FieldMatrix< RangeFieldType, dimDomain, dimDomain >, dimRange > HessianRangeType;

      /** \brief corresponding scalar function space */
      typedef typename FunctionSpaceTraits :: ScalarFunctionSpaceType
        ScalarFunctionSpaceType;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_FUNCTIONSPACEINTERFACE_HH
