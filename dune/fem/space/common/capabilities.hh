#ifndef DUNE_FEM_SPACE_COMMON_CAPABILITIES_HH
#define DUNE_FEM_SPACE_COMMON_CAPABILITIES_HH

#include <type_traits>

#include <dune/fem/quadrature/defaultquadratures.hh>

namespace Dune
{

  namespace Fem
  {

    namespace Capabilities
    {

      /** \class hasFixedPolynomialOrder
       *
       *  \brief specialize with \b true if polynomial order does
       *         not depend on the grid (part) entity
       */
      template< class DiscreteFunctionSpace >
      struct hasFixedPolynomialOrder
      {
        static const bool v = false;
      };



      /** \class hasStaticPolynomialOrder
       *
       *  \brief specialize with \b true if polynomial order fixed
       *         and compile time static
       */
      template< class DiscreteFunctionSpace >
      struct hasStaticPolynomialOrder
      {
        static const bool v = false;
        static const int order = -1;
      };



      /** \class isContinuous
       *
       *  \brief specialize with \b true if space is always continuous
       */
      template< class DiscreteFunctionSpace >
      struct isContinuous
      {
        static const bool v = false;
      };



      /** \class isLocalized
       *
       *  \brief specialize with \b true if the space is localized, *
       *  i.e., the basis function set is based on a shape function set.
       *
       *  We require, that a localized space has a method
\code
  ShapeFunctionSetType shapeFunctionSet ( const EntityType &entity );
\endcode
       */
      template< class DiscreteFunctionSpace >
      struct isLocalized
      {
        static const bool v = false;
      };



      /** \class isAdaptive
       *
       *  \brief specialize with \b true if space can be used with
       *         AdaptiveDiscreteFunction
       */
      template< class DiscreteFunctionSpace >
      struct isAdaptive
      {
        static const bool v = false;
      };



      /** \class threadSafe
       *
       *  \brief specialize with \b true if the space implementation
       *         is thread safe
       */
      template< class DiscreteFunctionSpace >
      struct threadSafe
      {
        static const bool v = false;
      };



      /** \class viewThreadSafe
       *
       * \brief specialize with \b true if the space implementation is
       *        thread safe, while it is not modified
       *
       */
      template< class DiscreteFunctionSpace >
      struct viewThreadSafe
      {
        static const bool v = false;
      };


      /** \class isHierarchic
       *
       *  \brief specialize with \b true if for a space the basis functions are
       *         sorted by the polynomial order, starting with the lowest order
       */
      template< class DiscreteFunctionSpace >
      struct isHierarchic
      {
        static const bool v = false;
      };


      /** \class DefaultQuadrature
       *
       *  \brief specialize when quadrature other than the standard quadrature
       *  should be used for volume and surface integral compution.
       */
      template< class DiscreteFunctionSpace >
      struct DefaultQuadrature
      {
        // traits specifying the quadrature points used for CachingQuadrature or ElementQuadrature.
        template <class F, int d>
        using DefaultQuadratureTraits = Dune::Fem::DefaultQuadratureTraits< F, d >;

        //! return quadrature order for volume quadratures for given polynomial order k
        static int volumeOrder( const int k ) {  return 2 * k; }
        //! return quadrature order for surface quadratures (i.e. over intersections) for given polynomial order k
        static int surfaceOrder( const int k )   {  return 2 * k + 1; }
      };


      namespace Impl
      {

        template< class DFS >
        std::true_type hasInterpolation ( const DFS &, decltype( std::declval< const DFS & >().interpolation() ) * = nullptr );

        std::false_type hasInterpolation ( ... );

      } // namespace Impl


      /**
       * \class hasInterpolation
       *
       * \brief determine whether a discrete function space provides a (local)
       *        interpolation
       *
       * \note This capability is generated automatically by SFINAE.
       */
      template< class DiscreteFunctionSpace >
      struct hasInterpolation
      {
        static const bool v = decltype( Impl::hasInterpolation( std::declval< const DiscreteFunctionSpace & >() ) )::value;
      };



      // const specialization
      // --------------------

      template< class DiscreteFunctionSpace >
      struct hasFixedPolynomialOrder< const DiscreteFunctionSpace >
      {
        static const bool v = hasFixedPolynomialOrder< DiscreteFunctionSpace >::v;
      };

      template< class DiscreteFunctionSpace >
      struct hasStaticPolynomialOrder< const DiscreteFunctionSpace >
      {
        static const bool v = hasStaticPolynomialOrder< DiscreteFunctionSpace >::v;
        static const int order = hasStaticPolynomialOrder< DiscreteFunctionSpace >::order;
      };

      template< class DiscreteFunctionSpace >
      struct isContinuous < const DiscreteFunctionSpace >
      {
        static const bool v = isContinuous< DiscreteFunctionSpace >::v;
      };

      template< class DiscreteFunctionSpace >
      struct isLocalized< const DiscreteFunctionSpace >
      {
        static const bool v = isLocalized< DiscreteFunctionSpace >::v;
      };

      template< class DiscreteFunctionSpace >
      struct isAdaptive< const DiscreteFunctionSpace >
      {
        static const bool v = isAdaptive< DiscreteFunctionSpace >::v;
      };

      template< class DiscreteFunctionSpace >
      struct threadSafe< const DiscreteFunctionSpace >
      {
        static const bool v = threadSafe< DiscreteFunctionSpace >::v;
      };

      template< class DiscreteFunctionSpace >
      struct viewThreadSafe< const DiscreteFunctionSpace >
      {
        static const bool v = viewThreadSafe< DiscreteFunctionSpace >::v;
      };

      template< class DiscreteFunctionSpace >
      struct isHierarchic< const DiscreteFunctionSpace >
      {
        static const bool v = isHierarchic< DiscreteFunctionSpace >::v;
      };

      template< class DiscreteFunctionSpace >
      struct hasInterpolation< const DiscreteFunctionSpace >
      {
        static const bool v = hasInterpolation< DiscreteFunctionSpace >::v;
      };

      template< class DiscreteFunctionSpace >
      struct DefaultQuadrature< const DiscreteFunctionSpace >
        : public DefaultQuadrature< DiscreteFunctionSpace >
      {};

    } // namespace Capabilities

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SPACE_COMMON_CAPABILITIES_HH
