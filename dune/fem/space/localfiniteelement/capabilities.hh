#ifndef DUNE_FEM_SPACE_LOCALFINITEELEMENT_CAPABILITIES_HH
#define DUNE_FEM_SPACE_LOCALFINITEELEMENT_CAPABILITIES_HH



#include <dune/fem/gridpart/common/capabilities.hh>
#include <dune/fem/space/common/capabilities.hh>
#include <dune/fem/space/shapefunctionset/selectcaching.hh>

#include <dune/fem/quadrature/interpolationquadrature.hh>

namespace Dune
{

  namespace Fem
  {

    // External Forward Declarations
    // -----------------------------

    template< class LFEMap, class FunctionSpace, class Storage = CachingStorage >
    class LocalFiniteElementSpace;

    template< class LFEMap, class FunctionSpace, class Storage = CachingStorage >
    class DiscontinuousLocalFiniteElementSpace;

    template< class FunctionSpace, class GridPart, unsigned int order,
              template< class, unsigned int > class PointSet>
    struct FixedOrderLagrangeFiniteElementMap;


    namespace Capabilities
    {

      template< class LFEMap, class FunctionSpace, class Storage >
      struct hasFixedPolynomialOrder< LocalFiniteElementSpace< LFEMap, FunctionSpace, Storage > >
      {
        static const bool v = false;
      };

      template< class LFEMap, class FunctionSpace, class Storage >
      struct hasFixedPolynomialOrder< DiscontinuousLocalFiniteElementSpace< LFEMap, FunctionSpace, Storage > >
      {
        static const bool v = false;
      };


      template< class LFEMap, class FunctionSpace, class Storage >
      struct hasStaticPolynomialOrder< LocalFiniteElementSpace< LFEMap, FunctionSpace, Storage > >
      {
        static const bool v = false;
        static const int order = 6; // default polynomial order if not specified otherwise
      };

      template< class LFEMap, class FunctionSpace, class Storage >
      struct hasStaticPolynomialOrder< DiscontinuousLocalFiniteElementSpace< LFEMap, FunctionSpace, Storage > >
      {
        static const bool v = false;
        static const int order = 6; // default polynomial order if not specified otherwise
      };


      template< class LFEMap, class FunctionSpace, class Storage >
      struct isContinuous< LocalFiniteElementSpace< LFEMap, FunctionSpace, Storage > >
      {
        static const bool v = false;
      };

      template< class LFEMap, class FunctionSpace, class Storage >
      struct isContinuous< DiscontinuousLocalFiniteElementSpace< LFEMap, FunctionSpace, Storage > >
      {
        static const bool v = false;
      };


      template< class LFEMap, class FunctionSpace, class Storage >
      struct isLocalized< LocalFiniteElementSpace< LFEMap, FunctionSpace, Storage > >
      {
        static const bool v = true;
      };

      template< class LFEMap, class FunctionSpace, class Storage >
      struct isLocalized< DiscontinuousLocalFiniteElementSpace< LFEMap, FunctionSpace, Storage > >
      {
        static const bool v = true;
      };


      template< class LFEMap, class FunctionSpace, class Storage >
      struct isAdaptive< LocalFiniteElementSpace< LFEMap, FunctionSpace, Storage > >
      {
        static const bool v = false;
      };

      template< class LFEMap, class FunctionSpace, class Storage >
      struct isAdaptive< DiscontinuousLocalFiniteElementSpace< LFEMap, FunctionSpace, Storage > >
      {
        static const bool v = false;
      };


      template< class LFEMap, class FunctionSpace, class Storage >
      struct threadSafe< LocalFiniteElementSpace< LFEMap, FunctionSpace, Storage > >
      {
        static const bool v = false;
      };

      template< class LFEMap, class FunctionSpace, class Storage >
      struct threadSafe< DiscontinuousLocalFiniteElementSpace< LFEMap, FunctionSpace, Storage > >
      {
        static const bool v = false;
      };


      template< class LFEMap, class FunctionSpace, class Storage >
      struct viewThreadSafe< LocalFiniteElementSpace< LFEMap, FunctionSpace, Storage > >
      {
        static const bool v = false;
      };

      template< class LFEMap, class FunctionSpace, class Storage >
      struct viewThreadSafe< DiscontinuousLocalFiniteElementSpace< LFEMap, FunctionSpace, Storage > >
      {
        static const bool v = false;
      };


#if HAVE_DUNE_LOCALFUNCTIONS
      namespace detail
      {

        struct DefaultQuadratureEquidistant
        {
          template <class F, int d>
          using DefaultQuadratureTraits = Dune::Fem::EquidistantQuadratureTraits< F, d >;

          // TODO: double check this
          static int volumeOrder ( const int k ) { return k; }
          static int surfaceOrder( const int k ) { return k; }
        };

        struct DefaultQuadratureGaussLobatto
        {
          template <class F, int d>
          using DefaultQuadratureTraits = Dune::Fem::GaussLobattoQuadratureTraits< F, d >;

          static int volumeOrder ( const int k ) { return (k > 0) ? (2 * k - 1) : 0; }
          static int surfaceOrder( const int k ) { return (k > 0) ? (2 * k - 1) : 0; }
        };

        struct DefaultQuadratureGaussLegendre
        {
          template <class F, int d>
          using DefaultQuadratureTraits = Dune::Fem::GaussLegendreQuadratureTraits< F, d >;

          static int volumeOrder ( const int k ) { return 2 * k + 1; }
          static int surfaceOrder( const int k ) { return 2 * k + 1; }
        };

        struct DefaultQuadratureCellCenters
        {
          template <class F, int d>
          using DefaultQuadratureTraits = Dune::Fem::CellCentersQuadratureTraits< F, d >;

          static int volumeOrder ( const int k ) { return k; }
          static int surfaceOrder( const int k ) { return k; }
        };

        // default uses the default values for all spaces (see space/common/capabilities.hh)
        template< class LFEMap >
        struct DefaultQuadratureSpec : public Dune::Fem::Capabilities::DefaultQuadrature< LFEMap >
        {};


        ///-----  Specialization for Equidistant  ----
        template < class FunctionSpace, class GridPart, unsigned int order >
        struct DefaultQuadratureSpec< Dune::Fem::FixedOrderLagrangeFiniteElementMap< FunctionSpace, GridPart, order, Dune::EquidistantPointSetDerived > >
          : public DefaultQuadratureEquidistant {};

        ///-----  Specialization for GaussLobatto ----
        template < class FunctionSpace, class GridPart, unsigned int order >
        struct DefaultQuadratureSpec< Dune::Fem::FixedOrderLagrangeFiniteElementMap< FunctionSpace, GridPart, order, Dune::GaussLobattoPointSet > >
          : public DefaultQuadratureGaussLobatto {};

        ///-----  Specialization for GaussLegendre ----
        template < class FunctionSpace, class GridPart, unsigned int order >
        struct DefaultQuadratureSpec< Dune::Fem::FixedOrderLagrangeFiniteElementMap< FunctionSpace, GridPart, order, Dune::GaussLegendrePointSet > >
          : public DefaultQuadratureGaussLegendre {};

        ///-----  Specialization for CellCenters  ----
        template < class FunctionSpace, class GridPart, unsigned int order >
        struct DefaultQuadratureSpec< Dune::Fem::FixedOrderLagrangeFiniteElementMap< FunctionSpace, GridPart, order, Dune::CellCentersPointSet > >
          : public DefaultQuadratureCellCenters {};
      } // end namespace detail

      template< class LFEMap, class FunctionSpace, class Storage >
      struct DefaultQuadrature< LocalFiniteElementSpace< LFEMap, FunctionSpace, Storage > >
        : public detail::DefaultQuadratureSpec< LFEMap >
      {
      };

      template< class LFEMap, class FunctionSpace, class Storage >
      struct DefaultQuadrature< DiscontinuousLocalFiniteElementSpace< LFEMap, FunctionSpace, Storage > >
        : public detail::DefaultQuadratureSpec< LFEMap >
      {
      };
#endif // #ifndef HAVE_DUNE_LOCALFUNCTIONS


    } // namespace Capabilities

  } // namespace Fem

} // namespace Dune


#endif // #ifndef DUNE_FEM_SPACE_LOCALFINITEELEMENT_CAPABILITIES_HH
