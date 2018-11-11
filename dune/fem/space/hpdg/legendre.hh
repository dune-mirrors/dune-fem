#ifndef DUNE_FEM_HPDG_SPACE_DISCONTINUOUSGALERKIN_LEGENDRE_HH
#define DUNE_FEM_HPDG_SPACE_DISCONTINUOUSGALERKIN_LEGENDRE_HH

#include <dune/grid/common/gridenums.hh>

#include <dune/fem/operator/matrix/istlmatrixadapter.hh>
#include <dune/fem/space/common/capabilities.hh>
#include <dune/fem/space/common/commoperations.hh>
#include <dune/fem/space/common/defaultcommhandler.hh>
#include <dune/fem/space/common/localrestrictprolong.hh>
#include <dune/fem/space/discontinuousgalerkin/localrestrictprolong.hh>

#include <dune/fem/space/basisfunctionset/hpdg/legendre.hh>

#include "blockmapper.hh"
#include "space.hh"

namespace Dune
{

  namespace Fem
  {

    namespace hpDG
    {

      // Internal forward declaration
      // ----------------------------

      template< class FunctionSpace, class GridPart, int order, bool caching = true >
      class LegendreDiscontinuousGalerkinSpace;

      template< class FunctionSpace, class GridPart, int order, bool caching = true >
      class HierarchicLegendreDiscontinuousGalerkinSpace;


#ifndef DOXYGEN

      // LegendreDiscontinuousGalerkinSpaceTraits
      // ----------------------------------------

      template< class FunctionSpace, class GridPart, int order, bool hierarchicalOrdering, bool caching >
      struct LegendreDiscontinuousGalerkinSpaceTraits
      {
        // select space implementation depending on basis function ordering
        typedef typename std::conditional< hierarchicalOrdering,
            HierarchicLegendreDiscontinuousGalerkinSpace< FunctionSpace, GridPart, order, caching >,
            LegendreDiscontinuousGalerkinSpace< FunctionSpace, GridPart, order, caching > >::type  DiscreteFunctionSpaceType;

        using FunctionSpaceType = FunctionSpace;

        using GridPartType = GridPart;

        using BasisFunctionSetsType = hpDG::LegendreBasisFunctionSets< FunctionSpaceType, GridPartType, order, hierarchicalOrdering, caching >;
        using BasisFunctionSetType = typename BasisFunctionSetsType::BasisFunctionSetType;

        static const int codimension = BasisFunctionSetType::EntityType::codimension;

        using BlockMapperType = hpDG::DiscontinuousGalerkinBlockMapper< GridPartType, BasisFunctionSetsType >;
        static const int localBlockSize = BasisFunctionSetsType::localBlockSize;

        typedef Hybrid::IndexRange< int, localBlockSize > LocalBlockIndices;

        template< class DiscreteFunction, class Operation = Dune::Fem::DFCommunicationOperation::Copy >
        struct CommDataHandle
        {
          using OperationType = Operation;
          using Type = Dune::Fem::DefaultCommunicationHandler< DiscreteFunction, Operation >;
        };
      };

#endif // #ifndef DOXYGEN



      // LegendreDiscontinuousGalerkinSpace
      // ----------------------------------

      /** \brief Implementation of an \f$hp\f$-adaptive discrete function space using product Legendre polynomials
       *
       *  \tparam FunctionSpace  a Dune::Fem::FunctionSpace
       *  \tparam GridPart  a Dune::Fem::GridPart
       *  \tparam order  maximum polynomial order per coordinate
       *  \tparam caching  enable/disable caching of quadratures
       *
       *  \ingroup DiscreteFunctionSpace_Implementation_Legendre
       */
      template< class FunctionSpace, class GridPart, int order, bool caching >
      class LegendreDiscontinuousGalerkinSpace
      : public hpDG::DiscontinuousGalerkinSpace< LegendreDiscontinuousGalerkinSpaceTraits< FunctionSpace, GridPart, order, false, caching > >
      {
        using BaseType = hpDG::DiscontinuousGalerkinSpace< LegendreDiscontinuousGalerkinSpaceTraits< FunctionSpace, GridPart, order, false, caching > >;

      public:
        using GridPartType = typename BaseType::GridPartType;
        using BasisFunctionSetsType = typename BaseType::BasisFunctionSetsType;

        explicit LegendreDiscontinuousGalerkinSpace ( GridPartType &gridPart,
                                                      const Dune::InterfaceType interface = Dune::InteriorBorder_All_Interface,
                                                      const Dune::CommunicationDirection direction = Dune::ForwardCommunication )
          : BaseType( gridPart, BasisFunctionSetsType{}, order, interface, direction )
        {}

        template< class Function >
        LegendreDiscontinuousGalerkinSpace ( GridPartType &gridPart, Function function,
                                             const Dune::InterfaceType interface = Dune::InteriorBorder_All_Interface,
                                             const Dune::CommunicationDirection direction = Dune::ForwardCommunication )
          : BaseType( gridPart, BasisFunctionSetsType{}, order, function, interface, direction )
        {}
      };

      // HierarchicLegendreDiscontinuousGalerkinSpace
      // --------------------------------------------

      /** \brief Implementation of an \f$hp\f$-adaptive discrete function space using product Legendre polynomials
       *
       *  \tparam FunctionSpace  a Dune::Fem::FunctionSpace
       *  \tparam GridPart  a Dune::Fem::GridPart
       *  \tparam order  maximum polynomial order per coordinate
       *  \tparam caching  enable/disable caching of quadratures
       *
       *  \ingroup DiscreteFunctionSpace_Implementation_Legendre
       */
      template< class FunctionSpace, class GridPart, int order, bool caching >
      class HierarchicLegendreDiscontinuousGalerkinSpace
      : public hpDG::DiscontinuousGalerkinSpace< LegendreDiscontinuousGalerkinSpaceTraits< FunctionSpace, GridPart, order, true, caching > >
      {
        using BaseType = hpDG::DiscontinuousGalerkinSpace< LegendreDiscontinuousGalerkinSpaceTraits< FunctionSpace, GridPart, order, true, caching > >;

      public:
        using GridPartType = typename BaseType::GridPartType;
        using BasisFunctionSetsType = typename BaseType::BasisFunctionSetsType;

        explicit HierarchicLegendreDiscontinuousGalerkinSpace ( GridPartType &gridPart,
                                                                const Dune::InterfaceType interface = Dune::InteriorBorder_All_Interface,
                                                                const Dune::CommunicationDirection direction = Dune::ForwardCommunication )
          : BaseType( gridPart, BasisFunctionSetsType{}, order, interface, direction )
        {}

        template< class Function >
        HierarchicLegendreDiscontinuousGalerkinSpace ( GridPartType &gridPart, Function function,
                                                       const Dune::InterfaceType interface = Dune::InteriorBorder_All_Interface,
                                                       const Dune::CommunicationDirection direction = Dune::ForwardCommunication )
          : BaseType( gridPart, BasisFunctionSetsType{}, order, function, interface, direction )
        {}
      };

     } // namespace hpDG



#ifndef DOXYGEN

     // DefaultLocalRestrictProlong
     // ---------------------------

    template< class FunctionSpace, class GridPart, int order, bool caching >
    class DefaultLocalRestrictProlong< hpDG::LegendreDiscontinuousGalerkinSpace< FunctionSpace, GridPart, order, caching > >
    : public DiscontinuousGalerkinLocalRestrictProlong< hpDG::LegendreDiscontinuousGalerkinSpace< FunctionSpace, GridPart, order, caching >, false >
    {
      using BaseType = DiscontinuousGalerkinLocalRestrictProlong< hpDG::LegendreDiscontinuousGalerkinSpace< FunctionSpace, GridPart, order, caching >, false >;

    public:
      explicit DefaultLocalRestrictProlong ( const typename BaseType::DiscreteFunctionSpaceType &space )
        : BaseType( space )
      {}
    };


     // DefaultLocalRestrictProlong
     // ---------------------------

    template< class FunctionSpace, class GridPart, int order, bool caching >
    class DefaultLocalRestrictProlong< hpDG::HierarchicLegendreDiscontinuousGalerkinSpace< FunctionSpace, GridPart, order, caching > >
    : public DiscontinuousGalerkinLocalRestrictProlong< hpDG::HierarchicLegendreDiscontinuousGalerkinSpace< FunctionSpace, GridPart, order, caching >, false >
    {
      using BaseType = DiscontinuousGalerkinLocalRestrictProlong< hpDG::HierarchicLegendreDiscontinuousGalerkinSpace< FunctionSpace, GridPart, order, caching >, false >;

    public:
      explicit DefaultLocalRestrictProlong ( const typename BaseType::DiscreteFunctionSpaceType &space )
        : BaseType( space )
      {}
    };


#if HAVE_DUNE_ISTL
    // ISTLParallelMatrixAdapter
    // -------------------------

    template< class Matrix, class FunctionSpace, class GridPart, int order, bool caching >
    struct ISTLParallelMatrixAdapter< Matrix, hpDG::LegendreDiscontinuousGalerkinSpace< FunctionSpace, GridPart, order, caching > >
    {
      using Type = DGParallelMatrixAdapter< Matrix >;
    };

    // ISTLParallelMatrixAdapter
    // -------------------------

    template< class Matrix, class FunctionSpace, class GridPart, int order, bool caching >
    struct ISTLParallelMatrixAdapter< Matrix, hpDG::HierarchicLegendreDiscontinuousGalerkinSpace< FunctionSpace, GridPart, order, caching > >
    {
      using Type = DGParallelMatrixAdapter< Matrix >;
    };
#endif // HAVE_DUNE_ISTL

#endif // #ifndef DOXYGEN


    namespace Capabilities
    {
      ////////////////////////////////////////////////////////////////////
      //  hpDG::LegendreDiscontinuousGalerkinSpace
      ////////////////////////////////////////////////////////////////////

      template< class FunctionSpace, class GridPart, int polOrder, bool caching >
      struct hasStaticPolynomialOrder< hpDG::LegendreDiscontinuousGalerkinSpace< FunctionSpace, GridPart, polOrder, caching > >
      {
        static const bool v = true;
        static const int order = polOrder;
      };

      template< class FunctionSpace, class GridPart, int polOrder, bool caching >
      struct isLocalized< hpDG::LegendreDiscontinuousGalerkinSpace< FunctionSpace, GridPart, polOrder, caching > >
      {
        static const bool v = true;
      };

      template< class FunctionSpace, class GridPart, int polOrder, bool caching >
      struct isAdaptive< hpDG::LegendreDiscontinuousGalerkinSpace< FunctionSpace, GridPart, polOrder, caching > >
      {
        static const bool v = true;
      };

      template< class FunctionSpace, class GridPart, int polOrder, bool caching >
      struct viewThreadSafe< hpDG::LegendreDiscontinuousGalerkinSpace< FunctionSpace, GridPart, polOrder, caching > >
      {
        static const bool v = true;
      };

      ////////////////////////////////////////////////////////////////////
      //  hpDG::HierarchicLegendreDiscontinuousGalerkinSpace
      ////////////////////////////////////////////////////////////////////

      template< class FunctionSpace, class GridPart, int polOrder, bool caching >
      struct hasStaticPolynomialOrder< hpDG::HierarchicLegendreDiscontinuousGalerkinSpace< FunctionSpace, GridPart, polOrder, caching > >
      {
        static const bool v = true;
        static const int order = polOrder;
      };

      template< class FunctionSpace, class GridPart, int polOrder, bool caching >
      struct isLocalized< hpDG::HierarchicLegendreDiscontinuousGalerkinSpace< FunctionSpace, GridPart, polOrder, caching > >
      {
        static const bool v = true;
      };

      template< class FunctionSpace, class GridPart, int polOrder, bool caching >
      struct isAdaptive< hpDG::HierarchicLegendreDiscontinuousGalerkinSpace< FunctionSpace, GridPart, polOrder, caching > >
      {
        static const bool v = true;
      };

      template< class FunctionSpace, class GridPart, int polOrder, bool caching >
      struct viewThreadSafe< hpDG::HierarchicLegendreDiscontinuousGalerkinSpace< FunctionSpace, GridPart, polOrder, caching > >
      {
        static const bool v = true;
      };

      template< class FunctionSpace, class GridPart, int polOrder, bool caching >
      struct isHierarchic< hpDG::HierarchicLegendreDiscontinuousGalerkinSpace< FunctionSpace, GridPart, polOrder, caching > >
      {
        static const bool v = true;
      };

    } // namespace Capabilities

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_HPDG_SPACE_DISCONTINUOUSGALERKIN_LEGENDRE_HH
