#ifndef DUNE_FEM_HPDG_SPACE_DISCONTINUOUSGALERKIN_LEGENDRE_HH
#define DUNE_FEM_HPDG_SPACE_DISCONTINUOUSGALERKIN_LEGENDRE_HH

#include <dune/grid/common/gridenums.hh>

#include <dune/fem/operator/matrix/istlmatrixadapter.hh>
#include <dune/fem/space/common/commoperations.hh>
#include <dune/fem/space/common/defaultcommhandler.hh>
#include <dune/fem/space/common/localrestrictprolong.hh>
#include <dune/fem/space/discontinuousgalerkin/localrestrictprolong.hh>

#include <dune/fem/hpdg/space/basisfunctionsets/legendre.hh>

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



#ifndef DOXYGEN

      // LegendreDiscontinuousGalerkinSpaceTraits
      // ----------------------------------------

      template< class FunctionSpace, class GridPart, int order, bool caching >
      struct LegendreDiscontinuousGalerkinSpaceTraits
      {
        using DiscreteFunctionSpaceType = hpDG::LegendreDiscontinuousGalerkinSpace< FunctionSpace, GridPart, order, caching >;

        using FunctionSpaceType = FunctionSpace;

        using GridPartType = GridPart;

        using BasisFunctionSetsType = hpDG::LegendreBasisFunctionSets< FunctionSpaceType, GridPartType, order, caching >;
        using BasisFunctionSetType = typename BasisFunctionSetsType::BasisFunctionSetType;

        static const int codimension = BasisFunctionSetType::EntityType::codimension;

        using BlockMapperType = hpDG::DiscontinuousGalerkinBlockMapper< GridPartType, BasisFunctionSetsType >;
        static const int localBlockSize = BasisFunctionSetsType::localBlockSize;

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
      : public hpDG::DiscontinuousGalerkinSpace< LegendreDiscontinuousGalerkinSpaceTraits< FunctionSpace, GridPart, order, caching > >
      {
        using BaseType = hpDG::DiscontinuousGalerkinSpace< LegendreDiscontinuousGalerkinSpaceTraits< FunctionSpace, GridPart, order, caching > >;

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



    // ISTLParallelMatrixAdapter
    // -------------------------

    template< class Matrix, class FunctionSpace, class GridPart, int order, bool caching >
    struct ISTLParallelMatrixAdapter< Matrix, hpDG::LegendreDiscontinuousGalerkinSpace< FunctionSpace, GridPart, order, caching > >
    {
      using Type = DGParallelMatrixAdapter< Matrix >;
    };

#endif // #ifndef DOXYGEN

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_HPDG_SPACE_DISCONTINUOUSGALERKIN_LEGENDRE_HH
