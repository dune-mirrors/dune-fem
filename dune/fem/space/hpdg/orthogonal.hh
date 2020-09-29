#ifndef DUNE_FEM_HPDG_SPACE_DISCONTINUOUSGALERKIN_ORTHOGONAL_HH
#define DUNE_FEM_HPDG_SPACE_DISCONTINUOUSGALERKIN_ORTHOGONAL_HH

#include <dune/grid/common/gridenums.hh>

#include <dune/fem/operator/matrix/istlmatrixadapter.hh>

#include <dune/fem/space/common/capabilities.hh>
#include <dune/fem/space/common/commoperations.hh>
#include <dune/fem/space/common/defaultcommhandler.hh>
#include <dune/fem/space/common/localrestrictprolong.hh>
#include <dune/fem/space/discontinuousgalerkin/localrestrictprolong.hh>

#include <dune/fem/space/shapefunctionset/selectcaching.hh>
#include <dune/fem/space/basisfunctionset/hpdg/orthogonal.hh>

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

      template< class FunctionSpace, class GridPart, int order, class Storage = Fem::CachingStorage >
      class OrthogonalDiscontinuousGalerkinSpace;



#ifndef DOXYGEN

      // OrthogonalDiscontinuousGalerkinSpaceTraits
      // ------------------------------------------

      template< class FunctionSpace, class GridPart, int order, class Storage >
      struct OrthogonalDiscontinuousGalerkinSpaceTraits
      {
        using DiscreteFunctionSpaceType = hpDG::OrthogonalDiscontinuousGalerkinSpace< FunctionSpace, GridPart, order, Storage >;

        using FunctionSpaceType = FunctionSpace;

        using GridPartType = GridPart;

        using BasisFunctionSetsType = hpDG::OrthogonalBasisFunctionSets< FunctionSpaceType, GridPartType, order, Storage >;
        using BasisFunctionSetType = typename BasisFunctionSetsType::BasisFunctionSetType;

        static const int codimension = BasisFunctionSetType::EntityType::codimension;

        using BlockMapperType = hpDG::DiscontinuousGalerkinBlockMapper< GridPartType, BasisFunctionSetsType >;
        static const int localBlockSize = BasisFunctionSetsType::localBlockSize;
        static_assert( localBlockSize == FunctionSpace::dimRange, " dimRange prob ");

        typedef Hybrid::IndexRange< int, localBlockSize > LocalBlockIndices;

        template< class DiscreteFunction, class Operation = Dune::Fem::DFCommunicationOperation::Copy >
        struct CommDataHandle
        {
          using OperationType = Operation;
          using Type = Dune::Fem::DefaultCommunicationHandler< DiscreteFunction, Operation >;
        };
      };

#endif // #ifndef DOXYGEN



      // OrthogonalDiscontinuousGalerkinSpace
      // ------------------------------------

      /** \brief Implementation of an \f$hp\f$-adaptive discrete function space using orthogonal polynomials
       *
       *  \tparam FunctionSpace  a Dune::Fem::FunctionSpace
       *  \tparam GridPart  a Dune::Fem::GridPart
       *  \tparam order  maximum polynomial order per coordinate
       *  \tparam Storage  for certain caching features
       *
       *  \ingroup DiscreteFunctionSpace_Implementation_Orthogonal
       */
      template< class FunctionSpace, class GridPart, int order, class Storage >
      class OrthogonalDiscontinuousGalerkinSpace
        : public hpDG::DiscontinuousGalerkinSpace< OrthogonalDiscontinuousGalerkinSpaceTraits< FunctionSpace, GridPart, order, Storage > >
      {
        using BaseType = hpDG::DiscontinuousGalerkinSpace< OrthogonalDiscontinuousGalerkinSpaceTraits< FunctionSpace, GridPart, order, Storage > >;

      public:

        static const int polynomialOrder = order ;

        using GridPartType = typename BaseType::GridPartType;
        using EntityType   = typename BaseType::EntityType;
        using BasisFunctionSetsType = typename BaseType::BasisFunctionSetsType;

        explicit OrthogonalDiscontinuousGalerkinSpace ( GridPartType &gridPart,
                                                         const Dune::InterfaceType interface = Dune::InteriorBorder_All_Interface,
                                                         const Dune::CommunicationDirection direction = Dune::ForwardCommunication )
          : BaseType( gridPart, BasisFunctionSetsType{}, order, interface, direction )
        {}

        template <class Function,
                  std::enable_if_t<
                    std::is_arithmetic<
                      decltype(Function(std::declval<const EntityType>()))>::value,int> i=0>
        OrthogonalDiscontinuousGalerkinSpace ( GridPartType &gridPart, Function function,
                                                const Dune::InterfaceType interface = Dune::InteriorBorder_All_Interface,
                                                const Dune::CommunicationDirection direction = Dune::ForwardCommunication )
          : BaseType( gridPart, BasisFunctionSetsType{}, order, function, interface, direction )
        {}
      };

     } // namespace hpDG

   /** \brief Local Mass Matrix for hierarchic Legendre space */
    template <class FunctionSpaceImp,
              class GridPartImp,
              int polOrd,
              class Storage,
              class VolumeQuadratureImp>
    class LocalMassMatrix<
      hpDG::OrthogonalDiscontinuousGalerkinSpace< FunctionSpaceImp, GridPartImp, polOrd, Storage >, VolumeQuadratureImp >
      : public LocalMassMatrixImplementationDgOrthoNormal<
          hpDG::OrthogonalDiscontinuousGalerkinSpace< FunctionSpaceImp, GridPartImp, polOrd, Storage >, VolumeQuadratureImp >
    {
      typedef hpDG::OrthogonalDiscontinuousGalerkinSpace<  FunctionSpaceImp, GridPartImp, polOrd, Storage > DiscreteFunctionSpaceImp;
      typedef LocalMassMatrixImplementationDgOrthoNormal< DiscreteFunctionSpaceImp, VolumeQuadratureImp > BaseType;
    public:
      using BaseType::BaseType;
    };


#ifndef DOXYGEN

     // DefaultLocalRestrictProlong
     // ---------------------------

    template< class FunctionSpace, class GridPart, int order, class Storage >
    class DefaultLocalRestrictProlong< hpDG::OrthogonalDiscontinuousGalerkinSpace< FunctionSpace, GridPart, order, Storage > >
      : public DiscontinuousGalerkinLocalRestrictProlong< hpDG::OrthogonalDiscontinuousGalerkinSpace< FunctionSpace, GridPart, order, Storage >, false >
    {
      using BaseType = DiscontinuousGalerkinLocalRestrictProlong< hpDG::OrthogonalDiscontinuousGalerkinSpace< FunctionSpace, GridPart, order, Storage >, false >;

    public:
      explicit DefaultLocalRestrictProlong ( const typename BaseType::DiscreteFunctionSpaceType &space )
        : BaseType( space )
      {}
    };


#endif // #ifndef DOXYGEN


    namespace Capabilities
    {
      ////////////////////////////////////////////////////////////////////
      //  hpDG::OrthogonalDiscontinuousGalerkinSpace
      ////////////////////////////////////////////////////////////////////

      template< class FunctionSpace, class GridPart, int polOrder, class Storage >
      struct hasStaticPolynomialOrder< hpDG::OrthogonalDiscontinuousGalerkinSpace< FunctionSpace, GridPart, polOrder, Storage > >
      {
        static const bool v = true;
        static const int order = polOrder;
      };

      template< class FunctionSpace, class GridPart, int polOrder, class Storage >
      struct isLocalized< hpDG::OrthogonalDiscontinuousGalerkinSpace< FunctionSpace, GridPart, polOrder, Storage > >
      {
        static const bool v = true;
      };

      template< class FunctionSpace, class GridPart, int polOrder, class Storage >
      struct isAdaptive< hpDG::OrthogonalDiscontinuousGalerkinSpace< FunctionSpace, GridPart, polOrder, Storage > >
      {
        static const bool v = true;
      };

      template< class FunctionSpace, class GridPart, int polOrder, class Storage >
      struct viewThreadSafe< hpDG::OrthogonalDiscontinuousGalerkinSpace< FunctionSpace, GridPart, polOrder, Storage > >
      {
        static const bool v = true;
      };

      template< class FunctionSpace, class GridPart, int polOrder, class Storage >
      struct isHierarchic< hpDG::OrthogonalDiscontinuousGalerkinSpace< FunctionSpace, GridPart, polOrder, Storage > >
      {
        static const bool v = true;
      };

    } // namespace Capabilities


  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_HPDG_SPACE_DISCONTINUOUSGALERKIN_ORTHOGONAL_HH
