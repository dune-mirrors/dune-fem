#ifndef DUNE_FEM_HPDG_SPACE_DISCONTINUOUSGALERKIN_ORTHOGONAL_HH
#define DUNE_FEM_HPDG_SPACE_DISCONTINUOUSGALERKIN_ORTHOGONAL_HH

#include <dune/grid/common/gridenums.hh>

#include <dune/fem/operator/matrix/istlmatrixadapter.hh>

#include <dune/fem/space/common/capabilities.hh>
#include <dune/fem/space/common/commoperations.hh>
#include <dune/fem/space/common/defaultcommhandler.hh>
#include <dune/fem/space/common/localrestrictprolong.hh>
#include <dune/fem/space/discontinuousgalerkin/localrestrictprolong.hh>

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

      template< class FunctionSpace, class GridPart, int order, bool caching = true >
      class OrthogonalDiscontinuousGalerkinSpace;



#ifndef DOXYGEN

      // OrthogonalDiscontinuousGalerkinSpaceTraits
      // ------------------------------------------

      template< class FunctionSpace, class GridPart, int order, bool caching >
      struct OrthogonalDiscontinuousGalerkinSpaceTraits
      {
        using DiscreteFunctionSpaceType = hpDG::OrthogonalDiscontinuousGalerkinSpace< FunctionSpace, GridPart, order, caching >;

        using FunctionSpaceType = FunctionSpace;

        using GridPartType = GridPart;

        using BasisFunctionSetsType = hpDG::OrthogonalBasisFunctionSets< FunctionSpaceType, GridPartType, order, caching >;
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
       *  \tparam caching  enable/disable caching of quadratures
       *
       *  \ingroup DiscreteFunctionSpace_Implementation_Orthogonal
       */
      template< class FunctionSpace, class GridPart, int order, bool caching >
      class OrthogonalDiscontinuousGalerkinSpace
        : public hpDG::DiscontinuousGalerkinSpace< OrthogonalDiscontinuousGalerkinSpaceTraits< FunctionSpace, GridPart, order, caching > >
      {
        using BaseType = hpDG::DiscontinuousGalerkinSpace< OrthogonalDiscontinuousGalerkinSpaceTraits< FunctionSpace, GridPart, order, caching > >;

      public:

        static const int polynomialOrder = order ;

        using GridPartType = typename BaseType::GridPartType;
        using BasisFunctionSetsType = typename BaseType::BasisFunctionSetsType;

        explicit OrthogonalDiscontinuousGalerkinSpace ( GridPartType &gridPart,
                                                         const Dune::InterfaceType interface = Dune::InteriorBorder_All_Interface,
                                                         const Dune::CommunicationDirection direction = Dune::ForwardCommunication )
          : BaseType( gridPart, BasisFunctionSetsType{}, order, interface, direction )
        {}

        template< class Function >
        OrthogonalDiscontinuousGalerkinSpace ( GridPartType &gridPart, Function function,
                                                const Dune::InterfaceType interface = Dune::InteriorBorder_All_Interface,
                                                const Dune::CommunicationDirection direction = Dune::ForwardCommunication )
          : BaseType( gridPart, BasisFunctionSetsType{}, order, function, interface, direction )
        {}
      };

     } // namespace hpDG

#if 0
   /** \brief Local Mass Matrix for hierarchic Legendre space */
    template <class FunctionSpaceImp,
              class GridPartImp,
              int polOrd,
              bool caching,
              class VolumeQuadratureImp>
    class LocalMassMatrix<
      hpDG::OrthogonalDiscontinuousGalerkinSpace< FunctionSpaceImp, GridPartImp, polOrd, caching >, VolumeQuadratureImp >
      : public LocalMassMatrixImplementationDgOrthoNormal<
          hpDG::OrthogonalDiscontinuousGalerkinSpace< FunctionSpaceImp, GridPartImp, polOrd, caching >, VolumeQuadratureImp >
    {
      typedef hpDG::OrthogonalDiscontinuousGalerkinSpace<  FunctionSpaceImp, GridPartImp, polOrd, caching > DiscreteFunctionSpaceImp;
      typedef LocalMassMatrixImplementationDgOrthoNormal< DiscreteFunctionSpaceImp, VolumeQuadratureImp > BaseType;
    public:
      LocalMassMatrix( const DiscreteFunctionSpaceImp& spc, const int volQuadOrd = -1 )
        : BaseType( spc, volQuadOrd )
      {}
    };
#endif



#ifndef DOXYGEN

     // DefaultLocalRestrictProlong
     // ---------------------------

    template< class FunctionSpace, class GridPart, int order, bool caching >
    class DefaultLocalRestrictProlong< hpDG::OrthogonalDiscontinuousGalerkinSpace< FunctionSpace, GridPart, order, caching > >
      : public DiscontinuousGalerkinLocalRestrictProlong< hpDG::OrthogonalDiscontinuousGalerkinSpace< FunctionSpace, GridPart, order, caching >, false >
    {
      using BaseType = DiscontinuousGalerkinLocalRestrictProlong< hpDG::OrthogonalDiscontinuousGalerkinSpace< FunctionSpace, GridPart, order, caching >, false >;

    public:
      explicit DefaultLocalRestrictProlong ( const typename BaseType::DiscreteFunctionSpaceType &space )
        : BaseType( space )
      {}
    };


#if HAVE_DUNE_ISTL
    // ISTLParallelMatrixAdapter
    // -------------------------

    template< class Matrix, class FunctionSpace, class GridPart, int order, bool caching >
    struct ISTLParallelMatrixAdapter< Matrix, hpDG::OrthogonalDiscontinuousGalerkinSpace< FunctionSpace, GridPart, order, caching > >
    {
      using Type = DGParallelMatrixAdapter< Matrix >;
    };
#endif // HAVE_DUNE_ISTL

#endif // #ifndef DOXYGEN


    namespace Capabilities
    {
      ////////////////////////////////////////////////////////////////////
      //  hpDG::OrthogonalDiscontinuousGalerkinSpace
      ////////////////////////////////////////////////////////////////////

      template< class FunctionSpace, class GridPart, int polOrder, bool caching >
      struct hasStaticPolynomialOrder< hpDG::OrthogonalDiscontinuousGalerkinSpace< FunctionSpace, GridPart, polOrder, caching > >
      {
        static const bool v = true;
        static const int order = polOrder;
      };

      template< class FunctionSpace, class GridPart, int polOrder, bool caching >
      struct isLocalized< hpDG::OrthogonalDiscontinuousGalerkinSpace< FunctionSpace, GridPart, polOrder, caching > >
      {
        static const bool v = true;
      };

      template< class FunctionSpace, class GridPart, int polOrder, bool caching >
      struct isAdaptive< hpDG::OrthogonalDiscontinuousGalerkinSpace< FunctionSpace, GridPart, polOrder, caching > >
      {
        static const bool v = true;
      };

      template< class FunctionSpace, class GridPart, int polOrder, bool caching >
      struct viewThreadSafe< hpDG::OrthogonalDiscontinuousGalerkinSpace< FunctionSpace, GridPart, polOrder, caching > >
      {
        static const bool v = true;
      };

      template< class FunctionSpace, class GridPart, int polOrder, bool caching >
      struct isHierarchic< hpDG::OrthogonalDiscontinuousGalerkinSpace< FunctionSpace, GridPart, polOrder, caching > >
      {
        static const bool v = true;
      };

    } // namespace Capabilities


  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_HPDG_SPACE_DISCONTINUOUSGALERKIN_ORTHOGONAL_HH
