#ifndef DUNE_FEM_HPDG_SPACE_DISCONTINUOUSGALERKIN_ANISOTROPIC_HH
#define DUNE_FEM_HPDG_SPACE_DISCONTINUOUSGALERKIN_ANISOTROPIC_HH

#include <algorithm>
#include <utility>

#include <dune/grid/common/gridenums.hh>

#include <dune/fem/operator/matrix/istlmatrixadapter.hh>
#include <dune/fem/space/common/commoperations.hh>
#include <dune/fem/space/common/defaultcommhandler.hh>
#include <dune/fem/space/common/localrestrictprolong.hh>
#include <dune/fem/space/discontinuousgalerkin/localrestrictprolong.hh>

#include <dune/fem/space/basisfunctionset/hpdg/anisotropic.hh>

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
      class AnisotropicDiscontinuousGalerkinSpace;



#ifndef DOXYGEN

      // AnisotropicDiscontinuousGalerkinSpaceTraits
      // -------------------------------------------

      template< class FunctionSpace, class GridPart, int order, bool caching >
      struct AnisotropicDiscontinuousGalerkinSpaceTraits
      {
        using DiscreteFunctionSpaceType = hpDG::AnisotropicDiscontinuousGalerkinSpace< FunctionSpace, GridPart, order >;

        using FunctionSpaceType = FunctionSpace;

        using GridPartType = GridPart;

        using BasisFunctionSetsType = hpDG::AnisotropicBasisFunctionSets< FunctionSpaceType, GridPartType, order, caching >;
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



      // AnisotropicDiscontinuousGalerkinSpace
      // --------------------------------

      /** \brief Implementation of an \f$hp\f$-adaptive discrete function space using anisotropic product Legendre polynomials
       *
       *  \tparam FunctionSpace  a Dune::Fem::FunctionSpace
       *  \tparam GridPart  a Dune::Fem::GridPart
       *  \tparam order  maximum polynomial order per coordinate
       *  \tparam caching  enable/disable caching of quadratures
       *
       *  \ingroup DiscreteFunctionSpace_Implementation_Anisotropic
       */
      template< class FunctionSpace, class GridPart, int order, bool caching >
      class AnisotropicDiscontinuousGalerkinSpace
      : public hpDG::DiscontinuousGalerkinSpace< AnisotropicDiscontinuousGalerkinSpaceTraits< FunctionSpace, GridPart, order, caching > >
      {
        using BaseType = hpDG::DiscontinuousGalerkinSpace< AnisotropicDiscontinuousGalerkinSpaceTraits< FunctionSpace, GridPart, order, caching > >;

      public:
        using GridPartType = typename BaseType::GridPartType;
        using BasisFunctionSetsType = typename BaseType::BasisFunctionSetsType;

        explicit AnisotropicDiscontinuousGalerkinSpace ( GridPartType &gridPart,
                                                         const Dune::InterfaceType interface = Dune::InteriorBorder_All_Interface,
                                                         const Dune::CommunicationDirection direction = Dune::ForwardCommunication )
          : BaseType( gridPart, BasisFunctionSetsType(), defaultKey(), interface, direction )
        {}

        template< class Function >
        AnisotropicDiscontinuousGalerkinSpace ( GridPartType &gridPart, Function function,
                                                const Dune::InterfaceType interface = Dune::InteriorBorder_All_Interface,
                                                const Dune::CommunicationDirection direction = Dune::ForwardCommunication )
          : BaseType( gridPart, BasisFunctionSetsType(), defaultKey(), function, interface, direction )
        {}

      private:
        static typename BaseType::KeyType defaultKey ()
        {
          typename BaseType::KeyType key;
          std::fill( key.begin(), key.end(), order );
          return std::move( key );
        }
      };

    } // namespace hpDG



#ifndef DOXYGEN

    // DefaultLocalRestrictProlong
    // ---------------------------

    template< class FunctionSpace, class GridPart, int order, bool caching >
    class DefaultLocalRestrictProlong< hpDG::AnisotropicDiscontinuousGalerkinSpace< FunctionSpace, GridPart, order, caching > >
    : public DiscontinuousGalerkinLocalRestrictProlong< hpDG::AnisotropicDiscontinuousGalerkinSpace< FunctionSpace, GridPart, order, caching >, false >
    {
      using BaseType = DiscontinuousGalerkinLocalRestrictProlong< hpDG::AnisotropicDiscontinuousGalerkinSpace< FunctionSpace, GridPart, order, caching >, false >;

    public:
      explicit DefaultLocalRestrictProlong ( const typename BaseType::DiscreteFunctionSpaceType &space )
        : BaseType( space )
      {}
    };


#if HAVE_DUNE_ISTL
    // ISTLParallelMatrixAdapter
    // -------------------------

    template< class Matrix, class Traits >
    struct ISTLParallelMatrixAdapter< Matrix, hpDG::DiscontinuousGalerkinSpace< Traits > >
    {
      using Type = DGParallelMatrixAdapter< Matrix >;
    };
#endif // HAVE_DUNE_ISTL

#endif //#ifndef DOXYGEN

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_HPDG_SPACE_DISCONTINUOUSGALERKIN_ANISOTROPIC_HH
