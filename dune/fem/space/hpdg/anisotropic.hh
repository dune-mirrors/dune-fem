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

#include <dune/fem/space/shapefunctionset/selectcaching.hh>
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

      template< class FunctionSpace, class GridPart, int order, class Storage = Fem::CachingStorage >
      class AnisotropicDiscontinuousGalerkinSpace;



#ifndef DOXYGEN

      // AnisotropicDiscontinuousGalerkinSpaceTraits
      // -------------------------------------------

      template< class FunctionSpace, class GridPart, int order, class Storage >
      struct AnisotropicDiscontinuousGalerkinSpaceTraits
      {
        using DiscreteFunctionSpaceType = hpDG::AnisotropicDiscontinuousGalerkinSpace< FunctionSpace, GridPart, order, Storage >;

        using FunctionSpaceType = FunctionSpace;

        using GridPartType = GridPart;

        using BasisFunctionSetsType = hpDG::AnisotropicBasisFunctionSets< FunctionSpaceType, GridPartType, order, Storage >;
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
       *  \tparam Storage  enable/disable Storage of quadratures
       *
       *  \ingroup DiscreteFunctionSpace_Implementation_Anisotropic
       */
      template< class FunctionSpace, class GridPart, int order, class Storage >
      class AnisotropicDiscontinuousGalerkinSpace
      : public hpDG::DiscontinuousGalerkinSpace< AnisotropicDiscontinuousGalerkinSpaceTraits< FunctionSpace, GridPart, order, Storage > >
      {
        using BaseType = hpDG::DiscontinuousGalerkinSpace< AnisotropicDiscontinuousGalerkinSpaceTraits< FunctionSpace, GridPart, order, Storage > >;

      public:
        using GridPartType = typename BaseType::GridPartType;
        using EntityType   = typename BaseType::EntityType;
        using BasisFunctionSetsType = typename BaseType::BasisFunctionSetsType;
        typedef typename BaseType::KeyType  KeyType;

        explicit AnisotropicDiscontinuousGalerkinSpace ( GridPartType &gridPart,
                                                         const Dune::InterfaceType interface = Dune::InteriorBorder_All_Interface,
                                                         const Dune::CommunicationDirection direction = Dune::ForwardCommunication )
          : BaseType( gridPart, BasisFunctionSetsType(), defaultKey(), interface, direction )
        {}

        explicit AnisotropicDiscontinuousGalerkinSpace ( GridPartType &gridPart,
                                                         const typename BaseType::KeyType key,
                                                         const Dune::InterfaceType interface = Dune::InteriorBorder_All_Interface,
                                                         const Dune::CommunicationDirection direction = Dune::ForwardCommunication )
          : BaseType( gridPart, BasisFunctionSetsType(), key, interface, direction )
        {}

        explicit AnisotropicDiscontinuousGalerkinSpace ( GridPartType &gridPart,
                                                         const std::vector<int>& key,
                                                         const Dune::InterfaceType interface = Dune::InteriorBorder_All_Interface,
                                                         const Dune::CommunicationDirection direction = Dune::ForwardCommunication )
          : BaseType( gridPart, BasisFunctionSetsType(), convert(key), interface, direction )
        {}

        template <class Function,
                  std::enable_if_t<
                    std::is_arithmetic<
                      decltype(Function(std::declval<const EntityType>()))>::value,int> i=0>
        AnisotropicDiscontinuousGalerkinSpace ( GridPartType &gridPart, Function function,
                                                const Dune::InterfaceType interface = Dune::InteriorBorder_All_Interface,
                                                const Dune::CommunicationDirection direction = Dune::ForwardCommunication )
          : BaseType( gridPart, BasisFunctionSetsType(), defaultKey(), function, interface, direction )
        {}

      private:
        KeyType convert( const std::vector<int>& v ) const
        {
          KeyType key;
          assert( key.size() == v.size() );
          for( unsigned int i=0; i<key.size(); ++i )
            key[ i ] = v[ i ];
          return key;
        }

        static KeyType defaultKey ()
        {
          KeyType key;
          std::fill( key.begin(), key.end(), order );
          return key;
        }
      };

    } // namespace hpDG



#ifndef DOXYGEN

    // DefaultLocalRestrictProlong
    // ---------------------------

    template< class FunctionSpace, class GridPart, int order, class Storage >
    class DefaultLocalRestrictProlong< hpDG::AnisotropicDiscontinuousGalerkinSpace< FunctionSpace, GridPart, order, Storage > >
    : public DiscontinuousGalerkinLocalRestrictProlong< hpDG::AnisotropicDiscontinuousGalerkinSpace< FunctionSpace, GridPart, order, Storage >, false >
    {
      using BaseType = DiscontinuousGalerkinLocalRestrictProlong< hpDG::AnisotropicDiscontinuousGalerkinSpace< FunctionSpace, GridPart, order, Storage >, false >;

    public:
      explicit DefaultLocalRestrictProlong ( const typename BaseType::DiscreteFunctionSpaceType &space )
        : BaseType( space )
      {}
    };


#endif //#ifndef DOXYGEN

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_HPDG_SPACE_DISCONTINUOUSGALERKIN_ANISOTROPIC_HH
