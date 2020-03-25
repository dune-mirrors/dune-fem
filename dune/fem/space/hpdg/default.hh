#ifndef DUNE_FEM_HPDG_SPACE_DISCONTINUOUSGALERKIN_DEFAULT_HH
#define DUNE_FEM_HPDG_SPACE_DISCONTINUOUSGALERKIN_DEFAULT_HH

#include <memory>

#include <dune/grid/common/gridenums.hh>

#include <dune/fem/operator/matrix/istlmatrixadapter.hh>
#include <dune/fem/space/common/commoperations.hh>
#include <dune/fem/space/common/defaultcommhandler.hh>
#include <dune/fem/space/common/localrestrictprolong.hh>
#include <dune/fem/space/discontinuousgalerkin/localrestrictprolong.hh>

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

      template< class BasisFunctionSets >
      class DefaultDiscontinuousGalerkinSpace;



#ifndef DOXYGEN

      // DefaultDiscontinuousGalerkinSpaceTraits
      // ---------------------------------------

      template< class BasisFunctionSets >
      struct DefaultDiscontinuousGalerkinSpaceTraits
      {
        using DiscreteFunctionSpaceType = hpDG::DefaultDiscontinuousGalerkinSpace< BasisFunctionSets >;

        using BasisFunctionSetsType = BasisFunctionSets;
        using BasisFunctionSetType = typename BasisFunctionSetsType::BasisFunctionSetType;

        using FunctionSpaceType = typename BasisFunctionSetType::FunctionSpaceType;
        using GridPartType = typename BasisFunctionSetsType::GridPartType;

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



      // DefaultDiscontinuousGalerkinSpace
      // ---------------------------------

      /** \brief Default implementation of an \f$hp\f$-adaptive discrete function space given a family of local basis function sets
       *
       *  \tparam BasisFunctionSets  a Dune::Fem::hpDG::BasisFunctionSets
       *
       *  \ingroup DiscreteFunctionSpace_Implementation_Default
       */
      template< class BasisFunctionSets >
      class DefaultDiscontinuousGalerkinSpace
        : public hpDG::DiscontinuousGalerkinSpace< DefaultDiscontinuousGalerkinSpaceTraits< BasisFunctionSets > >
      {
        using BaseType = hpDG::DiscontinuousGalerkinSpace< DefaultDiscontinuousGalerkinSpaceTraits< BasisFunctionSets > >;

      public:
        /** \copydoc Dune::Fem::hpDG::DiscontinuousGalerkinSpace::GridPartType */
        using GridPartType = typename BaseType::GridPartType;
        /** \copydoc Dune::Fem::hpDG::DiscontinuousGalerkinSpace::BasisFunctionSetsType */
        using BasisFunctionSetsType = typename BaseType::BasisFunctionSetsType;
        /** \copydoc Dune::Fem::hpDG::DiscontinuousGalerkinSpace::KeyType */
        using KeyType = typename BaseType::KeyType;

        /** \name Construction
         *  \{
         */

        explicit DefaultDiscontinuousGalerkinSpace ( GridPartType &gridPart,
                                                     const BasisFunctionSetsType &basisFunctionSets,
                                                     const KeyType &key,
                                                     const InterfaceType interface = InteriorBorder_All_Interface,
                                                     const CommunicationDirection direction = ForwardCommunication )
          : BaseType( gridPart, basisFunctionSets, key, interface, direction )
        {}

        /** \} */
      };



      // make_space
      // ----------

      /** \brief returns a new space instance for a given family of local basis function sets
       *
       *  \tparam BasisFunctionSets  a Dune::Fem::hpDG::BasisFunctionSets
       *
       *  \param[in] gridPart  a Dune::Fem::GridPart instance
       *  \param[in] basisFunctionSets  a Dune::Fem::hpDG::BasisFunctionSets instance
       *  \param[in] key  a default key
       *  \param[in] interface  interface type
       *  \param[in] direction  communication direction
       *
       *  \returns a std::unique_ptr with a new space instance
       *
       *  \ingroup DiscreteFunctionSpace_Implementation_Default
       */
      template< class BasisFunctionSets >
      std::unique_ptr< DefaultDiscontinuousGalerkinSpace< BasisFunctionSets > >
      make_space ( typename BasisFunctionSets::GridPartType &gridPart,
                   const BasisFunctionSets &basisFunctionSets,
                   const typename BasisFunctionSets::KeyType &key,
                   const InterfaceType interface = InteriorBorder_All_Interface,
                   const CommunicationDirection direction = ForwardCommunication )
      {
        using DiscreteFunctionSpaceType = DefaultDiscontinuousGalerkinSpace< BasisFunctionSets >;
        return std::unique_ptr< DiscreteFunctionSpaceType >( new DiscreteFunctionSpaceType( gridPart, basisFunctionSets, key, interface, direction ) );
      }

    } // namespace hpDG



#ifndef DOXYGEN

    // DefaultLocalRestrictProlong
    // ---------------------------

    template< class BasisFunctionSets >
    class DefaultLocalRestrictProlong< hpDG::DefaultDiscontinuousGalerkinSpace< BasisFunctionSets > >
      : public DiscontinuousGalerkinLocalRestrictProlong< hpDG::DefaultDiscontinuousGalerkinSpace< BasisFunctionSets >, !BasisFunctionSets::orthogonal() >
    {
      using BaseType = DiscontinuousGalerkinLocalRestrictProlong< hpDG::DefaultDiscontinuousGalerkinSpace< BasisFunctionSets >, !BasisFunctionSets::orthogonal() >;

    public:
      explicit DefaultLocalRestrictProlong ( const typename BaseType::DiscreteFunctionSpaceType &space )
        : BaseType( space )
      {}
    };



#endif // #ifndef DOXYGEN

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_HPDG_SPACE_DISCONTINUOUSGALERKIN_DEFAULT_HH
