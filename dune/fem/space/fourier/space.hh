#ifndef DUNE_FEM_SPACE_FOURIER_SPACE_HH
#define DUNE_FEM_SPACE_FOURIER_SPACE_HH

#include <cassert>
#include <limits>

#include <dune/fem/space/basisfunctionset/proxy.hh>
#include <dune/fem/space/common/defaultcommhandler.hh>
#include <dune/fem/space/common/discretefunctionspace.hh>
#include <dune/fem/space/mapper/codimensionmapper.hh>
#include <dune/fem/space/mapper/nonblockmapper.hh>
#include <dune/fem/version.hh>

#include "basisfunctionset.hh"
#include "basisfunctions.hh"
#include "declaration.hh"
#include "dofmapper.hh"


namespace Dune
{

  namespace Fem
  {

    // FourierDiscreteFunctionSpaceTraits
    // ----------------------------------

    template< class FunctionSpace, class GridPart, int order >
    struct FourierDiscreteFunctionSpaceTraits
    {
      typedef FourierDiscreteFunctionSpace< FunctionSpace, GridPart, order > DiscreteFunctionSpaceType; 

      typedef FunctionSpace FunctionSpaceType;
      typedef GridPart GridPartType;

      static const int codimension = 0;
      typedef typename GridPartType::template Codim< codimension >::EntityType EntityType;

      typedef typename FunctionSpaceType::ScalarFunctionSpaceType ScalarFunctionSpaceType;
      typedef FourierBasisFunctions< ScalarFunctionSpaceType, order > ScalarBasisFunctionsType;
      typedef FourierBasisFunctionSet< EntityType, ScalarBasisFunctionsType > ScalarBasisFunctionSetType;

      typedef ScalarBasisFunctionSetType BasisFunctionSetType;

      static const int localBlockSize = FunctionSpace::dimRange * NumFourierBasisFunctions< FunctionSpaceType::dimDomain, order >::v;

      typedef FourierDofMapper< GridPartType, order > BlockMapperType;
      typedef NonBlockMapper< BlockMapperType, localBlockSize > MapperType;

      template< class DiscreteFunction, class Operation = Dune::Fem::DFCommunicationOperation::Add >
      struct CommDataHandle
      {
        typedef Dune::Fem::DefaultCommunicationHandler< DiscreteFunction, Operation > Type;
        typedef Operation OperationType;
      };
    };



    // FourierDiscreteFunctionSpace
    // ----------------------------

    template< class FunctionSpace, class GridPart, int Order >
    class FourierDiscreteFunctionSpace
    : public DiscreteFunctionSpaceDefault< FourierDiscreteFunctionSpaceTraits< FunctionSpace, GridPart, Order > >
    {
      typedef FourierDiscreteFunctionSpaceTraits< FunctionSpace, GridPart, Order > ThisType;
      typedef DiscreteFunctionSpaceDefault< FourierDiscreteFunctionSpaceTraits< FunctionSpace, GridPart, Order > > BaseType;

    public:
      typedef typename BaseType::Traits Traits;

      static const int polynomialOrder = Order+1;

      typedef typename BaseType::FunctionSpaceType FunctionSpaceType;
      typedef typename BaseType::GridPartType GridPartType;
      typedef typename BaseType::EntityType EntityType;
      typedef typename BaseType::IntersectionType IntersectionType;

      typedef typename Traits::ScalarBasisFunctionsType ScalarBasisFunctionsType;

      typedef typename Traits::ScalarBasisFunctionSetType ScalarBasisFunctionSetType;
      typedef typename BaseType::BasisFunctionSetType BasisFunctionSetType;

      typedef typename BaseType::BlockMapperType BlockMapperType;
      typedef typename BaseType::MapperType MapperType;

    protected:
      static const InterfaceType defaultInterface = InteriorBorder_All_Interface;
      static const CommunicationDirection defaultDirection =  ForwardCommunication;

    public:
      using BaseType::order;

      explicit FourierDiscreteFunctionSpace ( GridPartType &gridPart,
                                              int order = std::numeric_limits< int >::max(),
                                              const InterfaceType commInterface = defaultInterface,
                                              const CommunicationDirection commDirection = defaultDirection )
      : BaseType( gridPart, commInterface, commDirection ),
        mapper_( blockMapper_ ),
        scalarBasisFunctions_( order )
      {}

      /** @copydoc Dune::Fem::DiscreteFunctionSpaceInterface::type */
      DFSpaceIdentifier type () const { return FourierSpace_id; }

      /** @copydoc Dune::Fem::DiscreteFunctionSpaceInterface::basisFunctionSet */
      BasisFunctionSetType basisFunctionSet ( const EntityType &entity ) const
      {
        return BasisFunctionSetType( entity, scalarBasisFunctions_ );
      }

      /** @copydoc Dune::Fem::DiscreteFunctionSpaceInterface::continuous */
      bool continuous () const { return true; }

      /** @copydoc Dune::Fem::DiscreteFunctionSpaceInterface::continuous */
      bool continuous ( const IntersectionType &intersection ) const { return true; }

      /** @copydoc Dune::Fem::DiscreteFunctionSpaceInterface::order */
      int order () const { return scalarBasisFunctions_.order(); }

      /** @copydoc Dune::Fem::DiscreteFunctionSpaceInterface::blockMapper */
      BlockMapperType &blockMapper () const { return blockMapper_; }

      /** @copydoc Dune::Fem::DiscreteFunctionSpaceInterface::mapper */
      DUNE_VERSION_DEPRECATED(1,4,remove)
      MapperType &mapper () const { return mapper_; }

    private:
      mutable BlockMapperType blockMapper_;
      mutable MapperType mapper_;
      ScalarBasisFunctionsType scalarBasisFunctions_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SPACE_FOURIER_SPACE_HH
