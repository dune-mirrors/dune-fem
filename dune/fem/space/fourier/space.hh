#ifndef DUNE_FEM_SPACE_FOURIER_SPACE_HH
#define DUNE_FEM_SPACE_FOURIER_SPACE_HH

#include <cassert>
#include <limits>

#include <dune/fem/common/hybrid.hh>
#include <dune/fem/function/common/functionset.hh>
#include <dune/fem/function/localfunction/localfunctionsetadapter.hh>
#include <dune/fem/space/basisfunctionset/simple.hh>
#include <dune/fem/space/basisfunctionset/vectorial.hh>
#include <dune/fem/space/common/defaultcommhandler.hh>
#include <dune/fem/space/common/discretefunctionspace.hh>
#include <dune/fem/space/mapper/codimensionmapper.hh>
#include <dune/fem/space/mapper/nonblockmapper.hh>

#include <dune/fem/space/fourier/capabilities.hh>
#include <dune/fem/space/fourier/declaration.hh>
#include <dune/fem/space/fourier/dofmapper.hh>
#include <dune/fem/space/fourier/functionset.hh>
#include <dune/fem/space/fourier/interpolate.hh>

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

      typedef GridPart GridPartType;
      typedef GridFunctionSpace< GridPartType, FunctionSpace > FunctionSpaceType;

      static const int codimension = 0;
      typedef typename GridPartType::template Codim< codimension >::EntityType EntityType;

      typedef typename FunctionSpaceType::ScalarFunctionSpaceType ScalarFunctionSpaceType;
      typedef FourierFunctionSet< ScalarFunctionSpaceType, order > FunctionSetType;
      typedef FunctionSetProxy< FunctionSetType > FunctionSetProxyType;
      typedef LocalFunctionSetAdapter< EntityType, FunctionSetProxyType > LocalFunctionSetType;
      typedef SimpleBasisFunctionSet< LocalFunctionSetType > ScalarBasisFunctionSetType;

      typedef VectorialBasisFunctionSet< ScalarBasisFunctionSetType, typename FunctionSpaceType::RangeType > BasisFunctionSetType;

      typedef Hybrid::IndexRange< int, FunctionSpaceType::dimRange * FourierFunctionSetSize< FunctionSpaceType::dimDomain, order >::v > LocalBlockIndices;

      typedef FourierDofMapper< GridPartType, order > BlockMapperType;

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

      typedef typename Traits::FunctionSetType FunctionSetType;

      typedef typename Traits::ScalarBasisFunctionSetType ScalarBasisFunctionSetType;
      typedef typename BaseType::BasisFunctionSetType BasisFunctionSetType;

      typedef typename BaseType::BlockMapperType BlockMapperType;

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
        functionSet_( order )
      {}

      /** @copydoc Dune::Fem::DiscreteFunctionSpaceInterface::type */
      DFSpaceIdentifier type () const { return FourierSpace_id; }

      /** @copydoc Dune::Fem::DiscreteFunctionSpaceInterface::basisFunctionSet */
      BasisFunctionSetType basisFunctionSet ( const EntityType &entity ) const
      {
        typedef typename Traits::LocalFunctionSetType LocalFunctionSetType;
        typedef typename Traits::ScalarBasisFunctionSetType ScalarBasisFunctionSetType;
        return BasisFunctionSetType( ScalarBasisFunctionSetType( LocalFunctionSetType( entity, &functionSet_ ) ) );
      }

      /** @copydoc Dune::Fem::DiscreteFunctionSpaceInterface::continuous */
      bool continuous () const { return true; }

      /** @copydoc Dune::Fem::DiscreteFunctionSpaceInterface::continuous */
      bool continuous ( const IntersectionType &intersection ) const { return true; }

      /** @copydoc Dune::Fem::DiscreteFunctionSpaceInterface::order */
      int order () const { return functionSet_.order(); }

      /** @copydoc Dune::Fem::DiscreteFunctionSpaceInterface::blockMapper */
      BlockMapperType &blockMapper () const { return blockMapper_; }

      const FunctionSetType &functionSet () const { return functionSet_; }

    private:
      mutable BlockMapperType blockMapper_;
      FunctionSetType functionSet_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SPACE_FOURIER_SPACE_HH
