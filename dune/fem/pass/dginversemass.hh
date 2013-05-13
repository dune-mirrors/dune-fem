#ifndef DUNE_FEM_PASS_DGINVERSEMASS_HH
#define DUNE_FEM_PASS_DGINVERSEMASS_HH

#include <cassert>
#include <iosfwd>
#include <sstream>

#include <dune/fem/common/tupleutility.hh>
#include <dune/fem/operator/1order/localmassmatrix.hh>
#include <dune/fem/pass/common/pass.hh>
#include <dune/fem/quadrature/cachingquadrature.hh>

namespace Dune
{

  namespace Fem
  {

    // DGInverseMassPassDiscreteModel
    // ------------------------------

    template< int functionalId, class PreviousPass >
    struct DGInverseMassPassDiscreteModel
    {
    private:
      typedef typename PreviousPass::PassIds PassIds;
      typedef typename PreviousPass::NextArgumentType LocalArgumentType;
      typedef typename PreviousPass::GlobalArgumentType GlobalArgumentType;
      typedef typename PushFrontTuple< LocalArgumentType, GlobalArgumentType * >::type TotalArgumentType;

    public:
      static const std::size_t functionalPosition
        = Dune::FirstTypeIndex< PassIds, Dune::integral_constant< int, functionalId > >::type::value;

      struct Traits
      {
        typedef typename Dune::tuple_element< functionalPosition, TotalArgumentType >::type DestinationType;
        typedef typename DestinationType::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
      };
    };



    // DGInverseMassPass
    // -----------------

    /** \brief Pass applying the local inverse mass matrix on each element
     *
     *  \tparam  functionalId   pass id of functional to convert
     *  \tparam  PreviousPass   type of previous pass
     *  \tparam  id             pass id
     */
    template< int functionalId, class PreviousPass, int id >
    class DGInverseMassPass
    : public Dune::Fem::Pass< DGInverseMassPassDiscreteModel< functionalId, PreviousPass >, PreviousPass, id >
    {
      typedef DGInverseMassPass< functionalId, PreviousPass, id > ThisType;
      typedef Dune::Fem::Pass< DGInverseMassPassDiscreteModel< functionalId, PreviousPass >, PreviousPass, id > BaseType;

    public:
      //! pass ids up to here (tuple of integral constants)
      typedef typename BaseType::PassIds PassIds;

      //! \brief argument type
      typedef typename BaseType::TotalArgumentType TotalArgumentType;
      //! \brief destination type
      typedef typename BaseType::DestinationType DestinationType;

      //! \brief discrete function space type
      typedef typename BaseType::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;

      //! \brief entity type
      typedef typename BaseType::EntityType EntityType;

    private:
      static const std::size_t functionalPosition
        = DGInverseMassPassDiscreteModel< functionalId, PreviousPass >::functionalPosition;

      typedef typename DiscreteFunctionSpaceType::GridPartType GridPartType;
      typedef LocalMassMatrix< DiscreteFunctionSpaceType, CachingQuadrature< GridPartType, 0 > > LocalMassMatrixType;

    public:
      using BaseType::passNumber;

      explicit DGInverseMassPass ( PreviousPass &previousPass,
                                   const DiscreteFunctionSpaceType &space )
      : BaseType( previousPass ),
        space_( space ),
        localMassMatrix_( space, 2*space.order() )
      {}

      void printTexInfo ( std::ostream &out ) const
      {
        BaseType::printTexInfo( out );
        out << "DGInverseMassPass: " << "\\\\" << std::endl;
      }

      void allocateLocalMemory ()
      {
        if( !destination_ )
        {
          std::ostringstream name;
          name << "DGInverseMassPass_" << passNumber();
          destination_ = new DestinationType( name.str(), space() );
          deleteHandler_ = &(BaseType::DeleteHandlerType::instance());
        }
      }

      const DiscreteFunctionSpaceType &space () const { return space_; }

    protected:
      void compute ( const TotalArgumentType &argument, DestinationType &destination ) const
      {
        typedef typename GridPartType::template Codim< 0 >::template Partition< All_Partition >::IteratorType IteratorType;
        typedef typename GridPartType::template Codim< 0 >::EntityType EntityType;

        const GridPartType &gridPart = space().gridPart();
        const IteratorType end = gridPart.template end< 0, All_Partition >();
        for( IteratorType it = gridPart.template begin< 0, All_Partition >(); it != end; ++it )
        {
          const EntityType &entity = *it;
          typename DestinationType::LocalFunctionType localDestination = destination.localFunction( entity );
          localDestination.assign( Dune::get< functionalPosition >( argument )->localFunction( entity ) );
          localMassMatrix_.applyInverse( entity, localDestination );
        }
      }

    protected:
      using BaseType::destination_;
      using BaseType::deleteHandler_;

    private:
      const DiscreteFunctionSpaceType &space_;
      LocalMassMatrixType localMassMatrix_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_PASS_APPLYLOCALOPERATOR_HH
