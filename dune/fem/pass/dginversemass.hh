#ifndef DUNE_FEM_PASS_DGINVERSEMASS_HH
#define DUNE_FEM_PASS_DGINVERSEMASS_HH

#if HAVE_DUNE_FEM_DG
#error "Outdated header, #include <dune/fem-dg/pass/dginversemass.hh> instead!"
#endif

#include <cassert>
#include <iosfwd>
#include <sstream>
#include <type_traits>

#include <dune/fem/common/tupleutility.hh>
#include <dune/fem/operator/1order/localmassmatrix.hh>
#include <dune/fem/pass/common/local.hh>
#include <dune/fem/quadrature/cachingquadrature.hh>
#include <dune/fem/pass/localdg/discretemodel.hh>

namespace Dune
{

  namespace Fem
  {

    // DGInverseMassPassDiscreteModel
    // ------------------------------

    template< int functionalId, class PreviousPass >
    struct DGInverseMassPassDiscreteModel : public DGAdaptiveDiscreteModel
    {
    private:
      typedef typename PreviousPass::PassIds PassIds;
      typedef typename PreviousPass::NextArgumentType LocalArgumentType;
      typedef typename PreviousPass::GlobalArgumentType GlobalArgumentType;
      typedef typename PushFrontTuple< LocalArgumentType, GlobalArgumentType * >::type TotalArgumentType;

    public:
      static const std::size_t functionalPosition
        = Dune::FirstTypeIndex< PassIds, std::integral_constant< int, functionalId > >::type::value;
      typedef typename std::tuple_element< functionalPosition, TotalArgumentType >::type DestinationPtrType;

      struct Traits
      {
        typedef typename std::remove_pointer< DestinationPtrType >::type DestinationType;
        typedef typename DestinationType::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
        typedef CachingQuadrature< typename DiscreteFunctionSpaceType::GridPartType, 0 > VolumeQuadratureType;
        typedef CachingQuadrature< typename DiscreteFunctionSpaceType::GridPartType, 1 > FaceQuadratureType;
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
    : public Dune::Fem::LocalPass< DGInverseMassPassDiscreteModel< functionalId, PreviousPass >, PreviousPass, id >
    {
      typedef DGInverseMassPass< functionalId, PreviousPass, id > ThisType;
      typedef Dune::Fem::LocalPass< DGInverseMassPassDiscreteModel< functionalId, PreviousPass >, PreviousPass, id > BaseType;

    public:
      //! type of the discrete model used
      typedef DGInverseMassPassDiscreteModel< functionalId, PreviousPass > DiscreteModelType;

      //! pass ids up to here (tuple of integral constants)
      typedef typename BaseType::PassIds PassIds;

      //! \brief argument type
      typedef typename BaseType::TotalArgumentType TotalArgumentType;
      //! \brief destination type
      typedef typename BaseType::DestinationType DestinationType;

      //! \brief discrete function space type
      typedef typename DiscreteModelType::Traits::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;

    private:
      static const std::size_t functionalPosition = DiscreteModelType::functionalPosition;

      typedef typename DiscreteFunctionSpaceType::GridPartType GridPartType;
      typedef typename DiscreteModelType::Traits::VolumeQuadratureType VolumeQuadratureType;
      typedef LocalMassMatrix< DiscreteFunctionSpaceType, VolumeQuadratureType > LocalMassMatrixType;

      typedef typename GridPartType::template Codim< 0 >::EntityType EntityType;
    public:
      using BaseType::passNumber;
      using BaseType::space;

      explicit DGInverseMassPass ( PreviousPass &previousPass,
                                   const DiscreteFunctionSpaceType &space )
      : BaseType( previousPass, space, "DGInverseMassPass" ),
        localMassMatrix_( space, 2*space.order() ),
        arg_( 0 ),
        dest_( 0 )
      {}

      //! constructor for use with thread pass
      explicit DGInverseMassPass ( const DiscreteModelType& discreteModel,
                                   PreviousPass &previousPass,
                                   const DiscreteFunctionSpaceType &space,
                                   const int volQuadOrd  = -1,
                                   const int faceQuadOrd = -1)
      : BaseType( previousPass, space, "DGInverseMassPass" ),
        localMassMatrix_( space, ( volQuadOrd == -1 ) ? (2*space.order()) : volQuadOrd ),
        arg_( 0 ),
        dest_( 0 )
      {}

      void printTexInfo ( std::ostream &out ) const
      {
        BaseType::printTexInfo( out );
        out << "DGInverseMassPass: " << "\\\\" << std::endl;
      }

      //! this pass needs no communication
      bool requireCommunication () const { return false ; }

      //! interface method
      void prepare( const TotalArgumentType &argument, DestinationType &destination ) const
      {
        arg_  = &argument;
        dest_ = &destination;
      }

      //! prepare for ThreadPass
      void prepare( const TotalArgumentType &argument, DestinationType &destination, const bool ) const
      {
        prepare( argument, destination );
      }

      //! finalize for ThreadPass
      void finalize( const TotalArgumentType &argument, DestinationType &destination, const bool ) const
      {
        finalize( argument, destination );
      }

      //! interface method
      void finalize( const TotalArgumentType &argument, DestinationType &destination ) const
      {
        arg_  = 0 ;
        dest_ = 0 ;
      }

      //! apply inverse mass matrix locally
      void applyLocal( const EntityType& entity ) const
      {
        assert( arg_ );
        assert( dest_ );
        typename DestinationType::LocalFunctionType localDestination = dest_->localFunction( entity );
        localDestination.assign( std::get< functionalPosition >( *arg_ )->localFunction( entity ) );
        localMassMatrix_.applyInverse( entity, localDestination );
      }

      //! apply local with neighbor checker (not needed here)
      template <class NBChecker>
      void applyLocal( const EntityType& entity, const NBChecker& ) const
      {
        applyLocal( entity );
      }

      /** \brief  apply local for all elements that do not need information from
       *          other processes (here all elements)
       */
      template <class NBChecker>
      void applyLocalInterior( const EntityType& entity, const NBChecker& ) const
      {
        applyLocal( entity );
      }

      /** \brief  apply local for all elements that need information from
       *          other processes (here no elements)
       */
      template <class NBChecker>
      void applyLocalProcessBoundary( const EntityType& entity, const NBChecker& ) const
      {
        DUNE_THROW(InvalidStateException,"DGInverseMassPass does not need a second phase for ThreadPass");
      }

    protected:
      void compute ( const TotalArgumentType &argument, DestinationType &destination ) const
      {
        // set pointer
        prepare( argument, destination );

        typedef typename GridPartType::template Codim< 0 >::template Partition< All_Partition >::IteratorType IteratorType;
        const GridPartType &gridPart = space().gridPart();
        const IteratorType end = gridPart.template end< 0, All_Partition >();
        for( IteratorType it = gridPart.template begin< 0, All_Partition >(); it != end; ++it )
        {
          applyLocal( *it );
        }

        // remove pointer
        finalize( argument, destination );
      }

    protected:
      using BaseType::destination_;
      using BaseType::deleteHandler_;

    private:
      LocalMassMatrixType localMassMatrix_;

      mutable const TotalArgumentType* arg_;
      mutable DestinationType*         dest_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_PASS_APPLYLOCALOPERATOR_HH
