#ifndef DUNE_FEM_RESTRICTPROLONGINTERFACE_HH
#define DUNE_FEM_RESTRICTPROLONGINTERFACE_HH

//- Dune includes
#include <dune/common/bartonnackmanifcheck.hh>
#include <dune/grid/common/capabilities.hh>

//- local includes
#include <dune/fem/function/localfunction/const.hh>
#include <dune/fem/function/localfunction/mutable.hh>
#include <dune/fem/space/common/localrestrictprolong.hh>
#include <dune/fem/gridpart/common/indexset.hh>

namespace Dune
{

  namespace Fem
  {

    /** @addtogroup RestrictProlongInterface

        Interface for restriction and prolongation operation of data
        on single elements.

        \remarks The Interface for a restriction and prolongation operation
        is defined by the class RestrictProlongInterface.


      @{
     */

    /*! @ingroup RestrictProlongInterface
        \brief Interface class defining the local behaviour of the
        restrict/prolong operation (using BN)

        \interfaceclass
     */
    template< class Traits >
    class RestrictProlongInterface
    {
      typedef RestrictProlongInterface< Traits > ThisType;

    public:
      //! \brief type of restrict-prolong operator implementation
      typedef typename Traits::RestProlImp RestProlImp;

      //! \brief field type of domain vector space
      typedef typename Traits::DomainFieldType DomainFieldType;

      /** \brief initialize restrict prolong object (if necessary) before adaptation takes place
       */
      void initialize ()
      {
        CHECK_AND_CALL_INTERFACE_IMPLEMENTATION( asImp().initialize() );
      }

      /** \brief finalize restrict prolong object (if necessary) after adaptation and dof compression was finished
       */
      void finalize ()
      {
        CHECK_AND_CALL_INTERFACE_IMPLEMENTATION( asImp().finalize() );
      }

      /** \brief explicit set volume ratio of son and father
       *
       *  \param[in]  weight  volume of son / volume of father
       *
       *  \note If this ratio is set, it is assume to be constant.
       */
      void setFatherChildWeight ( const DomainFieldType &weight ) const
      {
        CHECK_AND_CALL_INTERFACE_IMPLEMENTATION( asImp().setFatherChildWeight( weight ) );
      }

      //! restrict data to father
      template< class Entity >
      void restrictLocal ( const Entity &father, const Entity &son, bool initialize ) const
      {
        CHECK_AND_CALL_INTERFACE_IMPLEMENTATION( asImp().restrictLocal( father, son, initialize ) );
      }

      //! restrict data to father
      template< class Entity, class LocalGeometry >
      void restrictLocal ( const Entity &father, const Entity &son,
                           const LocalGeometry &geometryInFather,
                           bool initialize ) const
      {
        CHECK_AND_CALL_INTERFACE_IMPLEMENTATION( asImp().restrictLocal( father, son, geometryInFather, initialize ) );
      }

      //! finalize restriction on father
      template <class Entity>
      void restrictFinalize(const Entity &father) const
      {
        CHECK_AND_CALL_INTERFACE_IMPLEMENTATION( asImp().restrictFinalize( father ) );
      }

      //! prolong data to children
      template< class Entity >
      void prolongLocal ( const Entity &father, const Entity &son, bool initialize ) const
      {
        CHECK_AND_CALL_INTERFACE_IMPLEMENTATION( asImp().prolongLocal( father, son, initialize ) );
      }

      //! prolong data to children
      template< class Entity, class LocalGeometry >
      void prolongLocal ( const Entity &father, const Entity &son,
                          const LocalGeometry &geometryInFather,
                          bool initialize ) const
      {
        CHECK_AND_CALL_INTERFACE_IMPLEMENTATION( asImp().prolongLocal( father, son, geometryInFather, initialize ) );
      }

      /** \brief add discrete function to communicator
       *  \param[in]  comm  Communicator to add the discrete functions to
       */
      template< class Communicator >
      void addToList ( Communicator &comm )
      {
        CHECK_AND_CALL_INTERFACE_IMPLEMENTATION( asImp().addToList( comm ) );
      }

      /** \brief add discrete function to load balancer
       *  \param[in]  lb LoadBalancer to add the discrete functions to
       */
      template< class LoadBalancer >
      void addToLoadBalancer ( LoadBalancer &lb )
      {
        CHECK_AND_CALL_INTERFACE_IMPLEMENTATION( asImp().addToLoadBalancer( lb ) );
      }

    protected:
      /** \brief calculates the weight, i.e. (volume son)/(volume father)
          \param[in] father Father Entity
          \param[in] son Son Entity
          \return proportion between fahter and son volume
      */
      template< class Entity >
      DomainFieldType calcWeight ( const Entity &father, const Entity &son ) const
      {
        const DomainFieldType weight = son.geometry().volume() / father.geometry().volume();
        assert( weight > DomainFieldType( 0 ) );
        return weight;
      }

    protected:
      const RestProlImp &asImp () const { return static_cast< const RestProlImp & >( *this ); }
      RestProlImp &asImp () { return static_cast< RestProlImp & >( *this ); }
    };


    /** \brief Traits class for derivation from RestrictProlongInterface. */
    template< class Impl, class DomainField >
    struct RestrictProlongTraits
    {
      typedef Impl RestProlImp;
      typedef DomainField DomainFieldType;
    };



    /** \brief Interface default implementation for derived classes */
    template< class Traits >
    class RestrictProlongInterfaceDefault
    : public RestrictProlongInterface< Traits >
    {
      typedef RestrictProlongInterfaceDefault< Traits > ThisType;
      typedef RestrictProlongInterface< Traits > BaseType;

    public:
      typedef typename BaseType::DomainFieldType DomainFieldType;

    protected:
      //! return true if father and son have the same index
      template< class IndexSet, class Entity >
      bool entitiesAreCopies ( const IndexSet &indexSet,
                               const Entity &father, const Entity &son ) const
      {
        // adaptation and restriction/prolongation only works
        // for index set derived from Dune::Fem::IndexSet
        assert( Capabilities::isDuneFemIndexSet< IndexSet >::v );
        return (indexSet.index( father ) == indexSet.index( son ));
      }


    public:
      /** \copydoc RestrictProlongInterface::setFatherChildWeight(const DomainFieldType &weight) const*/
      void setFatherChildWeight ( const DomainFieldType &weight ) const {}

      void initialize () {}
      void finalize   () {}

      //! restrictFinalize is for some spaces where a local matrix inversion is carried out after restriction.
      template <class Entity>
      void restrictFinalize( const Entity& father ) const {}
    };



    /** \brief This is a wrapper for the default implemented
        restriction/prolongation operator, which only takes a discrete
        function template
     */
    template< class DiscreteFunction >
    class RestrictProlongDefault
    : public RestrictProlongInterfaceDefault< RestrictProlongTraits< RestrictProlongDefault< DiscreteFunction >, typename DiscreteFunction::DomainFieldType > >
    {
      typedef RestrictProlongDefault< DiscreteFunction > ThisType;
      typedef RestrictProlongInterfaceDefault< RestrictProlongTraits< ThisType, typename DiscreteFunction::DomainFieldType > > BaseType;

    public:
      typedef DiscreteFunction DiscreteFunctionType;

      typedef typename BaseType::DomainFieldType DomainFieldType;

      typedef typename DiscreteFunctionType::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;

      typedef ConstLocalFunction< DiscreteFunctionType > ConstLocalFunctionType;
      typedef MutableLocalFunction< DiscreteFunctionType > MutableLocalFunctionType;

      typedef typename DiscreteFunctionType::GridPartType GridPartType;

      typedef DefaultLocalRestrictProlong< DiscreteFunctionSpaceType > LocalRestrictProlongType;

      RestrictProlongDefault ( const RestrictProlongDefault& other )
        : discreteFunction_( const_cast< RestrictProlongDefault& > (other).discreteFunction_ ),
          constLf_( discreteFunction_ ),
          targetLf_( discreteFunction_ ),
          localRP_( discreteFunction_.space() )
      {
        // enable dof compression for this discrete function
        discreteFunction_.enableDofCompression();
      }

      explicit RestrictProlongDefault ( DiscreteFunctionType &discreteFunction )
      : discreteFunction_( discreteFunction ),
        constLf_( discreteFunction ),
        targetLf_( discreteFunction_ ),
        localRP_( discreteFunction_.space() )
      {
        // enable dof compression for this discrete function
        discreteFunction_.enableDofCompression();
      }

    protected:
      using BaseType::calcWeight;
      using BaseType::entitiesAreCopies;

    public:
      /** \brief explicit set volume ratio of son and father
       *
       *  \param[in]  weight  volume of son / volume of father
       *
       *  \note If this ratio is set, it is assume to be constant.
       */
      void setFatherChildWeight ( const DomainFieldType &weight ) const
      {
        localRP_.setFatherChildWeight( weight );
      }

      //! restrict data to father
      template< class Entity >
      void restrictLocal ( const Entity &father, const Entity &son, bool initialize ) const
      {
        assert( !father.isLeaf() );

        // convert from grid entities to grid part entities
        typedef typename GridPartType::template Codim< Entity::codimension >::EntityType GridPartEntityType;
        const GridPartType &gridPart = discreteFunction_.gridPart();
        const GridPartEntityType &gpFather = gridPart.convert( father );
        const GridPartEntityType &gpSon    = gridPart.convert( son );

        if( !entitiesAreCopies( gridPart.indexSet(), gpFather, gpSon ) )
          restrictLocal( gpFather, gpSon, son.geometryInFather(), initialize );
      }

      //! restrict data to father
      template< class Entity, class LocalGeometry >
      void restrictLocal ( const Entity &father, const Entity &son,
                           const LocalGeometry &geometryInFather,
                           bool initialize ) const
      {
        auto clf = bindGuard( constLf_, son );

        if( initialize || ! targetLf_.valid() )
        {
          targetLf_.bind( father );
          // unbind is done in restrictFinalize
        }

        assert( father == targetLf_.entity() );
        localRP_.restrictLocal( targetLf_, constLf_, geometryInFather, initialize );
      }

      //! finalize restriction on father
      template <class Entity>
      void restrictFinalize(const Entity &father) const
      {
        // make sure targetLf_ still is bound to father
        assert( discreteFunction_.gridPart().convert( father ) == targetLf_.entity() );

        localRP_.restrictFinalize( targetLf_ );

        // release target local function
        targetLf_.unbind();
      }


      //! prolong data to children
      template< class Entity >
      void prolongLocal ( const Entity &father, const Entity &son, bool initialize ) const
      {
        assert( !father.isLeaf() );

        // convert from grid entities to grid part entities
        typedef typename GridPartType::template Codim< Entity::codimension >::EntityType GridPartEntityType;
        const GridPartType &gridPart = discreteFunction_.gridPart();
        const GridPartEntityType &gpFather = gridPart.convert( father );
        const GridPartEntityType &gpSon    = gridPart.convert( son );

        if( !entitiesAreCopies( gridPart.indexSet(), gpFather, gpSon ) )
          prolongLocal( gpFather, gpSon, son.geometryInFather(), initialize );
      }

      //! prolong data to children
      template< class Entity, class LocalGeometry >
      void prolongLocal ( const Entity &father, const Entity &son,
                          const LocalGeometry &geometryInFather,
                          bool initialize ) const
      {
        // initialize father local function (only once)
        //if( initialize )
        {
          constLf_.unbind();
          constLf_.bind( father );
        }

        auto tlf = bindGuard( targetLf_, son );

        localRP_.prolongLocal( constLf_, targetLf_, geometryInFather, initialize );
        // write to discrete function
        discreteFunction_.setLocalDofs( son, targetLf_ );
      }

      //! add discrete function to communicator with given unpack operation
      template< class Communicator, class Operation >
      void addToList ( Communicator &comm, const Operation& op)
      {
        if( localRP_.needCommunication() )
          comm.addToList( discreteFunction_, op );
      }

      //! add discrete function to communicator
      template< class Communicator >
      void addToList ( Communicator &comm  )
      {
        if( localRP_.needCommunication() )
          comm.addToList( discreteFunction_ );
      }

      //! remove discrete function from communicator
      template< class Communicator >
      void removeFromList ( Communicator &comm )
      {
        if( localRP_.needCommunication() )
          comm.removeFromList( discreteFunction_ );
      }

      //! add discrete function to load balancer
      template< class LoadBalancer >
      void addToLoadBalancer ( LoadBalancer& lb )
      {
        lb.addToLoadBalancer( discreteFunction_ );
      }

    protected:
      DiscreteFunctionType &discreteFunction_;
      mutable ConstLocalFunctionType   constLf_;
      mutable MutableLocalFunctionType targetLf_;
      mutable LocalRestrictProlongType localRP_;
    };
    ///@}

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_RESTRICTPROLONGINTERFACE_HH
