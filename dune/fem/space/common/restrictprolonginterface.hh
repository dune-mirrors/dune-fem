#ifndef DUNE_FEM_RESTRICTPROLONGINTERFACE_HH
#define DUNE_FEM_RESTRICTPROLONGINTERFACE_HH

//- Dune includes 
#include <dune/common/bartonnackmanifcheck.hh>
#include <dune/grid/common/capabilities.hh>

//- local includes 
#include <dune/fem/gridpart/emptyindexset.hh>
#include <dune/fem/function/localfunction/temporarylocalfunction.hh>
#include <dune/fem/misc/combineinterface.hh>
#include <dune/fem/space/common/localrestrictprolong.hh>

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



  /*! \brief Allow the combination of two restrict/prolong instances
   */
  template <class I1,class I2>
  class RestrictProlongPair
  : public RestrictProlongInterface< RestrictProlongTraits
      < RestrictProlongPair<I1,I2>,
        typename PairOfInterfaces< I1, I2 >::T1Type::DomainFieldType > >,
    public PairOfInterfaces<I1,I2>
  {
    typedef PairOfInterfaces<I1,I2> BaseType;

  public:  
    typedef typename BaseType::T1Type::DomainFieldType DomainFieldType;
    dune_static_assert( (Conversion< DomainFieldType, typename BaseType::T2Type::DomainFieldType >::sameType),
                        "DomainFieldType doesn't match." );

    RestrictProlongPair(I1 i1, I2 i2)
    : PairOfInterfaces<I1,I2>(i1,i2)
    {}
    
    //! \copydoc RestrictProlongInterface::setFatherChildWeight 
    void setFatherChildWeight (const DomainFieldType& val) const {
      this->first().setFatherChildWeight(val);
      this->second().setFatherChildWeight(val);    
    }

    //! \copydoc RestrictProlongInterface::restrictLocal 
    template <class EntityType>
    void restrictLocal ( EntityType &father, EntityType &son, 
                         bool initialize ) const {
      this->first().restrictLocal(father,son,initialize);
      this->second().restrictLocal(father,son,initialize);    
    }

    //! \copydoc RestrictProlongInterface::prolongLocal 
    template <class EntityType>
    void prolongLocal ( EntityType &father, EntityType &son, 
                        bool initialize ) const {
      this->first().prolongLocal(father,son,initialize);
      this->second().prolongLocal(father,son,initialize);    
    }
    
    //! \copydoc RestrictProlongInterface::addToList 
    template <class CommunicatorImp>
    void addToList(CommunicatorImp& comm) 
    {
      this->first().addToList(comm); 
      this->second().addToList(comm);    
    }
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
      assert( indexSet.persistent() );
      return (indexSet.index( father ) == indexSet.index( son ));
    }

  public:
    /** \copydoc RestrictProlongInterface::setFatherChildWeight */
    void setFatherChildWeight ( const DomainFieldType &weight ) const {}
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
    typedef typename DiscreteFunctionType::LocalFunctionType LocalFunctionType;
    typedef typename DiscreteFunctionType::GridPartType GridPartType;

    typedef ConstLocalFunction< DiscreteFunctionType > ConstLocalFunctionType;
    typedef DefaultLocalRestrictProlong< DiscreteFunctionSpaceType > LocalRestrictProlongType;

    explicit RestrictProlongDefault ( DiscreteFunctionType &discreteFunction )
    : discreteFunction_( discreteFunction ),
      constLf_( discreteFunction ),
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
      constLf_.init( son );
      LocalFunctionType lfFather = discreteFunction_.localFunction( father );

      localRP_.restrictLocal( lfFather, constLf_, geometryInFather, initialize );
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
      constLf_.init( father );
      LocalFunctionType lfSon = discreteFunction_.localFunction( son );

      localRP_.prolongLocal( constLf_, lfSon, geometryInFather, initialize );
    }

    //! add discrete function to communicator 
    template< class Communicator >
    void addToList ( Communicator &comm )
    {
      if( localRP_.needCommunication() )
        comm.addToList( discreteFunction_ );
    }

  protected:
    DiscreteFunctionType &discreteFunction_;
    mutable ConstLocalFunctionType constLf_;
    mutable LocalRestrictProlongType localRP_;
  };
  ///@} 

  } // namespace Fem  

  // #if DUNE_FEM_COMPATIBILITY  
  // put this in next version 1.4 

  using Fem :: RestrictProlongDefault ;
  // #endif // DUNE_FEM_COMPATIBILITY

} // namespace Dune 

#endif // #ifndef DUNE_FEM_RESTRICTPROLONGINTERFACE_HH
