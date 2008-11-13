#ifndef DUNE_RESTRICTPROLONGINTERFACE_HH
#define DUNE_RESTRICTPROLONGINTERFACE_HH

//- Dune includes 
#include <dune/common/bartonnackmanifcheck.hh>
#include <dune/grid/common/capabilities.hh>

//- local includes 
#include <dune/fem/misc/combineinterface.hh>
#include <dune/fem/gridpart/emptyindexset.hh>

namespace Dune{
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
template <class RestProlImpTraits>
class RestrictProlongInterface {
public:  
  //! \brief type of restrict-prolong operator implementation 
  typedef typename RestProlImpTraits::RestProlImp RestProlImp;
  //! \brief field type of range vector space 
  typedef typename RestProlImpTraits::RangeFieldType RangeFieldType;
  
  /** \brief if weight is set, then its assumend 
      that we have always the same proportion between fahter and son volume 
      \param[in] weight proportion between fahter and son volume 
  */
  void setFatherChildWeight (const RangeFieldType& weight) const 
  {
    CHECK_AND_CALL_INTERFACE_IMPLEMENTATION(
        asImp().setFatherChildWeight(weight));
  }
  
  /** \brief restrict data to father 
      \param[in] father Father Entity 
      \param[in] son Son Entity 
      \param[in] initialize <b>true</b> if restrictLocal is called for first time for this father
  */
  template <class EntityType>
  void restrictLocal ( EntityType &father, 
                       EntityType &son, 
            		       bool initialize ) const 
  {
    CHECK_AND_CALL_INTERFACE_IMPLEMENTATION(
        asImp().restrictLocal(father,son,initialize));
  }
  
  /** \brief prolong data to children 
      \param[in] father Father Entity 
      \param[in] son Son Entity 
      \param[in] initialize <b>true</b> if restrictLocal is called for first time for this father
  */
  template <class EntityType>
  void prolongLocal ( EntityType &father, 
                      EntityType &son, 
            		      bool initialize ) const 
  {
    CHECK_AND_CALL_INTERFACE_IMPLEMENTATION(
        asImp().prolongLocal(father,son,initialize));
  }

  /** \brief add discrete function to communicator 
      \param[in] communicator Communicator to which internal discrete functions are added to 
  */
  template <class CommunicationManagerImp>
  void addToList(CommunicationManagerImp& communicator)
  {
    CHECK_AND_CALL_INTERFACE_IMPLEMENTATION(
        asImp().addToList(communicator));  
  }
  
protected:  
  /** \brief calculates the weight, i.e. (volume son)/(volume father)
      \param[in] father Father Entity 
      \param[in] son Son Entity 
      \return proportion between fahter and son volume
  */
  template< class EntityType >
  RangeFieldType calcWeight ( const EntityType &father, const EntityType &son ) const
  {
    const RangeFieldType weight = son.geometry().volume() / father.geometry().volume();
    assert( weight > RangeFieldType( 0 ) );
    return weight;
  }
  
private:
  RestProlImp& asImp() {
    return static_cast<RestProlImp&>(*this);
  }

  const RestProlImp& asImp() const {
    return static_cast<const RestProlImp&>(*this);
  }
};

/** \brief Traits class for derivation from RestrictProlongInterface. */
template <class Base>
struct RestrictProlongTraits {
  typedef Base RestProlImp;
  typedef double RangeFieldType;
};

/*! \brief Allow the combination of two restrict/prolong instances
 */
template <class I1,class I2>
class RestrictProlongPair : 
    public RestrictProlongInterface<RestrictProlongTraits<RestrictProlongPair<I1,I2> > >,
    public PairOfInterfaces<I1,I2> {
public:  
  typedef PairOfInterfaces<I1,I2> BaseType;
  typedef typename BaseType::T1Type::RangeFieldType RangeFieldType;
  typedef typename BaseType::T2Type::RangeFieldType RangeFieldType2;

  RestrictProlongPair(I1 i1, I2 i2) : PairOfInterfaces<I1,I2>(i1,i2) {
    //IsTrue<SameType<RangeFieldType,RangeFieldType2>::value>::yes();
  }
  
  //! if weight is set, then ists assumend that we have always the same
  //! proportion between fahter and son volume 
  void setFatherChildWeight (const RangeFieldType& val) const {
    this->first().setFatherChildWeight(val);
    this->second().setFatherChildWeight(val);    
  }
  //! restrict data to father 
  template <class EntityType>
  void restrictLocal ( EntityType &father, EntityType &son, 
		       bool initialize ) const {
    this->first().restrictLocal(father,son,initialize);
    this->second().restrictLocal(father,son,initialize);    
  }
  //! prolong data to children 
  template <class EntityType>
  void prolongLocal ( EntityType &father, EntityType &son, 
		      bool initialize ) const {
    this->first().prolongLocal(father,son,initialize);
    this->second().prolongLocal(father,son,initialize);    
  }
  
  //! prolong data to children 
  template <class CommunicatorImp>
  void addToList(CommunicatorImp& comm) 
  {
    this->first().addToList(comm); 
    this->second().addToList(comm);    
  }
};


/** \brief Interface default implementation for derived classes. 
    Note the difference to RestrictProlongDefaultImplementation, which
    represents the implementation for certain spaces. 
*/
template< class TraitsImp >
class RestrictProlongInterfaceDefault 
: public RestrictProlongInterface< TraitsImp >
{
  template< class IndexSetType, bool isFemIndexSet >
  struct Persistent
  {
    static bool check ( const IndexSetType &indexSet )
    {
      return indexSet.persistent();
    }
  };

  template< class IndexSetType >
  struct Persistent< IndexSetType, false >
  {
    static bool check ( const IndexSetType &indexSet )
    {
      return indexSet.adaptive();
    }
  };

  //! check persistence of index set (also includes backward compatibility) 
  template< class IndexSetType >
  bool checkPersistent ( const IndexSetType &indexSet ) const
  {
    const bool isFemIndexSet = Conversion< IndexSetType, EmptyIndexSet > :: exists;
    return Persistent< IndexSetType, isFemIndexSet > :: check( indexSet );
  }

protected:
  //! return true if father and son have the same index
  template< class IndexSetType, class EntityType >
  bool entitiesAreCopies ( const IndexSetType &indexSet,
                           const EntityType &father,
                           const EntityType &son ) const
  {
    assert( checkPersistent( indexSet ) );
    return (indexSet.index( father ) == indexSet.index( son ));
  }
};


/** \brief This is a general restriction/prolongation operator
    which is speciallized for some discreteFunctionSpaces (e.g. DG)
    ,i.e., over the second template argument.
    The RestrictProlongDefault class is then one which should be
    applied by the user.
 */
template <class DiscreteFunctionType,class DiscreteFunctionSpace>
class RestrictProlongDefaultImplementation;

/** \brief This is a wrapper for the default implemented
    restriction/prolongation operator, which only takes a discrete
    function template 
 */
template <class DiscreteFunctionType> 
class RestrictProlongDefault
: public RestrictProlongDefaultImplementation<
           DiscreteFunctionType,
           typename DiscreteFunctionType::DiscreteFunctionSpaceType >
{
  typedef RestrictProlongDefaultImplementation<
           DiscreteFunctionType,
           typename DiscreteFunctionType::DiscreteFunctionSpaceType >
  BaseType;
public:
  RestrictProlongDefault(DiscreteFunctionType& discreteFunction) :
    BaseType(discreteFunction) {
    discreteFunction.enableDofCompression();
  }
private:
  RestrictProlongDefault();
};



/** \brief This is a simple restriction/prolongation operator for
 piecewise constant data stored on elements. 
*/
template< class DiscreteFunctionType >
class RestrictProlongPieceWiseConstantData
: public RestrictProlongInterfaceDefault
  < RestrictProlongTraits< RestrictProlongPieceWiseConstantData< DiscreteFunctionType > > >
{
  typedef RestrictProlongInterfaceDefault
    < RestrictProlongTraits< RestrictProlongPieceWiseConstantData< DiscreteFunctionType > > >
    BaseType;

public:  
  typedef typename DiscreteFunctionType::LocalFunctionType LocalFunctionType;

  typedef typename DiscreteFunctionType::DiscreteFunctionSpaceType SpaceType; 
  typedef typename SpaceType :: GridPartType GridPartType;
  typedef typename SpaceType :: GridType GridType;

  typedef typename DiscreteFunctionType::RangeFieldType RangeFieldType;
  typedef typename DiscreteFunctionType::DomainType DomainType;

protected:
  using BaseType :: calcWeight;
  using BaseType :: entitiesAreCopies;

public:  
  //! Constructor
  explicit RestrictProlongPieceWiseConstantData( DiscreteFunctionType &df ) 
  : df_ (df), weight_( -1.0 )
  {}

  /** \brief explicit set volume ratio of son and father
   *
   *  \param[in]  weight  volume of son / volume of father
   *
   *  \note If this ratio is set, it is assume to be constant.
   */
  void setFatherChildWeight ( const RangeFieldType &weight ) const
  {
    weight_ = weight;
  }
  
  //! restrict data to father 
  template< class EntityType >
  void restrictLocal ( const EntityType &father, const EntityType &son, bool initialize ) const
  {
    // if father and son are copies, do nothing
    if( entitiesAreCopies( df_.space().indexSet(), father, son ) )
      return;

    assert( !father.isLeaf() );

    const RangeFieldType weight = (weight_ < 0.0) ? calcWeight( father, son ) : weight_; 

    assert( weight > 0.0 );
    
    LocalFunctionType lfFather = df_.localFunction( father );
    LocalFunctionType lfSon = df_.localFunction( son );

    const int numDofs = lfFather.numDofs();
    if( initialize )
    {
      for( int i = 0; i < numDofs; ++i )
        lfFather[ i ] = weight * lfSon[ i ];
    }
    else 
    {
      for( int i = 0; i < numDofs; ++i )
        lfFather[ i ] += weight * lfSon[ i ];
    }
  }

  //! prolong data to children 
  template< class EntityType >
  void prolongLocal ( const EntityType &father, const EntityType &son, bool initialize ) const
  {
    // if father and son are copies, do nothing
    if( entitiesAreCopies( df_.space().indexSet(), father, son ) )
      return;
    
    LocalFunctionType lfFather = df_.localFunction( father );
    LocalFunctionType lfSon = df_.localFunction( son );
    const int numDofs = lfFather.numDofs();
    for( int i = 0; i < numDofs; ++i )
      lfSon[ i ] = lfFather[ i ];
  }

  //! add discrete function to communicator 
  template <class CommunicatorImp> 
  void addToList(CommunicatorImp& comm) 
  {
    comm.addToList(df_);
  }

private:
  mutable DiscreteFunctionType & df_;
  mutable RangeFieldType weight_;
};

///@} 

} // end namespace Dune 

#endif
