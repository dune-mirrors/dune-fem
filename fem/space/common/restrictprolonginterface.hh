#ifndef DUNE_RESTRICTPROLONGINTERFACE_HH
#define DUNE_RESTRICTPROLONGINTERFACE_HH

//- Dune includes 
#include <dune/common/bartonnackmanifcheck.hh>

//- local includes 
#include <dune/fem/misc/combineinterface.hh>

namespace Dune{
/** @defgroup RestrictProlongInterface Restrict Prolong Interface
    @ingroup Adaptation
    Interface for restrict prolong operation on single elements.
  @{
 */

/*! @ingroup RestrictProlongInterface
    \brief Interface class defining the local behaviour of the
    restrict/prolong operation (using BN)
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
  template <class EntityType>
  RangeFieldType calcWeight (EntityType &father, EntityType &son) const
  {
    assert( (son.geometry().volume() / father.geometry().volume()) > 0.0 );
    return (son.geometry().volume() / father.geometry().volume());
  }
  
private:
  //! Barton-Nackman Trick
  RestProlImp& asImp() {
    return static_cast<RestProlImp&>(*this);
  }
  //! Barton-Nackman Trick
  const RestProlImp& asImp() const {
    return static_cast<const RestProlImp&>(*this);
  }
};

} // end namespace Dune 
#include <dune/fem/space/common/restrictprolongimpl.hh>
#endif
