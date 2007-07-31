#ifndef DUNE_ADAPTMANAGER_HH
#define DUNE_ADAPTMANAGER_HH

//- Dune includes 
#include <dune/common/array.hh>

//- local includes 
#include <dune/fem/space/common/dofmanager.hh>
#include <dune/fem/operator/common/objpointer.hh>

#include <dune/fem/space/common/communicationmanager.hh>
#include <dune/fem/space/common/loadbalancer.hh>

namespace Dune{

/** @defgroup Adaptation Adaptation 
    @{
 **/

/** \brief AdaptationManagerInterface class. 
 
 This Class is the result of a combination of different
 AdaptationOperators. It is the same principle as with Mapping. 
*/ 
class AdaptationManagerInterface : public LoadBalancerInterface 
{
public:
  //! \brief default constructor 
  AdaptationManagerInterface () : am_ (0) {}

  //! \brief destructor 
  virtual ~AdaptationManagerInterface () {}
  
  /** \brief on call of this method the internal adaptation operator is
      called. 
  */
  virtual void adapt ()  
  {
    //std::cout << "called AdaptationManagerInterface::adapt()" << std::endl;
    if(am_) am_->adapt();  
    else 
    {
      std::cerr << "WARNING: adapt! \n";
    }
  };

  /** \brief returns true if adaptation manager as adaptation method different to NONE
     \return <b>true</b> if adaptation method is not NONE, <b>false</b> otherwise
  */
  virtual bool adaptive () const  
  { 
    return (am_) ? (am_->adaptive()) : false; 
  } 

  /** \brief returns name of adaptation method 
     \return name of adaptation method 
  */
  virtual const char * methodName() const 
  {
    return (am_) ? (am_->methodName()) : "unknown method";
  }
    
  /** \brief Assignment operator, pointer to adaptation manager is stored 
      \return reference to this (i.e. *this) 
  */
  AdaptationManagerInterface & operator = (const AdaptationManagerInterface & am)
  {
      /** \todo This const-casting seems strange to me! */
    am_ = const_cast<AdaptationManagerInterface *> (&am);
    return (*this);
  }

  /** \brief @copydoc LoadBalancerInterface::loadBalance */
  virtual bool loadBalance () 
  { 
    return (am_) ? (am_->loadBalance()) : false; 
  }

  /** \brief @copydoc LoadBalancerInterface::balanceCounter */
  virtual int balanceCounter () const 
  { 
    return (am_) ? (am_->balanceCounter()) : 0; 
  }
 
private: 
  //! pointer to adaptation manager 
  AdaptationManagerInterface* am_;
};

//! deprecated typedef 
typedef AdaptationManagerInterface AdaptMapping;

} // end namespace Dune 

#include <dune/fem/space/common/adaptmanagerimpl.hh>
#endif
