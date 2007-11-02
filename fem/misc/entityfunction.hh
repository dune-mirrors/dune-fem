/**************************************************************************
**       Title: entityfunction class
**    $RCSfile$
**   $Revision$$Name$
**       $Date$
**   Copyright: GPL $Author$
** Description: Interface for functions (not a Dune::Function!), which 
**              have a quadrature-style evaluation possibilities
**
**************************************************************************/

#ifndef __ENTITYFUNCTION_HH__
#define __ENTITYFUNCTION_HH__

#warning "EntityFunction is deprecated. Use DiscreteFunctionAdapter instead"

#include <dune/common/bartonnackmanifcheck.hh>

namespace Dune
{  
  
/*======================================================================*/
/*!
 *  \class EntityFunction
 *  \brief The EntityFunction is an Interface for functions providing 
 *         evaluation functionality in a quadrature-style. 
 *  
 *  This type of function is more general than the Dune::Function: Functions, 
 *  which are depending on discretefunctions can be represented by this
 *  class and are then available for integration, interpolation, whatever...
 *  Of course Numerical components should be rewritten using this evaluation
 *  Syntax. A simple Wrapper around Dune::Function is available for use 
 *  of these functions as an EntityFunction
 */
/*======================================================================*/
  
  template <class EntityType, class EntityFunctionImp>
  class EntityFunction
  {
  public:

/*======================================================================*/
/*! 
 *  init: interface function, which must be provided by derived classes
 *        method should be called for setting the entity, on which (repeated)
 *        evaluations then can be performed by evaluate(...)
 *
 *   \param en the entity, on which subsequent evaluation will be performed
 */
/*======================================================================*/

    void init(const EntityType& en)
          {
            CHECK_AND_CALL_INTERFACE_IMPLEMENTATION( 
                (  asImp().template init<EntityType>( en ) ) );
            
            asImp().init( en );
          }

/*======================================================================*/
/*! 
 *  interface function, which must be provided by derived classes
 *  should perform local evaluation on the entity, which is 
 *  assumed to be set by init(...)
 *
 *   \param loc local coordinates of point to be evaluated
 *
 *   \param ret the return value of the local function evaluation
 */
/*======================================================================*/

    template <class DomainType, class RangeType>
    void evaluate(const DomainType& loc, RangeType& ret) const
          {
            CHECK_AND_CALL_INTERFACE_IMPLEMENTATION( 
                (  asImp().template evaluate<DomainType, 
                   RangeType>(loc, ret) ) );
            
            asImp().template evaluate<DomainType, 
                RangeType> ( loc, ret);
          }
    
  protected:
    //! Barton-Nackman method forwarding
    EntityFunctionImp &asImp() 
          { 
            return static_cast<EntityFunctionImp&>( *this ); 
          }
    
    //! Barton-Nackman method forwarding
    const EntityFunctionImp &asImp( ) const 
          { 
            return static_cast<const EntityFunctionImp&>( *this ); 
          }  
  }; // end of class EntityFunction

  
/*======================================================================*/
/*!
 *  \class EntityFunctionAdapter
 *
 *  \brief The EntityFunctionAdapter provides an entity evaluation method
 *         based on an instance of a Dune::Function
 */
/*======================================================================*/
  
  template <class EntityType, class FunctionType>
  class EntityFunctionAdapter : EntityFunction 
  < EntityType, EntityFunctionAdapter <EntityType, FunctionType> >
  {
  public:
    //! constructor with function as argument
    EntityFunctionAdapter(FunctionType& function)
            : function_(function), entityPtr_(NULL)
          {
          };

    //! implementation of the init-method as required by EntityFunction:
    //! does nothing except storing the adress of the entity

    inline void init(const EntityType& en)
          { 
            entityPtr_ = &en;
          };
    
    //! implemetation of the evaluation-method as required by EntityFunction
    //! as a simple mapping of the local to the global function
    //! Caution has to be taken in use: The entityPtr is assumed to be 
    //! valid! No check is performed here

    template <class DomainType, class RangeType>
    inline void evaluate(const DomainType& loc, RangeType& ret) const
          {
            const DomainType& glob = entityPtr_->geometry().global(loc);    
            function_.evaluate(glob, ret);
          };
    
  private:
    FunctionType& function_;    
    EntityType const * entityPtr_;
  }; // end class EntityFunctionAdapter
  
} // end namespace Dune

#endif
