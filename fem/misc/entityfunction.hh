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

#include <config.h>
#include <dune/grid/io/file/dgfparser/gridtype.hh>
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
  
  template <class EntityFunctionImp>
  class EntityFunction
  {
  public:
    //! interface function, which must be provided by derived classes
    template <class EntityType, class QuadratureType, class RangeType>
    void evaluate(EntityType& en, QuadratureType& quad, int p, RangeType& ret)
          {
            CHECK_AND_CALL_INTERFACE_IMPLEMENTATION( 
                (  asImp().template evaluate<EntityType, QuadratureType, 
                   RangeType>( en, quad, p, ret) ) );
            
            asImp().template evaluate<EntityType, QuadratureType, 
                RangeType> ( en, quad, p, ret);
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
 *  \class EntityFunctionWrapper
 *  \brief The EntityFunctionWrapper provides an entity evaluation method
 *         based on an instance of a Dune::Function
 */
/*======================================================================*/
  
  template <class FunctionType>
  class EntityFunctionWrapper : EntityFunction 
  < EntityFunctionWrapper <FunctionType> >
  {
  public:
    //! constructor with function as argument
    EntityFunctionWrapper(FunctionType& function)
            : function_(function)        
          {
          };
    
    //! implemetation of the evaluation method as required by EntityFunction
    template <class EntityType, class QuadratureType, class RangeType>
    inline void evaluate(EntityType& en, QuadratureType& quad, 
                         int p, RangeType& ret)
          {
            typedef typename EntityType::ctype CoordType;
            enum { dim = EntityType::dimension};

            const FieldVector< CoordType, dim > & 
                glob = en.geometry().global(quad.point(p));            
            function_.evaluate(glob, ret);
          };

  private:
    FunctionType& function_;    
  }; // end class EntityFunctionWrapper
  
} // end namespace Dune

#endif
