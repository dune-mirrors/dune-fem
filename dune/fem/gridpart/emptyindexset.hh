#ifndef DUNE_FEM_EMPTYINDEXSET_HH
#define DUNE_FEM_EMPTYINDEXSET_HH

//- system includes 
#include <iostream>
#include <string> 
#include <cassert>

#include <dune/common/exceptions.hh>

/** @file
 @brief Provides default empty index set class for persistent index sets. 
*/

namespace Dune
{

  namespace Fem 
  { 

    /*!
      The EmptyIndexSet implements all additional method of a DUNE fem index set with 
      an empty default implementation. 
    */
    class EmptyIndexSet
    {
      // dummy value 
      enum { myType = -1 };
    public:  
      //! return false mean the no memory has to be allocated 
      //! and no compress of date has to be done 
      bool compress () { 
        return false; 
      }

      //! return true if the index set is consecutive 
      bool consecutive () const
      {
        return false;
      }

      //! return true if the index set is persistent 
      bool persistent () const
      {
        return false;
      }

      //! do nothing here, because fathers index should already exist 
      template< class EntityType >
      void insertEntity ( const EntityType &entity )
      {}

      //! do nothing here, because fathers index should already exist 
      template< class EntityType >
      void removeEntity ( const EntityType &entity )
      {}

      //! nothing to do here 
      void resize ()
      {}

      //! no extra memory for restriction is needed
      int additionalSizeEstimate () const
      {
        return 0;
      }

      static int type ()
      {
        return myType;
      }

      //! we have no old size 
      int numberOfHoles ( const int codim ) const
      {
        return 0;
      }
      
      //! return old index, for dof manager only 
      int oldIndex ( const int hole, const int codim ) const
      {
        return 0;
      }
      
      //! return new index, for dof manager only 
      int newIndex ( const int hole, const int codim ) const
      {
        return 0;
      }
    };

  }  // namespace Fem 

} // namespace Dune 
#endif // #ifndef DUNE_FEM_EMPTYINDEXSET_HH
