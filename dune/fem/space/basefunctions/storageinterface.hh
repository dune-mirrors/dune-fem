#ifndef BASEFUNCTION_STORAGE_INTERFACE_HH
#define BASEFUNCTION_STORAGE_INTERFACE_HH

//- system includes 
#include <iostream>
#include <map>
#include <vector>
#include <list>

//- Dune includes 
#include <dune/common/geometrytype.hh>

namespace Dune {

  //! \brief Storage policy for base function sets.
  //! In a base function set, the base function values on quadrature points
  //! can either be cached or always recalculated. The storage policies
  //! CachingStorage and SimpleStorage do exactly that. The present class
  //! implements the common functionality and can be seen as a layer of
  //! abstraction in the access to basefunctions.
  template <int dimworld> 
  class StorageInterface
  {
      typedef StorageInterface<dimworld> ThisType;
      typedef std::list< ThisType* > StorageInterfaceListType;

      typedef StorageInterfaceListType* StorageInterfaceListPointer;
      typedef std::pair< size_t , int > QuadratureIdentifierType;
      typedef std::list< QuadratureIdentifierType > QuadratureListType;

      // return reference to list singleton pointer 
      static StorageInterfaceListPointer& storageListPtr () 
      {
        static StorageInterfaceListType* storageListObj = 0;
        return storageListObj;
      }

      // check if list is empty and if yes delete list 
      void checkAndDeleteStorageList() 
      {
        if( storageListPtr()->empty () )
        {
          delete storageListPtr();
          storageListPtr() = 0;
        }
      }

      // singelton implementation 
      static StorageInterfaceListType & storageList ()
      {
        // if list pointer is 0 then create new object 
        if( ! storageListPtr() )
        {
          storageListPtr() = new StorageInterfaceListType();
        }
        assert( storageListPtr() );
        return *(storageListPtr());
      }

      // singelton implementation 
      static QuadratureListType & quadratureList ()
      {
        static QuadratureListType quadratureListObj;
        return quadratureListObj;
      }
      
    public:
      //! Constructor, add me to the list of storages 
      StorageInterface()
      {
        storageList().push_back(this);
      }

      //! Destructor, remove me from the list of storages 
      virtual ~StorageInterface() 
      {
        typedef typename StorageInterfaceListType::iterator IteratorType;
        IteratorType endit = storageList().end();
        for(IteratorType it = storageList().begin(); it != endit; ++it)
        {
          if( (*it) == this )
          {
            storageList().erase(it);
            break;
          }
        }

        // if list is empty, then delete list 
        checkAndDeleteStorageList();
      }

      //! for a newly created storage cache all existing quadratures 
      template <class StorageImp>
      void cacheExistingQuadratures(StorageImp & storage)
      {
        typedef typename QuadratureListType::iterator IteratorType;
        IteratorType endit = quadratureList().end();
        for(IteratorType it = quadratureList().begin(); it != endit; ++it)
        {
          // get if and codim of quad 
          const size_t id = (*it).first;
          const int codim = (*it).second;
          storage.cacheQuadrature(id,codim);
        }
      }

      //! cache quadrature for given id and codim 
      virtual void cacheQuadrature(size_t id, int codim ) const = 0;

      //! return geometry type of base function set 
      virtual GeometryType geometryType() const = 0;

      template <class QuadratureType>
      static void registerQuadratureToStorages(const QuadratureType & quad)
      {
        const int codim = QuadratureType :: codimension;
        registerQuadratureToStorages(quad,codim);
      }

      //! register quadrature for all existing storages 
      template <class QuadratureType>
      static void registerQuadratureToStorages(const QuadratureType & quad, 
                                               const int codim)
      {
        const size_t id = quad.id();
        // store quadrature 
        QuadratureIdentifierType ident(id,codim);
        quadratureList().push_back(ident);

        typedef typename StorageInterfaceListType::iterator IteratorType;
        IteratorType endit = storageList().end();
        for(IteratorType it = storageList().begin(); it != endit; ++it)
        {
          (*it)->cacheQuadrature(id,codim);
        }
      }
  };

} // end namespace Dune                                                                                                                                                               
#endif
