#ifndef BASEFUNCTION_STORAGE_INTERFACE_HH
#define BASEFUNCTION_STORAGE_INTERFACE_HH

//- system includes 
#include <map>
#include <vector>
#include <list>

namespace Dune {

  //! \brief Storage policy for base function sets.
  //! In a base function set, the base function values on quadrature points
  //! can either be cached or always recalculated. The storage policies
  //! CachingStorage and SimpleStorage do exactly that. The present class
  //! implements the common functionality and can be seen as a layer of
  //! abstraction in the access to basefunctions.
  class StorageInterface
  {
      typedef std::list<StorageInterface *> StorageInterfaceListType;
      typedef std::pair< size_t , int > QuadratureIdentifierType;
      typedef std::list< QuadratureIdentifierType > QuadratureListType;

      // singelton implementation 
      static StorageInterfaceListType &  storageList ()
      {
        static StorageInterfaceListType storageListObj;
        return storageListObj;
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
      virtual ~StorageInterface() {
        typedef StorageInterfaceListType::iterator IteratorType;
        IteratorType endit = storageList().end();
        for(IteratorType it = storageList().begin(); it != endit; ++it)
        {
          if( (*it) == this )
          {
            storageList().erase(it);
            break;
          }
        }
      }

      //! for a newly created storage cache all existing quadratures 
      template <class StorageImp>
      void cacheExsistingQuadratures(StorageImp & storage)
      {
        typedef QuadratureListType::iterator IteratorType;
        IteratorType endit = quadratureList().end();
        for(IteratorType it = quadratureList().begin(); it != endit; ++it)
        {
          size_t id = (*it).first;
          int codim = (*it).second;
          storage.cacheQuadrature(id,codim);
        }
      }

      //! cache quadrature for given id and codim 
      virtual void cacheQuadrature(size_t id, int codim ) const = 0;

      template <class QuadratureType>
      static void registerQuadratureToStorages(const QuadratureType & quad)
      {
        int codim = QuadratureType :: codimension;
        registerQuadratureToStorages(quad.id(),codim);
      }

      //! register quadrature for all existing storages 
      static void registerQuadratureToStorages(size_t id, int codim)
      {
        // store quadrature 
        QuadratureIdentifierType ident(id,codim);
        quadratureList().push_back(ident);

        typedef StorageInterfaceListType::iterator IteratorType;
        IteratorType endit = storageList().end();
        for(IteratorType it = storageList().begin(); it != endit; ++it)
        {
          //std::cout << "add quad \n";
          (*it)->cacheQuadrature(id,codim);
        }
      }
  };

} // end namespace Dune                                                                                                                                                               
#endif
