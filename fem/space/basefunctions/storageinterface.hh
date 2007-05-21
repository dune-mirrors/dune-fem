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
  class StorageInterface
  {
      typedef std::list<StorageInterface *> StorageInterfaceListType;
      typedef std::pair< size_t , int > QuadIdPairType; 
      typedef std::pair< const GeometryType, QuadIdPairType > QuadratureIdentifierType;
      typedef std::list< QuadratureIdentifierType > QuadratureListType;

      // singelton implementation 
      static StorageInterfaceListType & storageList ()
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
      //! initialize singletons 
      static void initialize() 
      {
        storageList();
        quadratureList();
      }
      
      //! Constructor, add me to the list of storages 
      StorageInterface()
      {
        storageList().push_back(this);
      }

      //! Destructor, remove me from the list of storages 
      virtual ~StorageInterface() 
      {
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

      static bool geometryEqual(const GeometryType& storageGeom,
                                const GeometryType& quadGeom,
                                const int codim)
      {
        // only cache quad that have same geometry type 
        // if dim of quad > 1 compare also type of element 
        const bool sameType = (quadGeom.dim() > 1) ? 
              (storageGeom.basicType() == quadGeom.basicType()) : true;
        // for codim 1 the type are not necessarily the same
        return ( storageGeom == quadGeom 
            || ((quadGeom.dim() + codim) == storageGeom.dim() && sameType));
      }

      //! for a newly created storage cache all existing quadratures 
      template <class StorageImp>
      void cacheExsistingQuadratures(StorageImp & storage)
      {
        const GeometryType storageGeom = storage.geometryType(); 
        typedef QuadratureListType::iterator IteratorType;
        IteratorType endit = quadratureList().end();
        for(IteratorType it = quadratureList().begin(); it != endit; ++it)
        {
          const GeometryType& quadGeom = (*it).first;
          // get codim of quad 
          const int codim = (*it).second.second;
          // for codim 1 the type are not necessarily the same
          if( geometryEqual(storageGeom,quadGeom,codim) )
          {
            size_t id = (*it).second.first;
            storage.cacheQuadrature(id,codim);
          }
        }
      }

      //! cache quadrature for given id and codim 
      virtual void cacheQuadrature(size_t id, int codim ) const = 0;

      //! return geometry type of base function set 
      virtual GeometryType geometryType() const = 0;

      template <class QuadratureType>
      static void registerQuadratureToStorages(const QuadratureType & quad)
      {
        int codim = QuadratureType :: codimension;
        registerQuadratureToStorages(quad,codim);
      }

      //! register quadrature for all existing storages 
      template <class QuadratureType>
      static void registerQuadratureToStorages(const QuadratureType & quad, int codim)
      {
        const GeometryType quadGeom = quad.geometry();
        const size_t id = quad.id();
        // store quadrature 
        QuadratureIdentifierType ident(quadGeom,std::make_pair(id,codim));
        quadratureList().push_back(ident);

        typedef StorageInterfaceListType::iterator IteratorType;
        IteratorType endit = storageList().end();
        for(IteratorType it = storageList().begin(); it != endit; ++it)
        {
          // only cache quad that have same geometry type 
          const GeometryType storageGeom = (*it)->geometryType(); 
          // check if type are equal
          if( geometryEqual(storageGeom,quadGeom,codim) )
          {
            (*it)->cacheQuadrature(id,codim);
          }
        }
      }
  };

} // end namespace Dune                                                                                                                                                               
#endif
