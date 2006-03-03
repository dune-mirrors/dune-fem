#ifndef DUNE_BASEFUNCTIONSTORAGE_HH
#define DUNE_BASEFUNCTIONSTORAGE_HH

//- system includes 
#include <map>
#include <vector>
#include <list>

//- Dune includes 
#include <dune/common/array.hh>
#include <dune/fem/common/basefunctionfactory.hh>
#include <dune/grid/common/grid.hh>

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
}

//- local includes 
#include "../quadrature/cachequad.hh"

namespace Dune {

  template <class FunctionSpaceImp>
  class StorageBase : public StorageInterface 
  {
  public:
    typedef BaseFunctionFactory<FunctionSpaceImp> FactoryType;
    typedef typename FunctionSpaceImp::DomainType DomainType;
    typedef typename FunctionSpaceImp::RangeType RangeType;
    typedef typename FunctionSpaceImp::JacobianRangeType JacobianRangeType;

  public:
    //! Constructor
    StorageBase(const FactoryType& factory);
    //! Destructor (must be virtual)
    virtual ~StorageBase();

    //! Number of base functions
    inline
    int numBaseFunctions() const;

    //! Evaluation of base functions and their derivatives
    template <int diffOrd>
    inline
    void evaluate(int baseFunct,
                  const FieldVector<int, diffOrd>& diffVar,
                  const DomainType& xLocal, 
                  RangeType& result) const;

    //! Evaluation of jacobians
    inline
    void jacobian(int baseFunct, 
                  const DomainType& xLocal, 
                  JacobianRangeType& result) const;

    inline void cacheQuadrature(size_t id, int codim) const {}
  private:
    typedef typename FactoryType::BaseFunctionType BaseFunctionType;

  private:
    int storageSize_;
    BaseFunctionType** storage_;
    // std::vector<BaseFunctionType*> storage_;
    mutable FieldVector<int, 1> diffVar1_;
  };

  //! \brief A simple storage scheme which just forwards the calls to 
  //! the underlying base functions.
  template <class FunctionSpaceImp>
  class SimpleStorage : public StorageBase<FunctionSpaceImp> 
  {
  public:
    typedef BaseFunctionFactory<FunctionSpaceImp> FactoryType;
    typedef typename FunctionSpaceImp::DomainType DomainType;
    typedef typename FunctionSpaceImp::RangeType RangeType;
    typedef typename FunctionSpaceImp::JacobianRangeType JacobianRangeType;

  public:
    //! Constructor
    SimpleStorage(const FactoryType& factory);
    //! Destructor
    ~SimpleStorage();

    using StorageBase<FunctionSpaceImp>::evaluate;
    using StorageBase<FunctionSpaceImp>::jacobian;

    template <int diffOrd, class QuadratureType>
    inline
    void evaluate(int baseFunct, 
                  const FieldVector<int, diffOrd>& diffVar,
                  const QuadratureType& quad, int quadPoint, 
                  RangeType& result) const;
    
    template <class QuadratureType>
    inline
    void jacobian(int baseFunct, 
                  const QuadratureType& quad, int quadPoint, 
                  JacobianRangeType& result) const;

  };

  //! \brief Storage scheme which caches evaluations of base function values
  //! and derivatives.
  //! This storage scheme works in conjunction with the CacheQuadrature.
  //! \warning This works only for conforming grids so far!!!!!
  //! \todo Implement switch for non-conforming grids!!!!
  template <class FunctionSpaceImp>
  class CachingStorage : public StorageBase<FunctionSpaceImp>
  {
  public:
    typedef BaseFunctionFactory<FunctionSpaceImp> FactoryType;
    typedef typename FunctionSpaceImp::DomainType DomainType;
    typedef typename FunctionSpaceImp::RangeType RangeType;
    typedef typename FunctionSpaceImp::JacobianRangeType JacobianRangeType;

    friend class StorageInterface ;
  public:
    //! Constructor
    CachingStorage(const FactoryType& factory);
    //! Destructor
    ~CachingStorage();

    using StorageBase<FunctionSpaceImp>::evaluate;
    using StorageBase<FunctionSpaceImp>::jacobian;

    template <class QuadratureType>
    inline
    void evaluate(int baseFunct,
                  const FieldVector<int, 0>& diffVar,
                  const QuadratureType& quad, int quadPoint, 
                  RangeType& result) const;

    template <class GridType,int cdim>
    inline
    void evaluate(int baseFunct,
                  const FieldVector<int, 0>& diffVar,
                  const CachingQuadrature<GridType, cdim>& quad, int quadPoint,
                  RangeType& result) const;
    
    template <class QuadratureType>
    inline
    void evaluate(int baseFunct,
                  const FieldVector<int, 1>& diffVar,
                  const QuadratureType& quad, int quadPoint, 
                  RangeType& result) const;

    template <class GridType,int cdim>
    inline
    void evaluate(int baseFunct,
                  const FieldVector<int, 1>& diffVar,
                  const CachingQuadrature<GridType, cdim>& quad, int quadPoint,
                  RangeType& result) const;

    template <class QuadratureType>
    inline
    void jacobian(int baseFunct, 
                  const QuadratureType& quad, int quadPoint, 
                  JacobianRangeType& result) const;

    template <class GridType,int cdim>
    inline
    void jacobian(int baseFunct,
                  const CachingQuadrature<GridType, cdim>& quad, int quadPoint,
                  JacobianRangeType& result) const;
    
  private:
    typedef typename FunctionSpaceImp::RangeFieldType RealType;
    typedef Array<Array<RangeType> > RangeVectorType;
    typedef Array<Array<JacobianRangeType> > JacobianRangeVectorType;
    typedef std::map<size_t, bool> RangeStoredType;
    typedef std::map<size_t, bool> JacobianRangeStoredType;
    typedef typename RangeStoredType::iterator RangeIteratorType;
    typedef typename JacobianRangeStoredType::iterator JacobianRangeIteratorType;
    typedef std::pair<
      RangeIteratorType, JacobianRangeIteratorType> ReturnPairType;
    typedef Array<RangeVectorType> RangeContainerType;
    typedef Array<JacobianRangeVectorType> JacobianRangeContainerType;

  private:
    // caches the quadrature, see also addEntry.. 
    inline void cacheQuadrature(size_t id, int codim) const;
 
    // here a switch-case for codim is done and then addEntry called
    inline ReturnPairType addEntryInterface(size_t id , int codim) const;
    
    // here the real caching is done 
    template <int codim>
    inline ReturnPairType addEntry(size_t id) const;

  private:
    GeometryType elementGeometry_;

    mutable RangeContainerType ranges_;
    mutable JacobianRangeContainerType jacobians_;
    mutable RangeStoredType rangestored_;
    mutable JacobianRangeStoredType jacobianstored_;
  };

} // end namespace Dune

#include "basefunctionstorage.cc"

#endif
