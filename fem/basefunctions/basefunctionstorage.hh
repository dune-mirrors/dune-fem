#ifndef DUNE_BASEFUNCTIONSTORAGE_HH
#define DUNE_BASEFUNCTIONSTORAGE_HH

#include <map>
#include <vector>

#include <dune/common/array.hh>
#include <dune/fem/common/basefunctionfactory.hh>
#include <dune/grid/common/grid.hh>

#include <list>

namespace Dune {

  //! \brief Storage policy for base function sets.
  //! In a base function set, the base function values on quadrature points
  //! can either be cached or always recalculated. The storage policies
  //! CachingStorage and SimpleStorage do exactly that. The present class
  //! implements the common functionality and can be seen as a layer of
  //! abstraction in the access to basefunctions.
  class StorageInterface; 

  typedef std::list<StorageInterface *> StorageInterfaceListType;
  typedef std::pair< size_t , int > QuadratureIdentifierType; 
  typedef std::list< QuadratureIdentifierType > QuadratureListType; 
  static StorageInterfaceListType storageList_;
  static QuadratureListType quadratureList_; 
  
  class StorageInterface 
  {
    public: 
      StorageInterface() 
      { 
        //std::cout << "create storage " << this<< "\n";
        storageList_.push_back(this);
      }

      virtual ~StorageInterface() {
        typedef StorageInterfaceListType::iterator IteratorType;
        IteratorType endit = storageList_.end();
        for(IteratorType it = storageList_.begin(); it != endit; ++it)
        {
          if( (*it) == this )
          {
            storageList_.erase(it); 
            //std::cout << "remove storage from list " << this<< "\n";
            break;
          }
        }
      }

      template <class StorageImp>  
      void cacheExsistingQuadratures(StorageImp & storage) 
      {
        typedef QuadratureListType::iterator IteratorType;
        IteratorType endit = quadratureList_.end();
        for(IteratorType it = quadratureList_.begin(); it != endit; ++it)
        {
          size_t id = (*it).first;
          int codim = (*it).second;
          //std::cout << "update set for id " << id << " and cd " << codim << "\n";
          storage.addQuadrature(id,codim);
        }
      }

      virtual void addQuadrature(size_t id, int codim ) const = 0;

      template <class QuadratureType> 
      static void addQuadratureToList(const QuadratureType & quad) 
      {
        int codim = QuadratureType :: codimension;
        addQuadratureToList(quad.id(),codim); 
      }
      
      static void addQuadratureToList(size_t id, int codim) 
      {
        // store quadrature 
        QuadratureIdentifierType ident(id,codim); 
        quadratureList_.push_back(ident);
        
        typedef StorageInterfaceListType::iterator IteratorType;
        IteratorType endit = storageList_.end();
        for(IteratorType it = storageList_.begin(); it != endit; ++it)
        {
          //std::cout << "add quad \n";
          (*it)->addQuadrature(id,codim); 
        }
      }
  };
  
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

    inline void addQuadrature(size_t id, int codim) const {}
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

    template <class CacheQuadratureType>
    inline
    void evaluate(int baseFunct,
                  const FieldVector<int, 0>& diffVar,
                  const CacheQuadratureType& quad, int quadPoint, 
                  RangeType& result) const;

    template <class CacheQuadratureType>
    inline
    void evaluate(int baseFunct,
                  const FieldVector<int, 1>& diffVar,
                  const CacheQuadratureType& quad, int quadPoint, 
                  RangeType& result) const;

    template <class CacheQuadratureType>
    inline
    void jacobian(int baseFunct, 
                  const CacheQuadratureType& quad, int quadPoint, 
                  JacobianRangeType& result) const;

    
  private:
    typedef typename FunctionSpaceImp::RangeFieldType RealType;
    typedef Array<Array<RangeType> > RangeVectorType;
    typedef Array<Array<JacobianRangeType> > JacobianRangeVectorType;
    //typedef std::vector<std::vector<RangeType> > RangeVectorType;
    //typedef std::vector<std::vector<JacobianRangeType> > JacobianRangeVectorType;
    //typedef std::vector<RangeType>  RangeVectorType;
    //typedef std::vector<JacobianRangeType> JacobianRangeVectorType;
    typedef std::map<size_t, bool> RangeStoredType;
    typedef std::map<size_t, bool> JacobianRangeStoredType;
    typedef typename RangeStoredType::iterator RangeIteratorType;
    typedef typename JacobianRangeStoredType::iterator JacobianRangeIteratorType;
    typedef std::pair<
      RangeIteratorType, JacobianRangeIteratorType> ReturnPairType;
    //typedef std::vector<RangeVectorType> RangeContainerType;
    //typedef std::vector<JacobianRangeVectorType> JacobianRangeContainerType;
    typedef Array<RangeVectorType> RangeContainerType;
    typedef Array<JacobianRangeVectorType> JacobianRangeContainerType;

  private:
    // caches the quadrature, see also addEntry.. 
    inline void addQuadrature(size_t id, int codim) const;
 
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
