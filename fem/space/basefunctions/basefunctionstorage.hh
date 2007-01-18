#ifndef DUNE_BASEFUNCTIONSTORAGE_HH
#define DUNE_BASEFUNCTIONSTORAGE_HH

//- system includes 
#include <map>
#include <vector>
#include <list>

//- Dune includes 
#include <dune/common/array.hh>
#include <dune/grid/common/grid.hh>

//- local includes 
#include <dune/fem/space/common/basefunctionfactory.hh>
#include "storageinterface.hh"

namespace Dune {

  //! forward deklaration for CachingQuadrature 
  template <typename GridImp, int codim>
  class CachingQuadrature;
   
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
    const int storageSize_;
    BaseFunctionType** storage_;
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
    typedef CachingStorage<FunctionSpaceImp> ThisType;
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

    template <class GridPartType,int cdim>
    inline
    void evaluate(int baseFunct,
                  const FieldVector<int, 0>& diffVar,
                  const CachingQuadrature<GridPartType, cdim>& quad, int quadPoint,
                  RangeType& result) const;
    
    template <class QuadratureType>
    inline
    void evaluate(int baseFunct,
                  const FieldVector<int, 1>& diffVar,
                  const QuadratureType& quad, int quadPoint, 
                  RangeType& result) const;

    template <class GridPartType,int cdim>
    inline
    void evaluate(int baseFunct,
                  const FieldVector<int, 1>& diffVar,
                  const CachingQuadrature<GridPartType, cdim>& quad, int quadPoint,
                  RangeType& result) const;

    template <class QuadratureType>
    inline
    void jacobian(int baseFunct, 
                  const QuadratureType& quad, int quadPoint, 
                  JacobianRangeType& result) const;

    template <class GridPartType,int cdim>
    inline
    void jacobian(int baseFunct,
                  const CachingQuadrature<GridPartType, cdim>& quad, int quadPoint,
                  JacobianRangeType& result) const;
    
  private:
    typedef typename FunctionSpaceImp::DomainFieldType RealType;
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
