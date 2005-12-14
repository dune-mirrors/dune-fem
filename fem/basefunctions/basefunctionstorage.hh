#ifndef DUNE_BASEFUNCTIONSTORAGE_HH
#define DUNE_BASEFUNCTIONSTORAGE_HH

#include <map>
#include <vector>

#include <dune/fem/common/basefunctionfactory.hh>
#include <dune/grid/common/grid.hh>

namespace Dune {

  template <class FunctionSpaceImp>
  class StorageBase 
  {
  public:
    typedef BaseFunctionFactory<FunctionSpaceImp> FactoryType;
    typedef typename FunctionSpaceImp::DomainType DomainType;
    typedef typename FunctionSpaceImp::RangeType RangeType;
    typedef typename FunctionSpaceImp::JacobianRangeType JacobianRangeType;

  public:
    StorageBase(const FactoryType& factory);
    virtual ~StorageBase();

    int numBaseFunctions() const;

    template <int diffOrd>
    void evaluate(int baseFunct,
                  const FieldVector<int, diffOrd>& diffVar,
                  const DomainType& xLocal, 
                  RangeType& result) const;

    void jacobian(int baseFunct, 
                  const DomainType& xLocal, 
                  JacobianRangeType& result) const;

  private:
    typedef typename FactoryType::BaseFunctionType BaseFunctionType;

  private:
    std::vector<BaseFunctionType*> storage_;
    mutable FieldVector<int, 1> diffVar1_;
  };

  template <class FunctionSpaceImp>
  class SimpleStorage : public StorageBase<FunctionSpaceImp> 
  {
  public:
    typedef BaseFunctionFactory<FunctionSpaceImp> FactoryType;
    typedef typename FunctionSpaceImp::DomainType DomainType;
    typedef typename FunctionSpaceImp::RangeType RangeType;
    typedef typename FunctionSpaceImp::JacobianRangeType JacobianRangeType;

  public:
    SimpleStorage(const FactoryType& factory);
    ~SimpleStorage();

    using StorageBase<FunctionSpaceImp>::evaluate;
    using StorageBase<FunctionSpaceImp>::jacobian;

    template <int diffOrd, class QuadratureType>
    void evaluate(int baseFunct, 
                  const FieldVector<int, diffOrd>& diffVar,
                  const QuadratureType& quad, int quadPoint, 
                  RangeType& result) const;
    
    template <class QuadratureType>
    void jacobian(int baseFunct, 
                  const QuadratureType& quad, int quadPoint, 
                  JacobianRangeType& result) const;

  };

  template <class FunctionSpaceImp>
  class CachingStorage
  {
  public:
    typedef BaseFunctionFactory<FunctionSpaceImp> FactoryType;
    typedef typename FunctionSpaceImp::DomainType DomainType;
    typedef typename FunctionSpaceImp::RangeType RangeType;
    typedef typename FunctionSpaceImp::JacobianRangeType JacobianRangeType;

  public:
    CachingStorage(const FactoryType& factory);
    ~CachingStorage();

    using StorageBase<FunctionSpaceImp>::evaluate;
    using StorageBase<FunctionSpaceImp>::jacobian;

    template <class CacheQuadratureType>
    void evaluate(int baseFunct,
                  const FieldVector<int, 0>& diffVar,
                  const CacheQuadratureType& quad, int quadPoint, 
                  RangeType& result) const;

    template <class CacheQuadratureType>
    void evaluate(int baseFunct,
                  const FieldVector<int, 1>& diffVar,
                  const CacheQuadratureType& quad, int quadPoint, 
                  RangeType& result) const;
 
    template <class CacheQuadratureType>
    void jacobian(int baseFunct, 
                  const CacheQuadratureType& quad, int quadPoint, 
                  JacobianRangeType& result) const;

    
  private:
    typedef typename FunctionSpaceImp::RangeFieldType RealType;
    typedef std::vector<RangeType> RangeVectorType;
    typedef std::vector<JacobianRangeType> JacobianRangeVectorType;
    typedef std::map<size_t, RangeVectorType> RangeContainerType;
    typedef std::map<size_t, JacobianRangeVectorType> JacobianRangeContainerType;
    typedef typename RangeContainerType::iterator RangeIteratorType;
    typedef typename JacobianRangeContainerType::iterator JacobianRangeIteratorType;
    typedef std::pair<
      RangeIteratorType, JacobianRangeIteratorType> ReturnPairType;

  private:
    template <class CacheQuadratureType>
    ReturnPairType addEntry(const CacheQuadratureType& quad) const;

  private:
    GeometryType elementGeometry_;

    mutable RangeContainerType ranges_;
    mutable JacobianRangeContainerType jacobians_;
  };

} // end namespace Dune

#include "basefunctionstorage.cc"

#endif
