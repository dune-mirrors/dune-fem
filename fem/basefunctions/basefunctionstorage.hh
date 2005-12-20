#ifndef DUNE_BASEFUNCTIONSTORAGE_HH
#define DUNE_BASEFUNCTIONSTORAGE_HH

#include <map>
#include <vector>

#include <dune/fem/common/basefunctionfactory.hh>
#include <dune/grid/common/grid.hh>

namespace Dune {

  //! \brief Storage policy for base function sets.
  //! In a base function set, the base function values on quadrature points
  //! can either be cached or always recalculated. The storage policies
  //! CachingStorage and SimpleStorage do exactly that. The present class
  //! implements the common functionality and can be seen as a layer of
  //! abstraction in the access to basefunctions.
  template <class FunctionSpaceImp>
  class StorageBase 
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
    int numBaseFunctions() const;

    //! Evaluation of base functions and their derivatives
    template <int diffOrd>
    void evaluate(int baseFunct,
                  const FieldVector<int, diffOrd>& diffVar,
                  const DomainType& xLocal, 
                  RangeType& result) const;

    //! Evaluation of jacobians
    void jacobian(int baseFunct, 
                  const DomainType& xLocal, 
                  JacobianRangeType& result) const;

  private:
    typedef typename FactoryType::BaseFunctionType BaseFunctionType;

  private:
    std::vector<BaseFunctionType*> storage_;
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
    void evaluate(int baseFunct, 
                  const FieldVector<int, diffOrd>& diffVar,
                  const QuadratureType& quad, int quadPoint, 
                  RangeType& result) const;
    
    template <class QuadratureType>
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

  public:
    //! Constructor
    CachingStorage(const FactoryType& factory);
    //! Destructor
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
    typedef std::vector<std::vector<RangeType> > RangeVectorType;
    typedef std::vector<std::vector<JacobianRangeType> > JacobianRangeVectorType;
    //typedef std::vector<RangeType>  RangeVectorType;
    //typedef std::vector<JacobianRangeType> JacobianRangeVectorType;
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
