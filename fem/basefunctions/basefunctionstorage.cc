#include "../caching/pointprovider.hh"

using std::make_pair;

namespace Dune {

  //- class SimpleStorage
  template <class FunctionSpaceImp>
  StorageBase<FunctionSpaceImp>::StorageBase(const FactoryType& factory) :
    storage_(),
    diffVar1_(0)
  {
    for (int i = 0; i < factory.numBaseFunctions(); ++i) {
      storage_.push_back(factory.baseFunction(i));
    }
  }

  template <class FunctionSpaceImp>
  StorageBase<FunctionSpaceImp>::~StorageBase()
  {
    for (size_t i = 0; i < storage_.size(); ++i) {
      delete storage_[i];
      storage_[i] = 0;
    }
  }

  template <class FunctionSpaceImp>
  int StorageBase<FunctionSpaceImp>::numBaseFunctions() const 
  {
    return storage_.size();
  }

  template <class FunctionSpaceImp>
  template <int diffOrd>
  void StorageBase<FunctionSpaceImp>::
  evaluate(int baseFunct, const FieldVector<int, diffOrd>& diffVar, 
           const DomainType& xLocal, RangeType& result) const
  {
    assert(baseFunct >= 0 && baseFunct < numBaseFunctions());
    storage_[baseFunct]->evaluate(diffVar, xLocal, result);
  }

  template <class FunctionSpaceImp>
  void StorageBase<FunctionSpaceImp>::
  jacobian(int baseFunct, const DomainType& xLocal, 
           JacobianRangeType& result) const
  {
    assert(baseFunct >= 0 && baseFunct < numBaseFunctions());
    RangeType tmp;
    
    for (int i = 0; i < DomainType::size; ++i) {
      diffVar1_[0] = i;
      storage_[baseFunct]->evaluate(diffVar1_, xLocal, tmp);
      for (int j = 0; j < RangeType::size; ++j) {
        result[j][i] = tmp[j];
      }
    }
  }

  //- Simple storage
  template <class FunctionSpaceImp>
  SimpleStorage<FunctionSpaceImp>::SimpleStorage(const FactoryType& factory) :
    StorageBase<FunctionSpaceImp>(factory)
  {}

  template <class FunctionSpaceImp>
  SimpleStorage<FunctionSpaceImp>::~SimpleStorage() {}

  template <class FunctionSpaceImp>
  template <int diffOrd, class QuadratureType>
  void SimpleStorage<FunctionSpaceImp>::
  evaluate(int baseFunct, const FieldVector<int, diffOrd>& diffVar,
           const QuadratureType& quad, int quadPoint,
           RangeType& result) const 
  {
    this->evaluate(baseFunct, diffVar, quad.point(quadPoint), result);
  }

  template <class FunctionSpaceImp>
  template <class QuadratureType>
  void SimpleStorage<FunctionSpaceImp>::
  jacobian(int baseFunct, const QuadratureType& quad, int quadPoint,
           JacobianRangeType& result) const
  {
    this->jacobian(baseFunct, quad.point(quadPoint), result);
  }

  //- Caching storage
  template <class FunctionSpaceImp>
  CachingStorage<FunctionSpaceImp>::CachingStorage(const FactoryType& fac) :
    StorageBase<FunctionSpaceImp>(fac),
    elementGeometry_(fac.geo())
  {}

  template <class FunctionSpaceImp>
  CachingStorage<FunctionSpaceImp>::~CachingStorage() {}

  template <class FunctionSpaceImp>
  template <class CacheQuadratureType>
  void CachingStorage<FunctionSpaceImp>::
  evaluate(int baseFunct, const FieldVector<int, 0>& diffVar,
           const CacheQuadratureType& quad, int quadPoint,
           RangeType& result) const 
  {
    RangeIteratorType it = ranges_.find(quad.id());
    if (it == ranges_.end()) {
      it = addEntry(quad).first;
    }
    assert(it != ranges_.end());

    // calculate and return values
    result = it->second[quad.cachingPoint(quadPoint)];
  }

  template <class FunctionSpaceImp>
  template <class CacheQuadratureType>
  void CachingStorage<FunctionSpaceImp>::
  jacobian(int baseFunct, const CacheQuadratureType& quad, int quadPoint,
           JacobianRangeType& result) const
  {
    JacobianRangeIteratorType it = jacobians_.find(quad.id());
    if (it == jacobians_.end()) {
      it = addEntry(quad).second;
    }
    assert(it != jacobians_.end());

    // calculate and return values
    result = it->second[quad.cachingPoint(quadPoint)];
  }

  template <class FunctionSpaceImp>
  template <class CacheQuadratureType>
  void CachingStorage<FunctionSpaceImp>::
  evaluate(int baseFunct, const FieldVector<int, 1>& diffVar,
           const CacheQuadratureType& quad, int quadPoint,
           RangeType& result) const
  {
    JacobianRangeIteratorType it = jacobians_.find(quad.id());
    if (it == jacobians_.end()) {
      it = addEntry(quad).second;
    }
    assert(it != jacobians_.end());

    // calculate and return values
    result = it->second[quad.cachingPoint(quadPoint)][diffVar[0]];
  }

  template <class FunctionSpaceImp>
  template <class CacheQuadratureType>
  typename CachingStorage<FunctionSpaceImp>::ReturnPairType
  CachingStorage<FunctionSpaceImp>::
  addEntry(const CacheQuadratureType& quad) const 
  {
    typedef PointProvider<
      RealType,DomainType::size,CacheQuadratureType::codim> PointProviderType;
    typedef typename PointProviderType::GlobalPointVectorType PointVectorType;

    const PointVectorType& points = 
      PointProviderType::getPoints(quad.id(), elementGeometry_);

    RangeIteratorType rit =
      ranges_.insert(make_pair(quad.id(), 
                               RangeVectorType(this->numBaseFunctions()))).first;

    JacobianRangeIteratorType jit =
      jacobians_.insert(make_pair(quad.id(),
                                  JacobianRangeVectorType(this->numBaseFunctions()))).first;

    FieldVector<int, 0> diffVar;
    for (size_t i = 0; i < points.size(); ++i) {
      // evaluate value and jacobian and store it 
      this->evaluate(i, diffVar, points[i], rit->second[i]);
      this->jacobian(i, points[i], jit->second[i]);
    }

    return std::make_pair(rit, jit);
  }

} // end namespace Dune
