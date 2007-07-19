#include <dune/fem/quadrature/caching/pointprovider.hh>

using std::make_pair;

namespace Dune {

  //- class SimpleStorage
  template <class FunctionSpaceImp>
  StorageBase<FunctionSpaceImp>::StorageBase(const FactoryType& factory) :
    storageSize_(factory.numBaseFunctions()),
    storage_(new BaseFunctionType*[factory.numBaseFunctions()]),
    diffVar1_(0),
    elementGeometry_(factory.geometry())
  {
    for (int i = 0; i < factory.numBaseFunctions(); ++i) 
    {
      storage_[i]=factory.baseFunction(i);
    }
  }

  template <class FunctionSpaceImp>
  StorageBase<FunctionSpaceImp>::~StorageBase()
  {
    for (int i = 0; i < storageSize_; ++i) 
    {
      // delete base functions 
      delete storage_[i];
      storage_[i] = 0;
    }
    // delete storage vector 
    delete [] storage_;
  }

  template <class FunctionSpaceImp>
  int StorageBase<FunctionSpaceImp>::numBaseFunctions() const 
  {
    return storageSize_;
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
    
    for (int i = 0; i < DomainType::dimension; ++i) {
      diffVar1_[0] = i;
      storage_[baseFunct]->evaluate(diffVar1_, xLocal, tmp);
      for (int j = 0; j < RangeType::dimension; ++j) {
        result[j][i] = tmp[j];
      }
    }
  }

  //- Simple storage
  template <class FunctionSpaceImp>
  inline SimpleStorage<FunctionSpaceImp>::SimpleStorage(const FactoryType& factory) :
    StorageBase<FunctionSpaceImp>(factory)
  {}

  template <class FunctionSpaceImp>
  SimpleStorage<FunctionSpaceImp>::~SimpleStorage() {}

  template <class FunctionSpaceImp>
  template <int diffOrd, class QuadratureType>
  inline void SimpleStorage<FunctionSpaceImp>::
  evaluate(int baseFunct, const FieldVector<int, diffOrd>& diffVar,
           const QuadratureType& quad, int quadPoint,
           RangeType& result) const 
  {
    this->evaluate(baseFunct, diffVar, quad.point(quadPoint), result);
  }

  template <class FunctionSpaceImp>
  template <class QuadratureType>
  inline void SimpleStorage<FunctionSpaceImp>::
  jacobian(int baseFunct, const QuadratureType& quad, int quadPoint,
           JacobianRangeType& result) const
  {
    this->jacobian(baseFunct, quad.point(quadPoint), result);
  }

  //- Caching storage
  template <class FunctionSpaceImp>
  inline CachingStorage<FunctionSpaceImp>::CachingStorage(const FactoryType& fac) :
    StorageBase<FunctionSpaceImp>(fac)
  {
    this->cacheExsistingQuadratures(*this);
  }

  template <class FunctionSpaceImp>
  CachingStorage<FunctionSpaceImp>::~CachingStorage() {}

  template <class FunctionSpaceImp>
  template <class QuadratureType>
  inline void CachingStorage<FunctionSpaceImp>::
  evaluate(int baseFunct, const FieldVector<int, 0>& diffVar,
           const QuadratureType& quad, int quadPoint,
           RangeType& result) const 
  {
    assert( (Conversion<QuadratureType,CachingInterface>::exists) ?
             (rangestored_.find(quad.id()) != rangestored_.end()) : true );
    Evaluate<QuadratureType,
             Conversion<QuadratureType,CachingInterface>::exists >::
         evaluate(*this,baseFunct,diffVar,quad,quadPoint,ranges_,result);
  }

  template <class FunctionSpaceImp>
  template <class QuadratureType>
  inline void CachingStorage<FunctionSpaceImp>::
  jacobian(const int baseFunct, 
           const QuadratureType& quad, 
           const int quadPoint,
           JacobianRangeType& result) const
  {
    assert( (Conversion<QuadratureType,CachingInterface>::exists) ?
        ( jacobianstored_.find(quad.id()) != jacobianstored_.end()) : true);
    Evaluate<QuadratureType,
             Conversion<QuadratureType,CachingInterface>::exists >::
      jacobian(*this,baseFunct,quad,quadPoint,jacobians_,result);
  }

  template <class FunctionSpaceImp>
  template <class QuadratureType>
  inline void CachingStorage<FunctionSpaceImp>::
  evaluate(int baseFunct, const FieldVector<int, 1>& diffVar,
           const QuadratureType& quad, int quadPoint,
           RangeType& result) const
  {
    assert( (Conversion<QuadratureType,CachingInterface>::exists) ?
        ( jacobianstored_.find(quad.id()) != jacobianstored_.end()) : true);
    Evaluate<QuadratureType,
             Conversion<QuadratureType,CachingInterface>::exists >::
      evaluate(*this,baseFunct,diffVar,quad,quadPoint,jacobians_,result);
  }

  template <class FunctionSpaceImp>
  inline void CachingStorage<FunctionSpaceImp>::
  cacheQuadrature(size_t id, int codim) const 
  {
    RangeIteratorType it = rangestored_.find(id);
    if (it == rangestored_.end()) 
    {
      it = addEntryInterface(id,codim).first;
    }

    assert(rangestored_.find(id) != rangestored_.end());
    assert(jacobianstored_.find(id) != jacobianstored_.end());
  }

  template <class FunctionSpaceImp>
  inline typename CachingStorage<FunctionSpaceImp>::ReturnPairType
  CachingStorage<FunctionSpaceImp>::
  addEntryInterface(size_t id , int codim ) const 
  {
    enum { dimension = DomainType::dimension };
    switch (codim) 
    {
      case 0: return addEntry<0>(id);
      case 1: return addEntry<1>(id);
      default: assert(false); abort();
    }

    assert(false);
    abort();
    // only fake 
    return addEntry<0>(id);
  }

  template <class FunctionSpaceImp>
  template <int codimension>
  inline typename CachingStorage<FunctionSpaceImp>::ReturnPairType
  CachingStorage<FunctionSpaceImp>::
  addEntry(size_t id) const 
  {
    enum { dimension = DomainType::dimension };

    size_t quadId = id;

    typedef PointProvider<RealType, dimension, codimension> PointProviderType;
    typedef typename PointProviderType::GlobalPointVectorType PointVectorType;

    const PointVectorType& points = 
      PointProviderType::getPoints(quadId, this->elementGeometry_);

    assert(rangestored_.find(quadId) == rangestored_.end());
    RangeIteratorType rit =
      rangestored_.insert(make_pair(quadId,true)).first; 
    assert(rangestored_.find(quadId) != rangestored_.end());

    assert(jacobianstored_.find(quadId) == jacobianstored_.end());
    JacobianRangeIteratorType jit =
      jacobianstored_.insert(make_pair(quadId,true)).first;

    assert(jacobianstored_.find(quadId) != jacobianstored_.end());

    FieldVector<int, 0> diffVar;

    if ((size_t)ranges_.size()   <= quadId) ranges_.resize(quadId+10);
    if ((size_t)jacobians_.size()<= quadId) jacobians_.resize(quadId+10);

    ranges_[quadId].resize(points.size());
    jacobians_[quadId].resize(points.size());

    for (size_t i = 0; i < points.size(); ++i) {
      ranges_[quadId][i].resize(this->numBaseFunctions());
      jacobians_[quadId][i].resize(this->numBaseFunctions());

      for (int j = 0; j < this->numBaseFunctions(); ++j) 
      {
        // evaluate value and jacobian and store it 
        this->evaluate(j, diffVar, points[i], ranges_[quadId][i][j]);
        this->jacobian(j, points[i], jacobians_[quadId][i][j]);
      }
    }

    return std::make_pair(rit, jit);
  }

} // end namespace Dune
