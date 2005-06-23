#ifndef ADI_VECDISCFNC_CC
#define ADI_VECDISCFNC_CC

#include <config.h>
#include "vecdiscfnc.hh"

using namespace Dune;

namespace Adi {
  
  template<class CompositeFunctionSpaceImp,
           template <class> class ContainedFunctionImp>
  CompositeDiscreteFunction<CompositeFunctionSpaceImp,
                            ContainedFunctionImp>::
  CompositeDiscreteFunction(const CompositeFunctionSpaceImp& f) :
    DiscreteFunctionDefaultType(f) {
    for (int i = 0; i < CompositeFunctionSpaceImp::DimRange; ++i) {
      children_[i] = new ContainedFunctionType(f.containedFunctionSpace());
    }
  }

  template<class CompositeFunctionSpaceImp,
           template <class> class ContainedFunctionImp>
  CompositeDiscreteFunction<CompositeFunctionSpaceImp,
                            ContainedFunctionImp>::
  ~CompositeDiscreteFunction() {
    for (int i = 0; i < CompositeFunctionSpaceImp::DimRange; ++i) {
      delete children_[i];
      children_[i] = 0;
    }
  }

  template<class CompositeFunctionSpaceImp,
           template <class> class ContainedFunctionImp>
  CompositeDiscreteFunction<CompositeFunctionSpaceImp,
                            ContainedFunctionImp>::
  CompositeDiscreteFunction(const CompositeDiscreteFunction<CompositeFunctionSpaceImp,
                            ContainedFunctionImp>& other) :
    DiscreteFunctionDefaultType(other.functionSpace_)
  {
    // Copy objects behind pointers
    for (int i = 0; i < CompositeFunctionSpaceImp::DimRange; ++i) {
      children_[i] = new ContainedFunctionType(*(other.children_[i]));
    }
  }

  /* leads to virtual function override (definition in Vector<...>)
  template<class CompositeFunctionSpaceImp,
           template <class> class ContainedFunctionImp>
  CompositeDiscreteFunction<CompositeFunctionSpaceImp,
                            ContainedFunctionImp>&
  CompositeDiscreteFunction<CompositeFunctionSpaceImp,
                            ContainedFunctionImp>::
  operator= (const MyType& other) {
    for (int i = 0; i < DimRange; ++i) {
      // * this approach is too costly
      // delete children_[i];
      // children_[i] = new ContainedFunctionType(*(other.children_[i]));
      // call assignment in contained function
      *(children_[i]) = *(other.children_[i]);
    }
  }
  */

  template<class CompositeFunctionSpaceImp,
           template <class> class ContainedFunctionImp>
  ContainedFunctionImp<typename CompositeFunctionSpaceImp::ContainedFunctionSpaceType>&
  CompositeDiscreteFunction<CompositeFunctionSpaceImp,
                            ContainedFunctionImp>::
  getDiscreteFunction(int position) {
    assert(position < DimRange);
    return *(children_[position]);
  }

  template<class CompositeFunctionSpaceImp,
           template <class> class ContainedFunctionImp>
  LocalFunctionComposite<CompositeDiscreteFunction<CompositeFunctionSpaceImp,
                                                   ContainedFunctionImp> >
  CompositeDiscreteFunction<CompositeFunctionSpaceImp,
                            ContainedFunctionImp>::
  newLocalFunction() {
    return CompositeLocalFunctionType(this->functionSpace_, *this);
  }

  template<class CompositeFunctionSpaceImp,
           template <class> class ContainedFunctionImp>
  template<class EntityType>
  void
  CompositeDiscreteFunction<CompositeFunctionSpaceImp,
                            ContainedFunctionImp>::
  localFunction(EntityType& en,
                LocalFunctionComposite<MyType>& lf) {
    //lf.init(*this, en);
    for (int i = 0; i < DimRange; ++i) {
      children_[i]->localFunction(en, getChild(i, lf));
      //children_[i]->localFunction(en, *(lf.children_[i]));
    }
    lf.setNDof();
  }

  template<class CompositeFunctionSpaceImp,
           template <class> class ContainedFunctionImp>
  DofIteratorComposite<CompositeDiscreteFunction<CompositeFunctionSpaceImp,
                                                 ContainedFunctionImp> >
  CompositeDiscreteFunction<CompositeFunctionSpaceImp,
                            ContainedFunctionImp>::dbegin() {
    return DofIteratorType(*this);
  }

  template<class CompositeFunctionSpaceImp,
           template <class> class ContainedFunctionImp>
  DofIteratorComposite<CompositeDiscreteFunction<CompositeFunctionSpaceImp,
                                                 ContainedFunctionImp> >
  CompositeDiscreteFunction<CompositeFunctionSpaceImp,
                            ContainedFunctionImp>::dend() {
    return DofIteratorType(*this, 0);
  }

  template<class CompositeFunctionSpaceImp,
           template <class> class ContainedFunctionImp>
  ConstDofIteratorDefault<
    DofIteratorComposite<CompositeDiscreteFunction<CompositeFunctionSpaceImp,
                                                   ContainedFunctionImp> >
  >
  CompositeDiscreteFunction<CompositeFunctionSpaceImp,
                            ContainedFunctionImp>::dbegin() const {
     DofIteratorType tmp(const_cast<MyType&>(*this));
     return ConstDofIteratorType(tmp);
  }

  template<class CompositeFunctionSpaceImp,
           template <class> class ContainedFunctionImp>
  ConstDofIteratorDefault<
    DofIteratorComposite<CompositeDiscreteFunction<CompositeFunctionSpaceImp,
                                                   ContainedFunctionImp> >
  >
  CompositeDiscreteFunction<CompositeFunctionSpaceImp,
                            ContainedFunctionImp>::dend() const {
    DofIteratorType tmp(const_cast<MyType&>(*this), 0);
    return ConstDofIteratorType(tmp);
  }

  template<class CompositeFunctionSpaceImp,
           template <class> class ContainedFunctionImp>
  void
  CompositeDiscreteFunction<CompositeFunctionSpaceImp,
                            ContainedFunctionImp>::clear() {
    for (int i = 0; i < DimRange; ++i) {
      children_[i]->clear();
    }
  }

  template<class CompositeFunctionSpaceImp,
           template <class> class ContainedFunctionImp>
  void
  CompositeDiscreteFunction<CompositeFunctionSpaceImp,
                            ContainedFunctionImp>::set(RangeFieldType val) {
    for (int i = 0; i < DimRange; ++i) {
      children_[i]->set(val);
    }
  }

  //- class LocalFunctionComposite

  template <class CompositeFunctionImp>
  LocalFunctionComposite<CompositeFunctionImp>::
  LocalFunctionComposite(const CompositeFunctionSpaceType& spc, 
                         CompositeFunctionType& f) :
    fSpace_(spc)
  {
    typename CompositeFunctionType::ContainerType& containedFuncs = f.children_;
    for (int i = 0; i < DimRange; ++i) {
      const ContainedLocalFunctionType& tmp = containedFuncs[i]->newLocalFunction();
      children_[i] = new ContainedLocalFunctionType(tmp);
    }
  }

  template <class CompositeFunctionImp>
  LocalFunctionComposite<CompositeFunctionImp>::
  LocalFunctionComposite(const LocalFunctionComposite& other) :
    fSpace_(other.fSpace_) 
  {
    for (int i = 0; i < DimRange; ++i) {
      children_[i] = new ContainedLocalFunctionType(*(other.children_[i]));
    }
  }

  template <class CompositeFunctionImp>
  LocalFunctionComposite<CompositeFunctionImp>::
  ~LocalFunctionComposite() {
    for (int i = 0; i < DimRange; ++i) {
      delete children_[i];
      children_[i] = 0;
    }
  }

  template <class CompositeFunctionImp>
  typename LocalFunctionComposite<CompositeFunctionImp>::RangeFieldType&
  LocalFunctionComposite<CompositeFunctionImp>::
  operator[] (int num) {
    // Use fact that the number of dofs is the same for every contained func
    // * Cache that? -> Look which function can change the number in the contained local function to update correctly
    /* Uncached
      int nDofChildren = children_[0]->numberOfDofs();
      int funcNum = num/nDofChildren_;
      int position = num%nDofChildren_;
      
      return children_[funcNum]->operator[](position); 
    */
    return children_[num/nDofChildren_]->operator[](num%nDofChildren_);
  } 

  template <class CompositeFunctionImp>
  const typename LocalFunctionComposite<CompositeFunctionImp>::RangeFieldType&
  LocalFunctionComposite<CompositeFunctionImp>::
  operator[] (int num) const {
    /* Uncached
      int nDofChildren = children_[0]->numberOfDofs();
      int funcNum = num/nDofChildren_;
      int position = num%nDofChildren_;
 
      return children_[funcNum]->operator[](position);
    */
    return children_[num/nDofChildren_]->operator[](num%nDofChildren_);
  }

  template <class CompositeFunctionImp>
  int
  LocalFunctionComposite<CompositeFunctionImp>::
  numberOfDofs () const {
    // * More general, but slower variant
    /*
    int result = 0;
    for (int i = 0; i < DimRange; ++i) {
      result += children_[i]->numberOfDofs();
    }
    return result;
    */
    // * Faster variant, assuming identical contained spaces
      return DimRange*children_[0]->numberOfDofs();
  }

  template <class CompositeFunctionImp>
  template <class EntityType>
  void
  LocalFunctionComposite<CompositeFunctionImp>::
  evaluate (EntityType& en, 
            const DomainType& x, 
            RangeType& ret) const {
    ContainedRangeType tmp;
    for (int i = 0; i < DimRange; ++i) {
      children_[i]->evaluate(en, x, tmp);
      // * Assertion: RangeType is a scalar!
      ret[i] = tmp[0];
    }
  }

  template <class CompositeFunctionImp>
  template <class EntityType, class QuadratureType>
  void
  LocalFunctionComposite<CompositeFunctionImp>::
  evaluate (EntityType& en, 
            QuadratureType &quad, 
            int quadPoint, 
            RangeType& ret) const {
    ContainedRangeType tmp;
    for (int i = 0; i < DimRange; ++i) {
      children_[i]->evaluate(en, quad, quadPoint, tmp);
      // * Assertion: RangeType is a scalar!
      ret[i] = tmp[0];
    }
  }

  /*
  template <class CompositeFunctionImp>
  template <class EntityType>
  bool
  LocalFunctionComposite<CompositeFunctionImp>::
  init (CompositeFunctionImp& df, EntityType& en) const {
    //bool result = true;
    for (int i = 0; i < DimRange; ++i) {
      df.getDiscreteFunction(i).localFunction(en, *children_[i]);
      //result &= children_[i]->init(en);
    }
    nDofChildren_ = children_[0]->numberOfDofs();
    return true;
  }
  */

  //- class DofIteratorComposite
  template <class CompositeFunctionImp>
  DofIteratorComposite<CompositeFunctionImp>::
  DofIteratorComposite(const ThisType& other) :
    children_(other.children_),
    funcIter_(other.funcIter_),
    endFuncIter_(other.endFuncIter_),
    dofIter_(other.dofIter_),
    endDofIter_(other.endDofIter_) {}

  template <class CompositeFunctionImp>
  DofIteratorComposite<CompositeFunctionImp>&
  DofIteratorComposite<CompositeFunctionImp>::
  operator= (const ThisType& other) {
    if (*this != other) {
      children_ = other.children_;
      funcIter_ = other.funcIter_;
      endFuncIter_ = other.endFuncIter_;
      dofIter_ = other.dofIter_;
      endDofIter_ = other.endDofIter_;
    }
    return *this;
  }

  template <class CompositeFunctionImp>
  typename DofIteratorComposite<CompositeFunctionImp>::DofType&
  DofIteratorComposite<CompositeFunctionImp>::
  operator* () {
    return *dofIter_;
  }

  template <class CompositeFunctionImp>
  const  typename DofIteratorComposite<CompositeFunctionImp>::DofType & 
  DofIteratorComposite<CompositeFunctionImp>::
  operator* () const {
    return *dofIter_;
  }
  
  template <class CompositeFunctionImp>
  DofIteratorComposite<CompositeFunctionImp>&
  DofIteratorComposite<CompositeFunctionImp>::
  operator++ () {
    ++dofIter_;
    if (dofIter_ == endDofIter_) {
      ++funcIter_;
      if (funcIter_ != endFuncIter_) {
        dofIter_ = (**funcIter_).dbegin();
        endDofIter_ = (**funcIter_).dend();
      }
    }
    return *this;
  }
  
  template <class CompositeFunctionImp>
  bool 
  DofIteratorComposite<CompositeFunctionImp>::
  operator == (const ThisType& I ) const {
    return (dofIter_ == I.dofIter_) && (funcIter_ == I.funcIter_);
  }
  
  template <class CompositeFunctionImp>
  bool 
  DofIteratorComposite<CompositeFunctionImp>::
  operator != (const ThisType& I ) const {
    return !(*this == I);
  }
  
  template <class CompositeFunctionImp>
  void 
  DofIteratorComposite<CompositeFunctionImp>::
  reset () {
    assert(children_ != 0);
    funcIter_ = children_->begin();
    dofIter_ = (*funcIter_)->dbegin();
    endDofIter_ = (*funcIter_)->dend();
  }
  

} // end namespace Adi

#endif

