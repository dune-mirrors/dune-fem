#ifndef DUNE_NEWBASE_HH
#define DUNE_NEWBASE_HH

#include <dune/fem/common/basefunctions.hh>
#include <dune/fem/common/basefunctionfactory.hh>
#include <dune/fem/space/dofstorage.hh>
#include <dune/common/fvector.hh>


namespace Dune {

  template <class FunctionSpaceImp, template <class> class StorageImp>
  class StandardBaseFunctionSet;
  template <class FunctionSpaceImp, template <class> class StorageImp>
  class VectorialBaseFunctionSet;

  template <class FunctionSpaceImp, template <class> class StorageImp>
  struct StandardBaseFunctionSetTraits
  {
    typedef FunctionSpaceImp FunctionSpaceType;
    typedef StorageImp<FunctionSpaceType> StorageType;
    typedef StandardBaseFunctionSet<FunctionSpaceType, 
                                    StorageImp> BaseFunctionSetType;
    typedef BaseFunctionFactory<FunctionSpaceType> FactoryType;

  };

  template <class FunctionSpaceImp, template <class> class StorageImp>
  class StandardBaseFunctionSet : 
    public BaseFunctionSetDefault<StandardBaseFunctionSetTraits<FunctionSpaceImp, StorageImp> >
  {
  public:
    typedef StandardBaseFunctionSetTraits<FunctionSpaceImp, StorageImp> Traits;
    typedef typename FunctionSpaceImp::DomainType DomainType;
    typedef typename FunctionSpaceImp::RangeType RangeType;
    typedef typename FunctionSpaceImp::JacobianRangeType JacobianRangeType;

  public:
    StandardBaseFunctionSet(const typename Traits::FactoryType& factory) :
      storage_(factory)
    {}

    int numBaseFunctions() const;

    template <int diffOrd>
    void evaluate(int baseFunct,
                  const FieldVector<int, diffOrd>& diffVar,
                  const DomainType& xLocal,
                  RangeType& phi) const;

    template <int diffOrd, class QuadratureType>
    void evaluate (int baseFunct, 
                   const FieldVector<int, diffOrd> &diffVariable, 
                   QuadratureType & quad, 
                   int quadPoint, RangeType & phi) const;
      
  private:
    typename Traits::StorageType storage_;
  };

  //- VectorialBaseFunctionSet
  template <class FunctionSpaceImp, template <class> class StorageImp>
  struct VectorialBaseFunctionSetTraits 
  {
    typedef FunctionSpaceImp FunctionSpaceType;
    typedef VectorialBaseFunctionSet<
      FunctionSpaceType, StorageImp> BaseFunctionSetType;
  };

  template <class FunctionSpaceImp, template <class> class StorageImp>
  class VectorialBaseFunctionSet : 
    public BaseFunctionSetDefault<VectorialBaseFunctionSetTraits<FunctionSpaceImp, StorageImp> >
  {
  private:
    typedef BaseFunctionSetDefault<
      VectorialBaseFunctionSetTraits<FunctionSpaceImp, StorageImp> > BaseType;
    typedef FunctionSpaceImp FunctionSpaceType;
    typedef typename ToScalarFunctionSpace<
      FunctionSpaceImp>::Type ScalarFunctionSpaceType;
    typedef typename ScalarFunctionSpaceType::RangeType ScalarRangeType;
    typedef typename ScalarFunctionSpaceType::JacobianRangeType 
      ScalarJacobianRangeType;
    typedef BaseFunctionFactory<ScalarFunctionSpaceType> FactoryType;
    typedef typename FactoryType::BaseFunctionType BaseFunctionType;
    typedef StorageImp<ScalarFunctionSpaceType> StorageType;
    typedef VectorialBaseFunctionSetTraits<FunctionSpaceImp,StorageImp> Traits;
 
  public:
    typedef typename FunctionSpaceType::DomainType DomainType;
    typedef typename FunctionSpaceType::RangeType RangeType;
    typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;
    
    typedef typename FunctionSpaceType::RangeFieldType DofType;
        
  public:
    VectorialBaseFunctionSet(const FactoryType& factory) :
      storage_(factory),
      util_(FunctionSpaceType::DimRange),
      tmp_(0.),
      jTmp_(0.)
    {}

    ~VectorialBaseFunctionSet() {}

    int numBaseFunctions() const;

    template <int diffOrd>
    void evaluate(int baseFunct,
                  const FieldVector<int, diffOrd>& diffVar,
                  const DomainType& xLocal,
                  RangeType& phi) const;

    template <int diffOrd, class QuadratureType>
    void evaluate(int baseFunct, 
                  const FieldVector<deriType, diffOrd> &diffVariable, 
                  QuadratureType & quad, 
                  int quadPoint, RangeType & phi ) const;

    void jacobian(int baseFunct, const DomainType& xLocal, 
                  JacobianRangeType& gradPhi) const;

    template <class QuadratureImp>
    void jacobian(int baseFunct, QuadratureImp& quad, int quadPoint,
                  JacobianRangeType& gradPhi) const;

    // * add those methods with quadratures as well
    DofType evaluateSingle(int baseFunct, 
                           const DomainType& xLocal,
                           const RangeType& factor) const;
    
    template <class Entity>
    DofType evaluateGradientSingle(int baseFunct,
                                   Entity& en,
                                   const DomainType& xLocal,
                                   const JacobianRangeType& factor) const;

  private:
    StorageType storage_;
    DofConversionUtility<PointBased> util_;

    mutable FieldVector<int, 0> diffVar0_;
    mutable FieldVector<int, 1> diffVar1_;
    mutable ScalarRangeType tmp_;
    mutable ScalarJacobianRangeType jTmp_;
  };


} // end namespace Dune

#include "basefunctionsets.cc"

#endif
