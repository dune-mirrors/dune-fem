
namespace Dune {

  //- class StdBaseFunctionSet
  template <class FunctionSpaceImp, template <class> class StorageImp>
  int StdBaseFunctionSet<FunctionSpaceImp, StorageImp>::
  numBaseFunctions() const 
  {
    return storage_.numBaseFunctions();
  }

  template <class FunctionSpaceImp, template <class> class StorageImp>
  template <int diffOrd>
  void StdBaseFunctionSet<FunctionSpaceImp, StorageImp>::
  evaluate(int baseFunct,
           const FieldVector<int, diffOrd>& diffVar,
           const DomainType& xLocal,
           RangeType& phi) const 
  {
    storage_.evaluate(baseFunct, diffVar, xLocal, phi);
  }

  template <class FunctionSpaceImp, template <class> class StorageImp>
  template <int diffOrd, class QuadratureType>
  void StdBaseFunctionSet<FunctionSpaceImp, StorageImp>::
  evaluate(int baseFunct,
           const FieldVector<int, diffOrd>& diffVar,
           QuadratureType& quad, int quadPoint,
           RangeType& phi) const 
  {
    storage_.evaluate(baseFunct, diffVar, quad, quadPoint, phi);
  }

  template <class FunctionSpaceImp, template <class> class StorageImp>
  template <class QuadratureType>
  typename StdBaseFunctionSet<FunctionSpaceImp, StorageImp>::DofType
  StdBaseFunctionSet<FunctionSpaceImp, StorageImp>::
  evaluateSingle(int baseFunct,
                 const QuadratureType& quad, int quadPoint,
                 const RangeType& factor) const 
  {
    storage_.evaluate(baseFunct, diffVar0_, quad, quadPoint, tmp_);
    return tmp_*factor;
  }

  template <class FunctionSpaceImp, template <class> class StorageImp>
  template <class Entity, class QuadratureType>
  typename StdBaseFunctionSet<FunctionSpaceImp, StorageImp>::DofType
  StdBaseFunctionSet<FunctionSpaceImp, StorageImp>::
  evaluateGradientSingle(int baseFunct,
                         Entity& en,
                         const QuadratureType& quad, int quadPoint,
                         const JacobianRangeType& factor) const 
  {
    storage_.jacobian(baseFunct, quad, quadPoint, jTmp_);

    DofType result = 0;
    DomainType gradScaled(0.);
    for (int i = 0; i < FunctionSpaceImp::DimRange; ++i) {
      en.geometry().jacobianInverseTransposed(quad.point(quadPoint)).
        umv(jTmp_[i], gradScaled);
      result += gradScaled*factor[i];
    }
    return result;
  }
  
  //- class VecBaseFunctionSet
  template <class FunctionSpaceImp, template <class> class StorageImp>
  int VecBaseFunctionSet<FunctionSpaceImp, StorageImp>::
  numBaseFunctions() const 
  {
    return storage_.numBaseFunctions()*FunctionSpaceImp::DimRange;
  }

  template <class FunctionSpaceImp, template <class> class StorageImp>
  template <int diffOrd>
  void VecBaseFunctionSet<FunctionSpaceImp, StorageImp>::
  evaluate(int baseFunct,
           const FieldVector<int, diffOrd>& diffVar,
           const DomainType& xLocal,
           RangeType& phi) const 
  {
    ScalarRangeType tmp;
    storage_.evaluate(util_.containedDof(baseFunct), diffVar, xLocal, tmp);
    
    phi *= 0.0;
    phi[util_.component(baseFunct)] = tmp[0];
  }

  template <class FunctionSpaceImp, template <class> class StorageImp>
  template <int diffOrd, class QuadratureType>
  void VecBaseFunctionSet<FunctionSpaceImp, StorageImp>::
  evaluate(int baseFunct,
           const FieldVector<int, diffOrd>& diffVar,
           QuadratureType& quad, int quadPoint,
           RangeType& phi) const 
  {
    ScalarRangeType tmp;
    storage_.evaluate(util_.containedDof(baseFunct), diffVar, quad, quadPoint,
                      tmp);

    phi *= 0.;
    phi[util_.component(baseFunct)] = tmp[0];
  }

  template <class FunctionSpaceImp, template <class> class StorageImp>
  void VecBaseFunctionSet<FunctionSpaceImp, StorageImp>::
  jacobian(int baseFunct, const DomainType& xLocal, 
           JacobianRangeType& gradPhi) const 
  {
    gradPhi *= 0.;
    
    storage_.jacobian(util_.containedDof(baseFunct), xLocal, jTmp_); 
    gradPhi[util_.component(baseFunct)] = jTmp_[0];
  }

  template <class FunctionSpaceImp, template <class> class StorageImp>
  template <class QuadratureImp>
  void VecBaseFunctionSet<FunctionSpaceImp, StorageImp>::
  jacobian(int baseFunct, QuadratureImp& quad, int quadPoint, 
           JacobianRangeType& gradPhi) const 
  {
    gradPhi *= 0.;

    for (int i = 0; i < FunctionSpaceImp::DimDomain; ++i) {
      diffVar1_[0] = i;
      storage_.jacobian(util_.containedDof(baseFunct), quad, quadPoint, jTmp_); 
      gradPhi[util_.component(baseFunct)][i] = jTmp_[0];
    }
  }

  template <class FunctionSpaceImp, template <class> class StorageImp>
  typename VecBaseFunctionSet<FunctionSpaceImp, StorageImp>::DofType
  VecBaseFunctionSet<FunctionSpaceImp, StorageImp>::
  evaluateSingle(int baseFunct, 
                 const DomainType& xLocal,
                 const RangeType& factor) const 
  {
    //std::cout << "VecBaseFunctionSet::evaluateSingle" << std::endl;
    storage_.evaluate(util_.containedDof(baseFunct), diffVar0_, xLocal, tmp_);
    return factor[util_.component(baseFunct)]*tmp_[0];
  }
  
  template <class FunctionSpaceImp, template <class> class StorageImp>
  template <class QuadratureType>
  typename VecBaseFunctionSet<FunctionSpaceImp, StorageImp>::DofType
  VecBaseFunctionSet<FunctionSpaceImp, StorageImp>::
  evaluateSingle(int baseFunct,
                 const QuadratureType& quad, int quadPoint,
                 const RangeType& factor) const 
  {
    //std::cout << "VecBaseFunctionSet::evaluateSingle" << std::endl;
    storage_.evaluate(util_.containedDof(baseFunct), diffVar0_, quad, 
                      quadPoint, tmp_);
    return factor[util_.component(baseFunct)]*tmp_[0];
  }

  template <class FunctionSpaceImp, template <class> class StorageImp>
  template <class Entity>
  typename VecBaseFunctionSet<FunctionSpaceImp, StorageImp>::DofType
  VecBaseFunctionSet<FunctionSpaceImp, StorageImp>::
  evaluateGradientSingle(int baseFunct,
                         Entity& en,
                         const DomainType& xLocal,
                         const JacobianRangeType& factor) const 
  {
    DomainType gradScaled(0.);

    storage_.jacobian(util_.containedDof(baseFunct), xLocal, jTmp_);

    en.geometry().jacobianInverseTransposed(xLocal).umv(jTmp_[0], gradScaled);
    return gradScaled*factor[util_.component(baseFunct)];
  }

  template <class FunctionSpaceImp, template <class> class StorageImp>
  template <class Entity, class QuadratureType>
  typename VecBaseFunctionSet<FunctionSpaceImp, StorageImp>::DofType
  VecBaseFunctionSet<FunctionSpaceImp, StorageImp>::
  evaluateGradientSingle(int baseFunct,
                         Entity& en,
                         const QuadratureType& quad, int quadPoint,
                         const JacobianRangeType& factor) const 
  {
    storage_.jacobian(util_.containedDof(baseFunct), quad, quadPoint, jTmp_);

    DomainType gradScaled(0.);
    en.geometry().jacobianInverseTransposed(quad.point(quadPoint)).
      umv(jTmp_[0], gradScaled);
    return gradScaled*factor[util_.component(baseFunct)];
  }
} // end namespace Dune

