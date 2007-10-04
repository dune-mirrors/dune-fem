
namespace Dune {

  /*
  template <class A, int i> 
  struct SKPMeta
  {
    template <class B>
    static double skp(const A& a, const B& b)
    {
      return SKPMeta<A,i-1>::skp(a,b) + a[i]*b[i];
    }
  };
  
  template <class T> 
  struct SKPMeta<T,0>
  {
    template <class B>
    static double skp(const T& a, const B& b)
    {
      return a[0]*b[0];
    }
  };
  */
 
#if 0
  //- class StandardBaseFunctionSet
  template <class FunctionSpaceImp, template <class> class StorageImp>
  template <class QuadratureType>
  inline typename StandardBaseFunctionSet<FunctionSpaceImp, StorageImp>::DofType
  StandardBaseFunctionSet<FunctionSpaceImp, StorageImp>::
  evaluateSingle(int baseFunct,
                 const QuadratureType& quad, int quadPoint,
                 const RangeType& factor) const 
  {
    storage_.evaluate(baseFunct, diffVar0_, quad, quadPoint, tmp_);
    //return SKPMeta<RangeType,RangeType::dimension-1>::skp(tmp_,factor);
    return tmp_*factor;
  }
#endif

  template <class FunctionSpaceImp, template <class> class StorageImp>
  template <class Entity, class QuadratureType>
  inline typename StandardBaseFunctionSet<FunctionSpaceImp, StorageImp>::DofType
  StandardBaseFunctionSet<FunctionSpaceImp, StorageImp>::
  evaluateGradientSingle(const int baseFunct,
                         const Entity& en,
                         const QuadratureType& quad, 
                         const int quadPoint,
                         const JacobianRangeType& factor) const 
  {
    storage_.jacobian(baseFunct, quad, quadPoint, jTmp_);

    DofType result = 0.0;
    
    typedef FieldMatrix<DofType, FunctionSpaceImp::DimDomain,FunctionSpaceImp::DimDomain> JacobianInverseType;
    const JacobianInverseType& jti =
      en.geometry().jacobianInverseTransposed(quad.point(quadPoint));

    for (int i = 0; i < FunctionSpaceImp::DimRange; ++i) {
      DomainType gradScaled(0.);
      jti.umv(jTmp_[i], gradScaled);
      //result += SKPMeta<DomainType,DomainType::dimension-1>::skp(gradScaled,factor[i]);
      result += gradScaled*factor[i];
    }
    return result;
  }
  
  //- class VectorialBaseFunctionSet
  template <class FunctionSpaceImp, template <class> class StorageImp>
  inline int VectorialBaseFunctionSet<FunctionSpaceImp, StorageImp>::
  numBaseFunctions() const 
  {
    return storage_.numBaseFunctions()*FunctionSpaceImp::DimRange;
  }
  template <class FunctionSpaceImp, template <class> class StorageImp>
  inline int VectorialBaseFunctionSet<FunctionSpaceImp, StorageImp>::
  numDifferentBaseFunctions() const 
  {
    return storage_.numBaseFunctions();
  }

  template <class FunctionSpaceImp, template <class> class StorageImp>
  template <class QuadratureType>
  inline void VectorialBaseFunctionSet<FunctionSpaceImp, StorageImp>::
  evaluateScalar(int baseFunct, 
                 const QuadratureType & quad, int quadPoint, 
                 ScalarRangeType& phi) const
  {
    assert(baseFunct >= 0 && 
           baseFunct < numDifferentBaseFunctions());
    storage_.evaluate(baseFunct, diffVar0_, quad, quadPoint,
                      phi);
  }
  template <class FunctionSpaceImp, template <class> class StorageImp>
  inline void VectorialBaseFunctionSet<FunctionSpaceImp, StorageImp>::
  evaluateScalar(int baseFunct, 
                 const DomainType& xLocal, 
                 ScalarRangeType& phi) const
  {
    assert(baseFunct >= 0 && 
           baseFunct < numDifferentBaseFunctions());
    storage_.evaluate(baseFunct, diffVar0_, xLocal,
                      phi);
  }

  template <class FunctionSpaceImp, template <class> class StorageImp>
  template <int diffOrd>
  inline void VectorialBaseFunctionSet<FunctionSpaceImp, StorageImp>::
  evaluate(int baseFunct,
           const FieldVector<int, diffOrd>& diffVar,
           const DomainType& xLocal,
           RangeType& phi) const 
  {
    ScalarRangeType tmp;
    storage_.evaluate(util_.containedDof(baseFunct), diffVar, xLocal, tmp);
    
    phi = 0.0;
    phi[util_.component(baseFunct)] = tmp[0];
  }

  template <class FunctionSpaceImp, template <class> class StorageImp>
  template <int diffOrd, class QuadratureType>
  inline void VectorialBaseFunctionSet<FunctionSpaceImp, StorageImp>::
  evaluate(int baseFunct,
           const FieldVector<int, diffOrd>& diffVar,
           QuadratureType& quad, int quadPoint,
           RangeType& phi) const 
  {
    ScalarRangeType tmp;
    storage_.evaluate(util_.containedDof(baseFunct), diffVar, quad, quadPoint,
                      tmp);

    phi = 0.;
    phi[util_.component(baseFunct)] = tmp[0];
  }

  template <class FunctionSpaceImp, template <class> class StorageImp>
  inline void VectorialBaseFunctionSet<FunctionSpaceImp, StorageImp>::
  jacobianScalar(int baseFunct, const DomainType& xLocal, 
           ScalarJacobianRangeType& gradPhi) const 
  {
    assert(baseFunct >= 0 && 
           baseFunct < numDifferentBaseFunctions());
     storage_.jacobian(baseFunct, xLocal, gradPhi); 
  }

  template <class FunctionSpaceImp, template <class> class StorageImp>
  template <class QuadratureImp>
  inline void VectorialBaseFunctionSet<FunctionSpaceImp, StorageImp>::
  jacobianScalar(int baseFunct, QuadratureImp& quad, int quadPoint, 
                 ScalarJacobianRangeType& gradPhi) const 
  {
    assert(baseFunct >= 0 && 
           baseFunct < numDifferentBaseFunctions());
     storage_.jacobian(baseFunct, quad, quadPoint, gradPhi); 
  }
  template <class FunctionSpaceImp, template <class> class StorageImp>
  inline void VectorialBaseFunctionSet<FunctionSpaceImp, StorageImp>::
  jacobian(int baseFunct, const DomainType& xLocal, 
           JacobianRangeType& gradPhi) const 
  {
    gradPhi = 0.;
    
    storage_.jacobian(util_.containedDof(baseFunct), xLocal, jTmp_); 
    gradPhi[util_.component(baseFunct)] = jTmp_[0];
  }

  template <class FunctionSpaceImp, template <class> class StorageImp>
  template <class QuadratureImp>
  inline void VectorialBaseFunctionSet<FunctionSpaceImp, StorageImp>::
  jacobian(int baseFunct, QuadratureImp& quad, int quadPoint, 
           JacobianRangeType& gradPhi) const 
  {
    gradPhi = 0.;
    storage_.jacobian(util_.containedDof(baseFunct), quad, quadPoint, jTmp_); 
    gradPhi[util_.component(baseFunct)] = jTmp_[0];
  }

#if 0
  template <class FunctionSpaceImp, template <class> class StorageImp>
  inline typename VectorialBaseFunctionSet<FunctionSpaceImp, StorageImp>::DofType
  VectorialBaseFunctionSet<FunctionSpaceImp, StorageImp>::
  evaluateSingle(const int baseFunct, 
                 const DomainType& xLocal,
                 const RangeType& factor) const 
  {
    //std::cout << "VectorialBaseFunctionSet::evaluateSingle" << std::endl;
    storage_.evaluate(util_.containedDof(baseFunct), diffVar0_, xLocal, tmp_);
    return factor[util_.component(baseFunct)]*tmp_[0];
  }
  
  template <class FunctionSpaceImp, template <class> class StorageImp>
  template <class QuadratureType>
  inline typename VectorialBaseFunctionSet<FunctionSpaceImp, StorageImp>::DofType
  VectorialBaseFunctionSet<FunctionSpaceImp, StorageImp>::
  evaluateSingle(const int baseFunct,
                 const QuadratureType& quad, const int quadPoint,
                 const RangeType& factor) const 
  {
    //std::cout << "VectorialBaseFunctionSet::evaluateSingle" << std::endl;
    storage_.evaluate(util_.containedDof(baseFunct), diffVar0_, quad, 
                      quadPoint, tmp_);
    return factor[util_.component(baseFunct)]*tmp_[0];
  }
#endif

  template <class FunctionSpaceImp, template <class> class StorageImp>
  template <class Entity>
  inline typename VectorialBaseFunctionSet<FunctionSpaceImp, StorageImp>::DofType
  VectorialBaseFunctionSet<FunctionSpaceImp, StorageImp>::
  evaluateGradientSingle(const int baseFunct,
                         const Entity& en,
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
  inline typename VectorialBaseFunctionSet<FunctionSpaceImp, StorageImp>::DofType
  VectorialBaseFunctionSet<FunctionSpaceImp, StorageImp>::
  evaluateGradientSingle(const int baseFunct,
                         const Entity& en,
                         const QuadratureType& quad, int quadPoint,
                         const JacobianRangeType& factor) const 
  {
    storage_.jacobian(util_.containedDof(baseFunct), quad, quadPoint, jTmp_);

    DomainType gradScaled(0.);
    en.geometry().jacobianInverseTransposed(quad.point(quadPoint)).
      umv(jTmp_[0], gradScaled);
    //return SKPMeta<DomainType,DomainType::dimension-1>::skp(gradScaled,factor[util_.component(baseFunct)]);
    return gradScaled*factor[util_.component(baseFunct)];
  }
} // end namespace Dune

