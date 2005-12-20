#ifndef DUNE_NEWBASE_HH
#define DUNE_NEWBASE_HH

#include <dune/fem/common/basefunctions.hh>
#include <dune/fem/common/basefunctionfactory.hh>
#include <dune/fem/space/dofstorage.hh>
#include <dune/common/fvector.hh>


namespace Dune {

  // Forward declarations
  template <class FunctionSpaceImp, template <class> class StorageImp>
  class StandardBaseFunctionSet;
  template <class FunctionSpaceImp, template <class> class StorageImp>
  class VectorialBaseFunctionSet;

  //! Traits class for standard base function set
  template <class FunctionSpaceImp, template <class> class StorageImp>
  struct StandardBaseFunctionSetTraits
  {
    //! Export function space type
    typedef FunctionSpaceImp FunctionSpaceType;
    //! Type of the base function storage policy
    typedef StorageImp<FunctionSpaceType> StorageType;
    //! Exact type of the base function
    typedef StandardBaseFunctionSet<FunctionSpaceType, 
                                    StorageImp> BaseFunctionSetType;
    //! Factory type for the corresponding base functions (polymorphic)
    typedef BaseFunctionFactory<FunctionSpaceType> FactoryType;

  };

  //! \brief Standard base function set
  template <class FunctionSpaceImp, template <class> class StorageImp>
  class StandardBaseFunctionSet : 
    public BaseFunctionSetDefault<StandardBaseFunctionSetTraits<FunctionSpaceImp, StorageImp> >
  {
  public:
    typedef StandardBaseFunctionSetTraits<FunctionSpaceImp, StorageImp> Traits;
    typedef typename FunctionSpaceImp::DomainType DomainType;
    typedef typename FunctionSpaceImp::RangeType RangeType;
    typedef typename FunctionSpaceImp::RangeFieldType DofType;
    typedef typename FunctionSpaceImp::JacobianRangeType JacobianRangeType;

  public:
    //! Constructor
    StandardBaseFunctionSet(const typename Traits::FactoryType& factory) :
      storage_(factory),
      diffVar0_(0),
      tmp_(0.),
      jTmp_(0.)
    {}

    //! Total number of base functions
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

    template <class QuadratureType>
    DofType evaluateSingle(int baseFunct,
                           const QuadratureType& quad, int quadPoint,
                           const RangeType& factor) const;
      
    template <class Entity, class QuadratureType>
    DofType evaluateGradientSingle(int baseFunct,
                                   Entity& en,
                                   const QuadratureType& quad, int quadPoint,
                                   const JacobianRangeType& factor) const;
  private:
    typename Traits::StorageType storage_;
    
    mutable FieldVector<int, 0> diffVar0_;
    mutable RangeType tmp_;
    mutable JacobianRangeType jTmp_;
  };

  //- VectorialBaseFunctionSet
  template <class FunctionSpaceImp, template <class> class StorageImp>
  struct VectorialBaseFunctionSetTraits 
  {
    typedef FunctionSpaceImp FunctionSpaceType;
    typedef VectorialBaseFunctionSet<
      FunctionSpaceType, StorageImp> BaseFunctionSetType;
  };

  //! \brief Special base function implementation that takes advantage
  //! of the vectorial structure of the base functions.
  //! This base function can be used in conjunction with scalar basefunctions
  //! \phi_i which are extended to vectorial base functions like \Phi_j = 
  //! \phi_i e_k, where e_k = [ \kronecker_ik ]_i.
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
    typedef typename FunctionSpaceType::RangeFieldType RangeFieldType;
    typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;
    
    typedef typename FunctionSpaceType::RangeFieldType DofType;
        
  public:
    //! Constructor
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

    DofType evaluateSingle(int baseFunct, 
                           const DomainType& xLocal,
                           const RangeType& factor) const;
    
    template <class QuadratureType>
    DofType evaluateSingle(int baseFunct,
                           const QuadratureType& quad, int quadPoint,
                           const RangeType& factor) const;
      
    template <class Entity>
    DofType evaluateGradientSingle(int baseFunct,
                                   Entity& en,
                                   const DomainType& xLocal,
                                   const JacobianRangeType& factor) const;
    
    template <class Entity, class QuadratureType>
    DofType evaluateGradientSingle(int baseFunct,
                                   Entity& en,
                                   const QuadratureType& quad, int quadPoint,
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
