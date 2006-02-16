#ifndef DUNE_NEWBASE_HH
#define DUNE_NEWBASE_HH

#include <dune/fem/common/basefunctions.hh>
#include <dune/fem/common/basefunctionfactory.hh>
#include <dune/fem/space/dofstorage.hh>
#include <dune/common/fvector.hh>


namespace Dune {

  // Forward declarations
  template <class FunctionSpaceImp, template <class> class StorageImp>
  class StdBaseFunctionSet;
  template <class FunctionSpaceImp, template <class> class StorageImp>
  class VecBaseFunctionSet;

  //! Traits class for standard base function set
  template <class FunctionSpaceImp, template <class> class StorageImp>
  struct StdBaseFunctionSetTraits
  {
    //! Export function space type
    typedef FunctionSpaceImp FunctionSpaceType;
    //! Type of the base function storage policy
    typedef StorageImp<FunctionSpaceType> StorageType;
    //! Exact type of the base function
    typedef StdBaseFunctionSet<FunctionSpaceType, 
                                    StorageImp> BaseFunctionSetType;
    //! Factory type for the corresponding base functions (polymorphic)
    typedef BaseFunctionFactory<FunctionSpaceType> FactoryType;

  };

  //! \brief Standard base function set
  template <class FunctionSpaceImp, template <class> class StorageImp>
  class StdBaseFunctionSet : 
    public BaseFunctionSetDefault<StdBaseFunctionSetTraits<FunctionSpaceImp, StorageImp> >
  {
  public:
    typedef StdBaseFunctionSetTraits<FunctionSpaceImp, StorageImp> Traits;
    typedef typename FunctionSpaceImp::DomainType DomainType;
    typedef typename FunctionSpaceImp::RangeType RangeType;
    typedef typename FunctionSpaceImp::RangeFieldType DofType;
    typedef typename FunctionSpaceImp::JacobianRangeType JacobianRangeType;

  public:
    //! Constructor
    StdBaseFunctionSet(const typename Traits::FactoryType& factory) :
      storage_(factory),
      diffVar0_(0),
      tmp_(0.),
      jTmp_(0.)
    {}

    //! Total number of base functions
    inline
    int numBaseFunctions() const;

    template <int diffOrd>
    inline
    void evaluate(int baseFunct,
                  const FieldVector<int, diffOrd>& diffVar,
                  const DomainType& xLocal,
                  RangeType& phi) const;

    template <int diffOrd, class QuadratureType>
    inline
    void evaluate (int baseFunct, 
                   const FieldVector<int, diffOrd> &diffVariable, 
                   QuadratureType & quad, 
                   int quadPoint, RangeType & phi) const;

    template <class QuadratureType>
    inline
    DofType evaluateSingle(int baseFunct,
                           const QuadratureType& quad, int quadPoint,
                           const RangeType& factor) const;
      
    template <class Entity, class QuadratureType>
    inline
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

  //- VecBaseFunctionSet
  template <class FunctionSpaceImp, template <class> class StorageImp>
  struct VecBaseFunctionSetTraits 
  {
    typedef FunctionSpaceImp FunctionSpaceType;
    typedef VecBaseFunctionSet<
      FunctionSpaceType, StorageImp> BaseFunctionSetType;
  };

  //! \brief Special base function implementation that takes advantage
  //! of the vectorial structure of the base functions.
  //! This base function can be used in conjunction with scalar basefunctions
  //! \phi_i which are extended to vectorial base functions like \Phi_j = 
  //! \phi_i e_k, where e_k = [ \kronecker_ik ]_i.
  template <class FunctionSpaceImp, template <class> class StorageImp>
  class VecBaseFunctionSet : 
    public BaseFunctionSetDefault<VecBaseFunctionSetTraits<FunctionSpaceImp, StorageImp> >
  {
  private:
    typedef BaseFunctionSetDefault<
      VecBaseFunctionSetTraits<FunctionSpaceImp, StorageImp> > BaseType;
    typedef FunctionSpaceImp FunctionSpaceType;
    typedef typename ToScalarFunctionSpace<
      FunctionSpaceImp>::Type ScalarFunctionSpaceType;
    typedef typename ScalarFunctionSpaceType::RangeType ScalarRangeType;
    typedef typename ScalarFunctionSpaceType::JacobianRangeType 
      ScalarJacobianRangeType;
    typedef BaseFunctionFactory<ScalarFunctionSpaceType> FactoryType;
    typedef typename FactoryType::BaseFunctionType BaseFunctionType;
    typedef StorageImp<ScalarFunctionSpaceType> StorageType;
    typedef VecBaseFunctionSetTraits<FunctionSpaceImp,StorageImp> Traits;
 
  public:
    typedef typename FunctionSpaceType::DomainType DomainType;
    typedef typename FunctionSpaceType::RangeType RangeType;
    typedef typename FunctionSpaceType::RangeFieldType RangeFieldType;
    typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;
    
    typedef typename FunctionSpaceType::RangeFieldType DofType;
        
  public:
    //! Constructor
    VecBaseFunctionSet(const FactoryType& factory) :
      storage_(factory),
      util_(FunctionSpaceType::DimRange),
      tmp_(0.),
      jTmp_(0.)
    {}

    ~VecBaseFunctionSet() {}

    inline
    int numBaseFunctions() const;

    template <int diffOrd>
    inline
    void evaluate(int baseFunct,
                  const FieldVector<int, diffOrd>& diffVar,
                  const DomainType& xLocal,
                  RangeType& phi) const;

    template <int diffOrd, class QuadratureType>
    inline
    void evaluate(int baseFunct, 
                  const FieldVector<deriType, diffOrd> &diffVariable, 
                  QuadratureType & quad, 
                  int quadPoint, RangeType & phi ) const;

    inline
    void jacobian(int baseFunct, const DomainType& xLocal, 
                  JacobianRangeType& gradPhi) const;

    template <class QuadratureImp>
    inline
    void jacobian(int baseFunct, QuadratureImp& quad, int quadPoint,
                  JacobianRangeType& gradPhi) const;

    inline
    DofType evaluateSingle(int baseFunct, 
                           const DomainType& xLocal,
                           const RangeType& factor) const;
    
    template <class QuadratureType>
    inline
    DofType evaluateSingle(int baseFunct,
                           const QuadratureType& quad, int quadPoint,
                           const RangeType& factor) const;
      
    template <class Entity>
    inline
    DofType evaluateGradientSingle(int baseFunct,
                                   Entity& en,
                                   const DomainType& xLocal,
                                   const JacobianRangeType& factor) const;
    
    template <class Entity, class QuadratureType>
    inline
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
