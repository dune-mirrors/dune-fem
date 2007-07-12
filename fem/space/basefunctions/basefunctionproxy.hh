#ifndef DUNE_BASEFUNCTIONPROXY_HH
#define DUNE_BASEFUNCTIONPROXY_HH

namespace Dune{

template < class BaseFunctionSetImp > 
class SimpleBaseFunctionProxy; 

template <class BaseFunctionSetImp>
struct SimpleBaseFunctionProxyTraits  
{
  //! Export function space type
  typedef typename BaseFunctionSetImp :: Traits :: FunctionSpaceType FunctionSpaceType;

  //! Exact type of the base function
  typedef SimpleBaseFunctionProxy<BaseFunctionSetImp> BaseFunctionSetType;
};

//! BaseFunctionSetProxy Object 
template < class BaseFunctionSetImp > 
class SimpleBaseFunctionProxy
: public BaseFunctionSetDefault< SimpleBaseFunctionProxyTraits< BaseFunctionSetImp > >
{
public:
  typedef SimpleBaseFunctionProxyTraits< BaseFunctionSetImp > Traits;
  
private:
  typedef SimpleBaseFunctionProxy< BaseFunctionSetImp > ThisType;
  typedef BaseFunctionSetDefault< Traits > BaseType;

public:
  typedef typename Traits :: FunctionSpaceType FunctionSpaceType;

  enum { dimrange = FunctionSpaceType::DimRange };
  enum { dimRange = FunctionSpaceType::DimRange };
  enum { DimRange = FunctionSpaceType::DimRange };
  enum { dimDomain = FunctionSpaceType::DimDomain };
  enum { DimDomain = FunctionSpaceType::DimDomain };

  typedef typename FunctionSpaceType::RangeFieldType RangeFieldType;
  typedef typename FunctionSpaceType::DomainType DomainType;
  typedef typename FunctionSpaceType::RangeType RangeType;
  typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;

protected:
  // base function set 
  const BaseFunctionSetImp* baseSet_; 

public:
  inline SimpleBaseFunctionProxy ()
  : baseSet_( NULL )
  {
  }
    
  inline SimpleBaseFunctionProxy ( const BaseFunctionSetImp* baseSet )
  : baseSet_( baseSet ) 
  {
  }

  //! Constructor creating empty local function 
  inline SimpleBaseFunctionProxy( const ThisType& org )
  : baseSet_( org.baseSet_ )
  {
  }

  //! asignment operator 
  SimpleBaseFunctionProxy &operator= ( const ThisType& org )
  {
    baseSet_ = org.baseSet_; 
    return *this;
  }   

  //! destructor 
  ~SimpleBaseFunctionProxy() 
  { 
    // unset pointer to be sure its gone 
    baseSet_ = 0;
  }

  //! number of base functions
  inline
  int numBaseFunctions() const 
  {
    assert( this->baseSet_ );
    return baseSet_->numBaseFunctions(); 
  }
  
  //! @copydoc BaseFunctionSetInterface::evaluate 
  inline
  void evaluate(int baseFunct,
                const DomainType& xLocal,
                RangeType& phi) const
  {
    assert( this->baseSet_ );
    baseSet_->evaluate(baseFunct,
                               xLocal,
                               phi);
  }

  //! @copydoc BaseFunctionSetInterface::evaluate 
  template <int diffOrd>
  inline
  void evaluate(int baseFunct,
                const FieldVector<int, diffOrd>& diffVar,
                const DomainType& xLocal,
                RangeType& phi) const
  {
    assert( this->baseSet_ );
    baseSet_->evaluate(baseFunct,
                               diffVar,
                               xLocal,
                               phi);
  }

  //! @copydoc BaseFunctionSetInterface::evaluate 
  template <class QuadratureType>
  inline
  void evaluate (const int baseFunct, 
                 const QuadratureType & quad, 
                 const int quadPoint, 
                 RangeType & phi) const
  {
    assert( this->baseSet_ );
    baseSet_->evaluate(baseFunct,
                               quad,
                               quadPoint,
                               phi);
  }

  //! @copydoc BaseFunctionSetInterface::evaluate 
  template <int diffOrd, class QuadratureType>
  inline
  void evaluate (const int baseFunct, 
                 const FieldVector<int, diffOrd> &diffVariable, 
                 const QuadratureType & quad, 
                 const int quadPoint, 
                 RangeType & phi) const
  {
    assert( this->baseSet_ );
    baseSet_->evaluate(baseFunct,
                               diffVariable,
                               quad,
                               quadPoint,
                               phi);
  }

  //! @copydoc BaseFunctionSetDefault::evaluateSingle 
  template <class QuadratureType>
  inline
  RangeFieldType evaluateSingle(const int baseFunct,
                                const QuadratureType& quad, const int quadPoint,
                                const RangeType& factor) const
  {
    assert( this->baseSet_ );
    return baseSet_->evaluateSingle(baseFunct,
                                            quad,
                                            quadPoint,
                                            factor);
  }
    
  //! @copydoc BaseFunctionSetDefault::evaluateGradientSingle 
  template <class Entity, class QuadratureType>
  inline
  RangeFieldType evaluateGradientSingle(const int baseFunct,
                                        const Entity& en,
                                        const QuadratureType& quad, 
                                        const int quadPoint,
                                 const JacobianRangeType& factor) const
  {
    assert( this->baseSet_ );
    return baseSet_->
              evaluateGradientSingle(baseFunct,
                                     en,
                                     quad,
                                     quadPoint,
                                     factor);
  }


  //! @copydoc BaseFunctionSetDefault::jacobian
  inline
  void jacobian(int baseFunct, const DomainType& xLocal,
                JacobianRangeType& gradPhi) const
  {
    assert( this->baseSet_ );
    baseSet_->jacobian(baseFunct,
                               xLocal,
                               gradPhi);
  }

  //! @copydoc BaseFunctionSetDefault::jacobian 
  template <class QuadratureImp>
  inline
  void jacobian(int baseFunct, QuadratureImp& quad, int quadPoint,
                JacobianRangeType& gradPhi) const
  {
    assert( this->baseSet_ );
    baseSet_->jacobian(baseFunct,
                               quad,
                               quadPoint,
                               gradPhi);
  }
}; // end SimpleBaseFunctionProxy 

//! VectorialBaseFunctionSetProxy 
template < class BaseFunctionSetImp > 
class VectorialBaseFunctionProxy : public 
  SimpleBaseFunctionProxy< BaseFunctionSetImp > 
{
  typedef SimpleBaseFunctionProxy< BaseFunctionSetImp > BaseType;
    
public:  
  typedef typename BaseType :: Traits Traits;
  typedef typename Traits :: FunctionSpaceType FunctionSpaceType;

  typedef typename BaseFunctionSetImp :: ScalarRangeType ScalarRangeType;
  typedef typename BaseFunctionSetImp :: ScalarJacobianRangeType ScalarJacobianRangeType;

  enum { dimrange = FunctionSpaceType::DimRange };
  enum { dimRange = FunctionSpaceType::DimRange };
  enum { DimRange = FunctionSpaceType::DimRange };
  enum { dimDomain = FunctionSpaceType::DimDomain };
  enum { DimDomain = FunctionSpaceType::DimDomain };

  typedef typename FunctionSpaceType::RangeFieldType RangeFieldType;
  typedef typename FunctionSpaceType::DomainType DomainType;
  typedef typename FunctionSpaceType::RangeType RangeType;
  typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;

  //! constructor storing base function set pointer 
  VectorialBaseFunctionProxy(const BaseFunctionSetImp* b) 
    : BaseType(b) 
  {
  }

  //! Constructor creating empty local function 
  VectorialBaseFunctionProxy(const VectorialBaseFunctionProxy& org) 
    : BaseType(org) 
  {}

  //! asignment operator 
  VectorialBaseFunctionProxy& operator = (const VectorialBaseFunctionProxy& org) 
  {
    BaseType :: operator = (org);
    return *this;
  }   

  //! number of different base functions
  inline
  int numDifferentBaseFunctions() const 
  {
    assert( this->baseSet_ );
    return this->baseSet_->numDifferentBaseFunctions(); 
  }

  template <class QuadratureType>
  inline
  void evaluateScalar(const int baseFunct,
    const QuadratureType & quad, const int p,
    ScalarRangeType& phi) const
  {
    assert( this->baseSet_ );
    this->baseSet_->evaluateScalar(baseFunct,quad,p,phi);
  }
  
  inline
  void evaluateScalar(const int baseFunct,
    const DomainType& xLocal,
    ScalarRangeType& phi) const
  {
    assert( this->baseSet_ );
    this->baseSet_->evaluateScalar(baseFunct,xLocal,phi);
  }

  inline
  void jacobianScalar(const int baseFunct, const DomainType& xLocal,
                      ScalarJacobianRangeType& gradPhi) const
  {
    assert( this->baseSet_ );
    this->baseSet_->jacobianScalar(baseFunct,xLocal,gradPhi);
  }

  template <class QuadratureImp>
  inline
  void jacobianScalar(const int baseFunct, 
                      const QuadratureImp& quad, 
                      const int quadPoint,
                      ScalarJacobianRangeType& gradPhi) const
  {
    assert( this->baseSet_ );
    this->baseSet_->jacobianScalar(baseFunct,quad,quadPoint,gradPhi);
  } 
};

} // end namespace Dune
#endif
