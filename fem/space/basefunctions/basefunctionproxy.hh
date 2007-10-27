#ifndef DUNE_BASEFUNCTIONPROXY_HH
#define DUNE_BASEFUNCTIONPROXY_HH

namespace Dune
{

  /** \addtogroup BaseFunction
   *  \{
   */
    
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

  /** \class SimpleBaseFunctionProxy
   *  \brief proxy object for base function sets
   */
  template< class BaseFunctionSetImp > 
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
    const BaseFunctionSetImp *baseSet_; 

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

    /** \copydoc Dune::BaseFunctionSetInterface::numBaseFunctions */
    inline int numBaseFunctions() const 
    {
      return baseFunctionSet().numBaseFunctions(); 
    }
    
    /** \çopydoc Dune::BaseFunctionSetInterface::geometryType */
    inline GeometryType geometryType () const
    {
      return baseFunctionSet().geometryType();
    }
   
    /** \copydoc Dune::BaseFunctionSetInterface::evaluate(const int baseFunction,const FieldVector<deriType,diffOrd> &diffVariable,const DomainType &x,RangeType &phi) const */
    template< int diffOrd, class PointType >
    inline void evaluate ( const int baseFunction,
                           const FieldVector< deriType, diffOrd > &diffVariable,
                           const PointType &x,
                           RangeType &phi ) const
    {
      baseFunctionSet().evaluate( baseFunction, diffVariable, x, phi );
    }
    
    /** \copydoc Dune::BaseFunctionSetInterface::evaluate(const int baseFunction,const FieldVector<deriType,diffOrd> &diffVariable,const QuadratureType &quadrature,const int quadPoint,RangeType &phi) const */
    template< int diffOrd, class QuadratureType >
    inline void evaluate ( const int baseFunction,
                           const FieldVector< deriType, diffOrd > &diffVariable,
                           const QuadratureType &quadrature,
                           const int quadPoint,
                           RangeType &phi ) const
    {
      baseFunctionSet().evaluate( baseFunction, diffVariable, quadrature, quadPoint, phi );
    }

    /** \copydoc Dune::BaseFunctionSetInterface::evaluate(const int baseFunction,const DomainType &x,RangeType &phi) const */
    template< class PointType >
    inline void evaluate ( const int baseFunction,
                           const PointType &x,
                           RangeType &phi ) const
    {
      baseFunctionSet().evaluate( baseFunction, x, phi );
    }

    /** \copydoc Dune::BaseFunctionSetInterface::evaluate(const int baseFunction,const QuadratureType &quadrature,const int quadPoint,RangeType &phi) const */
    template< class QuadratureType >
    inline void evaluate ( const int baseFunction,
                           const QuadratureType &quadrature,
                           const int quadPoint,
                           RangeType &phi ) const
    {
      baseFunctionSet().evaluate( baseFunction, quadrature, quadPoint, phi );
    }
    
    /** \copydoc Dune::BaseFunctionSetDefault::jacobian(const int baseFunction,const DomainType &x,JacobianRangeType &phi) const */
    template< class PointType >
    inline void jacobian( const int baseFunction,
                          const PointType &x,
                          JacobianRangeType &phi ) const
    {
      baseFunctionSet().jacobian( baseFunction, x, phi );
    }

    /** \copydoc Dune::BaseFunctionSetDefault::jacobian(const int baseFunction,const QuadratureType &quadrature,const int quadPoint,JacobianRangeType &phi) const */
    template< class QuadratureType >
    inline void jacobian ( const int baseFunction,
                           const QuadratureType &quadrature,
                           const int quadPoint,
                           JacobianRangeType &phi ) const
    {
      baseFunctionSet().jacobian( baseFunction, quadrature, quadPoint, phi );
    }

    /** \copydoc Dune::BaseFunctionSetInterface::evaluateSingle(const int baseFunction,const DomainType &x,const RangeType &psi) const */
    template< class PointType >
    inline RangeFieldType evaluateSingle ( const int baseFunction,
                                           const PointType &x,
                                           const RangeType &psi ) const
    {
      return baseFunctionSet().evaluateSingle( baseFunction, x, psi );
    }

    /** \copydoc Dune::BaseFunctionSetInterface::evaluateSingle(const int baseFunction,const QuadratureType &quadrature,const int quadPoint,const RangeType &psi) const */
    template< class QuadratureType >
    inline RangeFieldType evaluateSingle ( const int baseFunction,
                                           const QuadratureType &quadrature,
                                           const int quadPoint,
                                           const RangeType &psi ) const
    {
      return baseFunctionSet().evaluateSingle
        ( baseFunction, quadrature, quadPoint, psi );
    }
    
     /** \copydoc Dune::BaseFunctionSetInterface::evaluateGradientSingle(const int baseFunction,const EntityType &entity,const DomainType &x,const JacobianRangeType &psi) const */
    template< class EntityType, class PointType >
    inline RangeFieldType evaluateGradientSingle ( const int baseFunction,
                                                   const EntityType &entity,
                                                   const PointType &x,
                                                   const JacobianRangeType &psi ) const
    {
      return baseFunctionSet().evaluateGradientSingle( baseFunction, entity, x, psi );
    }
     
    /** \copydoc Dune::BaseFunctionSetInterface::evaluateGradientSingle(const int baseFunction,const EntityType &entity,const QuadratureType &quadrature,const int quadPoint,const JacobianRangeType &psi) const */
    template< class EntityType, class QuadratureType >
    inline RangeFieldType evaluateGradientSingle ( const int baseFunction,
                                                   const EntityType &entity,
                                                   const QuadratureType &quadrature,
                                                   const int quadPoint,
                                                   const JacobianRangeType &psi ) const
    {
      return baseFunctionSet().evaluateGradientSingle
        ( baseFunction, entity, quadrature, quadPoint, psi );
    }
   
  protected:
    inline const BaseFunctionSetImp &baseFunctionSet () const
    {
      assert( this->baseSet_ );
      return *baseSet_;
    }
  }; // end SimpleBaseFunctionProxy 



  /** \class VectorialBaseFunctionSetProxy
   */
  template< class BaseFunctionSetImp > 
  class VectorialBaseFunctionProxy
  : public SimpleBaseFunctionProxy< BaseFunctionSetImp >
  {
  private:
    typedef VectorialBaseFunctionProxy< BaseFunctionSetImp > ThisType;
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

  protected:
    using BaseType :: baseFunctionSet;

  public:
    //! default constructor 
    inline VectorialBaseFunctionProxy()
    : BaseType () 
    {
    }
      
    //! constructor storing base function set pointer 
    inline VectorialBaseFunctionProxy( const BaseFunctionSetImp *baseSet )
    : BaseType( baseSet ) 
    {
    }

    //! Constructor creating empty local function 
    inline VectorialBaseFunctionProxy( const VectorialBaseFunctionProxy &org )
    : BaseType( org )
    {
    }

    //! asignment operator 
    inline VectorialBaseFunctionProxy& operator = (const VectorialBaseFunctionProxy& org) 
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

    template< class PointType >
    inline void evaluateScalar ( const int baseFunction,
                                 const PointType &x,
                                 ScalarRangeType &phi ) const
    {
      baseFunctionSet().evaluateScalar( baseFunction, x, phi );
    }

    template< class PointType >
    inline void jacobianScalar ( const int baseFunction,
                                 const PointType &x,
                                 ScalarJacobianRangeType &gradPhi ) const
    {
      baseFunctionSet().jacobianScalar( baseFunction, x, gradPhi );
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

  /** \} */

} // end namespace Dune

#endif
