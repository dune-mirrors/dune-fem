#ifndef DUNE_BASEFUNCTIONPROXY_HH
#define DUNE_BASEFUNCTIONPROXY_HH

#include <dune/fem/space/common/basefunctionsetinterface.hh>

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

    enum { dimRange = FunctionSpaceType :: dimRange };
    enum { dimDomain = FunctionSpaceType :: dimDomain };

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
    
    /** \copydoc Dune::BaseFunctionSetInterface::geometryType */
    inline GeometryType geometryType () const
    {
      return baseFunctionSet().geometryType();
    }
   
    /** \copydoc Dune::BaseFunctionSetInterface::evaluate(const int baseFunction,const FieldVector<deriType,diffOrd> &diffVariable,const PointType &x,RangeType &phi) const */
    template< int diffOrd, class PointType >
    inline void evaluate ( const int baseFunction,
                           const FieldVector< deriType, diffOrd > &diffVariable,
                           const PointType &x,
                           RangeType &phi ) const
    {
      baseFunctionSet().evaluate( baseFunction, diffVariable, x, phi );
    }

    /** \copydoc Dune::BaseFunctionSetInterface::evaluate(const int baseFunction,const PointType &x,RangeType &phi) const */
    template< class PointType >
    inline void evaluate ( const int baseFunction,
                           const PointType &x,
                           RangeType &phi ) const
    {
      baseFunctionSet().evaluate( baseFunction, x, phi );
    }

    /** \copydoc Dune::BaseFunctionSetDefault::jacobian(const int baseFunction,const PointType &x,JacobianRangeType &phi) const */
    template< class PointType >
    inline void jacobian( const int baseFunction,
                          const PointType &x,
                          JacobianRangeType &phi ) const
    {
      baseFunctionSet().jacobian( baseFunction, x, phi );
    }

    /** \copydoc Dune::BaseFunctionSetInterface::evaluateSingle(const int baseFunction,const PointType &x,const RangeType &psi) const */
    template< class PointType >
    inline RangeFieldType evaluateSingle ( const int baseFunction,
                                           const PointType &x,
                                           const RangeType &psi ) const
    {
      return baseFunctionSet().evaluateSingle( baseFunction, x, psi );
    }

     /** \copydoc Dune::BaseFunctionSetInterface::evaluateGradientSingle(const int baseFunction,const EntityType &entity,const PointType &x,const JacobianRangeType &psi) const */
    template< class EntityType, class PointType >
    inline RangeFieldType evaluateGradientSingle ( const int baseFunction,
                                                   const EntityType &entity,
                                                   const PointType &x,
                                                   const JacobianRangeType &psi ) const
    {
      return baseFunctionSet().evaluateGradientSingle( baseFunction, entity, x, psi );
    }
    
  protected:
    inline const BaseFunctionSetImp &baseFunctionSet () const
    {
      assert( this->baseSet_ );
      return *baseSet_;
    }
  }; // end SimpleBaseFunctionProxy 



  /** \class VectorialBaseFunctionSetProxy */
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

    enum { dimRange = FunctionSpaceType :: dimRange };
    enum { dimDomain = FunctionSpaceType :: dimDomain };

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

    /** \brief copy constructor
     *
     *  \param[in]  other  VectorialBaseFunctionProxy to copy
     */
    inline VectorialBaseFunctionProxy ( const ThisType &other )
    : BaseType( other )
    {
    }

    /** \brief assign another VectorialBaseFunctionProxy to this one
     *
     *  \param[in]  other  VectorialBaseFunctionProxy to copy
     */
    inline ThisType &operator= ( const ThisType &other ) 
    {
      BaseType :: operator= ( other );
      return *this;
    }   

    //! number of different base functions
    inline int numDifferentBaseFunctions () const 
    {
      return baseFunctionSet().numDifferentBaseFunctions(); 
    }
    
    template< int diffOrd, class PointType >
    inline void evaluateScalar ( const int baseFunction,
                                 const FieldVector< int, diffOrd > &diffVariable,
                                 const PointType &x,
                                 ScalarRangeType &phi ) const
    {
      baseFunctionSet().evaluateScalar( baseFunction, diffVariable, x, phi );
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
  };

  /** \} */

} // end namespace Dune

#endif
