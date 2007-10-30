#ifndef DUNE_FUNCTION_HH
#define DUNE_FUNCTION_HH

//- Dune includes 
#include <dune/common/bartonnackmanifcheck.hh>
#include <dune/common/fvector.hh>

//- local includes 
#include <dune/fem/operator/common/mapping.hh>

namespace Dune
{
  /** @addtogroup Functions
      Functions are mapping from one finite dimensional
      vector space into another, e.g.,
      \f$K^n\f$ into \f$L^m\f$.
      They are element of a 
      \ref FunctionSpaceInterface "function space".

      \remark The interface is given by Function.
      @{
  **/
    
  /*! \brief 
      Abstract class representing a function


      Template parameters are:
      -  FunctionSpaceImp      type of the function space where the function 
                               belongs to.
      -  FunctionImp           type of the implemented function (Barton-Nackman)
      
      @interfaceclass
  **/
  template< class FunctionSpaceImp, class FunctionImp >
  class Function
  :  public Mapping < typename FunctionSpaceImp :: DomainFieldType,
                      typename FunctionSpaceImp :: RangeFieldType,
                      typename FunctionSpaceImp :: DomainType, 
                      typename FunctionSpaceImp :: RangeType > 
  {
  public:
    //! type of function space this function belongs to 
    typedef FunctionSpaceImp FunctionSpaceType;

  private:
    typedef Function< FunctionSpaceType, FunctionImp > ThisType;
      
    // FunctionImp should be exported as FunctionType
    typedef Function<FunctionSpaceImp,FunctionImp> FunctionType;

  public:
    //! domain type (from function space)
    typedef typename FunctionSpaceType::DomainType DomainType ;
    //! range type (from function space)
    typedef typename FunctionSpaceType::RangeType RangeType ;
    //! jacobian type (from function space)
    typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;
    //! hessian type (from function space)
    typedef typename FunctionSpaceType::HessianRangeType  HessianRangeType;


  protected:
    const FunctionSpaceType &functionSpace_;

  public:
    //! constructor
    inline explicit Function ( const FunctionSpaceType &fSpace )
    : functionSpace_( fSpace )
    {
    }   

    //! copy constructor
    inline Function ( const ThisType &org )
    : functionSpace_( org.functionSpace_ )
    {
    }   

    //! destructor 
    virtual ~Function ()
    {
    }

  private:
    // Disallow copying
    ThisType &operator= ( const ThisType &other );
    
  public:
    /** \brief application operator call evaluate 
        \param[in] arg argument 
        \param[out] dest destination, i.e. f(arg) 
    */
    virtual void operator()(const DomainType & arg, RangeType & dest) const 
    {
      evaluate(arg,dest);
    }

    /** \brief evaluate the function
     *
     *  \param[in]  x    evaluation point 
     *  \param[out] ret  value of the function in x
     */
    inline void evaluate ( const DomainType &x,
                           RangeType &ret ) const
    {
      CHECK_AND_CALL_INTERFACE_IMPLEMENTATION( asImp().evaluate( x, ret ) );
    }

    /** \brief evaluate a derivative of the function
     * 
     *  \param[in]  diffVariable  vector describing the partial derivative to
     *                            evaluate
     *  \param[in]  x             evaluation point
     *  \param[out] ret           value of the derivative in x 
     */    
    template< int diffOrder >
    inline void evaluate ( const FieldVector< deriType, diffOrder > &diffVariable,
                           const DomainType &x,
                           RangeType &ret ) const
    {
      CHECK_AND_CALL_INTERFACE_IMPLEMENTATION
        ( asImp().evaluate( diffVariable, x, ret ) );
    }

    /** \brief Obtain the related function space
     * 
     *  \returns a reference to the function space 
     */ 
    inline const FunctionSpaceType &space () const
    {
      return functionSpace_;
    }

  protected:
    const FunctionImp &asImp () const
    {
      return static_cast< const FunctionImp & >( *this );
    }
   
    FunctionImp &asImp ()
    {
      return static_cast< FunctionImp & >( *this );
    }
   
  private:
    //! Helper function for Mapping
    //! With this function, a combined mapping can choose the right application
    //! operator (i.e. the one from Mapping itself, or from Function/Operator)
    //! \note: Do not override this definition
    virtual void apply (const DomainType& arg, RangeType& dest) const 
    {
      operator()(arg, dest);
    }
  };

  ///@} 
}
#endif
