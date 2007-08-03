#ifndef DUNE_LOCALFUNCTION_HH
#define DUNE_LOCALFUNCTION_HH

#include <dune/common/bartonnackmanifcheck.hh>

namespace Dune
{

  /** @addtogroup LocalFunction  

  On every element from a discrete function the local funtion can be accessed.
  With the local function one has access to the dof and on the other hand to 
  the base function set of this actual element. Therefore this is called a 
  local function. 

  \remarks 
  The interface for using a LocalFunction is defined by the class
  LocalFunctionInterface.

  @{
 */
  
//****************************************************************************
//
//  --LocalFunctionInterface 
//
//****************************************************************************
/** \brief 
   The LocalFunctionInterface is the Interface to local function which
   form the discrete Function 
   */
template < class DiscreteFunctionSpaceType, class LocalFunctionImp > 
class LocalFunctionInterface 
{
public:
  //- this are the types for the derived classes 
  //! type of domain field, i.e. type of coordinate component
  typedef typename DiscreteFunctionSpaceType :: DomainFieldType DomainFieldType;
  //! type of range field, i.e. dof type 
  typedef typename DiscreteFunctionSpaceType::RangeFieldType RangeFieldType;
  //! type of domain, i.e. type of coordinates 
  typedef typename DiscreteFunctionSpaceType::DomainType DomainType;
  //! type of range, i.e. result of evaluation 
  typedef typename DiscreteFunctionSpaceType::RangeType RangeType;
  //! type of jacobian, i.e. type of evaluated gradient 
  typedef typename DiscreteFunctionSpaceType::JacobianRangeType JacobianRangeType;
  
  //! type of base function set  
  typedef typename DiscreteFunctionSpaceType::BaseFunctionSetType BaseFunctionSetType; 

  /** \brief access to local dofs (read-write)
      \param[in] num local dof number 
      \return reference to dof 
  */
  RangeFieldType& operator [] (const int num) 
  {
    CHECK_INTERFACE_IMPLEMENTATION(asImp().operator [] ( num ));
    return asImp().operator [] ( num );
  }

  /** \brief access to local dofs (read-only)
      \param[in] num local dof number 
      \return reference to dof 
  */
  const RangeFieldType& operator [] (const int num) const 
  {
    CHECK_INTERFACE_IMPLEMENTATION(asImp().operator [] ( num ));
    return asImp().operator [] ( num );
  }

  /** \brief return the number of local dofs for this local function 
   *  \return number of local dofs 
   */
  int numDofs() const 
  {
    CHECK_INTERFACE_IMPLEMENTATION(asImp().numDofs());
    return asImp().numDofs();
  }
  
  /** \brief evaluate local function.
      \param[in] x local coordinate of evaluation point 
      \param[out] ret return value 
  */
  void evaluate( const DomainType &x, RangeType &ret )
  {
    CHECK_AND_CALL_INTERFACE_IMPLEMENTATION(
       asImp().evaluate( x, ret ));
  }

  /** \brief evaluate local function.
      \param[in] quad Quadrature 
      \param[in] quadPoint number of quadrature point for evaluation 
      \param[out] ret return value 
  */
  template <class QuadratureType> 
  void evaluate( const QuadratureType &quad, int quadPoint, 
                 RangeType &ret )
  {
    CHECK_AND_CALL_INTERFACE_IMPLEMENTATION(
        asImp().evaluate( quad, quadPoint, ret ));
  }

  /** \brief evaluate jacobian on reference element
      \param[in] x local coordinate 
      \param[out] grad return value 
  */
  void jacobian ( const DomainType &x, 
                  JacobianRangeType &grad ) 
  {
    CHECK_AND_CALL_INTERFACE_IMPLEMENTATION(
      asImp().jacobian( x, grad ));
  }

  /** \brief evaluate jacobian on reference element
      \param[in] quad Quadrature 
      \param[in] quadPoint number of quadrature point for evaluation 
      \param[out] grad return value 
  */
  template <class QuadratureType>
  void jacobian ( const QuadratureType &quad, int quadPoint , 
                  JacobianRangeType &grad ) 
  {
    CHECK_AND_CALL_INTERFACE_IMPLEMENTATION(
      asImp().jacobian( quad, quadPoint, grad ));
  }

  /** \brief axpy operation for local function 
      \param[in] quad Quadrature
      \param[in] quadPoint number of quadrature point 
      \param[in] factor axpy factor  
  */
  template <class QuadratureType>
  inline void axpy(const QuadratureType& quad,
                   const int quadPoint, 
                   const RangeType& factor)
  {
    CHECK_AND_CALL_INTERFACE_IMPLEMENTATION( 
      asImp().axpy( quad, quadPoint , factor ));
  }

  /** \brief axpy operation for local function 
      \param[in] quad Quadrature
      \param[in] quadPoint number of quadrature point 
      \param[in] factor axpy gradient factor  
  */
  template <class QuadratureType> 
  inline void axpy(const QuadratureType& quad, 
                   const int quadPoint, 
                   const JacobianRangeType& factor)
  {
    CHECK_AND_CALL_INTERFACE_IMPLEMENTATION( 
      asImp().axpy( quad, quadPoint , factor ));
  }

  /** \brief axpy operation for local function 
      \param[in] quad Quadrature
      \param[in] quadPoint number of quadrature point 
      \param[in] factor1 axpy factor 
      \param[in] factor2 gradient axpy factor 
  */
  template <class QuadratureType>
  inline void axpy(const QuadratureType &quad,
                   const int quadPoint ,
                   const RangeType& factor1,
                   const JacobianRangeType & factor2)
  {
    CHECK_AND_CALL_INTERFACE_IMPLEMENTATION( 
      asImp().axpy( quad, quadPoint, factor1, factor2 ));
  }

  //! \brief return reference to corresponding base function set of local function
  const BaseFunctionSetType& baseFunctionSet() const 
  {
    CHECK_INTERFACE_IMPLEMENTATION(asImp().baseFunctionSet());
    return asImp().baseFunctionSet();
  }
  
protected:
  //! Barton-Nackman trick 
  LocalFunctionImp& asImp() 
  { 
    return static_cast<LocalFunctionImp&>(*this);
  }
  //! Barton-Nackman trick 
  const LocalFunctionImp& asImp() const  
  {
    return static_cast<const LocalFunctionImp&>(*this);
  }
}; // end LocalFunctionInterface



  //************************************************************************
  //
  //  --LocalFunctionDefault 
  //
  /** \brief The Interface to the dune programmer, use this class to derive 
     the own implementation. But only the methods declared in the interface
     class must be implemented. 
  */
  //************************************************************************
  template< class DiscreteFunctionSpaceImp, class LocalFunctionImp > 
  class LocalFunctionDefault
  : public LocalFunctionInterface< DiscreteFunctionSpaceImp, LocalFunctionImp >
  {
  public:
    typedef DiscreteFunctionSpaceImp DiscreteFunctionSpaceType;

  private:
    typedef LocalFunctionDefault< DiscreteFunctionSpaceType, LocalFunctionImp > ThisType;
    typedef LocalFunctionInterface< DiscreteFunctionSpaceType, LocalFunctionImp > BaseType;

    using BaseType :: asImp;
    using BaseType :: numDofs;
    
  public:
    typedef typename DiscreteFunctionSpaceType :: DomainType DomainType;
    typedef typename DiscreteFunctionSpaceType :: RangeType RangeType;

    typedef typename DiscreteFunctionSpaceType :: JacobianRangeType JacobianRangeType;
   
    typedef typename DiscreteFunctionSpaceType :: DomainFieldType DomainFieldType;
    typedef typename DiscreteFunctionSpaceType :: RangeFieldType RangeFieldType;
    
  private:
    mutable DomainType xLoc_;

  public:
    //! Constructor
    LocalFunctionDefault ()
    : xLoc_( 0 )
    {
    }

    /** \brief += operation 
        \param[in] v vector that is added to local function (can also be local function)
        \return reference to local function (i.e. *this)
    */
    template< class VectorType >
    LocalFunctionImp &operator+= ( const VectorType &v )
    {
      unsigned int numDofs = this->numDofs(); 
      
      assert( numDofs == v.size() );
      for( unsigned int i = 0; i < numDofs; ++i )
        (*this)[ i ] += v[ i ];

      return asImp();
    }
    
    /** \brief -= operation 
        \param[in] v vector that is substracted from local function (can also be local function)
        \return reference to local function (i.e. *this)
    */
    template< class VectorType >
    LocalFunctionImp &operator-= ( const VectorType &v )
    {
      unsigned int numDofs = this->numDofs(); 
      
      assert( numDofs == v.size() );
      for( unsigned int i = 0; i < numDofs; ++i )
        (*this)[ i ] -= v[ i ];
      
      return asImp();
    }

    /** \brief evaluate the local function on 
               real world coordinate x and return ret (calls local method of entitys geometry)
         \param[in] entity Entity x is belonging to 
         \param[in] x global evaluation coordinate 
         \param[out] ret return value 
    */
    template< class EntityType >
    void evaluateGlobal ( const EntityType &entity, 
                          const DomainType &x, 
                          RangeType &ret )
    {
      xLoc_ = entity.geometry().local( x );
      evaluate( xLoc_, ret );
    }

    /** \brief default implementation of evaluation 
               of local function using a quadrature 
              (calls evaluate method with local coordiante)
        \param[in] quad Quadrature
        \param[in] quadPoint number of quadrature point      
        \param[out] ret return value 
    */
    template< class QuadratureType >
    void evaluate ( const QuadratureType &quad,
                    const int quadPoint,
                    RangeType &ret )
    {
      evaluate( quad.point( quadPoint ), ret );
    }
    
    /** \brief evaluate jacobian of the local function on 
               real world coordinate x and return ret 
               (calls local method of entitys geometry if not overloaded)
         \param[in] entity Entity x is belonging to 
         \param[in] x global evaluation coordinate 
         \param[out] grad  return value 
    */
    template< class EntityType >
    void jacobianGlobal ( const EntityType&  entity,
                          const DomainType&  x, 
                          JacobianRangeType& grad ) 
    {
      xLoc_ = entity.geometry().local( x );
      jacobian( xLoc_, grad );
    }

    /** \brief default implementation of evaluation 
               of jacobian using a quadrature 
              (calls jacobian method with local coordiante if not
              overloaded)
        \param[in] quad Quadrature
        \param[in] quadPoint number of quadrature point      
        \param[out] grad return value 
    */
    template< class QuadratureType >
    void jacobian ( const QuadratureType &quad,
                    const int quadPoint,
                    JacobianRangeType &grad )
    {
      jacobian( quad.point( quadPoint ), grad );
    }

    /** \brief size method to make local 
        function conform to a vector-like interface
        \return number of local dofs 
    */
    inline size_t size () const
    {
      return numDofs();
    }
  }; // end LocalFunctionDefault

///@} 
} // end namespace Dune 
#endif
