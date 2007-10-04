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
/** \class LocalFunctionInterface
 *  \brief interface for local functions
 *  
 *  Local functions are used to represend a discrete function on one entity.
 *  The LocalFunctionInterface defines the functionality that can be expected
 *  from such a local function.
 */
template< class DiscreteFunctionSpaceImp, class LocalFunctionImp > 
class LocalFunctionInterface 
{
public:
  //! type of the discrete function space, this local function belongs to
  typedef DiscreteFunctionSpaceImp DiscreteFunctionSpaceType;

private:
  typedef LocalFunctionInterface< DiscreteFunctionSpaceType, LocalFunctionImp >
    ThisType;

public:
  //- this are the types for the derived classes 
  //! type of domain field, i.e. type of coordinate component
  typedef typename DiscreteFunctionSpaceType :: DomainFieldType DomainFieldType;
  //! type of range field, i.e. dof type 
  typedef typename DiscreteFunctionSpaceType :: RangeFieldType RangeFieldType;
  //! type of domain, i.e. type of coordinates 
  typedef typename DiscreteFunctionSpaceType :: DomainType DomainType;
  //! type of range, i.e. result of evaluation 
  typedef typename DiscreteFunctionSpaceType :: RangeType RangeType;
  //! type of jacobian, i.e. type of evaluated gradient 
  typedef typename DiscreteFunctionSpaceType :: JacobianRangeType JacobianRangeType;

  //! type of base function set  
  typedef typename DiscreteFunctionSpaceType :: BaseFunctionSetType BaseFunctionSetType; 
  
  /** \brief access to local dofs (read-only)
   *
   *  \param[in]  num  local dof number 
   *  \return reference to dof 
  */
  const RangeFieldType &operator[] ( const int num ) const
  {
    CHECK_INTERFACE_IMPLEMENTATION( asImp().operator[]( num ) );
    return asImp().operator[]( num );
  }

  /** \brief access to local dofs (read-write)
   *
   *  \param[in]  num  local DoF number
   *  \return reference to DoF
   */
  inline RangeFieldType &operator[] ( const int num )
  {
    CHECK_INTERFACE_IMPLEMENTATION( asImp().operator[]( num ) );
    return asImp().operator[]( num );
  }

  /** \brief obtain the number of local DoFs
   *
   *  Obtain the number of local DoFs of this local function. The value is
   *  identical to the number of base functons on the entity.
   *  
   *  \returns number of local DoFs
   */
  int numDofs() const 
  {
    CHECK_INTERFACE_IMPLEMENTATION(asImp().numDofs());
    return asImp().numDofs();
  }
  
  /** \brief evaluate the local function
   *
   *  \param[in]   x    evaluation point in local coordinates 
   *  \param[out]  ret  value of the function in the given point
  */
  inline void evaluate ( const DomainType &x,
                         RangeType &ret ) const
  {
    CHECK_AND_CALL_INTERFACE_IMPLEMENTATION
      ( asImp().evaluate( x, ret ) );
  }

  /** \brief evaluate the local function in a quadrature point
   *
   *  \param[in]   quadrature  quadrature to use
   *  \param[in]   quadPoint   number of the quadrature point within the
   *                           quadrature
   *  \param[out]  ret         value of the function in the quadrature point
   */
  template< class QuadratureType >
  inline void evaluate( const QuadratureType &quadrature,
                        const int quadPoint, 
                        RangeType &ret ) const
  {
    CHECK_AND_CALL_INTERFACE_IMPLEMENTATION
      ( asImp().evaluate( quadrature, quadPoint, ret ) );
  }

  /** \brief evaluate Jacobian of the local function
   *
   *  \note Though the Jacobian is evaluated on the reference element, the
   *        return value is the Jacobian with respect to the actual entity.
   *
   *  \param[in]   x    evaluation point in local coordinates
   *  \param[out]  ret  Jacobian of the function in the evaluation point
   */
  inline void jacobian ( const DomainType &x, 
                         JacobianRangeType &ret ) const
  {
    CHECK_AND_CALL_INTERFACE_IMPLEMENTATION
      ( asImp().jacobian( x, ret ) );
  }

  /** \brief evaluate Jacobian of the local function in a quadrature point
   *
   *  \note Though the Jacobian is evaluated on the reference element, the
   *        return value is the Jacobian with respect to the actual entity.
   *
   *  \param[in]   quadrature  quadrature to use
   *  \param[in]   quadPoint   number of the quadrature point within the
   *                           quadrature
   *  \param[out]  ret         Jacobian of the function in the quadrature point
   */
  template< class QuadratureType >
  void jacobian ( const QuadratureType &quadrature,
                  const int quadPoint,
                  JacobianRangeType &ret ) const
  {
    CHECK_AND_CALL_INTERFACE_IMPLEMENTATION
      ( asImp().jacobian( quadrature, quadPoint, ret ) );
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

  /** \brief obtain the base function set for this local function
   *
   *  \returns reference to the base function set
   */
  const BaseFunctionSetType &baseFunctionSet () const 
  {
    CHECK_INTERFACE_IMPLEMENTATION( asImp().baseFunctionSet() );
    return asImp().baseFunctionSet();
  }
  
protected:
  // Barton-Nackman trick 
  inline const LocalFunctionImp &asImp() const
  {
    return static_cast< const LocalFunctionImp & >( *this );
  }

  // Barton-Nackman trick 
  inline LocalFunctionImp &asImp()
  { 
    return static_cast< LocalFunctionImp & >( *this );
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
                          RangeType &ret ) const
    {
      xLoc_ = entity.geometry().local( x );
      evaluate( xLoc_, ret );
    }

    /** \copydoc Dune::LocalFunctionInterface::evaluate(const QuadratureType &quadrature,const int quadPoint,RangeType &ret) const
     *
     *  \note The default implementation just calls
     *  \code
     *  evaluate( quadrature.point( quadPoint ), ret );
     *  \endcode
     */
    template< class QuadratureType >
    void evaluate ( const QuadratureType &quadrature,
                    const int quadPoint,
                    RangeType &ret ) const
    {
      evaluate( quadrature.point( quadPoint ), ret );
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

    /** \copydoc Dune::LocalFunctionInterface::jacobian(const QuadratureType &quadrature,const int quadPoint,JacobianRangeType &ret) const
     *
     *  \note The default implementation just calls
     *  \code
     *  jacobian( quadrature.point( quadPoint ), ret );
     *  \endcode
     */
    template< class QuadratureType >
    void jacobian ( const QuadratureType &quadrature,
                    const int quadPoint,
                    JacobianRangeType &ret ) const
    {
      jacobian( quadrature.point( quadPoint ), ret );
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
