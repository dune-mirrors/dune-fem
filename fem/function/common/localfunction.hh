#ifndef DUNE_LOCALFUNCTION_HH
#define DUNE_LOCALFUNCTION_HH

#include <dune/common/bartonnackmanifcheck.hh>

namespace Dune
{

  /** @defgroup LocalFunction LocalFunction 
     @ingroup DiscreteFunction

  On every element from a discrete function the local funtion can be accessed.
  With the local function one has access to the dof and on the other hand to 
  the base function set of this actual element. Therefore this is called a 
  local function. 

  @{
 */
  
//****************************************************************************
//
//  --LocalFunctionInterface 
//
//! The LocalFunctionInterface is the Interface to local function which
//! form the discrete Function 
//
//****************************************************************************
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
      \param num local dof number 
  */
  RangeFieldType& operator [] (const int num) 
  {
    CHECK_INTERFACE_IMPLEMENTATION(asImp().operator [] ( num ));
    return asImp().operator [] ( num );
  }

  /** \brief access to local dofs (read-only)
      \param num local dof number 
  */
  const RangeFieldType& operator [] (const int num) const 
  {
    CHECK_INTERFACE_IMPLEMENTATION(asImp().operator [] ( num ));
    return asImp().operator [] ( num );
  }

  //! \brief return the number of local dofs for this local function 
  int numDofs() const 
  {
    CHECK_INTERFACE_IMPLEMENTATION(asImp().numDofs());
    return asImp().numDofs();
  }
  
  /** \brief evaluate local function.
      \param x local coordinate of evaluation point 
      \param ret return value 
  */
  void evaluate( const DomainType &x, RangeType &ret )
  {
    CHECK_AND_CALL_INTERFACE_IMPLEMENTATION(
       asImp().evaluate( x, ret ));
  }

  /** \brief evaluate local function.
      \param quad Quadrature 
      \param quadPoint number of quadrature point for evaluation 
      \param ret return value 
  */
  template <class QuadratureType> 
  void evaluate( const QuadratureType &quad, int quadPoint, 
                 RangeType &ret )
  {
    CHECK_AND_CALL_INTERFACE_IMPLEMENTATION(
        asImp().evaluate( quad, quadPoint, ret ));
  }

  /** \brief evaluate jacobian on reference element
      \param x local coordinate 
      \param grad return value 
  */
  void jacobian ( const DomainType &x, 
                  JacobianRangeType &grad ) 
  {
    CHECK_AND_CALL_INTERFACE_IMPLEMENTATION(
      asImp().jacobian( x, grad ));
  }

  /** \brief evaluate jacobian on reference element
      \param quad Quadrature 
      \param quadPoint number of quadrature point for evaluation 
      \param grad return value 
  */
  template <class QuadratureType>
  void jacobian ( const QuadratureType &quad, int quadPoint , 
                  JacobianRangeType &grad ) 
  {
    CHECK_AND_CALL_INTERFACE_IMPLEMENTATION(
      asImp().jacobian( quad, quadPoint, grad ));
  }

  /** \brief axpy operation for local function 
      \param quad Quadrature
      \param quadPoint number of quadrature point 
      \param factor axpy factor  
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
      \param quad Quadrature
      \param quadPoint number of quadrature point 
      \param factor axpy gradient factor  
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
      \param quad Quadrature
      \param quadPoint number of quadrature point 
      \param factor1 axpy factor 
      \param factor2 gradient axpy factor 
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
  //! The Interface to the dune programmer, use this class to derive 
  //! the own implementation. But only the methods declared in the interface
  //! class must be implemented. 
  //!
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
        \param v vector that is added to local function (can also be local function)
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
        \param v vector that is substracted from local function (can also be local function)
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
         \param entity Entity x is belonging to 
         \param x global evaluation coordinate 
         \param ret return value 
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
        \param quad Quadrature
        \param quadPoint number of quadrature point      
        \param ret return value 
    */
    template< class QuadratureType >
    void evaluate ( const QuadratureType &quad,
                    const int quadPoint,
                    RangeType &ret )
    {
      evaluate( quad.point( quadPoint ), ret );
    }
    
    /** \brief evaluate jacobian of the local function on 
               real world coordinate x and return ret (calls local method of entitys geometry)
         \param entity Entity x is belonging to 
         \param x global evaluation coordinate 
         \param grad  return value 
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
              (calls jacobian method with local coordiante)
        \param quad Quadrature
        \param quadPoint number of quadrature point      
        \param grad return value 
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
    */
    inline size_t size () const
    {
      return numDofs();
    }
  }; // end LocalFunctionDefault

} // end namespace Dune 
#endif
