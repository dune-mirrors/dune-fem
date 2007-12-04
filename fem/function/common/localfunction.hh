#ifndef DUNE_LOCALFUNCTION_HH
#define DUNE_LOCALFUNCTION_HH

#include <dune/fem/misc/bartonnackmaninterface.hh>

#ifndef DUNE_FEM_COMPATIBILITY
#define DUNE_FEM_COMPATIBILITY 1
#endif

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



  /** \class LocalFunctionInterface
   *  \brief interface for local functions
   *  
   *  Local functions are used to represend a discrete function on one entity.
   *  The LocalFunctionInterface defines the functionality that can be expected
   *  from such a local function.
   */
  template< class DiscreteFunctionSpaceImp, class LocalFunctionImp > 
  class LocalFunctionInterface 
  : public BartonNackmanInterface
    < LocalFunctionInterface< DiscreteFunctionSpaceImp, LocalFunctionImp >,
      LocalFunctionImp >
  {
  public:
    //! type of the discrete function space, this local function belongs to
    typedef DiscreteFunctionSpaceImp DiscreteFunctionSpaceType;

    //! type of implementation (Barton-Nackman)
    typedef LocalFunctionImp LocalFunctionType;

  private:
    typedef LocalFunctionInterface< DiscreteFunctionSpaceType, LocalFunctionType >
      ThisType;
    typedef BartonNackmanInterface< ThisType, LocalFunctionType > BaseType;

  public:
    //! field type of the domain 
    typedef typename DiscreteFunctionSpaceType :: DomainFieldType DomainFieldType;
    //! field type of the range
    typedef typename DiscreteFunctionSpaceType :: RangeFieldType RangeFieldType;
    //! type of domain vectors, i.e., type of coordinates
    typedef typename DiscreteFunctionSpaceType :: DomainType DomainType;
    //! type of range vectors, i.e., type of function values
    typedef typename DiscreteFunctionSpaceType :: RangeType RangeType;
    //! type of Jacobian, i.e., type of evaluated Jacobian matrix
    typedef typename DiscreteFunctionSpaceType :: JacobianRangeType JacobianRangeType;

    //! type of base function set  
    typedef typename DiscreteFunctionSpaceType :: BaseFunctionSetType
      BaseFunctionSetType; 

  protected:
    using BaseType :: asImp;

  public:
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
    template< class PointType >
    inline void evaluate ( const PointType &x,
                           RangeType &ret ) const
    {
      CHECK_AND_CALL_INTERFACE_IMPLEMENTATION
        ( asImp().evaluate( x, ret ) );
    }

#if DUNE_FEM_COMPATIBILITY
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
#endif

    /** \brief evaluate Jacobian of the local function
     *
     *  \note Though the Jacobian is evaluated on the reference element, the
     *        return value is the Jacobian with respect to the actual entity.
     *
     *  \param[in]   x    evaluation point in local coordinates
     *  \param[out]  ret  Jacobian of the function in the evaluation point
     */
    template< class PointType >
    inline void jacobian ( const PointType &x, 
                           JacobianRangeType &ret ) const
    {
      CHECK_AND_CALL_INTERFACE_IMPLEMENTATION
        ( asImp().jacobian( x, ret ) );
    }

#if DUNE_FEM_COMPATIBILITY
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
#endif
    
    /** \brief axpy operation for local function
     *
     *  Denoting the DoFs of the local function by \f$u_i\f$ and the base
     *  functions by \f$\varphi_i\f$, this function performs the following
     *  operation:
     *  \f[
     *  u_i = u_i + factor \cdot \varphi_i( x )
     *  \f]
     *
     *  \param[in]  x       point to evaluate base functions in
     *  \param[in]  factor  axpy factor
     */
    template< class PointType >
    inline void axpy ( const PointType &x,
                       const RangeType &factor )
    {
      CHECK_AND_CALL_INTERFACE_IMPLEMENTATION
        ( asImp().axpy( x, factor ) );
    }
  
#if DUNE_FEM_COMPATIBILITY
    /** \brief axpy operation for local function
     *
     *  Denoting the DoFs of the local function by \f$u_i\f$ and the base
     *  functions by \f$\varphi_i\f$, this function performs the following
     *  operation:
     *  \f[
     *  u_i = u_i + factor \cdot \varphi_i( x )
     *  \f]
     *
     *  \param[in]  quadrature  quadrature to use
     *  \param[in]  quadPoint   number of the quadrature point wihin the
     *                          quadrature
     *  \param[in]  factor      axpy factor
     */
    template< class QuadratureType >
    inline void axpy ( const QuadratureType &quadrature,
                       const int quadPoint, 
                       const RangeType &factor )
    {
      CHECK_AND_CALL_INTERFACE_IMPLEMENTATION
        ( asImp().axpy( quadrature, quadPoint, factor ) );
    }
#endif
    
    /** \brief axpy operation for local function
     *
     *  Denoting the DoFs of the local function by \f$u_i\f$ and the base
     *  functions by \f$\varphi_i\f$, this function performs the following
     *  operation:
     *  \f[
     *  u_i = u_i + factor \cdot \nabla\varphi_i( x )
     *  \f]
     *
     *  \param[in]  x       point to evaluate jacobian of base functions in
     *  \param[in]  factor  axpy factor
     */
    template< class PointType >
    inline void axpy ( const PointType &x,
                       const JacobianRangeType &factor)
    {
      CHECK_AND_CALL_INTERFACE_IMPLEMENTATION( asImp().axpy( x, factor ) );
    }

#if DUNE_FEM_COMPATIBILITY
    /** \brief axpy operation for local function
     *
     *  Denoting the DoFs of the local function by \f$u_i\f$ and the base
     *  functions by \f$\varphi_i\f$, this function performs the following
     *  operation:
     *  \f[
     *  u_i = u_i + factor \cdot \nabla\varphi_i( x )
     *  \f]
     *
     *  \param[in]  quadrature  quadrature to use
     *  \param[in]  quadPoint   number of the quadrature point wihin the
     *                          quadrature
     *  \param[in]  factor      axpy factor
     */
    template< class QuadratureType >
    inline void axpy ( const QuadratureType &quadrature,
                       const int quadPoint, 
                       const JacobianRangeType &factor)
    {
      CHECK_AND_CALL_INTERFACE_IMPLEMENTATION
        ( asImp().axpy( quadrature, quadPoint, factor ) );
    }
#endif
    
    /** \brief axpy operation for local function
     *
     *  Denoting the DoFs of the local function by \f$u_i\f$ and the base
     *  functions by \f$\varphi_i\f$, this function performs the following
     *  operation:
     *  \f[
     *  u_i = u_i + factor1 \cdot \varphi_i( x ) + factor2 \cdot \nabla\varphi_i( x )
     *  \f]
     *
     *  \param[in]  x        point to evaluate base functions in
     *  \param[in]  factor1  axpy factor for \f$\varphi( x )\f$
     *  \param[in]  factor2  axpy factor for \f$\nabla\varphi( x )\f$
     */
    template< class PointType >
    inline void axpy ( const PointType &x,
                       const RangeType &factor1,
                       const JacobianRangeType &factor2 )
    {
      CHECK_AND_CALL_INTERFACE_IMPLEMENTATION
        ( asImp().axpy( x, factor1, factor2 ) );
    }
   
#if DUNE_FEM_COMPATIBILITY
    /** \brief axpy operation for local function
     *
     *  Denoting the DoFs of the local function by \f$u_i\f$ and the base
     *  functions by \f$\varphi_i\f$, this function performs the following
     *  operation:
     *  \f[
     *  u_i = u_i + factor1 \cdot \varphi_i( x ) + factor2 \cdot \nabla\varphi_i( x )
     *  \f]
     *
     *  \param[in]  quadrature  quadrature to use
     *  \param[in]  quadPoint   number of the quadrature point wihin the
     *                          quadrature
     *  \param[in]  factor1     axpy factor for \f$\varphi( x )\f$
     *  \param[in]  factor2     axpy factor for \f$\nabla\varphi( x )\f$
     */
    template< class QuadratureType >
    inline void axpy ( const QuadratureType &quadrature,
                       const int quadPoint,
                       const RangeType &factor1,
                       const JacobianRangeType &factor2 )
    {
      CHECK_AND_CALL_INTERFACE_IMPLEMENTATION
        ( asImp().axpy( quadrature, quadPoint, factor1, factor2 ) );
    }
#endif

    /** \brief obtain the base function set for this local function
     *
     *  \returns reference to the base function set
     */
    const BaseFunctionSetType &baseFunctionSet () const 
    {
      CHECK_INTERFACE_IMPLEMENTATION( asImp().baseFunctionSet() );
      return asImp().baseFunctionSet();
    }
  }; // end LocalFunctionInterface



  /** \brief default implementation of LocalFunctionInterface */
  template< class DiscreteFunctionSpaceImp, class LocalFunctionImp > 
  class LocalFunctionDefault
  : public LocalFunctionInterface< DiscreteFunctionSpaceImp, LocalFunctionImp >
  {
  public:
    //! type of the discrete function space, this local function belongs to
    typedef DiscreteFunctionSpaceImp DiscreteFunctionSpaceType;

    //! type of the implementation (Barton-Nackman)
    typedef LocalFunctionImp LocalFunctionType;

  private:
    typedef LocalFunctionDefault< DiscreteFunctionSpaceType, LocalFunctionType >
      ThisType;
    typedef LocalFunctionInterface< DiscreteFunctionSpaceType, LocalFunctionType >
      BaseType;
   
  public:
    typedef typename DiscreteFunctionSpaceType :: DomainFieldType DomainFieldType;
    typedef typename DiscreteFunctionSpaceType :: RangeFieldType RangeFieldType;
 
    typedef typename DiscreteFunctionSpaceType :: DomainType DomainType;
    typedef typename DiscreteFunctionSpaceType :: RangeType RangeType;

    typedef typename DiscreteFunctionSpaceType :: JacobianRangeType JacobianRangeType;
   
  protected:
    using BaseType :: asImp;

  public:
    using BaseType :: numDofs;
 
  public:
    //! Constructor
    LocalFunctionDefault ()
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

#if DUNE_FEM_COMPATIBILITY
    /** \copydoc Dune::LocalFunctionInterface::evaluate(const QuadratureType &quadrature,const int quadPoint,RangeType &ret) const
     *
     *  \note The default implementation just calls
     *  \code
     *  evaluate( quadrature.point( quadPoint ), ret );
     *  \endcode
     */
    template< class QuadratureType >
    inline void evaluate ( const QuadratureType &quadrature,
                           const int quadPoint,
                           RangeType &ret ) const
    {
      asImp().evaluate( quadrature[ quadPoint ], ret );
    }
#endif

#if DUNE_FEM_COMPATIBILITY
    /** \copydoc Dune::LocalFunctionInterface::jacobian(const QuadratureType &quadrature,const int quadPoint,JacobianRangeType &ret) const
     *
     *  \note The default implementation just calls
     *  \code
     *  jacobian( quadrature.point( quadPoint ), ret );
     *  \endcode
     */
    template< class QuadratureType >
    inline void jacobian ( const QuadratureType &quadrature,
                           const int quadPoint,
                           JacobianRangeType &ret ) const
    {
      asImp().jacobian( quadrature[ quadPoint ], ret );
    }
#endif
    
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
