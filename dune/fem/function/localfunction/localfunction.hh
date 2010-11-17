#ifndef DUNE_FEM_LOCALFUNCTION_HH
#define DUNE_FEM_LOCALFUNCTION_HH

#include <dune/fem/misc/engineconcept.hh>
#include <dune/fem/space/common/dofstorage.hh>
#include <dune/fem/function/common/function.hh>

namespace Dune
{

  /** \addtogroup LocalFunction
   *
   * On every element from a discrete function the local funtion can be
   * accessed. With the local function one has access to the dof and on the
   * other hand to the base function set of this actual element. Therefore
   * this is called a local function.
   *
   * \remarks The interface for using a LocalFunction is defined by the class
   *          LocalFunction.
   *
   * \{
   */



  // LocalFunction
  // -------------

  /** \class LocalFunction
   *  \brief interface for local functions
   *  
   *  Local functions are used to represend a discrete function on one entity.
   *  The LocalFunctionInterface defines the functionality that can be expected
   *  from such a local function.
   */
  template< class LFTraits >
  class LocalFunction
  : public EngineWrapper< typename LFTraits::LocalFunctionImpType, typename LFTraits::LocalFunctionUserType >
  {
    typedef EngineWrapper< typename LFTraits::LocalFunctionImpType, typename LFTraits::LocalFunctionUserType >
      BaseType;

  public:
    //! type of the traits
    typedef LFTraits Traits;

    //! type of the local function (this type!)
    typedef LocalFunction< Traits > LocalFunctionType;

    //! type of the discrete function space, the local function belongs to
    typedef typename Traits::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;

    //! type of the local function implementation (engine concept)
    typedef typename Traits::LocalFunctionImpType LocalFunctionImpType;

    //! type of the user implementation (Barton-Nackman)
    typedef typename Traits::LocalFunctionUserType LocalFunctionUserType;

  private:
    typedef typename DiscreteFunctionSpaceType::GridType GridType;

  public:
    //! type of the entity, the local function lives on
    typedef typename GridType::template Codim< 0 >::Entity EntityType;

    //! field type of the domain
    typedef typename DiscreteFunctionSpaceType::DomainFieldType DomainFieldType;
    //! field type of the range
    typedef typename DiscreteFunctionSpaceType::RangeFieldType RangeFieldType;
    //! type of domain vectors, i.e., type of coordinates
    typedef typename DiscreteFunctionSpaceType::DomainType DomainType;
    //! type of range vectors, i.e., type of function values
    typedef typename DiscreteFunctionSpaceType::RangeType RangeType;
    //! type of the Jacobian, i.e., type of evaluated Jacobian matrix
    typedef typename DiscreteFunctionSpaceType::JacobianRangeType JacobianRangeType;
    //! type of the Hessian
    typedef typename DiscreteFunctionSpaceType::HessianRangeType HessianRangeType;

    //! type of local coordinates 
    typedef typename EntityType::Geometry::LocalCoordinate LocalCoordinateType;

    //! dimension of the domain
    static const int dimDomain = DiscreteFunctionSpaceType::dimDomain;
    //! dimension of the range
    static const int dimRange = DiscreteFunctionSpaceType::dimRange;

    //! type of base function set  
    typedef typename DiscreteFunctionSpaceType::BaseFunctionSetType BaseFunctionSetType;

    //! type of the Jacobian of a base function
    typedef typename BaseFunctionSetType::JacobianRangeType LocalJacobianRangeType;
    // type of the Hessian of a base function
    typedef typename BaseFunctionSetType::HessianRangeType LocalHessianRangeType;

    /** \brief access to local dofs (read-only)
     *
     *  \param[in]  num  local dof number 
     *  \return reference to dof 
     */
    const RangeFieldType &operator[] ( const int num ) const
    {
      return asImp()[ num ];
    }

    /** \brief access to local dofs (read-write)
     *
     *  \param[in]  num  local DoF number
     *  \return reference to DoF
     */
    RangeFieldType &operator[] ( const int num )
    {
      return asImp()[ num ];
    }

    /** \brief add another local function to this one
     *
     *  \note The local function to add may be any implementation of a
     *        LocalFunction.
     *
     *  \param[in]  lf  local function to add
     *
     *  \returns a reference to this local function (i.e., *this)
     */
    template< class T >
    LocalFunctionType &operator+= ( const LocalFunction< T > &lf )
    {
      asImp() += lf;
      return *this;
    }

    /** \brief assign all DoFs of this local function 
     *
     *  \param[in]  lf  local function to assign DoFs from 
     */
    template< class T >
    void assign ( const LocalFunction< T > &lf )
    {
      asImp().assign(lf);
    }

    /** \brief set all DoFs to zero */
    void clear ()
    {
      asImp().clear();
    }

    /** \brief subtract another local function to this one
     *
     *  \note The local function to suctract may be any implementation of a
     *        LocalFunction.
     *
     *  \param[in]  lf  local function to subtract
     *
     *  \returns a reference to this local function (i.e., *this)
     */
    template< class T >
    LocalFunctionType &operator-= ( const LocalFunction< T > &lf )
    {
      asImp() -= lf;
      return *this;
    }

    /** \brief add a multiple of another local function to this one
     *
     *  \note The local function to add may be any implementation of a
     *        LocalFunction.
     *
     *  \param[in]  s   scalar factor to scale lf with
     *  \param[in]  lf  local function to add
     *
     *  \returns a reference to this local function (i.e., *this)
     */
    template< class T >
    LocalFunctionType &axpy ( const RangeFieldType s, const LocalFunction< T > &lf )
    {
      asImp().axpy( s, lf );
      return *this;
    }

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
    void axpy ( const PointType &x, const RangeType &factor )
    {
      asImp().axpy( x, factor );
    }
  
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
    void axpy ( const PointType &x, const JacobianRangeType &factor)
    {
      asImp().axpy( x, factor );
    }

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
    void axpy ( const PointType &x, const RangeType &factor1, const JacobianRangeType &factor2 )
    {
      asImp().axpy( x, factor1, factor2 );
    }

    /** \brief obtain the order of this local function
     *
     *  The order of a local function refers to the polynomial order required
     *  to integrate it exactly.
     *
     *  \note It is not completely clear what this value should be, e.g., for
     *        bilinear base functions.
     *
     *  \returns order of the local function
     */
    int order () const
    {
      return asImp().order();
    }

    /** \brief obtain the base function set for this local function
     *
     *  \returns reference to the base function set
     */
    const BaseFunctionSetType &baseFunctionSet () const 
    {
      return asImp().baseFunctionSet();
    }

    /** \brief obtain the entity, this local function lives on
     *
     *  \returns reference to the entity
     */
    const EntityType &entity () const
    {
      return asImp().entity();
    }
    
    /** \brief evaluate a partial deriviaive of the local function
     *
     *  \param[in]   diffVariable  vector describing the desired partial
     *                             derivative
     *  \param[in]   x             evaluation point in local coordinates 
     *  \param[out]  ret           value of the function in the given point
     */
    template< int diffOrder, class PointType >
    void evaluate ( const FieldVector< int, diffOrder > &diffVariable,
                    const PointType &x, RangeType &ret ) const
    {
      asImp().evaluate( diffVariable, x, ret );
    }
   
    /** \brief evaluate the local function
     *
     *  \param[in]   x    evaluation point in local coordinates 
     *  \param[out]  ret  value of the function in the given point
     */
    template< class PointType >
    void evaluate ( const PointType &x, RangeType &ret ) const
    {
      asImp().evaluate( x, ret );
    }

    void init ( const EntityType &entity )
    {
      asImp().init( entity );
    }

    /** \brief evaluate Jacobian of the local function
     *
     *  \note Though the Jacobian is evaluated on the reference element, the
     *        return value is the Jacobian with respect to the actual entity.
     *
     *  \param[in]   x    evaluation point in local coordinates
     *  \param[out]  ret  Jacobian of the function in the evaluation point
     */
    template< class PointType >
    void jacobian ( const PointType &x, JacobianRangeType &ret ) const
    {
      asImp().jacobian( x, ret );
    }

    /** \brief evaluate Hessian of the local function
     *
     *  \note Though the Hessian is evaluated on the reference element, the
     *        return value is the Hessian with respect to the actual entity.
     *
     *  \param[in]   x        evaluation point in local coordinates
     *  \param[out]  hessian  Hessian of the function in the evaluation point
     */
    template< class PointType >
    void hessian ( const PointType &x, HessianRangeType &hessian ) const
    {
      asImp().hessian( x, hessian );
    }

    /** \brief obtain the number of local DoFs
     *
     *  Obtain the number of local DoFs of this local function. The value is
     *  identical to the number of base functons on the entity.
     *  
     *  \returns number of local DoFs
     */
    int numDofs () const 
    {
      return asImp().numDofs();
    }

    /** \brief obtain the number of local DoFs in the scalar case 
     *
     *  Obtain the number of local DoFs of the scalar case 
     *  of this local function. The value is
     *  identical to the number of base functons on the entity.
     *  
     *  \returns number of local DoFs, scalar case 
     */
    int numScalarDofs () const 
    {
      return asImp().numScalarDofs();
    }

    template< class QuadratureType, class VectorType >
    void axpyQuadrature ( const QuadratureType &quad, const VectorType &factorVec )
    {
      assert( factorVec.size() > 0 );
      axpyQuadrature( quad, factorVec, factorVec[ 0 ] );
    }

    template< class QuadratureType, class RangeVectorType, class JacobianRangeVectorType >
    void axpyQuadrature ( const QuadratureType &quad,
                          const RangeVectorType& rangeVector,
                          const JacobianRangeVectorType& jacobianVector )
    {
      asImp().baseFunctionSet().axpyRanges( quad, rangeVector, asImp() );
      asImp().baseFunctionSet().axpyJacobians( quad, entity().geometry(), jacobianVector, asImp() );
    }
    
    template< class QuadratureType, class VectorType  >
    void evaluateQuadrature( const QuadratureType &quad, VectorType &result ) const
    {
      assert( result.size() > 0 );
      evaluateQuadrature( quad, result, result[ 0 ] );
    }

  protected:  
    template< class QuadratureType, class VectorType  >
    void axpyQuadrature ( const QuadratureType &quad, const VectorType &factorVec, const RangeType & )
    {
      asImp().baseFunctionSet().axpyRanges( quad, factorVec, asImp() );
    }

    template< class QuadratureType, class VectorType  >
    void axpyQuadrature ( const QuadratureType &quad, const VectorType& factorVec, const JacobianRangeType & )
    {
      asImp().baseFunctionSet().axpyJacobians( quad, entity().geometry(), factorVec, asImp() );
    }

    // evaluate local function and store results in vector of RangeTypes 
    template< class QuadratureType, class VectorType  >
    void evaluateQuadrature( const QuadratureType &quad, VectorType &result, const RangeType & ) const
    {
      asImp().baseFunctionSet().evaluateRanges( quad, asImp(), result );
    }

    // evaluate jacobian of local function and store result in vector of
    // JacobianRangeTypes 
    template< class QuadratureType, class VectorType >
    void evaluateQuadrature( const QuadratureType &quad, VectorType &result, const JacobianRangeType & ) const
    {
      asImp().baseFunctionSet().evaluateJacobians( quad, entity().geometry(), asImp(), result );
    }

    using BaseType::asImp;
  };



  /** \class LocalFunctionDefault
   *  \brief default implementation of a LocalFunction
   */
  template< class DiscreteFunctionSpace, class LocalFunctionImp >
  class LocalFunctionDefault
  : public EngineDefault< LocalFunctionImp >
  {
    typedef EngineDefault< LocalFunctionImp > BaseType;

  public:
    //! type of  discrete function space the local function belongs to
    typedef DiscreteFunctionSpace DiscreteFunctionSpaceType;

  private:
    typedef typename DiscreteFunctionSpaceType::GridType GridType;

  public:
    //! type of the entity, the local function lives on
    typedef typename GridType::template Codim< 0 >::Entity EntityType;

    //! field type of the domain
    typedef typename DiscreteFunctionSpaceType::DomainFieldType DomainFieldType;
    //! field type of the range
    typedef typename DiscreteFunctionSpaceType::RangeFieldType RangeFieldType;
    //! type of domain vectors, i.e., type of coordinates
    typedef typename DiscreteFunctionSpaceType::DomainType DomainType;
    //! type of range vectors, i.e., type of function values
    typedef typename DiscreteFunctionSpaceType::RangeType RangeType;
    //! type of Jacobian, i.e., type of evaluated Jacobian matrix
    typedef typename DiscreteFunctionSpaceType::JacobianRangeType JacobianRangeType;
    //! type of the Hessian
    typedef typename DiscreteFunctionSpaceType::HessianRangeType HessianRangeType;

    //! type of base function set  
    typedef typename DiscreteFunctionSpaceType::BaseFunctionSetType
      BaseFunctionSetType;

    //! dimension of the domain
    enum { dimDomain = DiscreteFunctionSpaceType::dimDomain };
    //! dimension of the range
    enum { dimRange = DiscreteFunctionSpaceType::dimRange };

  protected:
    //! type of entity's geometry
    typedef typename EntityType::Geometry GeometryType;

  public:
    /** \copydoc Dune::LocalFunction::operator+=(const LocalFunction<T> &lf) */
    template< class T >
    void operator+= ( const LocalFunction< T > &lf );

    /** \copydoc Dune::LocalFunction::operator-=(const LocalFunction<T> &lf) */
    template< class T >
    void operator-= ( const LocalFunction< T > &lf );

    template< class T >
    void assign ( const LocalFunction< T > &lf );

    void clear ( );
    
    /** \copydoc Dune::LocalFunction::axpy(const RangeFieldType s,const LocalFunction<T> &lf) */
    template< class T >
    void axpy ( const RangeFieldType s, const LocalFunction< T > &lf );

    /** \copydoc Dune::LocalFunction::evaluate(const FieldVector<int,diffOrder> &diffVariable,const PointType &x,RangeType &ret) const */
    template< int diffOrder, class PointType >
    void evaluate ( const FieldVector< int, diffOrder > &diffVariable,
                    const PointType &x, RangeType &ret ) const;

    /** \copydoc Dune::LocalFunction::evaluate(const PointType &x,RangeType &ret) const */
    template< class PointType >
    void evaluate ( const PointType &x, RangeType &ret ) const;

    /** \copydoc Dune::LocalFunction::jacobian(const PointType &x,JacobianRangeType &ret) const */
    template< class PointType >
    void jacobian ( const PointType &x, JacobianRangeType &ret ) const;

    template< class PointType >
    void hessian ( const PointType &x, HessianRangeType &hessian ) const;

    /** \copydoc Dune::LocalFunction::axpy(const PointType &x,const RangeType &factor) */
    template< class PointType >
    void axpy( const PointType &x, const RangeType &factor );

    /** \copydoc Dune::LocalFunction::axpy(const PointType &x,const JacobianRangeType &factor) */
    template< class PointType >
    void axpy( const PointType &x, const JacobianRangeType &factor );

    /** \copydoc Dune::LocalFunction::axpy(const PointType &x,const RangeType &factor1,const JacobianRangeType &factor2) */
    template< class PointType >
    void axpy( const PointType &x, const RangeType &factor1, const JacobianRangeType &factor2 );

    int numScalarDofs () const
    {
      const int numDofs = asImp().numDofs();
      assert( numDofs % dimRange == 0 );
      return numDofs / dimRange;
    }

  protected:
    using BaseType::asImp;
  };

  /** \} */



  // Implementation of LocalFunctionDefault
  // --------------------------------------

  template< class DiscreteFunctionSpace, class LocalFunctionImp >
  template< class T >
  inline void LocalFunctionDefault< DiscreteFunctionSpace, LocalFunctionImp >
    ::operator+= ( const LocalFunction< T > &lf )
  {
    asImp().axpy( RangeFieldType(1), lf );
  }



  template< class DiscreteFunctionSpace, class LocalFunctionImp >
  template< class T >
  inline void LocalFunctionDefault< DiscreteFunctionSpace, LocalFunctionImp >
    ::operator-= ( const LocalFunction< T > &lf )
  {
    asImp().axpy( RangeFieldType(-1), lf );
  }


  template< class DiscreteFunctionSpace, class LocalFunctionImp >
  template< class T >
  inline void LocalFunctionDefault< DiscreteFunctionSpace, LocalFunctionImp >
    ::assign ( const LocalFunction< T > &lf )
  {
    const int numDofs = asImp().numDofs();
    assert( numDofs == lf.numDofs() );

    for( int i = 0; i < numDofs; ++i )
      asImp()[ i ] = lf[ i ];
  }


  template< class DiscreteFunctionSpace, class LocalFunctionImp >
  inline void LocalFunctionDefault< DiscreteFunctionSpace, LocalFunctionImp >
    ::clear ()
  {
    const int numDofs = asImp().numDofs();
    for( int i = 0; i < numDofs; ++i )
      asImp()[ i ] = 0.0;
  }
  
  
  template< class DiscreteFunctionSpace, class LocalFunctionImp >
  template< class T >
  inline void LocalFunctionDefault< DiscreteFunctionSpace, LocalFunctionImp >
    ::axpy ( const RangeFieldType s, const LocalFunction< T > &lf )
  {
    const int numDofs = asImp().numDofs();
    assert( numDofs == lf.numDofs() );

    for( int i = 0; i < numDofs; ++i )
      asImp()[ i ] += s * lf[ i ];
  }


  template< class DiscreteFunctionSpace, class LocalFunctionImp >
  template< int diffOrder, class PointType >
  inline void LocalFunctionDefault< DiscreteFunctionSpace, LocalFunctionImp >
    ::evaluate ( const FieldVector< int, diffOrder > &diffVariable,
                 const PointType &x, RangeType &ret ) const
  {
    asImp().baseFunctionSet().evaluateAll( diffVariable, x, asImp(), ret);
  }

  
  template< class DiscreteFunctionSpace, class LocalFunctionImp >
  template< class PointType >
  inline void LocalFunctionDefault< DiscreteFunctionSpace, LocalFunctionImp >
    :: evaluate ( const PointType &x,
                  RangeType &ret ) const
  {
    asImp().baseFunctionSet().evaluateAll( x, asImp(), ret);
  }


  template< class DiscreteFunctionSpace, class LocalFunctionImp >
  template< class PointType >
  inline void LocalFunctionDefault< DiscreteFunctionSpace, LocalFunctionImp >
    ::jacobian ( const PointType &x, JacobianRangeType &ret ) const
  {
    asImp().baseFunctionSet().jacobianAll( 
        x, 
        asImp().entity().geometry().jacobianInverseTransposed( coordinate( x ) ),
        asImp(), ret);
  }


  template< class DiscreteFunctionSpace, class LocalFunctionImp >
  template< class PointType >
  inline void LocalFunctionDefault< DiscreteFunctionSpace, LocalFunctionImp >
    ::hessian ( const PointType &x, HessianRangeType &hessian ) const
  {
    asImp().baseFunctionSet().hessianAll( x, asImp().entity().geometry(), asImp(), hessian );
  }
 
  
  template< class DiscreteFunctionSpace, class LocalFunctionImp >
  template< class PointType >
  inline void LocalFunctionDefault< DiscreteFunctionSpace, LocalFunctionImp >
    ::axpy ( const PointType &x, const RangeType &factor )
  {
    asImp().baseFunctionSet().axpy( x, factor, asImp() );
  }


  template< class DiscreteFunctionSpace, class LocalFunctionImp >
  template< class PointType >
  inline void LocalFunctionDefault< DiscreteFunctionSpace, LocalFunctionImp >
    ::axpy ( const PointType &x, const JacobianRangeType &factor )
  {
    asImp().baseFunctionSet().axpy( 
            x, 
            asImp().entity().geometry().jacobianInverseTransposed( coordinate( x ) ),
            factor, 
            asImp() );
  }

  
  template< class DiscreteFunctionSpace, class LocalFunctionImp >
  template< class PointType >
  inline void LocalFunctionDefault< DiscreteFunctionSpace, LocalFunctionImp >
    ::axpy ( const PointType &x, const RangeType &factor1, const JacobianRangeType &factor2 )
  {
    asImp().baseFunctionSet().axpy( 
          x, 
          asImp().entity().geometry().jacobianInverseTransposed( coordinate( x ) ),
          factor1, factor2, asImp() );
  }


} // end namespace Dune 

#endif // #ifndef DUNE_FEM_LOCALFUNCTION_HH
