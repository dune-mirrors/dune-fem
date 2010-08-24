#ifndef DUNE_BASEFUNCTIONSETINTERFACE_HH
#define DUNE_BASEFUNCTIONSETINTERFACE_HH

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/geometrytype.hh>

#include <dune/fem/version.hh>
#include <dune/fem/misc/bartonnackmaninterface.hh>
#include <dune/fem/misc/deritype.hh>

namespace Dune
{

  /** \addtogroup BaseFunction
   *  \{
   */

  //------------------------------------------------------------------------
  //-
  //-  --BaseFunctionSetInterface
  //-
  //------------------------------------------------------------------------

  /** \class BaseFunctionSetInterface
   *  \brief interface class for base function sets
   */  
  template< class TraitsImp > 
  class BaseFunctionSetInterface
  : public BartonNackmanInterface< BaseFunctionSetInterface< TraitsImp >,
                                   typename TraitsImp :: BaseFunctionSetType >
  {
  public:
    //! type of the traits
    typedef TraitsImp Traits;

    //! type of base function set implementation (Barton-Nackman) 
    typedef typename Traits :: BaseFunctionSetType BaseFunctionSetType;

  private:
    typedef BaseFunctionSetInterface< Traits > ThisType;
    typedef BartonNackmanInterface< ThisType, BaseFunctionSetType > BaseType;
    
  public:
    //! type of function space 
    typedef typename Traits :: FunctionSpaceType FunctionSpaceType;
    //! type of domain vector, i.e. element of R^dimDomain 
    typedef typename FunctionSpaceType::DomainType DomainType;
    //! type of range vector, i.e. element of R^dimRange 
    typedef typename FunctionSpaceType::RangeType RangeType;
    //! type of domain vector components, i.e. double 
    typedef typename FunctionSpaceType::DomainFieldType DomainFieldType;
    //! type of range vector components, i.e. double 
    typedef typename FunctionSpaceType::RangeFieldType RangeFieldType;
    //! range type of gradient 
    typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;
    //! range type of second derivate 
    typedef typename FunctionSpaceType::HessianRangeType  HessianRangeType;
    //! dimension of domain 
    enum { dimDomain = FunctionSpaceType :: dimDomain };
    //! dimension of range 
    enum { dimRange  = FunctionSpaceType :: dimRange };

  protected:
    using BaseType :: asImp;

  public:
    /** \brief number of base functions 
        \return number of base functions 
    */
    inline int numBaseFunctions () const 
    {
      CHECK_INTERFACE_IMPLEMENTATION(asImp().numBaseFunctions());
      return asImp().numBaseFunctions();
    }

    /** \brief obtain type of geometry
     * 
     *  \returns GeometryType of the base function set
     */
    inline GeometryType geometryType () const
    {
      CHECK_INTERFACE_IMPLEMENTATION( asImp().geometryType() );
      return asImp().geometryType();
    }

    /** \brief evaluate a derivative of the base function
     *
     *  Evaluates a partial derivative of the base function. If the derivative is
     *  of order 0, the base function itself is evaluated.
     *  
     *  \param[in]   baseFunction  number of base function to evaluate
     *  \param[in]   diffVariable  vector describing the partial derivative to
     *                             evaluate
     *  \param[in]   x             point within reference element to evaluate the
     *                             base function in
     *  \param[out]  phi           return value, i.e., value of the base function
     */
    template< int diffOrd, class PointType >
    inline void evaluate ( const int baseFunction,
                           const FieldVector< deriType, diffOrd > &diffVariable,
                           const PointType &x,
                           RangeType &phi ) const
    {
      CHECK_AND_CALL_INTERFACE_IMPLEMENTATION
        ( asImp().evaluate( baseFunction, diffVariable, x, phi ) );
    }

    /** \brief evaluate the base function
     *
     *  Evaluates the base function in a point (specified in local coordinates).
     *  
     *  \param[in]   baseFunction  number of base function to evaluate
     *  \param[in]   x             point within reference element to evaluate the
     *                             base function in
     *  \param[out]  phi           return value, i.e., value of the base function
     */
    template< class PointType >
    inline void evaluate ( const int baseFunction,
                           const PointType &x,
                           RangeType &phi ) const
    {
      CHECK_AND_CALL_INTERFACE_IMPLEMENTATION
        ( asImp().evaluate( baseFunction, x , phi ) );
    }

    /** \brief evaluate the Jacobian of the base function
     *
     *  Evaluates the Jacobian of base function in a point (specified in local
     *  coordinates.
     *
     *  \note The Jacobian is calculated with respect to the reference element.
     *        To obtain the jacobian with respect to a real entity, it must be
     *        multiplied (from the left) by
     *  \code
     *  entity.geometry().jacobianInverseTransposed()
     *  \endcode
     *
     *  \note While the Jacobian of base functions is returned with respect to
     *        the reference element, local functions return the Jacobian with
     *        respect to the real entity.
     *  
     *  \param[in]   baseFunction  number of base function to evaluate
     *  \param[in]   x             point within reference element to evaluate the
     *                             Jacobian in
     *  \param[out]  phi           return value, i.e., value of the Jacobian
     */
    template< class PointType >
    inline void jacobian ( const int baseFunction,
                           const PointType &x,
                           JacobianRangeType &phi ) const 
    {
      CHECK_AND_CALL_INTERFACE_IMPLEMENTATION
        ( asImp().jacobian( baseFunction, x, phi ) );
    }

    
    /** \brief evaluate the base function and multiply the result by a vector
     *
     *  It is quite common, that the base functions have the special structure
     *  \f{displaymath}
     *  \varphi( x ) = \sum_{k=1}^n \omega_k( x ) e_k,
     *  \f}
     *  where \f$e_k\f$ denotes the k-th unit vector. In this case, scalar
     *  products by \f$\psi\f$ can be evaluated very efficiently, since
     *  \f{displaymath}
     *  \varphi( x ) \cdot \psi = \sum_{k=1}^n \omega_k( x ) e_k \cdot \psi
     *                          = \sum_{k=1}^n \omega_k( x ) \psi_k.
     *  \f}
     *  This function provides the possibility to numerically use this
     *  information.
     *  
     *  \param[in]  baseFunction  number of the base function to evaluate
     *  \param[in]  x             point within reference element to evaluate
     *                            the base function in
     *  \param[in]  psi           vector to multiply the function value with
     *
     *  \returns the scalar product between the value of the base function and
     *           the specified vector
     */
    template< class PointType >
    DUNE_VERSION_DEPRECATED(1,2,remove)
    RangeFieldType evaluateSingle ( const int baseFunction,
                                    const PointType &x,
                                    const RangeType &psi ) const
    {
      CHECK_INTERFACE_IMPLEMENTATION
        ( asImp().evaluateSingle( baseFunction, x, psi ) );
      return asImp().evaluateSingle( baseFunction, x, psi );
    }

    /** \brief evaluate gradient of basefunction on given entity (uses
        jacobianInverseTransposed) and multiply with factor, return is RangeFieldType 
        \param[in]  baseFunction  number of base functions to evaluate jacobian  
        \param[in]  entity        entity gradient of base function is evaluated on 
        \param[in]  x             local point in reference element 
        \param[in]  psi           factor to multiply with 
        \return return scalar product between gradient of base function and factor 
    */
    template< class EntityType, class PointType >
    DUNE_VERSION_DEPRECATED(1,2,remove)
    inline RangeFieldType evaluateGradientSingle( const int baseFunction,
                                                  const EntityType &entity,
                                                  const PointType &x,
                                                  const JacobianRangeType &psi ) const
    {
      CHECK_INTERFACE_IMPLEMENTATION
        ( asImp().evaluateGradientSingle( baseFunction, entity, x, psi ) );
      return asImp().evaluateGradientSingle( baseFunction, entity, x, psi );
    }
  };



  //------------------------------------------------------------------------
  //-  
  //-  --BaseFunctionSetDefault
  //-
  //------------------------------------------------------------------------
  /** \brief default implementation of the BaseFunctionSetInterface */
  template< class TraitsImp >
  class BaseFunctionSetDefault
  : public BaseFunctionSetInterface< TraitsImp >
  {
  public:
    //! type of the traits
    typedef TraitsImp Traits;

  private:
    typedef BaseFunctionSetDefault< Traits > ThisType;
    typedef BaseFunctionSetInterface< Traits > BaseType;

  public:
    typedef typename Traits :: BaseFunctionSetType BaseFunctionSetType;
    typedef typename Traits :: FunctionSpaceType FunctionSpaceType;
    typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;
    typedef typename FunctionSpaceType::HessianRangeType  HessianRangeType;

    static const int dimDomain = BaseType::dimDomain;
    static const int dimRange = BaseType::dimRange;

    //! number of rows of jacobian range type 
    enum { dimRow = JacobianRangeType::rows };
    //! number of columns of jacobian range type 
    enum { dimCol = JacobianRangeType::cols };
   
    typedef typename FunctionSpaceType::DomainFieldType DomainFieldType;
    typedef typename FunctionSpaceType::RangeFieldType RangeFieldType;

    typedef typename FunctionSpaceType::DomainType DomainType ;
    typedef typename FunctionSpaceType::RangeType RangeType ;

  protected:
    using BaseType::asImp;

  public:
    using BaseType::numBaseFunctions;

    /** \copydoc Dune::BaseFunctionSetInterface::evaluate(const int baseFunction,const FieldVector<deriType,diffOrd> &diffVariable,const PointType &x,RangeType &phi) const */
    template< int diffOrd, class PointType >
    inline void evaluate ( const int baseFunction,
                           const FieldVector< deriType, diffOrd > &diffVariable,
                           const PointType &x,
                           RangeType &phi ) const
    {
      asImp().evaluate( baseFunction, diffVariable, coordinate( x ), phi );
    }
    
    /** \copydoc Dune::BaseFunctionSetInterface::evaluate(const int baseFunction,const PointType &x,RangeType &phi) const */
    template< class PointType >
    void evaluate ( const int baseFunction, const PointType &x, RangeType &phi ) const
    {
      FieldVector< deriType, 0 > diffVar;
      asImp().evaluate( baseFunction, diffVar, x, phi );
    }

    // This declaration is only to avoid a problem with specialization
    template< class PointType >
    void evaluate ( const unsigned int baseFunction, const PointType &x, RangeType &phi ) const
    {
      evaluate( int( baseFunction ), x, phi );
    }

    template< int diffOrder, class PointType, class LocalDofVectorType >
    void evaluateAll ( const FieldVector< int, diffOrder > &diffVariable,
                       const PointType &x, const LocalDofVectorType& dofs, RangeType &ret ) const
    {
      ret = 0;

      const int numBase = numBaseFunctions();
      for( int i = 0; i < numBase; ++i )
      {
        RangeType phi;
        evaluate( i, diffVariable, x, phi );
        ret.axpy( dofs[ i ], phi );
      }
    }

    template< class PointType, class LocalDofVectorType >
    void evaluateAll ( const PointType &x, const LocalDofVectorType &dofs, RangeType &ret ) const
    {
      ret = 0;

      const int numBase = numBaseFunctions();
      for( int i = 0; i < numBase; ++i )
      {
        RangeType phi;
        evaluate( i, x, phi );
        ret.axpy( dofs[ i ], phi );
      }
    }

    template< class PointType, class RangeVectorType >
    void evaluateAll ( const PointType &x, RangeVectorType &ret ) const
    {
      const int numBase = numBaseFunctions();
      ret.resize( numBase );
      for( int i = 0; i < numBase; ++i )
        evaluate( i, x, ret[ i ] );
    }

    /** \copydoc Dune::BaseFunctionSetInterface::jacobian(const int baseFunction,const PointType &x, JacobianRangeType &phi) const */
    template< class PointType >
    void jacobian ( const int baseFunction, const PointType &x, JacobianRangeType &phi ) const;

    template< class PointType, class GeometryJacobianInverseType,
              class LocalDofVectorType, class GlobalJacobianRangeType >
    void jacobianAll ( const PointType &x,
                       const GeometryJacobianInverseType& gjit, 
                       const LocalDofVectorType& dofs, 
                       GlobalJacobianRangeType &ret ) const;

    template< class PointType, class GeometryJacobianInverseType,
              class GlobalJacobianRangeVectorType >
    void jacobianAll ( const PointType &x,
                       const GeometryJacobianInverseType& gjit, 
                       GlobalJacobianRangeVectorType &ret ) const;

    template< class PointType >
    void hessian ( const int baseFunction, const PointType &x, HessianRangeType &hessian ) const;

    template< class PointType, class Geometry, class LocalDofVectorType, class GlobalHessianRangeType >
    void hessianAll ( const PointType &x, const Geometry &geometry,
                      const LocalDofVectorType &dofs, GlobalHessianRangeType &ret ) const;

    /** \copydoc Dune::BaseFunctionSetInterface::evaluateSingle(const int baseFunction,const PointType &x,const RangeType &psi) const */
    template< class PointType >
    DUNE_VERSION_DEPRECATED(1,2,remove)
    inline RangeFieldType evaluateSingle ( const int baseFunction,
                                           const PointType &x,
                                           const RangeType &psi ) const
    {
      RangeType phi;
      asImp().evaluate( baseFunction, x, phi );
      return phi * psi;
    }

    /** \copydoc Dune::BaseFunctionSetInterface::evaluateGradientSingle(const int baseFunction,const EntityType &entity,const PointType &x,const JacobianRangeType &psi) const */
    template< class EntityType, class PointType >
    DUNE_VERSION_DEPRECATED(1,2,remove)
    inline RangeFieldType evaluateGradientSingle( const int baseFunction,
                                                  const EntityType &entity,
                                                  const PointType &x,
                                                  const JacobianRangeType &psi ) const
    {
      typedef typename EntityType :: Geometry GeometryType;
      typedef FieldMatrix< typename GeometryType :: ctype,
                           GeometryType :: mydimension,
                           GeometryType :: mydimension >
        GeometryJacobianType;

      const GeometryType &geometry = entity.geometry();
      const GeometryJacobianType &jacobianInverseTransposed
        = geometry.jacobianInverseTransposed( coordinate( x ) );
 
      JacobianRangeType gradPhi;
      asImp().jacobian( baseFunction, x, gradPhi );

      RangeFieldType result = 0;
      for( int i = 0; i < FunctionSpaceType :: dimRange; ++i )
      {
        DomainType gradScaled( 0 );
        jacobianInverseTransposed.umv( gradPhi[ i ], gradScaled );
        result += gradScaled * psi[ i ];
      }
      return result;
    }

    template< class PointType, class LocalDofVectorType >
    void axpy ( const PointType &x,
                const RangeType &rangeFactor,
                LocalDofVectorType& dofs ) const
    {
      const int numBase = numBaseFunctions();
      for( int i = 0; i < numBase; ++i )
      {
        RangeType phi;
        evaluate( i, x, phi );
        dofs[ i ] += phi * rangeFactor;
      }
    }

    template< class PointType, class GeometryJacobianInverseType, 
              class GlobalJacobianRangeType, class LocalDofVectorType >
    void axpy ( const PointType &x,
                const GeometryJacobianInverseType& gjit,
                const GlobalJacobianRangeType &jacFactor,
                LocalDofVectorType& dofs ) const
    {
      GlobalJacobianRangeType jacFactorInv;
      for( int r = 0; r < dimRange; ++r )
        gjit.mtv( jacFactor[ r ], jacFactorInv[ r ] );

      const int numBase = numBaseFunctions();
      for( int i = 0; i < numBase; ++i )
      {
        JacobianRangeType grad;
        jacobian( i, x, grad );
        for( int r = 0; r < dimRange; ++r )
          dofs[ i ] += grad[ r ] * jacFactorInv[ r ];
      }
    }
 
    template< class PointType, class GeometryJacobianInverseType, 
              class GlobalJacobianRangeType, class LocalDofVectorType >
    void axpy ( const PointType &x,
                const GeometryJacobianInverseType& gjit, 
                const RangeType &rangeFactor,
                const GlobalJacobianRangeType &jacFactor,
                LocalDofVectorType& dofs ) const
    {
      GlobalJacobianRangeType jacFactorInv;
      for( int r = 0; r < dimRange; ++r )
        gjit.mtv( jacFactor[ r ], jacFactorInv[ r ] );

      const int numBase = numBaseFunctions();
      for( int i = 0; i < numBase; ++i )
      {
        RangeType phi;
        evaluate( i, x, phi );
        dofs[ i ] += phi * rangeFactor;

        JacobianRangeType grad;
        jacobian( i, x, grad );
        for( int r = 0; r < dimRange; ++r )
          dofs[ i ] += grad[ r ] * jacFactorInv[ r ];
      }
    }

    template< class QuadratureType,
              class LocalDofVectorType,
              class RangeVectorType>
    inline void
    evaluateRanges ( const QuadratureType& quad,
                     const LocalDofVectorType& dofs,
                     RangeVectorType &rangeVector) const
    {
      const size_t quadNop = quad.nop();
      for( size_t qp = 0; qp < quadNop; ++qp )
      {
        evaluate( quad[ qp ], dofs, rangeVector[ qp ] );
      }
    }

    template< class QuadratureType,
              class Geometry,
              class LocalDofVectorType,
              class JacobianRangeVectorType>
    inline void
    evaluateJacobians ( const QuadratureType& quad,
                        const Geometry& geometry,
                        const LocalDofVectorType& dofs,
                        JacobianRangeVectorType &jacVector ) const
    {
      const size_t quadNop = quad.nop();
      for( size_t qp = 0; qp < quadNop; ++qp )
      {
        jacobian( quad[ qp ], 
                  geometry.jacobianInverseTransposed( quad.point( qp ) ),
                  dofs, 
                  jacVector[ qp ] );
      }
    }

    ////////////////////////////////////////////////////
    //  axpyRanges and axpyJacobians 
    ////////////////////////////////////////////////////
    template< class QuadratureType,
              class RangeVectorType,
              class LocalDofVectorType >
    inline void axpyRanges ( const QuadratureType& quad,
                             const RangeVectorType &rangeFactors,
                             LocalDofVectorType& dofs ) const
    {
      const size_t quadNop = quad.nop();
      for( size_t qp = 0; qp < quadNop; ++qp )
      {
        axpy( quad[ qp ], rangeFactors[ qp ], dofs );
      }
    }

    template< class QuadratureType,
              class Geometry,
              class JacobianVectorType,
              class LocalDofVectorType >
    inline void axpyJacobians ( const QuadratureType& quad,
                                const Geometry& geometry,
                                const JacobianVectorType &jacVector,
                                LocalDofVectorType& dofs) const
    {
      const size_t quadNop = quad.nop();
      for( size_t qp = 0; qp < quadNop; ++qp )
      {
        axpy( quad[ qp ], 
              geometry.jacobianInverseTransposed( quad.point( qp ) ),
              jacVector[ qp ],
              dofs );
      }
    }
  };

  /** \} */



  // Implementation of BaseFunctionSetDefault
  // ----------------------------------------

  template< class TraitsImp >
  template< class PointType >
  inline void BaseFunctionSetDefault< TraitsImp >
    ::jacobian ( const int baseFunction, const PointType &x, JacobianRangeType &phi ) const
  {
    FieldVector< deriType, 1 > diffVar;
    deriType &i = diffVar[ 0 ];
    // create temporary variable here 
    RangeType tmp;
    for( i = 0; i < dimCol; ++i )
    {
      asImp().evaluate( baseFunction, diffVar, x, tmp );
      for( int j = 0; j < dimRow; ++j )
        phi[ j ][ i ] = tmp[ j ];
    }
  }


  template< class TraitsImp >
  template< class PointType, class GeometryJacobianInverseType,
            class LocalDofVectorType, class GlobalJacobianRangeType >
  inline void BaseFunctionSetDefault< TraitsImp >
    ::jacobianAll ( const PointType &x,
                    const GeometryJacobianInverseType &gjit, 
                    const LocalDofVectorType &dofs, 
                    GlobalJacobianRangeType &ret ) const
  {
    JacobianRangeType refJacobian( 0 );
    const int numBase = numBaseFunctions();
    for( int i = 0; i < numBase; ++i )
    {
      JacobianRangeType grad;
      jacobian( i, x, grad );
      
      for( int r = 0; r < dimRange; ++r )
        refJacobian[ r ].axpy( dofs[ i ], grad[ r ] );
    }

    for( int r = 0; r < dimRange; ++r )
      gjit.mv( refJacobian[ r ], ret[ r ] );
  }


  template< class TraitsImp >
  template< class PointType, class GeometryJacobianInverseType,
            class GlobalJacobianRangeVectorType >
  inline void BaseFunctionSetDefault< TraitsImp >
    ::jacobianAll ( const PointType &x,
                    const GeometryJacobianInverseType& gjit, 
                    GlobalJacobianRangeVectorType &ret ) const
  {
    const int numBase = numBaseFunctions();
    ret.resize( numBase );
    for( int i = 0; i < numBase; ++i )
    {
      JacobianRangeType grad;
      jacobian( i, x, grad );
      for( int r = 0; r < dimRange; ++r )
        gjit.mv( grad[ r ], ret[ i ][ r ] );
    }
  }


  template< class TraitsImp >
  template< class PointType >
  inline void BaseFunctionSetDefault< TraitsImp >
    ::hessian ( const int baseFunction, const PointType &x, HessianRangeType &hessian ) const
  {
    FieldVector< deriType, 2 > diffVar;
    deriType &i = diffVar[ 0 ];
    deriType &j = diffVar[ 1 ];
    RangeType tmp;
    for( i = 0; i < dimDomain; ++i )
    {
      // We use symmetrized evaluation of the hessian, since calling
      // evaluate is in general quite expensive
      for( j = 0; j < i; ++j )
      {
        asImp().evaluate( baseFunction, diffVar, x, tmp );
        for( int k = 0; k < dimRange; ++k )
          hessian[ k ][ i ][ j ] = hessian[ k ][ j ][ i ] = tmp[ k ];
      }
      assert( j == i );
      asImp().evaluate( baseFunction, diffVar, x, tmp );
      for( int k = 0; k < dimRange; ++k )
        hessian[ k ][ i ][ i ] = tmp[ k ];
    }
  }


  template< class TraitsImp >
  template< class PointType, class Geometry, class LocalDofVectorType, class GlobalHessianRangeType >
  inline void BaseFunctionSetDefault< TraitsImp >
    ::hessianAll ( const PointType &x, const Geometry &geometry,
                   const LocalDofVectorType &dofs, GlobalHessianRangeType &ret ) const
  {
    const int numBase = numBaseFunctions();

    HessianRangeType refHessian( typename HessianRangeType::field_type( 0 ) );
    for( int b = 0; b < numBase; ++b )
    {
      HessianRangeType H;
      hessian( b, x, H );
      for( int r = 0; r < dimRange; ++r )
        refHessian[ r ].axpy( dofs[ b ], H[ r ] );
    }

    const typename Geometry::Jacobian &Jgeo
      = geometry.jacobianInverseTransposed( coordinate( x ) );
    ret = typename GlobalHessianRangeType::field_type( 0 );
    for( int i = 0; i < Geometry::coorddimension; ++i )
    {
      for( int j = 0; j < dimDomain; ++j )
      {
        for( int r = 0; r < dimRange; ++r )
        {
          FieldVector< RangeFieldType, Geometry::coorddimension > tmp;
          Jgeo.mv( refHessian[ r ][ j ], tmp );
          ret[ r ][ i ].axpy( Jgeo[ i ][ j ], tmp );
        }
      }
    }

    // for nonaffine geometries we have to add Du D^2 F^-1
    if( !geometry.affine() )
      DUNE_THROW( NotImplemented, "hessianAll: not implemented for nonaffine geometries." );
  }

} // end namespace Dune 

#endif // #ifndef DUNE_BASEFUNCTIONSETINTERFACE_HH
