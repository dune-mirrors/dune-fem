#ifndef DUNE_BASEFUNCTIONSETINTERFACE_HH
#define DUNE_BASEFUNCTIONSETINTERFACE_HH

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

#if HAVE_DUNE_GEOMETRY
#include <dune/geometry/type.hh>
#else
#include <dune/common/geometrytype.hh>
#endif

#include <dune/fem/version.hh>
#include <dune/fem/misc/bartonnackmaninterface.hh>

namespace Dune
{
  namespace Fem 
  {

  // BaseFunctionSetInterface
  // ------------------------

  /** \class BaseFunctionSetInterface
   *  \ingroup BaseFunction
   *  \brief interface class for base function sets
   */ 
  template< class TraitsImp > 
  class BaseFunctionSetInterface
  : public BartonNackmanInterface< BaseFunctionSetInterface< TraitsImp >, typename TraitsImp::BaseFunctionSetType >
  {
    typedef BaseFunctionSetInterface< TraitsImp > ThisType;
    typedef BartonNackmanInterface< BaseFunctionSetInterface< TraitsImp >, typename TraitsImp::BaseFunctionSetType > BaseType;

  public:
    //! type of the traits
    typedef TraitsImp Traits;

    //! type of base function set implementation (Barton-Nackman) 
    typedef typename Traits::BaseFunctionSetType BaseFunctionSetType;

    //! type of function space 
    typedef typename Traits::FunctionSpaceType FunctionSpaceType;
    //! type of domain vector, i.e. element of R^dimDomain 
    typedef typename FunctionSpaceType::DomainType DomainType;
    //! type of range vector, i.e. element of R^dimRange 
    typedef typename FunctionSpaceType::RangeType RangeType;
    //! type of domain vector components, e.g., double 
    typedef typename FunctionSpaceType::DomainFieldType DomainFieldType;
    //! type of range vector components, e.g., double 
    typedef typename FunctionSpaceType::RangeFieldType RangeFieldType;
    //! range type of gradient 
    typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;
    //! range type of second derivate 
    typedef typename FunctionSpaceType::HessianRangeType  HessianRangeType;

    //! dimension of domain
    static const int dimDomain = FunctionSpaceType::dimDomain;
    //! dimension of range
    static const int dimRange  = FunctionSpaceType::dimRange;

  protected:
    using BaseType::asImp;

  public:
    /** \todo please doc me */
    template< class Point, class DofVector >
    void axpy ( const Point &x, const RangeType &valueFactor, DofVector &dofs ) const
    {
      CHECK_AND_CALL_INTERFACE_IMPLEMENTATION( asImp().axpy( x, valueFactor, dofs ) );
    }

    /** \todo please doc me */
    template< class Point, class GeometryJacobianInverse, class GlobalJacobianRange, class DofVector >
    void axpy ( const Point &x, const GeometryJacobianInverse &gjit,
                const GlobalJacobianRange &jacobianFactor, DofVector &dofs ) const
    {
      CHECK_AND_CALL_INTERFACE_IMPLEMENTATION( asImp().axpy( x, gjit, jacobianFactor, dofs ) );
    }

    /** \todo please doc me */
    template< class Point, class GeometryJacobianInverse, class GlobalJacobianRange, class DofVector >
    void axpy ( const Point &x, const GeometryJacobianInverse &gjit,
                const RangeType &valueFactor, const GlobalJacobianRange &jacobianFactor,
                DofVector &dofs ) const
    {
      CHECK_AND_CALL_INTERFACE_IMPLEMENTATION( asImp().axpy( x, gjit, valueFactor, jacobianFactor, dofs ) );
    }

    /** \brief evaluate all base function
     *
     *  Evaluates all base functions in a point (specified in local
     *  coordinates).
     *  
     *  \param[in]   x             point within reference element to evaluate the
     *                             base function in
     *  \param[out]  values        values of all base functions
     */
    template< class Point, class RangeArray >
    void evaluateAll ( const Point &x, RangeArray &values ) const
    {
      CHECK_AND_CALL_INTERFACE_IMPLEMENTATION( asImp().evaluateAll( x, values ) );
    }

    /** \brief evaluate derivative of a linear combination
     *
     *  Evaluates a partial derivative of a linear combination of the base
     *  functions in a point (specified in local coordinates).
     *  If the derivative is of order 0, the base functions themselves are
     *  evaluated.
     *
     *  \note This method only makes sense for shape function sets.
     *
     *  \param[in]   diffVariable  vector describing the partial derivative to
     *                             evaluate
     *  \param[in]   x             point within reference element to evaluate the
     *                             base function in
     *  \param[in]   dofs          local degree of freedom
     *  \param[out]  value         value of the derivative of the linear
     *                             combination
     */
    template< int diffOrder, class Point, class DofVector >
    void evaluateAll ( const FieldVector< int, diffOrder > &diffVariable,
                       const Point &x, const DofVector &dofs, RangeType &value ) const
    {
      CHECK_AND_CALL_INTERFACE_IMPLEMENTATION
        ( asImp().evaluateAll( diffVariable, x, dofs, value ) );
    }

    /** \brief evaluate a linear combination
     *
     *  Evaluates a linear combination of the base functions in a point
     *  (specified in local coordinates).
     *  
     *  \param[in]   x             point within reference element to evaluate the
     *                             base function in
     *  \param[in]   dofs          local degree of freedom
     *  \param[out]  value         value of the linear combination
     */
    template< class Point, class DofVector >
    void evaluateAll ( const Point &x, const DofVector &dofs, RangeType &value ) const
    {
      CHECK_AND_CALL_INTERFACE_IMPLEMENTATION
        ( asImp().evaluateAll( x, dofs, value ) );
    }

    /** \brief obtain type of geometry
     * 
     *  \returns GeometryType of the base function set
     */
    GeometryType geometryType () const
    {
      CHECK_INTERFACE_IMPLEMENTATION( asImp().geometryType() );
      return asImp().geometryType();
    }

    /** \brief \todo please doc me */
    template< class Point, class Geometry, class DofVector, class GlobalHessianRange >
    void hessianAll ( const Point &x, const Geometry &geometry,
                      const DofVector &dofs, GlobalHessianRange &hessian ) const
    {
      CHECK_AND_CALL_INTERFACE_IMPLEMENTATION
        ( asImp().hessianAll( x, geometry, dofs, hessian ) );
    }

    /** \brief \todo please doc me */
    template< class Point, class Geometry, class GlobalHessianRangeArray >
    void hessianAll ( const Point &x, const Geometry &geometry, GlobalHessianRangeArray &hessians ) const
    {
      CHECK_AND_CALL_INTERFACE_IMPLEMENTATION
        ( asImp().hessianAll( x, geometry, hessians ) );
    }

    /** \brief \todo please doc me */
    template< class Point, class GeometryJacobianInverse, class DofVector, class GlobalJacobianRange >
    void jacobianAll ( const Point &x, const GeometryJacobianInverse &gjit, 
                       const DofVector &dofs, GlobalJacobianRange &jacobian ) const
    {
      CHECK_AND_CALL_INTERFACE_IMPLEMENTATION
        ( asImp().evaluateAll( x, gjit, dofs, jacobian ) );
    }

    /** \brief \todo please doc me */
    template< class Point, class GeometryJacobianInverse, class GlobalJacobianRangeArray >
    void jacobianAll ( const Point &x, const GeometryJacobianInverse &gjit,
                       GlobalJacobianRangeArray &jacobians ) const
    {
      CHECK_AND_CALL_INTERFACE_IMPLEMENTATION
        ( asImp().evaluateAll( x, gjit, jacobians ) );
    }

    /** \brief size of this base function set
     *  \returns number of base functions
     */
    size_t size () const
    {
      CHECK_INTERFACE_IMPLEMENTATION( asImp().size() );
      return asImp().size();
    }


    // deprecated methods

    /** \brief evaluate a derivative of a base function
     *
     *  Evaluates a partial derivative of a base function in a point (specified
     *  in local coordinates).
     *  If the derivative is  of order 0, the base function itself is
     *  evaluated.
     *  
     *  \note This method only makes sense for shape function sets.
     *
     *  \param[in]   baseFunction  number of base function to evaluate
     *  \param[in]   diffVariable  vector describing the partial derivative to
     *                             evaluate
     *  \param[in]   x             point within reference element to evaluate the
     *                             base function in
     *  \param[out]  value         return value, i.e., value of the base function
     */
    template< int diffOrder, class Point >
    DUNE_DEPRECATED
    void evaluate ( const int baseFunction,
                    const FieldVector< int, diffOrder > &diffVariable,
                    const Point &x,
                    RangeType &value ) const
    {
      CHECK_AND_CALL_INTERFACE_IMPLEMENTATION
        ( asImp().evaluate( baseFunction, diffVariable, x, value ) );
    }

    /** \brief evaluate the base function
     *
     *  Evaluates the base function in a point (specified in local coordinates).
     *  
     *  \param[in]   baseFunction  number of base function to evaluate
     *  \param[in]   x             point within reference element to evaluate the
     *                             base function in
     *  \param[out]  value         return value, i.e., value of the base function
     */
    template< class Point >
    DUNE_DEPRECATED
    void evaluate ( const int baseFunction, const Point &x, RangeType &value ) const
    {
      CHECK_AND_CALL_INTERFACE_IMPLEMENTATION
        ( asImp().evaluate( baseFunction, x , value ) );
    }

    /** \brief \todo please doc me */
    template< class Point >
    void hessian ( const int baseFunction, const Point &x, HessianRangeType &hessian ) const
    {
      CHECK_AND_CALL_INTERFACE_IMPLEMENTATION
        ( asImp().hessian( baseFunction, x, hessian ) );
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
     *  \param[out]  jacobian      return value, i.e., value of the Jacobian
     */
    template< class Point >
    DUNE_DEPRECATED
    void jacobian ( const int baseFunction, const Point &x, JacobianRangeType &jacobian ) const
    {
      CHECK_AND_CALL_INTERFACE_IMPLEMENTATION
        ( asImp().jacobian( baseFunction, x, jacobian ) );
    }

    /** \brief number of base functions 
        \return number of base functions 
    */
    DUNE_DEPRECATED
    int numBaseFunctions () const
    {
      return size();
    }
  };



  // BaseFunctionSetDefault
  // ----------------------

  /** \class BaseFunctionSetDefault
   *  \ingroup BaseFunction
   *  \brief default implementation of the BaseFunctionSetInterface
   */
  template< class TraitsImp >
  class BaseFunctionSetDefault
  : public BaseFunctionSetInterface< TraitsImp >
  {
    typedef BaseFunctionSetDefault< TraitsImp > ThisType;
    typedef BaseFunctionSetInterface< TraitsImp > BaseType;

  public:
    //! type of the traits
    typedef TraitsImp Traits;

    typedef typename BaseType::BaseFunctionSetType BaseFunctionSetType;
    typedef typename BaseType::FunctionSpaceType FunctionSpaceType;

    typedef typename BaseType::JacobianRangeType JacobianRangeType;
    typedef typename BaseType::HessianRangeType  HessianRangeType;

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
    using BaseType::size;

    /** \copydoc Dune::Fem::BaseFunctionSetInterface::evaluate(const int baseFunction,const FieldVector<int,diffOrder> &diffVariable,const Point &x,RangeType &value) const */
    template< int diffOrder, class Point >
    DUNE_DEPRECTED
    void evaluate ( const int baseFunction,
                    const FieldVector< int, diffOrder > &diffVariable,
                    const Point &x,
                    RangeType &value ) const
    {
      asImp().evaluate( baseFunction, diffVariable, coordinate( x ), value );
    }
    
    /** \copydoc Dune::Fem::BaseFunctionSetInterface::evaluate(const int baseFunction,const Point &x,RangeType &value) const */
    template< class Point >
    DUNE_DEPRECTED
    void evaluate ( const int baseFunction, const Point &x, RangeType &value ) const
    {
      FieldVector< int, 0 > diffVariable;
      asImp().evaluate( baseFunction, diffVariable, x, value );
    }

#if 0
    // This declaration is only to avoid a problem with specialization
    template< class Point >
    void evaluate ( const unsigned int baseFunction, const Point &x, RangeType &value ) const DUNE_DEPRECATED
    {
      evaluate( int( baseFunction ), x, value );
    }
#endif

    template< int diffOrder, class Point, class DofVector >
    void evaluateAll ( const FieldVector< int, diffOrder > &diffVariable,
                       const Point &x, const DofVector &dofs, RangeType &value ) const
    {
      value = RangeType( 0 );
      const size_t numBaseFunctions = size();
      for( size_t i = 0; i < numBaseFunctions; ++i )
      {
        RangeType phi;
        asImp().evaluate( i, diffVariable, x, phi );
        value.axpy( dofs[ i ], phi );
      }
    }

    template< class Point, class DofVector >
    void evaluateAll ( const Point &x, const DofVector &dofs, RangeType &value ) const
    {
      value = RangeType( 0 );
      const size_t numBaseFunctions = size();
      for( size_t i = 0; i < numBaseFunctions; ++i )
      {
        RangeType phi;
        asImp().evaluate( i, x, phi );
        value.axpy( dofs[ i ], phi );
      }
    }

    template< class Point, class RangeArray >
    void evaluateAll ( const Point &x, RangeArray &values ) const
    {
      const size_t numBaseFunctions = size();
      values.resize( numBaseFunctions );
      for( size_t i = 0; i < numBaseFunctions; ++i )
        asImp().evaluate( i, x, values[ i ] );
    }

    /** \copydoc Dune::Fem::BaseFunctionSetInterface::jacobian(const int baseFunction,const Point &x,JacobianRangeType &jacobian) const */
    template< class Point >
    DUNE_DEPRECATED
    void jacobian ( const int baseFunction, const Point &x, JacobianRangeType &jacobian ) const;

    template< class Point, class GeometryJacobianInverse, class DofVector, class GlobalJacobianRange >
    void jacobianAll ( const Point &x, const GeometryJacobianInverse &gjit, 
                       const DofVector &dofs, GlobalJacobianRange &jacobian ) const;

    template< class Point, class GeometryJacobianInverse, class GlobalJacobianRangeArray >
    void jacobianAll ( const Point &x, const GeometryJacobianInverse &gjit,
                       GlobalJacobianRangeArray &jacobians ) const;

    template< class Point >
    DUNE_DEPRECTED
    void hessian ( const int baseFunction, const Point &x, HessianRangeType &hessian ) const;

    template< class Point, class Geometry, class DofVector, class GlobalHessianRange >
    void hessianAll ( const Point &x, const Geometry &geometry,
                      const DofVector &dofs, GlobalHessianRange &hessian ) const;

    template< class Point, class Geometry, class GlobalHessianRangeArray >
    void hessianAll ( const Point &x, const Geometry &geometry, GlobalHessianRangeArray &hessians ) const;

    template< class Point, class DofVector >
    void axpy ( const Point &x, const RangeType &valueFactor, DofVector &dofs ) const
    {
      const size_t numBaseFunctions = size();
      for( size_t i = 0; i < numBaseFunctions; ++i )
      {
        RangeType tmpValue;
        asImp().evaluate( i, x, tmpValue );
        dofs[ i ] += tmpValue * valueFactor;
      }
    }

    template< class Point, class GeometryJacobianInverse, class GlobalJacobianRange, class DofVector >
    void axpy ( const Point &x, const GeometryJacobianInverse &gjit,
                const GlobalJacobianRange &jacobianFactor, DofVector &dofs ) const
    {
      JacobianRangeType tmpJacobianFactor;
      for( int r = 0; r < dimRange; ++r )
        gjit.mtv( jacobianFactor[ r ], tmpJacobianFactor[ r ] );

      const size_t numBaseFunctions = size();
      for( size_t i = 0; i < numBaseFunctions; ++i )
      {
        JacobianRangeType tmpJacobian;
        asImp().jacobian( i, x, tmpJacobian );
        for( int r = 0; r < dimRange; ++r )
          dofs[ i ] += tmpJacobian[ r ] * tmpJacobianFactor[ r ];
      }
    }
 
    template< class Point, class GeometryJacobianInverse, class GlobalJacobianRange, class DofVector >
    void axpy ( const Point &x, const GeometryJacobianInverse &gjit,
                const RangeType &valueFactor, const GlobalJacobianRange &jacobianFactor,
                DofVector &dofs ) const
    {
      JacobianRangeType tmpJacobianFactor;
      for( int r = 0; r < dimRange; ++r )
        gjit.mtv( jacobianFactor[ r ], tmpJacobianFactor[ r ] );

      const size_t numBaseFunctions = size();
      for( size_t i = 0; i < numBaseFunctions; ++i )
      {
        RangeType tmpValue;
        asImp().evaluate( i, x, tmpValue );
        dofs[ i ] += tmpValue * valueFactor;

        JacobianRangeType tmpJacobian;
        asImp().jacobian( i, x, tmpJacobian );
        for( int r = 0; r < dimRange; ++r )
          dofs[ i ] += tmpJacobian[ r ] * tmpJacobianFactor[ r ];
      }
    }

    template< class QuadratureType,
              class LocalDofVectorType,
              class RangeVectorType>
    void
    evaluateRanges ( const QuadratureType& quad,
                     const LocalDofVectorType& dofs,
                     RangeVectorType &rangeVector) const
    {
      const size_t quadNop = quad.nop();
      for( size_t qp = 0; qp < quadNop; ++qp )
      {
        asImp().evaluateAll( quad[ qp ], dofs, rangeVector[ qp ] );
      }
    }

    template< class QuadratureType,
              class Geometry,
              class LocalDofVectorType,
              class JacobianRangeVectorType>
    void
    evaluateJacobians ( const QuadratureType& quad,
                        const Geometry& geometry,
                        const LocalDofVectorType& dofs,
                        JacobianRangeVectorType &jacVector ) const
    {
      const size_t quadNop = quad.nop();
      for( size_t qp = 0; qp < quadNop; ++qp )
      {
        asImp().jacobianAll( quad[ qp ], 
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
    void axpyRanges ( const QuadratureType& quad,
                             const RangeVectorType &rangeFactors,
                             LocalDofVectorType& dofs ) const
    {
      const size_t quadNop = quad.nop();
      for( size_t qp = 0; qp < quadNop; ++qp )
        asImp().axpy( quad[ qp ], rangeFactors[ qp ], dofs );
    }

    template< class QuadratureType,
              class Geometry,
              class JacobianVectorType,
              class LocalDofVectorType >
    void axpyJacobians ( const QuadratureType& quad,
                                const Geometry& geometry,
                                const JacobianVectorType &jacVector,
                                LocalDofVectorType& dofs) const
    {
      const size_t quadNop = quad.nop();
      for( size_t qp = 0; qp < quadNop; ++qp )
      {
        asImp().axpy( quad[ qp ], 
                      geometry.jacobianInverseTransposed( quad.point( qp ) ),
                      jacVector[ qp ],
                      dofs );
      }
    }
  };



  // Implementation of BaseFunctionSetDefault
  // ----------------------------------------

  template< class TraitsImp >
  template< class Point >
  inline void BaseFunctionSetDefault< TraitsImp >
    ::jacobian ( const int baseFunction, const Point &x, JacobianRangeType &jacobian ) const
  {
    FieldVector< int, 1 > diffVariable;
    int &i = diffVariable[ 0 ];
    for( i = 0; i < dimCol; ++i )
    {
      RangeType tmp;
      asImp().evaluate( baseFunction, diffVariable, x, tmp );
      for( int j = 0; j < dimRange; ++j )
        asImp().jacobian[ j ][ i ] = tmp[ j ];
    }
  }


  template< class TraitsImp >
  template< class Point, class GeometryJacobianInverse, class DofVector, class GlobalJacobianRange >
  inline void BaseFunctionSetDefault< TraitsImp >
    ::jacobianAll ( const Point &x, const GeometryJacobianInverse &gjit, 
                    const DofVector &dofs, GlobalJacobianRange &jacobian ) const
  {
    JacobianRangeType refJacobian( 0 );
    const size_t numBaseFunctions = size();
    for( size_t i = 0; i < numBaseFunctions; ++i )
    {
      JacobianRangeType tmp;
      asImp().jacobian( i, x, tmp );
      
      for( int r = 0; r < dimRange; ++r )
        refJacobian[ r ].axpy( dofs[ i ], tmp[ r ] );
    }

    for( int r = 0; r < dimRange; ++r )
      gjit.mv( refJacobian[ r ], jacobian[ r ] );
  }


  template< class TraitsImp >
  template< class Point, class GeometryJacobianInverse, class GlobalJacobianRangeArray >
  inline void BaseFunctionSetDefault< TraitsImp >
    ::jacobianAll ( const Point &x, const GeometryJacobianInverse &gjit,
                    GlobalJacobianRangeArray &jacobians ) const
  {
    const int numBaseFunctions = size();
    jacobians.resize( numBaseFunctions );
    for( size_t i = 0; i < numBaseFunctions; ++i )
    {
      JacobianRangeType tmp;
      asImp().jacobian( i, x, tmp );
      for( int r = 0; r < dimRange; ++r )
        gjit.mv( tmp[ r ], jacobians[ i ][ r ] );
    }
  }


  template< class TraitsImp >
  template< class Point >
  inline void BaseFunctionSetDefault< TraitsImp >
    ::hessian ( const int baseFunction, const Point &x, HessianRangeType &hessian ) const
  {
    FieldVector< int, 2 > diffVar;
    int &i = diffVar[ 0 ];
    int &j = diffVar[ 1 ];
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
  template< class Point, class Geometry, class DofVector, class GlobalHessianRange >
  inline void BaseFunctionSetDefault< TraitsImp >
    ::hessianAll ( const Point &x, const Geometry &geometry,
                   const DofVector &dofs, GlobalHessianRange &hessian ) const
  {
    const size_t numBaseFunctions = size();

    HessianRangeType refHessian( typename HessianRangeType::field_type( 0 ) );
    for( size_t b = 0; b < numBaseFunctions; ++b )
    {
      HessianRangeType tmp;
      asImp().hessian( b, x, tmp );
      for( int r = 0; r < dimRange; ++r )
        refHessian[ r ].axpy( dofs[ b ], tmp[ r ] );
    }

    const typename Geometry::Jacobian &gjit
      = geometry.jacobianInverseTransposed( coordinate( x ) );
    hessian = typename HessianRangeType::field_type( 0 );
    for( int j = 0; j < Geometry::coorddimension; ++j )
    {
      for( int k = 0; k < dimDomain; ++k )
      {
        for( int r = 0; r < dimRange; ++r )
        {
          FieldVector< RangeFieldType, Geometry::coorddimension > tmp;
          gjit.mv( refHessian[ r ][ k ], tmp );
          hessian[ r ][ j ].axpy( gjit[ j ][ k ], tmp );
        }
      }
    }

    // for nonaffine geometries we have to add Du D^2 F^-1
    if( !geometry.affine() )
      DUNE_THROW( NotImplemented, "hessianAll: not implemented for nonaffine geometries." );
  }

  template< class TraitsImp >
  template< class Point, class Geometry, class GlobalHessianRangeArray >
  inline void BaseFunctionSetDefault< TraitsImp >
   ::hessianAll ( const Point &x, const Geometry &geometry, GlobalHessianRangeArray &hessians ) const
  {
    const int numBaseFunctions = size();
    hessians.resize( numBaseFunctions );
    const typename Geometry::Jacobian &gjit
      = geometry.jacobianInverseTransposed( coordinate( x ) );
    for( size_t i = 0; i < numBaseFunctions; ++i )
    {
      HessianRangeType refHessian;
      asImp().hessian( i, x, refHessian );
      hessians[ i ] = typename HessianRangeType::field_type( 0 );
      for( int j = 0; j < Geometry::coorddimension; ++j )
      {
        for( int k = 0; k < dimDomain; ++k )
        {
          for( int r = 0; r < dimRange; ++r )
          {
            FieldVector< RangeFieldType, Geometry::coorddimension > tmp;
            gjit.mv( refHessian[ r ][ k ], tmp );
            hessians[ i ][ r ][ j ].axpy( gjit[ j ][ k ], tmp );
          }
        }
      }
    }
  }

  } // end namespace Fem

} // end namespace Dune 

#endif // #ifndef DUNE_BASEFUNCTIONSETINTERFACE_HH
