#ifndef DUNE_BASEFUNCTIONSETINTERFACE_HH
#define DUNE_BASEFUNCTIONSETINTERFACE_HH

//- local includes 
#include <dune/fem/space/common/basefunctioninterface.hh>

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
  template< class BaseFunctionSetTraits > 
  class BaseFunctionSetInterface 
  {
  private:
    typedef BaseFunctionSetInterface< BaseFunctionSetTraits > ThisType;
    
  public:
    //! type of function space 
    typedef typename BaseFunctionSetTraits::FunctionSpaceType FunctionSpaceType;
    //! type of base function set implementation 
    typedef typename BaseFunctionSetTraits::BaseFunctionSetType BaseFunctionSetType;
    //! type of domain vector, i.e. element of R^DimDomain 
    typedef typename FunctionSpaceType::DomainType DomainType;
    //! type of range vector, i.e. element of R^DimRange 
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
    enum { DimDomain = FunctionSpaceType::DimDomain };
    //! dimension of range 
    enum { DimRange  = FunctionSpaceType::DimRange  };

    //! type of BaseFunctions 
    typedef BaseFunctionInterface<FunctionSpaceType> BaseFunctionType;

  public:
    //! \brief empty constructor
    inline BaseFunctionSetInterface ()
    {
    }

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
    template< int diffOrd >
    inline void evaluate ( const int baseFunction,
                           const FieldVector< deriType, diffOrd > &diffVariable,
                           const DomainType &x,
                           RangeType &phi ) const
    {
      CHECK_AND_CALL_INTERFACE_IMPLEMENTATION
        ( asImp().evaluate( baseFunction, diffVariable, x, phi ) );
    }

    /** \brief evaluate a derivative of the base function in a quadrature point
     *
     *  Evaluates a partial derivative of the base function in a quadrature
     *  point. If the derivative is of order 0, the base function itself is
     *  evaluated.
     *  
     *  \param[in]   baseFunction  number of base function to evaluate
     *  \param[in]   diffVariable  vector describing the partial derivative to
     *                             evaluate
     *  \param[in]   quadrature    Quadrature to use
     *  \param[in]   quadPoint     number of the evaluation point within the
     *                             quadrature
     *  \param[out]  phi           return value, i.e., value of the base function
     */
    template< int diffOrd, class QuadratureType >
    inline void evaluate ( const int baseFunction,
                           const FieldVector< deriType, diffOrd > &diffVariable,
                           const QuadratureType &quadrature,
                           const int quadPoint,
                           RangeType &phi ) const
    {
      CHECK_AND_CALL_INTERFACE_IMPLEMENTATION
        ( asImp().evaluate( baseFunction, diffVariable, quadrature, quadPoint, phi ) );
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
    inline void evaluate ( const int baseFunction,
                           const DomainType &x,
                           RangeType &phi ) const
    {
      CHECK_AND_CALL_INTERFACE_IMPLEMENTATION
        ( asImp().evaluate( baseFunction, x , phi ) );
    }

    /** \brief evaluate a derivative of the base function in a quadrature point
     *
     *  Evaluates the base function in a quadrature point.
     *  
     *  \param[in]   baseFunction  number of base function to evaluate
     *  \param[in]   quadrature    Quadrature to use
     *  \param[in]   quadPoint     number of the evaluation point within the
     *                             quadrature
     *  \param[out]  phi           return value, i.e., value of the base function
     */
    template< class QuadratureType >
    inline void evaluate ( const int baseFunction,
                           const QuadratureType &quadrature,
                           const int quadPoint,
                           RangeType &phi ) const
    {
      CHECK_AND_CALL_INTERFACE_IMPLEMENTATION
        ( asImp().evaluate ( baseFunction, quadrature, quadPoint, phi ) );
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
    inline void jacobian ( const int baseFunction,
                           const DomainType &x,
                           JacobianRangeType &phi ) const 
    {
      CHECK_AND_CALL_INTERFACE_IMPLEMENTATION
        ( asImp().jacobian( baseFunction, x, phi ) );
    }
    

    /** \brief evaluate the Jacobian of the base function
     *
     *  Evaluates the Jacobian of base function in a quadrature point.
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
     *  \param[in]   quadrature    Quadrature to use
     *  \param[in]   quadPoint     number of the evaluation point within the
     *                             quadrature
     *  \param[out]  phi           return value, i.e., value of the Jacobian
     */
    template< class QuadratureType >
    inline void jacobian ( const int baseFunction, 
                           const QuadratureType &quadrature,
                           const int quadPoint,
                           JacobianRangeType &phi ) const
    {
      CHECK_AND_CALL_INTERFACE_IMPLEMENTATION
        ( asImp().jacobian( baseFunction, quadrature, quadPoint, phi ) );
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
    RangeFieldType evaluateSingle ( const int baseFunction,
                                    const DomainType &x,
                                    const RangeType &psi ) const
    {
      CHECK_INTERFACE_IMPLEMENTATION
        ( asImp().evaluateSingle( baseFunction, x, psi ) );
      return asImp().evaluateSingle( baseFunction, x, psi );
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
     *  \param[in]  quadrature    Quadrature to use
     *  \param[in]  quadPoint     number of the evaluation point within the
     *                            quadrature
     *  \param[in]  psi           vector to multiply the function value with
     *
     *  \returns the scalar product between the value of the base function and
     *           the specified vector
     */
    template< class QuadratureType >
    inline RangeFieldType evaluateSingle ( const int baseFunction,
                                           const QuadratureType &quadrature,
                                           const int quadPoint,
                                           const RangeType &psi ) const
    {
      CHECK_INTERFACE_IMPLEMENTATION
        ( asImp().evaluateSingle( baseFunction, quadrature, quadPoint, psi ) );
      return asImp().evaluateSingle( baseFunction, quadrature, quadPoint, psi );
    }

    /** \brief evaluate gradient of basefunction on given entity (uses
        jacobianInverseTransposed) and multiply with factor, return is RangeFieldType 
        \param[in] baseFunct number of base functions to evaluate jacobian  
        \param[in] entity Entity gradient of base function is evaluated on 
        \param[in] x local point in reference element 
        \param[in] factor factor to multiply with 
        \return return scalar product between gradient of base function and factor 
    */
    template <class Entity>
    RangeFieldType evaluateGradientSingle(const int baseFunct,
                                   const Entity& en,
                                   const DomainType& x,
                                   const JacobianRangeType& factor) const
    {
      CHECK_INTERFACE_IMPLEMENTATION(
        asImp().evaluateGradientSingle( baseFunct, en, x, factor ));
      return asImp().evaluateGradientSingle( baseFunct, en, x, factor );
    }

    /** \brief evaluate gradient of basefunction on given entity (uses
        jacobianInverseTransposed) and multiply with factor, return is RangeFieldType 
        \param[in] baseFunct number of base functions to evaluate jacobian  
        \param[in] entity Entity gradient of base function is evaluated on 
        \param[in] quad Quadrature 
        \param[in] quadPoint number of quadrature point 
        \param[in] factor factor to multiply with 
        \return return scalar product between gradient of base function and factor 
    */
    template <class Entity, class QuadratureType>
    RangeFieldType evaluateGradientSingle(const int baseFunct,
                                   const Entity& en,
                                   const QuadratureType& quad, 
                                   const int quadPoint,
                                   const JacobianRangeType& factor) const
    {
      CHECK_INTERFACE_IMPLEMENTATION(
        asImp().evaluateGradientSingle( baseFunct, en, quad, quadPoint, factor ));
      return asImp().evaluateGradientSingle( baseFunct, en, quad, quadPoint, factor );
    }

  protected:
    // Barton-Nackman trick 
    inline const BaseFunctionSetType &asImp() const
    {
      return static_cast< const BaseFunctionSetType & >( *this );
    }

    // Barton-Nackman trick 
    inline BaseFunctionSetType &asImp()
    {
      return static_cast< BaseFunctionSetType & >( *this );
    }
  };



  //------------------------------------------------------------------------
  //-  
  //-  --BaseFunctionSetDefault
  //-
  //------------------------------------------------------------------------
  /** \brief The BaseFunctionSetDefault implements some methods of the interface BaseFunctionSetInterface
      using the interface of the BaseFunctionSetSlimInterface.
  */
  template< class BaseFunctionSetTraits >
  class BaseFunctionSetDefault
  : public BaseFunctionSetInterface< BaseFunctionSetTraits >
  {
  private:
    typedef BaseFunctionSetDefault< BaseFunctionSetTraits > ThisType;
    typedef BaseFunctionSetInterface< BaseFunctionSetTraits > BaseType;

  public:
    typedef typename BaseFunctionSetTraits::BaseFunctionSetType BaseFunctionSetType;
    typedef typename BaseFunctionSetTraits::FunctionSpaceType FunctionSpaceType;
    typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;

    //! number of rows of jacobian range type 
    enum { dimRow = JacobianRangeType::rows };
    //! number of columns of jacobian range type 
    enum { dimCol = JacobianRangeType::cols };
    
    typedef typename FunctionSpaceType::DomainType DomainType ;
    typedef typename FunctionSpaceType::RangeType RangeType ;

    typedef typename FunctionSpaceType::HessianRangeType  HessianRangeType;
    typedef typename FunctionSpaceType::RangeFieldType RangeFieldType;

  protected:
    using BaseType :: asImp;

  public:
    //! set the default diffVar Types 
    inline BaseFunctionSetDefault ()
    : BaseType()
    {
#if 0
      for( int i = 0; i < dimCol; ++i )
        jacobianDiffVar_[i] = i;
#endif
    }
   
    /** \copydoc Dune::BaseFunctionSetInterface::evaluate(const int baseFunction,const FieldVector<deriType,diffOrd> &diffVariable,const QuadratureType &quadrature,const int quadPoint,RangeType &phi) const */
    template< int diffOrd, class QuadratureType >
    inline void evaluate ( const int baseFunction,
                           const FieldVector< deriType, diffOrd > &diffVariable,
                           const QuadratureType &quadrature,
                           const int quadPoint,
                           RangeType &phi ) const
    {
      asImp().evaluate( baseFunction, diffVariable, quadrature.point( quadPoint ), phi );
    }

    /** \copydoc Dune::BaseFunctionSetInterface::evaluate(const int baseFunction,const DomainType &x,RangeType &phi) const */
    inline void evaluate ( const int baseFunction,
                           const DomainType &x,
                           RangeType &phi ) const
    {
      FieldVector< deriType, 0 > diffVar;
      asImp().evaluate( baseFunction, diffVar, x, phi );
    }

    /** \copydoc Dune::BaseFunctionSetInterface::evaluate(const int baseFunction,const QuadratureType &quadrature,const int quadPoint,RangeType &phi) const */
    template< class QuadratureType >
    inline void evaluate ( const int baseFunction,
                           const QuadratureType &quadrature,
                           const int quadPoint,
                           RangeType &phi ) const
    {
      FieldVector< deriType, 0 > diffVar;
      asImp().evaluate( baseFunction, diffVar, quadrature, quadPoint, phi );
    }

    /** \copydoc Dune::BaseFunctionSetInterface::jacobian(const int baseFunction,const DomainType &x, JacobianRangeType &phi) const */
    inline void jacobian ( const int baseFunction,
                           const DomainType &x,
                           JacobianRangeType &phi ) const
    {
      FieldVector< deriType, 1 > diffVar;
      deriType &i = diffVar[ 0 ];
      for( i = 0; i < dimCol; ++i )
      {
        RangeType tmp;
        asImp().evaluate( baseFunction, diffVar, x, tmp );
        for( int j = 0; j < dimRow; ++j )
          phi[ j ][ i ] = tmp[ j ];
      }
    }

    /** \copydoc Dune::BaseFunctionSetInterface::jacobian(const int baseFunction,const QuadratureType &quadrature,const int quadPoint,JacobianRangeType &phi) const */
    template< class QuadratureType >
    inline void jacobian ( const int baseFunction,
                           const QuadratureType &quadrature,
                           const int quadPoint,
                           JacobianRangeType &phi ) const
    {
      FieldVector< deriType, 1 > diffVar;
      deriType &i = diffVar[ 0 ];
      for( i = 0; i < dimCol; ++i )
      {
        RangeType tmp;
        asImp().evaluate( baseFunction, diffVar, quadrature, quadPoint, tmp );
        for( int j = 0; j < dimRow; ++j )
          phi[ j ][ i ] = tmp[ j ];
      }
    }

    /** \copydoc Dune::BaseFunctionSetInterface::evaluateSingle(const int baseFunction,const DomainType &x,const RangeType &psi) const */
    inline RangeFieldType evaluateSingle ( const int baseFunction,
                                           const DomainType &x,
                                           const RangeType &psi ) const
    {
      RangeType phi;
      evaluate( baseFunction, x, phi );
      return phi * psi;
    }

    /** \copydoc Dune::BaseFunctionSetInterface::evaluateSingle(const int baseFunction,const QuadratureType &quadrature,const int quadPoint,const RangeType &psi) const */
    template< class QuadratureType >
    inline RangeFieldType evaluateSingle ( const int baseFunction,
                                           const QuadratureType &quadrature,
                                           const int quadPoint,
                                           const RangeType &psi ) const
    {
      RangeType phi;
      evaluate( baseFunction, quadrature, quadPoint, phi );
      return phi * psi;
    }

    /** \brief @copydoc BaseFunctionSetInterface::evaluateGradientSingle */
    template <class Entity>
    RangeFieldType evaluateGradientSingle(const int baseFunct,
                                   const Entity& en,
                                   const DomainType& xLocal,
                                   const JacobianRangeType& factor) const
    {
      JacobianRangeType gradPhi(0.);
      jacobian(baseFunct, xLocal, gradPhi);

      RangeFieldType result = 0;
      for (int i = 0; i < FunctionSpaceType::DimRange; ++i) 
      {
        DomainType gradScaled(0.);
        en.geometry().jacobianInverseTransposed(xLocal).
          umv(gradPhi[i], gradScaled);
        result += gradScaled*factor[i];
      }
      return result;
    }

    /** \brief @copydoc BaseFunctionSetInterface::evaluateGradientSingle */
    template <class Entity, class QuadratureType>
    RangeFieldType evaluateGradientSingle(const int baseFunct,
                                   const Entity& en,
                                   const QuadratureType& quad, 
                                   const int quadPoint,
                                   const JacobianRangeType& factor) const
    {
      return evaluateGradientSingle(baseFunct, en, quad.point(quadPoint),factor);
    }

#if 0
  private: 
    // just diffVariable for evaluation of the functions 
    const FieldVector<deriType, 0> diffVariable_;

    // diff variable for jacobian evaluation 
    FieldVector<deriType, 1> jacobianDiffVar_[dimCol];
#endif
  };

  /** \} */

} // end namespace Dune 

#endif
