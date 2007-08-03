#ifndef DUNE_BASEFUNCTIONSETINTERFACE_HH
#define DUNE_BASEFUNCTIONSETINTERFACE_HH

//- local includes 
#include <dune/fem/space/common/basefunctioninterface.hh>

namespace Dune{

/** @addtogroup BaseFunction
    @{
*/

//------------------------------------------------------------------------
//-
//-  --BaseFunctionSetInterface
//-
//------------------------------------------------------------------------

/** 
   \brief 
    
   The BaseFunctionSetInterface describes the interface for
   BaseFunctionSets. 
  
   Why has the BaseFunctionInterface class virtual methods?
  
   Because the we want to have different base functions in 
   at the same time and we havent found a way to do this with 
   Barton-Nackman. But there is a solution which dosent cost to much
   efficiency, for example in FastBaseFunctionSet all evaluations of the
   BaseFunction are cached in a vector, though the virtual method evaluate
   of the BaseFunctionImp has only to be called once for the same quadrature
   rule. If you change the quadrature rule then on the first call the
   evaluations are stored again.
   This method brings us flexebility and effeciency. 
*/  
template<class BaseFunctionSetTraits> 
class BaseFunctionSetInterface 
{
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
  BaseFunctionSetInterface () {}

  //! \brief empty destructor 
  virtual ~BaseFunctionSetInterface() {}
  
  /** \brief number of base functions 
      \return number of base functions 
  */
  int numBaseFunctions() const 
  {
    CHECK_INTERFACE_IMPLEMENTATION(asImp().numBaseFunctions());
    return asImp().numBaseFunctions();
  }

  /** \brief interface method for evaluation, evaluation of jacobian, and hessian of base functions 
      \param[in] baseFunct number of base function to evaluate 
      \param[in] diffVariable length determines the derivative (i.e. 0 is
             evaluate, 1 is gradient, 2 hessian, ... ) and the value
             determines the component 
      \param[in] x point in reference element for evaluation
      \param[out] phi return value 
  */
  template <int diffOrd>
  void evaluate (const int baseFunct, 
                 const FieldVector<deriType, diffOrd> &diffVariable, 
                 const DomainType & x, 
                 RangeType & phi ) const 
  {
    CHECK_AND_CALL_INTERFACE_IMPLEMENTATION(
        asImp().evaluate( baseFunct, diffVariable, x, phi ));
  }

  /** \brief interface method for evaluation, evaluation of jacobian, and hessian of base functions 
      \param[in] baseFunct number of base function to evaluate 
      \param[in] diffVariable length determines the derivative (i.e. 0 is
             evaluate, 1 is gradient, 2 hessian, ... ) and the value
             determines the component 
      \param[in] quad Quadrature 
      \param[in] quadPoint number of quadrature point 
      \param[out] phi return value 
  */
  template <int diffOrd, class QuadratureImp>
  void evaluate (const int baseFunct, 
                 const FieldVector<deriType, diffOrd> &diffVariable, 
                 const QuadratureImp& quad, 
                 const int quadPoint, 
                 RangeType & phi ) const 
  {
    CHECK_AND_CALL_INTERFACE_IMPLEMENTATION(
      asImp().evaluate( baseFunct, diffVariable, quad, quadPoint, phi ));
  }

  /** \brief return type of geometry
      \return GeometryType of base function set
  */
  inline GeometryType geometryType () const
  {
    CHECK_INTERFACE_IMPLEMENTATION( asImp().geometryType() );
    return asImp().geometryType();
  }

  /** \brief evaluation of base functions 
      \param[in] baseFunct number of base function to evaluate 
      \param[in] x coordiante in reference element to evaluate base function on 
      \param[out] phi return value, i.e. value of base function 
  */ 
  void evaluate(const int baseFunct, 
                const DomainType & x, 
                RangeType & phi) const
  {
    CHECK_AND_CALL_INTERFACE_IMPLEMENTATION(
      asImp().evaluate(baseFunct, x , phi));
  }

  /** \brief evaluation of base functions 
      \param[in] baseFunct number of base function to evaluate 
      \param[in] quad Quadrature 
      \param[in] quadPoint number of quadrature point 
      \param[out] phi return value, i.e. value of base function 
  */ 
  template <class QuadratureImp>
  void evaluate(const int baseFunct, 
                const QuadratureImp & quad, 
                const int quadPoint, 
                RangeType & phi) const 
  {
    CHECK_AND_CALL_INTERFACE_IMPLEMENTATION(
      asImp().evaluate( baseFunct, quad, quadPoint, phi ));
  }

  /** \brief evaluation of jacobian of base function
      \param[in] baseFunct number of base function to evaluate jacobian 
      \param[in] x coordiante in reference element to evaluate base function on 
      \param[out] phi return value, i.e. gradient on reference elememnt (multiply with jacobianInverseTransposed)
  */ 
  void jacobian(const int baseFunct, 
                const DomainType& x, 
                JacobianRangeType & phi) const 
  {
    CHECK_AND_CALL_INTERFACE_IMPLEMENTATION(
      asImp().jacobian( baseFunct, x, phi ));
  }

  /** \brief default jacobian using the evaluate interface 
      \param[in] baseFunct number of base function to evaluate jacobian 
      \param[in] quad Quadrature 
      \param[in] quadPoint number of quadrature point 
      \param[out] phi return value, i.e. gradient on reference elememnt (multiply with jacobianInverseTransposed)
  */ 
  template <class QuadratureImp>
  void jacobian ( const int baseFunct, 
                  const QuadratureImp & quad, 
                  const int quadPoint, 
                  JacobianRangeType & phi ) const 
  {
    CHECK_AND_CALL_INTERFACE_IMPLEMENTATION(
      asImp().jacobian( baseFunct, quad, quadPoint, phi ));
  }

  /** \brief evaluate basefunction and multiply with factor
      \param[in] baseFunct number of base functions to evaluate 
      \param[in] x local point in reference element 
      \param[in] factor factor to multiply with 
      \return return scalar product between base function and factor 
    */
  RangeFieldType evaluateSingle(const int baseFunct, 
                         const DomainType& x,
                         const RangeType& factor) const
  {
    CHECK_INTERFACE_IMPLEMENTATION(
      asImp().evaluateSingle( baseFunct, x, factor ));
    return asImp().evaluateSingle( baseFunct, x, factor );
  }

  /** \brief evaluate basefunction and multiply with factor
      \param[in] baseFunct number of base functions to evaluate 
      \param[in] quad Quadrature 
      \param[in] quadPoint number of quadrature point 
      \param[in] factor factor to multiply with 
      \return return scalar product between base function and factor 
  */
  template <class QuadratureType>
  RangeFieldType evaluateSingle(const int baseFunct,
                         const QuadratureType& quad, 
                         const int quadPoint,
                         const RangeType& factor) const
  {
    CHECK_INTERFACE_IMPLEMENTATION(
      asImp().evaluateSingle( baseFunct, quad, quadPoint, factor ));
    return asImp().evaluateSingle( baseFunct, quad, quadPoint, factor );
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

private:
  //! Barton-Nackman trick 
  BaseFunctionSetType &asImp() { 
    return static_cast<BaseFunctionSetType&>(*this); 
  }
  
  //! Barton-Nackman trick 
  const BaseFunctionSetType &asImp() const { 
    return static_cast<const BaseFunctionSetType&>(*this); 
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
template<class BaseFunctionSetTraits> 
class BaseFunctionSetDefault  : 
  public BaseFunctionSetInterface <BaseFunctionSetTraits> 
{
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

public:
  //! set the default diffVar Types 
  BaseFunctionSetDefault () : 
    BaseFunctionSetInterface<BaseFunctionSetTraits> ()
  {
    for(int i=0; i<dimCol; i++)
      jacobianDiffVar_[i] = i;
  };
 
  //! destructor 
  virtual ~BaseFunctionSetDefault() {}

  /** \brief @copydoc BaseFunctionSetInterface::evaluate */
  void evaluate(const int baseFunct, 
                const DomainType & x, 
                RangeType & phi) const
  {
    CHECK_AND_CALL_INTERFACE_IMPLEMENTATION(
      asImp().evaluate(baseFunct, diffVariable_ , x , phi));
  }

  /** \brief @copydoc BaseFunctionSetInterface::evaluate */
  template <class QuadratureImp>
  void evaluate(const int baseFunct, 
                const QuadratureImp & quad, 
                const int quadPoint, 
                RangeType & phi) const 
  {
    CHECK_AND_CALL_INTERFACE_IMPLEMENTATION(
      asImp().evaluate( baseFunct, diffVariable_ , quad, quadPoint, phi ));
  }

  /** \brief @copydoc BaseFunctionSetInterface::jacobian */
  void jacobian(const int baseFunct, 
                const DomainType& x, 
                JacobianRangeType & phi) const 
  {
    RangeType tmp;
    for(int i=0; i<dimCol; i++)
    {
      CHECK_AND_CALL_INTERFACE_IMPLEMENTATION(
        asImp().evaluate( baseFunct, jacobianDiffVar_[i] , x , tmp ));
      for(int j=0; j<dimRow; j++)
        phi[j][i] = tmp[j];
    }
  }

  /** \brief @copydoc BaseFunctionSetInterface::jacobian */
  template <class QuadratureImp>
  void jacobian ( const int baseFunct, 
                  const QuadratureImp & quad, 
                  const int quadPoint, 
                  JacobianRangeType & phi ) const 
  {
    RangeType tmp;
  
    for(int i=0; i<dimCol; i++)
    {
      asImp().evaluate( baseFunct, jacobianDiffVar_[i] , quad, quadPoint, tmp );
      for(int j=0; j<dimRow; j++)
        phi[j][i] = tmp[j];
    }
  }

  /** \brief @copydoc BaseFunctionSetInterface::evaluateSingle */
  RangeFieldType evaluateSingle(const int baseFunct, 
                         const DomainType& xLocal,
                         const RangeType& factor) const
  {
    RangeType phi(0.);
    evaluate(baseFunct, xLocal, phi);
    return phi*factor;
  }

  /** \brief @copydoc BaseFunctionSetInterface::evaluateSingle */
  template <class QuadratureType>
  RangeFieldType evaluateSingle(const int baseFunct,
                         const QuadratureType& quad, 
                         const int quadPoint,
                         const RangeType& factor) const
  {
    return evaluateSingle(baseFunct, quad.point(quadPoint), factor);
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

private: 
  //! just diffVariable for evaluation of the functions 
  const FieldVector<deriType, 0> diffVariable_;

  //! diff variable for jacobian evaluation 
  FieldVector<deriType, 1> jacobianDiffVar_[dimCol];

  //! Barton-Nackman trick 
  BaseFunctionSetType &asImp() { return static_cast<BaseFunctionSetType&>(*this); }
  //! Barton-Nackman trick 
  const BaseFunctionSetType &asImp() const 
  { 
    return static_cast<const BaseFunctionSetType&>(*this); 
  }
};
///@}
} // end namespace Dune 
#endif
