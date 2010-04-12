#ifndef DUNE_FEM_OPERATOR_HH
#define DUNE_FEM_OPERATOR_HH

#include <dune/fem/operator/common/mapping.hh>

namespace Dune
{

  namespace Fem
  {

    /** \class Operator
     *  \brief abstract operator
     *
     *  Operators map a discrete function onto another discrete function.
     *  Their interface is described by the abstract class Operator.
     *
     *  \tparam  DomainFunction  type of discrete function for the domain
     *  \tparam  RangeFunction   type of discrete function for the range
     *
     *  \interfaceclass
     */
    template< class DomainFunction, class RangeFunction >
    struct Operator
    {
      /** \brief type of discrete function in the operator's domain */
      typedef DomainFunction DomainFunctionType;
      /** \brief type of discrete function in the operator's range */
      typedef RangeFunction RangeFunctionType;

      /** \brief field type of the operator's domain */
      typedef typename DomainFunctionType::RangeFieldType DomainFieldType;
      /** \brief field type of the operator's range */
      typedef typename RangeFunctionType::RangeFieldType RangeFieldType;

      /** \brief application operator
       *
       *  \param[in]   u  argument discrete function
       *  \param[out]  w  destination discrete function
       *
       *  \note This method has to be implemented by all derived classes.
       */
      virtual void operator() ( const DomainFunctionType &u, RangeFunctionType &w ) const = 0;
    };

  } // end namespace Dune::Fem



/** @addtogroup OperatorCommon 
    Operators are mappings from function spaces into function spaces.

    \remarks 
    The most general interface for a mapping is defined by Mapping. From
    Mapping the Function is derived and also Operator is derived. 
    Operator defines the interface for operations on discrete functions
    while Function is an interface for functions like sinus.
    
    @{
 */

/** \brief An abstract operator
 Interface class for Operators. Operators are applied to Functions and
 the result is a Function again. 

 \interfaceclass
*/
template <typename DFieldType, typename RFieldType,
          typename DType , typename RType>
class Operator
: public Mapping <DFieldType,RFieldType,DType,RType>,
  public Fem::Operator< DType, RType >
{
  typedef Fem::Operator< DType, RType > BaseType;

protected: 
  //! \brief type of mapping base class 
  typedef Mapping <DFieldType,RFieldType,DType,RType> MappingType;
  
public:
  //- remember template parameters for derived classes  
  typedef DType DomainType;
  typedef RType  RangeType;
  typedef DFieldType DomainFieldType;
  typedef RFieldType RangeFieldType;

  using BaseType::operator();

#if 0
  /** \brief Application operator 
      \param[in] arg argument 
      \param[out] dest destination 
      \note This method has to be implemented by all derived classes. 
  */
  virtual void operator() (const DomainType& arg, RangeType& dest) const = 0;
#endif
 
protected:
  /** \brief The method apply calls the application operator. The method 
      has to be implemented here, because this method called when a mapping list 
      is evaluated. 
      \param[in] arg argument 
      \param[out] dest destination 
  */
  virtual void apply (const DomainType& arg, RangeType& dest) const 
  {
    this->operator() (arg, dest); 
  }
}; // end class Operator 

///@}

} // end namespace Dune 

#endif // #ifndef DUNE_FEM_OPERATOR_HH
