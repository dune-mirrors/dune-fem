#ifndef DUNE_OPERATOR_HH
#define DUNE_OPERATOR_HH

#include "mapping.hh"

namespace Dune
{
/** @defgroup OperatorCommon Operators
    Operators are mappings from function spaces into function spaces.
    @{
 */

/** \brief An abstract operator
 Interface class for Operators. Operators are applied to Functions and
 the result is a Function again. 
*/
template <typename DFieldType, typename RFieldType,
          typename DType , typename RType>
class Operator : public Mapping <DFieldType,RFieldType,DType,RType>
{
protected: 
  //! \brief type of mapping base class 
  typedef Mapping <DFieldType,RFieldType,DType,RType> MappingType;
  
public:
  //- remember template parameters for derived classes  
  typedef DType DomainType;
  typedef RType  RangeType;
  typedef DFieldType DomainFieldType;
  typedef RFieldType RangeFieldType;

  /** \brief Application operator 
      \param[in] arg argument 
      \param[out] dest destination 
      \note This method has to be implemented by all derived classes. 
  */
  virtual void operator() (const DomainType& arg, RangeType& dest) const = 0;
 
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

} // end namespace Dune 
#endif
