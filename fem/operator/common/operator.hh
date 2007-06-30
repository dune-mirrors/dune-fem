#ifndef DUNE_OPERATOR_HH
#define DUNE_OPERATOR_HH

#include "mapping.hh"

namespace Dune
{
/** @defgroup OperatorCommon Operators
  Operators are mappings from function spaces into function spaces.
 */

/** @defgroup Operator Operator Interface
    @ingroup OperatorCommon

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
  typedef Mapping <DFieldType,RFieldType,DType,RType> MappingType;
  
public:
  //! remember template parameters for derived classes  
  typedef DType DomainType;
  typedef RType  RangeType;
  typedef DFieldType DomainFieldType;
  typedef RFieldType RangeFieldType;

  /** \brief Application operator 
     \note This method has to be implemented by all derived classes. 
  */
  virtual void operator() (const DomainType& arg, RangeType& dest) const = 0;
 
protected:
  /** \brief The method apply is the virtual version of the 
      application operator. This method must be implemented by all 
      derived classes. 
  */
  virtual void apply (const DomainType& arg, RangeType& dest) const 
  {
    this->operator() (arg, dest); 
  }
}; // end class Operator 

/** @} end documentation group */

} // end namespace Dune 
#endif
