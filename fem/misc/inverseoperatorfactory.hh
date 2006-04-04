#ifndef DUNE_INVERSEOPERATORFACTORY_HH
#define DUNE_INVERSEOPERATORFACTORY_HH

#include "../operator/common/mapping.hh"
#include "../operator/inverseoperators.hh"

namespace Dune {

/** \brief Abstract factory class for inverse operators
    Concrete factory classes derive from this class and overwrite method
    createOperator. For each concrete inverse operator class, a factory class
    is derived.
 */
template <typename DiscreteFunctionType>
class InverseOperatorFactory {
public:
  //! Definition of a generic mapping with identical domain and range field
  typedef Mapping<typename DiscreteFunctionType::DomainFieldType,
                  typename DiscreteFunctionType::RangeFieldType,
                  DiscreteFunctionType,DiscreteFunctionType> MappingType;

  virtual ~InverseOperatorFactory() {}

  /** \brief Method creating inverse operator from the given operator
      \note The client is responsible for deleting the return inverse operator!
   */
  virtual MappingType* createOperator(const MappingType& op) = 0;

private:

};

/** \brief Concrete factory class for CGInverseOperator
*/
template <typename DiscreteFunctionType>
class CGInverseOperatorFactory : 
  public InverseOperatorFactory<DiscreteFunctionType> {
public:
  //! The mapping type from the base class
  typedef typename InverseOperatorFactory<DiscreteFunctionType>::MappingType MappingType;
  /** \brief Constructor
      Constructs a factory producing CGInverseOperator objects with the
      tolerance specifications here.
  */
  CGInverseOperatorFactory(double redEps,
                           double absLimit,
                           int maxIter,
                           int verbose) :
    redEps_(redEps),
    absLimit_(absLimit),
    maxIter_(maxIter),
    verbose_(verbose) {}

  //* Method returning CGInverseOperator
  virtual MappingType* createOperator(const MappingType& op) {
    return new CGInverseOperator<DiscreteFunctionType>(op, 
                                                       redEps_, 
                                                       absLimit_, 
                                                       maxIter_, 
                                                       verbose_);
  }

private:
  //* Reduce error each step by
  double redEps_;

  //* Minimal error to reach
  double absLimit_;

  //* Number of maximal iterations
  int maxIter_;

  //* Level of output
  int verbose_;
};

} // end namepsace Dune
#endif
