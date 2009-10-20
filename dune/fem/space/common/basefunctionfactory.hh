#ifndef DUNE_BASEFUNCTIONFACTORY_HH
#define DUNE_BASEFUNCTIONFACTORY_HH

#include <dune/common/geometrytype.hh>

#include "basefunctioninterface.hh"

namespace Dune {

  //! \brief Interface class for the generation of base functions.
  //! For every concrete set of base functions, derive your own concrete
  //! base function factory.
  template <class FunctionSpaceImp>
  class BaseFunctionFactory 
  {
  public:
    typedef BaseFunctionInterface<FunctionSpaceImp> BaseFunctionType;
  public:
    //! constructor storing geometry type 
    BaseFunctionFactory(GeometryType geo) :
      geo_(geo)
    {}

    virtual ~BaseFunctionFactory() {}

    //! return pointer to num-th base function 
    virtual BaseFunctionType* baseFunction(int num) const = 0;

    //! return number of base function 
    virtual int numBaseFunctions() const = 0;

    //! return geometry type to which base functions belong
    GeometryType geometry() const { return geo_; }
  private:
    GeometryType geo_;
  };

} // end namespace Dune
#endif
