#ifndef DUNE_FEM_BASEFUNCTIONFACTORY_HH
#define DUNE_FEM_BASEFUNCTIONFACTORY_HH

#include <dune/common/geometrytype.hh>

#include <dune/fem/space/basefunctions/basefunctioninterface.hh>

namespace Dune
{

  //! \brief Interface class for the generation of base functions.
  //! For every concrete set of base functions, derive your own concrete
  //! base function factory.
  template< class FunctionSpace >
  struct BaseFunctionFactory 
  {
    typedef FunctionSpace FunctionSpaceType;

    typedef BaseFunctionInterface< FunctionSpaceType > BaseFunctionType;

    //! constructor storing geometry type 
    explicit BaseFunctionFactory ( const GeometryType &geo )
    : geo_( geo )
    {}

    virtual ~BaseFunctionFactory() {}

    //! return pointer to num-th base function 
    virtual BaseFunctionType *baseFunction ( int num ) const = 0;

    //! return number of base function 
    virtual int numBaseFunctions () const = 0;

    //! return geometry type to which base functions belong
    const GeometryType &geometry () const { return geo_; }

  private:
    GeometryType geo_;
  };

} // namespace Dune

#endif // #ifndef DUNE_FEM_BASEFUNCTIONFACTORY_HH
