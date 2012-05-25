#ifndef DUNE_FVSPACEBASEFUNCTIONS_HH
#define DUNE_FVSPACEBASEFUNCTIONS_HH

//- Dune includes 
#include <dune/common/misc.hh>
#include <dune/grid/common/grid.hh>

//- Dune-Fem includes 
#include <dune/fem/space/basefunctions/basefunctioninterface.hh>
#include <dune/fem/space/basefunctions/basefunctionfactory.hh>

namespace Dune
{

  //! definition of FVBaseFunction, implementation via specialization 
  template< class FunctionSpaceType, int polOrd >
  class FVBaseFunction;
           
  //! Piecewise const base functions for all types of elements 
  template< class FunctionSpaceType >
  class FVBaseFunction< FunctionSpaceType, 0 >
  : public BaseFunctionInterface< FunctionSpaceType >
  {
    typedef BaseFunctionInterface< FunctionSpaceType > BaseType;

    enum { dimRange = FunctionSpaceType::dimRange };
    int baseNum_;

    typedef typename FunctionSpaceType::DomainType DomainType;
    typedef typename FunctionSpaceType::RangeType RangeType;
    typedef typename FunctionSpaceType::RangeFieldType RangeFieldType;
    
  public:
    FVBaseFunction ( const int baseNum )
    : baseNum_ ( baseNum ) 
    { 
      assert( (baseNum_ >= 0) && (baseNum_ < dimRange) );
    }
    
    virtual void evaluate ( const FieldVector<int, 0> &diffVariable, 
                            const DomainType & x, RangeType & phi) const 
    {
      phi = 0;
      phi[baseNum_] = 1;
    }

    virtual void evaluate ( const FieldVector<int, 1> &diffVariable, 
                            const DomainType & x, RangeType & phi) const 
    {
      phi = 0;
    }

    virtual void evaluate ( const FieldVector<int, 2> &diffVariable, 
                            const DomainType & x, RangeType & phi) const 
    {
      phi = 0;
    }
  };



  //! Factory class for base functions
  template< class FunctionSpace, int polOrd >
  class FVBaseFunctionFactory
  : public BaseFunctionFactory< FunctionSpace >
  {
    dune_static_assert( polOrd == 0, "Only implemented for PolOrd=0." );
      
  public:
    typedef FunctionSpace FunctionSpaceType;
    typedef BaseFunctionInterface< FunctionSpaceType > BaseFunctionType;

    static const int dimRange = FunctionSpaceType::dimRange;

    FVBaseFunctionFactory ( const GeometryType &type )
    : BaseFunctionFactory< FunctionSpaceType >( type )
    {}

    virtual BaseFunctionType *baseFunction ( const int i ) const
    {
      return new FVBaseFunction< FunctionSpaceType, polOrd >( i );
    }
    
    virtual int numBaseFunctions() const 
    {
      return dimRange;
    }
  };

} // end namespace Dune

#endif // #ifndef DUNE_FVSPACEBASEFUNCTIONS_HH
