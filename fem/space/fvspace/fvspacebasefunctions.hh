#ifndef DUNE_FVSPACEBASEFUNCTIONS_HH
#define DUNE_FVSPACEBASEFUNCTIONS_HH

//- Dune includes 
#include <dune/common/misc.hh>
#include <dune/grid/common/grid.hh>

//- Dune-Fem includes 
#include <dune/fem/space/common/basefunctionfactory.hh>
#include <dune/fem/space/common/geometryconversion.hh>

namespace Dune {

//! definition of FVBaseFunction, implementation via specialization 
template<class FunctionSpaceType, GeometryIdentifier::IdentifierType ElType, int polOrd> 
class FVBaseFunction;
         
//! Piecewise const base functions for all types of elements 
template<class FunctionSpaceType, GeometryIdentifier::IdentifierType ElType>
class FVBaseFunction < FunctionSpaceType, ElType, 0 >  
: public BaseFunctionInterface<FunctionSpaceType> 
{
  enum { dimRange = FunctionSpaceType::dimRange };
  int baseNum_;

  typedef typename FunctionSpaceType::DomainType DomainType;
  typedef typename FunctionSpaceType::RangeType RangeType;
  typedef typename FunctionSpaceType::RangeFieldType RangeFieldType;
  
public:
  FVBaseFunction (int baseNum) : 
    BaseFunctionInterface<FunctionSpaceType> (),
    baseNum_ ( baseNum ) 
  { 
    assert((baseNum_ >= 0) || (baseNum_ < dimRange));
  }
  
  virtual void evaluate ( const FieldVector<deriType, 0> &diffVariable, 
                          const DomainType & x, RangeType & phi) const 
  {
    phi = 0;
    phi[baseNum_] = 1;
  }

  virtual void evaluate ( const FieldVector<deriType, 1> &diffVariable, 
                          const DomainType & x, RangeType & phi) const 
  {
    phi = 0;
  }

  virtual void evaluate ( const FieldVector<deriType, 2> &diffVariable, 
                          const DomainType & x, RangeType & phi) const 
  {
    phi = 0;
  }
};

//! default definition stays empty because implementation via
//! specialization
template <GeometryIdentifier::IdentifierType ElType, int polOrd ,int dimrange > 
struct FVDefinition;

//! FV Space Definition for all types 
template <GeometryIdentifier::IdentifierType ElType, int dimrange > 
struct FVDefinition< ElType, 0, dimrange > 
{
  //! number of base functions for polOrd = 0
  enum { numOfBaseFct = dimrange }; 
};

//! Factory class for base functions
template <class ScalarFunctionSpaceImp, int polOrd>
class FVBaseFunctionFactory :
  public BaseFunctionFactory<ScalarFunctionSpaceImp>
{
  // at the moment only for polOrd 0
  CompileTimeChecker<polOrd == 0> only_implemented_for_polOrd_0;
    
public:
  typedef ScalarFunctionSpaceImp FunctionSpaceType;
  typedef BaseFunctionInterface<FunctionSpaceType> BaseFunctionType;
public:
  FVBaseFunctionFactory(GeometryType geo) :
    BaseFunctionFactory<FunctionSpaceType>(geo)
  {}

  virtual ~FVBaseFunctionFactory() {}

  virtual BaseFunctionType* baseFunction(int i) const 
  {
    switch (GeometryIdentifier::fromGeo(this->geometry())) 
    {
      case GeometryIdentifier::Line:
        return new 
          FVBaseFunction<FunctionSpaceType, GeometryIdentifier::Line, polOrd>(i);
      case GeometryIdentifier::Triangle:
        return new 
          FVBaseFunction<FunctionSpaceType, GeometryIdentifier::Triangle, polOrd>(i);
      case GeometryIdentifier::Tetrahedron:
        return new 
          FVBaseFunction<FunctionSpaceType, GeometryIdentifier::Tetrahedron, polOrd>(i);
      case GeometryIdentifier::Quadrilateral:
        return new 
          FVBaseFunction<FunctionSpaceType, GeometryIdentifier::Quadrilateral, polOrd>(i);
      case GeometryIdentifier::Hexahedron:
        return new 
          FVBaseFunction<FunctionSpaceType, GeometryIdentifier::Hexahedron, polOrd>(i);
      case GeometryIdentifier::Prism:
        return new 
          FVBaseFunction<FunctionSpaceType, GeometryIdentifier::Prism, polOrd>(i);
      case GeometryIdentifier::Pyramid:
        return new 
          FVBaseFunction<FunctionSpaceType, GeometryIdentifier::Pyramid, polOrd>(i);
      default:
        DUNE_THROW(NotImplemented, 
                   "The chosen geometry type is not implemented");
    }
    return 0;
  }
  
  virtual int numBaseFunctions() const 
  {
    const int dimRange = FunctionSpaceType::dimRange;
    switch (GeometryIdentifier::fromGeo(this->geometry())) 
    {
      case GeometryIdentifier::Line:
        return  
          FVDefinition<GeometryIdentifier::Line, polOrd, dimRange>::numOfBaseFct;
      case GeometryIdentifier::Triangle:
        return  
          FVDefinition<GeometryIdentifier::Triangle, polOrd, dimRange>::numOfBaseFct;
      case GeometryIdentifier::Tetrahedron:
        return  
          FVDefinition<GeometryIdentifier::Tetrahedron, polOrd, dimRange>::numOfBaseFct;
      case GeometryIdentifier::Quadrilateral:
        return  
          FVDefinition<GeometryIdentifier::Quadrilateral, polOrd, dimRange>::numOfBaseFct;
      case GeometryIdentifier::Hexahedron:
        return  
          FVDefinition<GeometryIdentifier::Hexahedron, polOrd, dimRange>::numOfBaseFct;
      case GeometryIdentifier::Prism:
        return  
          FVDefinition<GeometryIdentifier::Prism, polOrd, dimRange>::numOfBaseFct;
      case GeometryIdentifier::Pyramid:
        return 
          FVDefinition<GeometryIdentifier::Pyramid, polOrd, dimRange>::numOfBaseFct;
      default:
        DUNE_THROW(NotImplemented, 
                   "The chosen geometry type is not implemented");
    }
    return 0;
  }
};

} // end namespace Dune
#endif
