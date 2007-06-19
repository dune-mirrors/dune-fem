#ifndef VTKIO_HH
#define VTKIO_HH
#include <dune/grid/io/file/vtk/vtkwriter.hh>
namespace Dune {
template <class DF>
class VTKFunctionWrapper : 
  public VTKWriter<typename DF::FunctionSpaceType::GridPartType::GridType,
                   typename DF::FunctionSpaceType::GridPartType::IndexSetType>::VTKFunction {
public:
  typedef DF DiscreteFunctionType;
  typedef typename DF::LocalFunctionType LocalFunctionType;
  typedef typename DF::FunctionSpaceType FunctionSpaceType;
  typedef typename FunctionSpaceType::DomainType DomainType;
  typedef typename FunctionSpaceType::RangeType RangeType;
  typedef typename FunctionSpaceType::GridPartType GridPartType;
  typedef typename GridPartType::EntityCodim0Type EntityType;
  enum {DimRange=FunctionSpaceType::DimRange,
        DimDomain=FunctionSpaceType::DimDomain};
  //! constructor taking discrete function 
  VTKFunctionWrapper(const DiscreteFunctionType& df,
                     bool vector)
    : discFunc_(df),
      vector_(vector) {}
  //! return number of components
  virtual int ncomps () const {
    return DimRange;
  }
  //! evaluate single component comp in
  //! the entity
  virtual double evaluate (int comp, const EntityType& e, const DomainType& xi) const {
    const LocalFunctionType lf = discFunc_.localFunction(e);
    RangeType val;
    lf.evaluate(xi,val);
    return val[comp];
  }
  //! get name
  virtual std::string name () const {
    return discFunc_.name();
  }
  //! virtual destructor
  virtual ~VTKFunctionWrapper () {}
private:
  const DF& discFunc_;
  const bool vector_;
};
template <class GridPartType>
class VTKIO : 
  public VTKWriter<typename GridPartType::GridType,
                   typename GridPartType::IndexSetType> {
  typedef VTKWriter<typename GridPartType::GridType,
                    typename GridPartType::IndexSetType> BaseType;
  public:
    VTKIO(GridPartType& part) : 
      BaseType(part.grid(),part.indexSet()) {
      }
    template <class DF>
    void addCellData(DF& df,bool vector=false) {
      BaseType::addCellData(new VTKFunctionWrapper<DF>(df,vector)); 
    }
    template <class DF>
    void addVertexData(DF& df,bool vector=false) {
      BaseType::addVertexData(new VTKFunctionWrapper<DF>(df,vector)); 
    }
};
}
#endif
