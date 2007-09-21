#ifndef VTKIO_HH
#define VTKIO_HH

#include <dune/grid/io/file/vtk/vtkwriter.hh>

namespace Dune
{

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
    typedef typename FunctionSpaceType::IteratorType::Entity EntityType;
    enum {DimRange=FunctionSpaceType::DimRange,
          DimDomain=FunctionSpaceType::DimDomain};
    //! constructor taking discrete function 
    VTKFunctionWrapper(const DiscreteFunctionType& df,
                       int component)
      : discFunc_(df),
        component_(component) {}
    //! return number of components
    virtual int ncomps () const {
      return (component_>=0) ? 1 : DimRange;
    }
    //! evaluate single component comp in
    //! the entity
    virtual double evaluate (int comp, const EntityType& e, const DomainType& xi) const {
      const LocalFunctionType lf = discFunc_.localFunction(e);
      RangeType val;
      lf.evaluate(xi,val);
      if (component_<0)
        return val[comp];
      else 
        return val[component_];
    }
    //! get name
    virtual std::string name () const {
      if (component_<0) 
        return discFunc_.name();
      std::stringstream ret;
      ret << discFunc_.name() << "-" << component_;
      return ret.str();
    }
    //! virtual destructor
    virtual ~VTKFunctionWrapper () {}
  private:
    const DF& discFunc_;
    const int component_;
  };



  template< class GridPartImp >
  class VTKIO
  : public VTKWriter< typename GridPartImp :: GridType,
                      typename GridPartImp :: IndexSetType >
  {
  public:
    typedef GridPartImp GridPartType;

    typedef typename GridPartType :: GridType GridType;
    typedef typename GridPartType :: IndexSetType IndexSetType;

  private:
    typedef VTKIO< GridPartType > ThisType;
    typedef VTKWriter< GridType, IndexSetType > BaseType;
    
    const GridPartType& gridPart_;
  public:
    //! constructor  
    VTKIO( const GridPartType &gridPart, VTKOptions::DataMode dm = VTKOptions::conforming )
    : BaseType( gridPart.grid(), gridPart.indexSet() , dm )
    , gridPart_( gridPart )
    {
    }

    //! return grid part 
    const GridPartType& gridPart() const { return gridPart_; }

    template< class DF >
    void addCellData( DF &df,
                      bool vector=false )
    {
      typedef typename DF::FunctionSpaceType FunctionSpaceType;
      enum {DimRange=FunctionSpaceType::DimRange};
      if (vector) {
        BaseType::addCellData(new VTKFunctionWrapper<DF>(df,-1)); 
      } else {
        for (int i=0;i<DimRange;++i) 
          BaseType::addCellData(new VTKFunctionWrapper<DF>(df,i)); 
      }
    }
    
    template< class DF >
    void addVertexData( DF &df,
                        bool vector=false )
    {
      typedef typename DF::FunctionSpaceType FunctionSpaceType;
      enum {DimRange=FunctionSpaceType::DimRange};
      if (vector) {
        BaseType::addVertexData(new VTKFunctionWrapper<DF>(df,-1)); 
      } else {
        for (int i=0;i<DimRange;++i) 
          BaseType::addVertexData(new VTKFunctionWrapper<DF>(df,i)); 
      }
    }
  };
  
}

#endif
