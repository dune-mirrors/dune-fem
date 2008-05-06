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
    enum {dimRange=FunctionSpaceType::dimRange,
          dimDomain=FunctionSpaceType::dimDomain};
    //! constructor taking discrete function 
    VTKFunctionWrapper(const DiscreteFunctionType& df,
                       int component,bool vector)
      : discFunc_(df),
        vector_(vector),
        component_(component) {}
    //! return number of components
    virtual int ncomps () const {
      return (!vector_) ? 1 : dimDomain;
    }
    //! evaluate single component comp in
    //! the entity
    virtual double evaluate (int comp, const EntityType& e, const DomainType& xi) const {
      const LocalFunctionType lf = discFunc_.localFunction(e);
      RangeType val;
      lf.evaluate(xi,val);
      if (vector_)
        return val[comp + component_];
      else 
        return val[component_];
    }
    //! get name
    virtual std::string name () const {
      if (vector_) {
        std::stringstream ret_vec;
        ret_vec << discFunc_.name() << "-vec-" << component_;
        return ret_vec.str(); 
      }
      std::stringstream ret;
      ret << discFunc_.name() << "-" << component_;
      return ret.str();
    }
    //! virtual destructor
    virtual ~VTKFunctionWrapper () {}
  private:
    const DF& discFunc_;
    const bool vector_;
    const int component_;
  };

  //! /brief Output using VTK
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
    void addCellData( DF &df)
    {
      typedef typename DF::FunctionSpaceType FunctionSpaceType;
      enum {dimRange=FunctionSpaceType::dimRange};
      enum {dimDomain=FunctionSpaceType::dimDomain};
      for (int i=0;i<dimRange;++i) 
         BaseType::addCellData(new VTKFunctionWrapper<DF>(df,i,false)); 
    }

    template< class DF >
    void addVectorCellData( DF &df,
                      int startPoint = 0 )
    {
      typedef typename DF::FunctionSpaceType FunctionSpaceType;
      enum {dimRange=FunctionSpaceType::dimRange};
      enum {dimDomain=FunctionSpaceType::dimDomain};
      BaseType::addCellData(new VTKFunctionWrapper<DF>(df,startPoint,true)); 
    }

    template< class DF >
    void addVertexData( DF &df )
    {
      typedef typename DF::FunctionSpaceType FunctionSpaceType;
      enum {dimRange=FunctionSpaceType::dimRange};
      enum {dimDomain=FunctionSpaceType::dimDomain};
      for (int i=0;i<dimRange;++i) 
        BaseType::addVertexData(new VTKFunctionWrapper<DF>(df,i,false)); 
    }

    template< class DF >
    void addVectorVertexData( DF &df,
                      int startPoint = 0 )
    {
      typedef typename DF::FunctionSpaceType FunctionSpaceType;
      enum {dimRange=FunctionSpaceType::dimRange};
      enum {dimDomain=FunctionSpaceType::dimDomain};
      BaseType::addVertexData(new VTKFunctionWrapper<DF>(df,startPoint,true)); 
    }

  };
  
}

#endif
