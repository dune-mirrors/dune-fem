#ifndef VTKIO_HH
#define VTKIO_HH

#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>
#include <dune/fem/gridpart/gridpartview.hh>

namespace Dune
{

  template< class DF >
  class VTKFunctionWrapper
  : public VTKWriter< typename DF::DiscreteFunctionSpaceType::GridPartType::GridViewType >::VTKFunction
  {
    typedef VTKFunctionWrapper< DF > ThisType;

  public:
    typedef DF DiscreteFunctionType;

    typedef typename DiscreteFunctionType::LocalFunctionType LocalFunctionType;
    typedef typename DiscreteFunctionType::FunctionSpaceType FunctionSpaceType;
    typedef typename FunctionSpaceType::DomainType DomainType;
    typedef typename FunctionSpaceType::RangeType RangeType;
    typedef typename FunctionSpaceType::GridPartType GridPartType;
    typedef typename FunctionSpaceType::IteratorType::Entity EntityType;

    static const int dimRange = FunctionSpaceType::dimRange;
    static const int dimDomain = FunctionSpaceType::dimDomain;

    //! constructor taking discrete function 
    VTKFunctionWrapper ( const DiscreteFunctionType& df,
                         int component, bool vector )
    : discFunc_( df ),
      vector_( vector ),
      component_( component )
    {}

    //! virtual destructor
    virtual ~VTKFunctionWrapper ()
    {}

    //! return number of components
    virtual int ncomps () const
    {
      return (!vector_) ? 1 : 3; // dimDomain;
    }

    //! evaluate single component comp in
    //! the entity
    virtual double evaluate (int comp, const EntityType& e, const DomainType& xi) const {
      const LocalFunctionType lf = discFunc_.localFunction(e);
      RangeType val;
      lf.evaluate(xi,val);
      if (vector_) {
        if (comp>dimDomain)
          return 0;
        else;
          return val[comp + component_];
      }
      else 
        return val[component_];
    }

    //! get name
    virtual std::string name () const
    {
      if (vector_) {
        std::stringstream ret_vec;
        ret_vec << discFunc_.name() << "-vec-" << component_;
        return ret_vec.str(); 
      }
      std::stringstream ret;
      ret << discFunc_.name() << "-" << component_;
      return ret.str();
    }

  private:
    const DiscreteFunctionType &discFunc_;
    const bool vector_;
    const int component_;
  };



  //! /brief Output using VTK
  template< class GridPart >
  class VTKIO
  {
    typedef VTKIO< GridPart > ThisType;

    typedef VTKWriter< typename GridPart::GridViewType > VTKWriterType;

  public:
    typedef GridPart GridPartType;

    typedef typename GridPartType::GridType GridType;
    typedef typename GridPartType::IndexSetType IndexSetType;

  protected:
    VTKIO ( const GridPartType &gridPart, VTKWriterType *vtkWriter )
    : gridPart_( gridPart ),
      vtkWriter_( vtkWriter )
    {}

  public:
    //! constructor  
    explicit VTKIO ( const GridPartType &gridPart,
                     VTKOptions::DataMode dm = VTKOptions::conforming )
    : gridPart_( gridPart ),
      vtkWriter_( new VTKWriterType( gridPart.gridView(), dm ) )
    {}

    ~VTKIO ()
    {
      delete vtkWriter_;
    }

    //! return grid part 
    const GridPartType &gridPart () const
    {
      return gridPart_;
    }

    template< class DF >
    void addCellData( DF &df)
    {
      typedef typename DF::FunctionSpaceType FunctionSpaceType;
      enum {dimRange=FunctionSpaceType::dimRange};
      enum {dimDomain=FunctionSpaceType::dimDomain};
      for (int i=0;i<dimRange;++i) 
         vtkWriter_->addCellData(new VTKFunctionWrapper<DF>(df,i,false)); 
    }

    template< class DF >
    void addVectorCellData( DF &df,
                      int startPoint = 0 )
    {
      typedef typename DF::FunctionSpaceType FunctionSpaceType;
      enum {dimRange=FunctionSpaceType::dimRange};
      enum {dimDomain=FunctionSpaceType::dimDomain};
      vtkWriter_->addCellData(new VTKFunctionWrapper<DF>(df,startPoint,true)); 
    }

    template< class DF >
    void addVertexData( DF &df )
    {
      typedef typename DF::FunctionSpaceType FunctionSpaceType;
      enum {dimRange=FunctionSpaceType::dimRange};
      enum {dimDomain=FunctionSpaceType::dimDomain};
      for (int i=0;i<dimRange;++i) 
        vtkWriter_->addVertexData(new VTKFunctionWrapper<DF>(df,i,false)); 
    }

    template< class DF >
    void addVectorVertexData( DF &df, int startPoint = 0 )
    {
      typedef typename DF::FunctionSpaceType FunctionSpaceType;
      enum {dimRange=FunctionSpaceType::dimRange};
      enum {dimDomain=FunctionSpaceType::dimDomain};
      vtkWriter_->addVertexData(new VTKFunctionWrapper<DF>(df,startPoint,true)); 
    }

    void clear ()
    {
      vtkWriter_->clear();
    }

    std::string write ( const std::string &name,
                        VTKOptions::OutputType type = VTKOptions::ascii )
    {
      return vtkWriter_->write( name, type );
    }

    std::string pwrite ( const char *name, const char *path, const char *extendpath,
                         VTKOptions::OutputType type = VTKOptions::ascii )
    {
      return vtkWriter_->pwrite( name, path, extendpath, type );
    }

  private:
    const GridPartType& gridPart_;
    VTKWriterType *vtkWriter_;
  };



  template< class GridPart >
  class SubsamplingVTKIO
  : public VTKIO< GridPart >
  {
    typedef SubsamplingVTKIO< GridPart > ThisType;
    typedef VTKIO< GridPart > BaseType;

    typedef SubsamplingVTKWriter< typename GridPart::GridViewType > VTKWriterType;

  public:
    typedef GridPart GridPartType;

    explicit SubsamplingVTKIO ( const GridPartType &gridPart,
                                unsigned int level = 0,
                                bool coerceToSimplex = false )
    : BaseType( gridPart, new VTKWriterType( gridPart.gridView(), level, coerceToSimplex ) )
    {}
  };
  
}

#endif // #ifndef VTKIO_HH
