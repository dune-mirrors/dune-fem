#ifndef DUNE_FEM_VTKIO_HH
#define DUNE_FEM_VTKIO_HH

#include <dune/common/typetraits.hh>

#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/fem/version.hh>

#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>

#include <dune/fem/gridpart/gridpartview.hh>

namespace Dune
{

  // Internal Forward Declarations
  // -----------------------------

  template< class GridPart, bool subsampling = false >
  class VTKIO;



  // Compatibility definiton
  // -----------------------

#if ! DUNE_VERSION_NEWER(DUNE_GRID,2,1,0)
typedef VTKOptions VTK;
#endif



  // VTKFunctionWrapper
  // ------------------

  template< class DF >
  class VTKFunctionWrapper
  : public VTKWriter< typename DF::DiscreteFunctionSpaceType::GridPartType::GridViewType >::VTKFunction
  {
    typedef VTKFunctionWrapper< DF > ThisType;

  public:
    typedef DF DiscreteFunctionType;

    typedef typename DiscreteFunctionType::LocalFunctionType LocalFunctionType;
    typedef typename DiscreteFunctionType::DiscreteFunctionSpaceType
      DiscreteFunctionSpaceType;

    static const int dimRange = DiscreteFunctionSpaceType::dimRange;
    static const int dimDomain = DiscreteFunctionSpaceType::dimDomain;

    typedef typename DiscreteFunctionSpaceType::DomainType DomainType;
    typedef typename DiscreteFunctionSpaceType::RangeType RangeType;

    typedef typename DiscreteFunctionSpaceType::GridPartType GridPartType;
    typedef typename DiscreteFunctionSpaceType::IteratorType::Entity EntityType;

    typedef typename EntityType::Geometry::LocalCoordinate LocalCoordinateType;


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
    virtual double evaluate (int comp, const EntityType& e, const LocalCoordinateType& xi) const {
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



  // VTKIOBase
  // ---------

  //! /brief Output using VTK
  template< class GridPart, bool subsampling >
  class VTKIOBase
  {
    typedef VTKIOBase< GridPart, subsampling > ThisType;

  protected:
    typedef typename GridPart::GridViewType GridViewType;

    class VTKWriter;
    class SubsamplingVTKWriter;

    typedef typename SelectType< subsampling, SubsamplingVTKWriter, VTKWriter >::Type
      VTKWriterType;

  public:
    typedef GridPart GridPartType;

    typedef typename GridPartType::GridType GridType;
    typedef typename GridPartType::IndexSetType IndexSetType;

  protected :
    VTKIOBase ( const GridPartType &gridPart, VTKWriterType *vtkWriter )
    : gridPart_( gridPart ),
      vtkWriter_( vtkWriter )
    {}

  public:
    ~VTKIOBase ()
    {
      delete vtkWriter_;
    }

    //! return grid part 
    const GridPartType &gridPart () const
    {
      return gridPart_;
    }

    template< class DF >
    void addCellData( DF &df )
    {
      static const int dimRange = DF::FunctionSpaceType::dimRange;
      for( int i = 0;i < dimRange; ++i )
        vtkWriter_->addCellData( new VTKFunctionWrapper< DF >( df, i, false ) );
    }

    template< class DF >
    void addVectorCellData( DF &df, int startPoint = 0 )
    {
      vtkWriter_->addCellData( new VTKFunctionWrapper< DF >( df, startPoint, true ) );
    }

    template< class DF >
    void addVertexData( DF &df )
    {
      static const int dimRange = DF::FunctionSpaceType::dimRange;
      for( int i = 0;i < dimRange; ++i )
        vtkWriter_->addVertexData( new VTKFunctionWrapper< DF >( df, i, false ) );
    }

    template< class DF >
    void addVectorVertexData( DF &df, int startPoint = 0 )
    {
      vtkWriter_->addVertexData( new VTKFunctionWrapper< DF >( df, startPoint, true ) );
    }

    void clear ()
    {
      vtkWriter_->clear();
    }

    std::string write ( const std::string &name,
                        VTK::OutputType type = VTK::ascii )
    {
      return vtkWriter_->write( name, type );
    }

    std::string pwrite ( const std::string &name,
                         const std::string &path,
                         const std::string &extendpath,
                         VTK::OutputType type = VTK::ascii )
    {
      return vtkWriter_->pwrite( name.c_str(), path.c_str(), extendpath.c_str(), type );
    }

    std::string write ( const std::string &name,
                        VTK::OutputType type,
                        const int rank, 
                        const int size )
    {
      return vtkWriter_->write( name, type, rank, size );
    }

  private:
    const GridPartType &gridPart_;
    VTKWriterType *vtkWriter_;
  };



  // VTKIOBase::VTKWriter
  // -----------------------

  template< class GridPart, bool subsampling >
  class VTKIOBase< GridPart, subsampling >::VTKWriter
  : public Dune::VTKWriter< GridViewType >
  {
    typedef Dune::VTKWriter< GridViewType > BaseType;

  public:
    // make all write methods public for data convert
    using BaseType::write;
    using BaseType::pwrite;

    //! constructor  
    VTKWriter( const GridPartType &gridPart, 
               VTK::DataMode dm = VTK::conforming )
    : BaseType( gridPart.gridView(), dm )
    {}
  };



  // VTKIOBase::SubSamplingVTKWriter
  // --------------------------

  template< class GridPart, bool subsampling >
  class VTKIOBase< GridPart, subsampling >::SubsamplingVTKWriter
  : public Dune::SubsamplingVTKWriter< GridViewType >
  {
    typedef Dune::SubsamplingVTKWriter< GridViewType > BaseType;

  public:
    // make all write methods public for data convert
    using BaseType::write;
    using BaseType::pwrite;

    //! constructor  
    SubsamplingVTKWriter( const GridPartType &gridPart, 
                          const int level, 
                          bool coerceToSimplex = false )
    : BaseType( gridPart.gridView(), level, coerceToSimplex )
    {}
  };



  // VTKIO (without subsampling)
  // ---------------------------

  template< class GridPart >
  class VTKIO< GridPart, false >
  : public VTKIOBase< GridPart, false >
  {
    typedef VTKIO< GridPart > ThisType;
    typedef VTKIOBase< GridPart, false > BaseType;

    typedef typename BaseType::VTKWriterType VTKWriterType;

  public:
    typedef typename BaseType::GridPartType GridPartType;

    explicit VTKIO ( const GridPartType &gridPart,
                     VTK::DataMode dm = VTK::conforming )
    : BaseType( gridPart, new VTKWriterType( gridPart, dm ) )
    {}
  };



  // VTKIO (with subsampling)
  // ------------------------

  template< class GridPart >
  class VTKIO< GridPart, true >
  : public VTKIOBase< GridPart , true >
  {
    typedef VTKIO< GridPart, true > ThisType;
    typedef VTKIOBase< GridPart, true > BaseType;

     typedef typename BaseType::VTKWriterType VTKWriterType;

  public:
    typedef typename BaseType::GridPartType GridPartType;

    explicit VTKIO ( const GridPartType &gridPart,
                     unsigned int level = 0,
                     bool coerceToSimplex = false )
    : BaseType( gridPart, new VTKWriterType( gridPart, level, coerceToSimplex ) )
    {}
  };



  // SubsamplingVTKIO
  // ----------------

  template< class GridPart >
  class SubsamplingVTKIO
  : public VTKIO< GridPart, true >
  {
    typedef SubsamplingVTKIO< GridPart > ThisType;
    typedef VTKIO< GridPart, true > BaseType;

  public:
    typedef typename BaseType::GridPartType GridPartType;

    explicit SubsamplingVTKIO ( const GridPartType &gridPart,
                                unsigned int level = 0,
                                bool coerceToSimplex = false )
    : BaseType( gridPart, level, coerceToSimplex )
    {}
  };

}

#endif // #ifndef DUNE_FEM_VTKIO_HH
