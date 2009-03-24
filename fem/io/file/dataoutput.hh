#ifndef DUNE_DATAOUTPUT_HH
#define DUNE_DATAOUTPUT_HH

#ifndef USE_VTKWRITER
#define USE_VTKWRITER 1
#endif

//- Dune includes 
#include <dune/fem/io/file/iointerface.hh>
#include <dune/fem/io/file/iotuple.hh>
#include <dune/fem/io/parameter.hh>
#include <dune/fem/solver/timeprovider.hh>
#include <dune/fem/space/common/loadbalancer.hh>
#include <dune/fem/gridpart/adaptiveleafgridpart.hh>
#include <dune/fem/space/lagrangespace.hh>
#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/quadrature/cachequad.hh>
#include <dune/fem/space/dgspace/dgadaptiveleafgridpart.hh>

#if USE_VTKWRITER
#include <dune/fem/io/file/vtkio.hh>
#include <dune/fem/operator/projection/vtxprojection.hh>
#endif

#ifndef USE_GRAPE 
// define whether to use grape of not 
#define USE_GRAPE HAVE_GRAPE
#endif

#if USE_GRAPE
#include <dune/grid/io/visual/grapedatadisplay.hh>
#endif

namespace Dune {
template <class GridImp, 
          class DataImp> 
class DataOutput;

/** \brief Parameter class for Dune::DataWriter

    Structure providing the main parameters used to setup the Dune::DataWriter.
    By default these parameters are set through the Dune::Parameter class, i.e.,
    can be set in a parameter file. To override this a user can derive from this
    structure and modify any subset of the parameters. An instance of this modified
    class can then be passed in the construction of the Dune::DataWriter.
 */
struct DataOutputParameters : public LocalParameter<DataOutputParameters,DataOutputParameters> {
  //! path where the data is stored (path are always relative to fem.commonOutputPath)
  virtual std::string path() const {
    return Parameter::getValue<std::string>("fem.io.path","");
  }
  //! base of file name for data file (fem.io.datafileprefix)
  virtual std::string prefix() const {
    return Parameter::getValue<std::string>("fem.io.datafileprefix");
  }
  //! format of output (fem.io.outputformat)
  virtual int outputformat() const {
    return Parameter::getValidValue<int>("fem.io.outputformat",2,
           ValidateInterval<int,true,true>(1,3));
  }
  //! use online grape display (fem.io.grapedisplay)
  virtual bool grapedisplay() const {
    return (Parameter::getValue("fem.io.grapedisplay",0) == 1);
  }
  //! save data every savestep interval (fem.io.savestep)
  virtual double savestep() const {
    return Parameter::getValue<double>("fem.io.savestep",0);
  }
  //! save data every savecount calls to write method (fem.io.savecount)
  virtual int savecount() const {
    return Parameter::getValue<int>("fem.io.savecount",0);
  }
  //! number for first data file (no parameter available)
  virtual int startcounter() const {
    return 0;
  }
  //! method used for conditional data output - default
  //! value passed as argument.
  virtual bool willWrite(bool write) const {
    return write;
  }
};

/** @ingroup DiscFuncIO 
   \brief Implementation of the Dune::IOInterface. 
   This class manages data output.
   Available output formats are GRAPE, VTK and VTK Vertex projected
   using the VtxProjection operator. Details can be
   found in \ref DiscFuncIO.
*/    
template <class GridImp, 
          class DataImp> 
class DataOutput : public IOInterface 
{

protected:  
#if USE_VTKWRITER
  template <class VTKIOType>
  class VTKListEntry 
  {
  protected: 
    VTKListEntry () {}
  public:
    virtual ~VTKListEntry () {}
    virtual void add(VTKIOType &) const = 0;
  };

  template <class VTKIOType, class DFType> 
  class VTKFunc : public VTKListEntry<VTKIOType>
  {
    typedef typename VTKIOType :: GridPartType GridPartType ;
    typedef typename DFType :: DiscreteFunctionSpaceType :: FunctionSpaceType FunctionSpaceType;

    typedef LagrangeDiscreteFunctionSpace< FunctionSpaceType,
            GridPartType, 1> LagrangeSpaceType;
    typedef AdaptiveDiscreteFunction< LagrangeSpaceType >  NewFunctionType;

    const DFType& df_;
    LagrangeSpaceType space_;
    mutable NewFunctionType* func_;
    
    VTKFunc( const VTKFunc&) ; 
  public:
    VTKFunc( const GridPartType& gridPart, const DFType& df ) 
      : df_ (df) 
      , space_( const_cast<GridPartType&> (gridPart) )
      , func_(0)
    {
      // space_.setDescription(df.space().getDescription());
    }

    ~VTKFunc () 
    {
      delete func_;
    }

    virtual void add(VTKIOType& vtkio) const 
    {
      func_ = new NewFunctionType ( df_.name(), space_ );
      WeightDefault<GridPartType> weight;
      VtxProjectionImpl::project( df_, *func_, weight );
      vtkio.addVertexData( *func_ );
      // vtkio.addVectorVertexData( *func_ );
    }
  };

  template <class Entity>
  class GnuplotOutputer {
    std::ostream& out_;
    const Entity& en_;
  public:
    //! Constructor
    GnuplotOutputer(std::ostream& out,
                    const Entity& en) : 
    out_(out),
    en_(en)
    {
    }

    //! Applies the setting on every DiscreteFunction/LocalFunction pair.
    template <class DFType>
    void visit(DFType* df) 
    {
      typedef typename DFType :: Traits Traits;
      typedef typename Traits :: LocalFunctionType LocalFunctionType;
      typedef typename Traits :: DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
      typedef typename DiscreteFunctionSpaceType :: IteratorType IteratorType;
      typedef typename DiscreteFunctionSpaceType :: GridPartType GridPartType;

      typedef typename DiscreteFunctionSpaceType :: DomainType DomainType;
      typedef typename DiscreteFunctionSpaceType :: RangeType RangeType;
 
      enum{ dimDomain = DiscreteFunctionSpaceType :: dimDomain };
      enum{ dimRange = DiscreteFunctionSpaceType :: dimRange };

      CachingQuadrature<GridPartType,0> quad(en_,df->space().order());
      LocalFunctionType lf = df->localFunction(en_);
      for (int i=0;i<quad.nop();++i) {
        RangeType u;
        DomainType x = en_.geometry().global(quad.point(i));
        lf.evaluate(quad[i],u);
        for (int i = 0; i < dimDomain; ++i) 
          out_ << x[i] << " ";
        for (int i = 0; i < dimRange; ++i) 
          out_ << u[i] << "   ";
        out_ << std::endl;
      }
    }
  };

  template <class VTKOut>
  class VTKOutputerDG {
  public:
    //! Constructor
    VTKOutputerDG(VTKOut& vtkOut) :
      vtkOut_(vtkOut)
    {
    }

    //! Applies the setting on every DiscreteFunction/LocalFunction pair.
    template <class DFType>
    void visit(DFType* df) 
    {
      if(df->space().order() > 0)
      {
        vtkOut_.addVertexData( *df );
        // vtkOut_.addVectorVertexData( *df );
      }
      else 
      {
        vtkOut_.addCellData( *df ); 
        // vtkOut_.addVectorCellData( *df );
      }
    }

  private:
    VTKOut & vtkOut_;
  };

  template <class VTKOut>
  class VTKOutputerLagrange {
  public:
    //! Constructor
    VTKOutputerLagrange( VTKOut& vtkOut ) :
      vtkOut_(vtkOut), vec_()
    {
    }

    //! destructor, delete entries of list 
    ~VTKOutputerLagrange() 
    {
      for(size_t i=0; i<vec_.size(); ++i)
      {
        delete vec_[i];
        vec_[i] = 0;
      }
    }

    //! Applies the setting on every DiscreteFunction/LocalFunction pair.
    template <class DFType>
    void visit(DFType* df) 
    {
      typedef VTKFunc<VTKOut,DFType> EntryType; 
      EntryType* entry = new EntryType(vtkOut_.gridPart(), *df);
      entry->add( vtkOut_ );
      vec_.push_back( entry );
    }

  private:
    VTKOut & vtkOut_;
    typedef VTKListEntry<VTKOut> VTKListEntryType;
    std::vector<VTKListEntryType *> vec_;
  };
#endif

protected:  
  enum OutputFormat { vtk = 1 , vtkvtx = 2 , gnuplot = 3 };

  //! \brief type of grid used 
  typedef GridImp GridType;
  //! \brief type of data tuple 
  typedef DataImp OutPutDataType; 

  //! \brief type of this class 
  typedef DataOutput< GridImp, DataImp > ThisType;
  const GridType& grid_;
  mutable OutPutDataType& data_; 
  const int myRank_;
  bool verbose_;

  // name for data output 
  std::string path_;
  std::string datapref_;
  // use grape display
  bool grapeDisplay_;
  
  // counter for sequence of files
  mutable int writeStep_;
  // counter for call to write
  mutable int writeCalls_;
  // next point in time to save data
  mutable double saveTime_;
  // time interval to save files
  double saveStep_;
  // number of write call between writting files
  int saveCount_;
  // grape, vtk or ...
  OutputFormat outputFormat_;
  mutable std::ofstream sequence_;
  const DataOutputParameters* param_;
public: 
 /** \brief Constructor creating data output class 
    \param grid corresponding grid 
    \param data Tuple containing discrete functions to write 
    \param parameter structure for tuning the behavior of the Dune::DataOutput 
                     defaults to Dune::DataOutputParameters
  */
  DataOutput(const GridType & grid, 
             OutPutDataType& data,
             const DataOutputParameters& parameter=DataOutputParameters())
    : grid_( grid ),
      data_( data ), 
      myRank_(grid_.comm().rank()),
      writeStep_(0),
      writeCalls_(0),
      saveTime_(0.0),  // why 0.0?
      saveStep_(-1),
      saveCount_(-1),
      outputFormat_(vtkvtx),
      param_(parameter.clone())
  {
    // initialize class 
    init(parameter);
  }

  /** \brief Constructor creating data writer 
    \param grid corresponding grid 
    \param data Tuple containing discrete functions to write 
    \param tp   a time provider to set time
    \param parameter structure for tuning the behavior of the Dune::DataOutput 
                     defaults to Dune::DataOutputParameters
  */ 
  DataOutput(const GridType & grid, 
             OutPutDataType& data, 
             const TimeProviderBase& tp,
             const DataOutputParameters& parameter=DataOutputParameters())
    : grid_(grid),
      data_(data),
      myRank_(grid_.comm().rank()),
      writeStep_(0),
      writeCalls_(0),
      saveTime_(0),
      saveStep_(-1),
      saveCount_(-1),
      outputFormat_(vtkvtx),
      param_(parameter.clone())
  {
    init(parameter);

    // set old values according to new time 
    if (saveStep_>0) {
      while(saveTime_ < tp.time()) 
      {
        ++writeStep_;
        saveTime_ += saveStep_;
      }
    }
  }

  //! destructor 
  virtual ~DataOutput() {}

protected:  
  //! initialize data writer 
  void init(const DataOutputParameters& parameter) 
  {
    // read verbose parameter 
    verbose_ = Parameter::verbose(); 

    IOInterface :: createGlobalPath ( grid_.comm(), Parameter::commonOutputPath() );
    path_ = Parameter::commonOutputPath()+"/"+parameter.path();
    // create path if not already exists 
    IOInterface :: createGlobalPath ( grid_.comm(), path_ );

    // write parameter file 
    
    // add prefix for data file
    datapref_ += parameter.prefix();

    int outputFormat = parameter.outputformat();
    switch( outputFormat ) 
    {
      case 1: outputFormat_ = vtk; break;
      case 2: outputFormat_ = vtkvtx; break;
      case 3: outputFormat_ = gnuplot; break;
      default:
        DUNE_THROW(NotImplemented,"DataWriter::init: wrong output format");
    }

    grapeDisplay_ = parameter.grapedisplay();

    // get parameters for data writing
    saveStep_ = parameter.savestep();
    saveCount_ = parameter.savecount();

    writeStep_ = parameter.startcounter();

    if (myRank_==0) {
      std::string name = path_ + "/" + datapref_;
      name += ".series";
      std::cout << "opening file: " << name << std::endl;
      sequence_.open(name.c_str());
      if (!sequence_)
        std::cout << "could not write sequence file" << std::endl;
    }
    Parameter::write("parameter.log");
  }
private:
  bool willWrite(const TimeProviderBase& tp) const 
  {
    return ( (saveStep_>0 && tp.time() >= saveTime_ ) || 
             // tp.end() ||
           (saveCount_>0 && writeCalls_%saveCount_ == 0) );
  }
  bool willWrite() const 
  {
    return (saveCount_>0 && writeCalls_%saveCount_ == 0) ;
  }
public:
  /** \brief write given data to disc
  */
  void write() const 
  {
    if (param_->willWrite(willWrite())) {
      if (sequence_) 
        sequence_ << writeStep_ << " " << writeStep_ << std::endl;
      writeData();
    }
  }
  void write(const TimeProviderBase& tp) const
  {
    if (param_->willWrite(willWrite(tp))) {
      if (sequence_)
        sequence_ << writeStep_ << " " << tp.time() << std::endl;
      writeData();
    }
  }
  void writeData() const
  {
    // check online display 
    display(); 
    switch (outputFormat_)
    {
#if USE_VTKWRITER
    case vtk : 
    case vtkvtx :
      // write data in vtk output format 
      writeVTKOutput();
      break;
#endif
    case gnuplot :
      writeGnuPlotOutput();
      break;
    otherwise:
      DUNE_THROW(NotImplemented,"DataWriter::write: wrong output format");
    }
    // only write info for proc 0, otherwise on large number of procs
    // this is to much output 
    if(myRank_ <= 0)
    {
      std::cout << "DataOutput["<<myRank_<<"]::write data with step number " << writeStep_ << std::endl;
    }
    saveTime_ += saveStep_;
    ++writeCalls_;
    ++writeStep_;
  }

  /** Return output path name */
  const std::string& path() const {
    return path_;
  }

protected:  
#if USE_VTKWRITER
  void writeVTKOutput() const 
  {
    // check whether to use vertex data of discontinuous data 
    const bool vertexData = (outputFormat_ == vtkvtx);

    // check whether we have parallel run  
    const bool parallel = (grid_.comm().size() > 1);

    // generate filename, with path only for serial run  
    std::string name = genFilename( (parallel) ? "" : path_, datapref_, writeStep_ );

#ifdef YASPGRID 
    if( vertexData )
    {
      static bool called = false; 
      if( ! called ) 
      {
        std::cerr << "WARNING: vertexData output not working with YaspGrid! \n"; 
        called = true;
      }
    }

#else 
    if( vertexData ) 
    {
      // generate adaptive leaf grid part 
      // do not use leaf grid part since this will 
      // create the grids leaf index set, which might not be wanted. 
      typedef AdaptiveLeafGridPart< GridType > GridPartType; 
      GridPartType gridPart( const_cast<GridType&> (grid_) );

      // create vtk output handler 
      typedef VTKIO < GridPartType > VTKIOType; 
      VTKIOType vtkio ( gridPart, VTKOptions::conforming );

      // add all functions 
      ForEachValue<OutPutDataType> forEach(data_); 
      VTKOutputerLagrange< VTKIOType > io( vtkio );
      forEach.apply( io );

      if( parallel )
      {
        // write all data for parallel runs  
        vtkio.pwrite( name.c_str(), path_.c_str(), "." , Dune::VTKOptions::binaryappended );
      }
      else
      {
        // write all data serial 
        vtkio.write( name.c_str(), Dune::VTKOptions::binaryappended );
      }
    }
    else
#endif 
    {
      // generate adaptive leaf grid part 
      // do not use leaf grid part since this will 
      // create the grids leaf index set, which might not be wanted. 
      typedef AdaptiveLeafGridPart< GridType > GridPartType; 
      GridPartType gridPart( const_cast<GridType&> (grid_) );
      {
        // create vtk output handler 
        typedef VTKIO < GridPartType > VTKIOType; 
        VTKIOType vtkio ( gridPart, VTKOptions::nonconforming );

        // add all functions 
        ForEachValue<OutPutDataType> forEach(data_); 
        VTKOutputerDG< VTKIOType > io( vtkio );
        forEach.apply( io );

        // write all data 
        if( parallel )
        {
          // write all data for parallel runs  
          vtkio.pwrite( name.c_str(), path_.c_str(), "." , Dune::VTKOptions::binaryappended );
        }
        else
        {
          // write all data serial 
          vtkio.write( name.c_str(), Dune::VTKOptions::binaryappended );
        }
      }
    }
  }
#endif

  // write to gnuplot file format
  void writeGnuPlotOutput() const
  {
    // generate filename
    std::string name = genFilename( path_, datapref_, writeStep_ );
    name += ".gnu";
    std::ofstream gnuout(name.c_str());

    // generate adaptive leaf grid part 
    // do not use leaf grid part since this will 
    // create the grids leaf index set, which might not be wanted. 
    typedef DGAdaptiveLeafGridPart< GridType > GridPartType; 
    GridPartType gridPart( const_cast<GridType&> (grid_) );
    // start iteration
    typedef typename GridPartType::template Codim<0>::IteratorType IteratorType;
    IteratorType endit = gridPart.template end<0>();
    for (IteratorType it = gridPart.template begin<0>(); it != endit; ++it) {
      ForEachValue<OutPutDataType> forEach(data_); 
      GnuplotOutputer< typename GridPartType::EntityCodim0Type > 
                     io( gnuout,*it );
      forEach.apply( io );
    }
  }

  //! display data with grape 
  virtual void display() const 
  {
    if( grapeDisplay_ )
    {
#if USE_GRAPE
      GrapeDataDisplay<GridType> grape(grid_);
      IOTuple<OutPutDataType>::addToDisplay(grape,data_);
      grape.display();
#else 
      std::cerr <<"WARNING: HAVE_GRAPE == 0, but grapeDisplay == true, recompile with grape support for online display!" << std::endl;
#endif
    }
  }

}; // end class DataOutput

} // end namespace DataIO 
#endif
