#ifndef DUNE_DATAWRITER_HH
#define DUNE_DATAWRITER_HH

#define USE_VTKWRITER

//- Dune includes 
#include <dune/fem/io/file/iointerface.hh>
#include <dune/fem/io/file/iotuple.hh>

#include <dune/fem/space/common/adaptiveleafgridpart.hh>
#include <dune/fem/space/lagrangespace.hh>
#include <dune/fem/function/adaptivefunction.hh>

#ifdef USE_VTKWRITER
#include <dune/fem/io/file/vtkio.hh>
#include <dune/fem/operator/projection/vtxprojection.hh>
#endif

// define whether to use grape of not 
#define USE_GRAPE HAVE_GRAPE

#if USE_GRAPE
#include <dune/grid/io/visual/grapedatadisplay.hh>
#endif

namespace Dune {

/** @ingroup Checkpointing 
   \brief Implementation of the IOInterface. 
   This class is managing our data output. 
   Available output formats are GRAPE, VTK and VTK Vertex projected
   using the VtxProjection operator. 
*/    
template <class GridImp, 
          class DataImp> 
class DataWriter : public IOInterface 
{

protected:  
#ifdef USE_VTKWRITER
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
      }
      else 
      {
        vtkOut_.addCellData( *df ); 
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
  enum OutputFormat { grape = 0 , vtk = 1 , vtkvtx = 2 };

  typedef GridImp GridType;
  typedef DataImp OutPutDataType; 
  
  const GridType& grid_;
  mutable OutPutDataType& data_; 

  // name for data output 
  std::string path_;
  std::string datapref_;

  bool grapeDisplay_;
  
  mutable int writeStep_;
  mutable double saveTime_;

  double saveStep_;
  double endTime_;
  const int myRank_;

  bool verbose_;

  // grape of vtk 
  OutputFormat outputFormat_;

public: 
  /** \brief Constructor creating data writer 
    \param grid corresponding grid 
    \param data Tuple containing discrete functions to write 
    \param paramfile parameter file containing parameters for data writer
    \note Possible values are (copy this to your parameter file):

    # verbosity (0 = no, 1 = yes)
    verbose: 0
    
    # OutputPath 
    OutputPath: ./
    
    # name for data set  
    OutputPrefix: solution
    
    # format of output: 0 = GRAPE, 1 = VTK, 2 = VTK vertex data
    OutputFormat: 0
    
    # GrapeDisplay (0 = no, 1 = yes)
    GrapeDisplay: 0 
    
    # SaveStep (write data every `saveStep' time period 
    SaveStep: 0.01

    # StartTime of simulation
    StartTime: 0.0

    # EndTime of simulation
    EndTime: 1.0
  */
  DataWriter(const GridType & grid, OutPutDataType& data, 
             const std::string paramfile) DUNE_DEPRECATED 
    : grid_(grid), data_(data) 
    , writeStep_(0)
    , saveTime_(0.0)
    , saveStep_(0.1)
    , endTime_(1.0)
    , myRank_(grid_.comm().rank())
    , outputFormat_(grape)
  {
    init(paramfile);
  }

  /** \brief Constructor creating data writer 
    \param grid corresponding grid 
    \param gridName corresponding macro grid name (needed for structured grids)
    \param data Tuple containing discrete functions to write 
    \param paramfile parameter file containing parameters for data writer
    \note Possible values are (copy this to your parameter file):

    # verbosity (0 = no, 1 = yes)
    verbose: 0
    
    # OutputPath 
    OutputPath: ./
    
    # name for data set  
    OutputPrefix: solution
    
    # format of output: 0 = GRAPE, 1 = VTK, 2 = VTK vertex data
    OutputFormat: 0
    
    # GrapeDisplay (0 = no, 1 = yes)
    GrapeDisplay: 0 
    
    # SaveStep (write data every `saveStep' time period 
    SaveStep: 0.01

    # StartTime of simulation
    StartTime: 0.0

    # EndTime of simulation
    EndTime: 1.0
  */
  DataWriter(const GridType & grid, 
             const std::string& gridName, 
             OutPutDataType& data, 
             const std::string paramfile)
    : grid_(grid), data_(data) 
    , writeStep_(0)
    , saveTime_(0.0)
    , saveStep_(0.1)
    , endTime_(1.0)
    , myRank_(grid_.comm().rank())
    , outputFormat_(grape)
  {
    saveMacroGrid( gridName );
    init(paramfile);
  }

  /** \brief Constructor creating data writer 
    \param grid corresponding grid 
    \param data Tuple containing discrete functions to write 
    \param paramfile parameter file containing parameters for data writer
    \param time for restarted jobs get time from outside 
  */ 
  DataWriter(const GridType & grid, OutPutDataType& data, 
             const std::string paramfile,
             double time)
    : grid_(grid), data_(data) 
    , writeStep_(0) 
    , saveTime_(0.0)
    , saveStep_(0.1)
    , endTime_(1.0)
    , myRank_(grid_.comm().rank())
    , outputFormat_(grape)
  {
    init(paramfile);

    // set old values according to new time 
    while(saveTime_ < time) 
    {
      ++writeStep_;
      saveTime_ += saveStep_;
    }
  }

  //! destructor 
  virtual ~DataWriter() {}

protected:  
  void init(const std::string& paramfile) 
  {
    const bool verboseOutput = (myRank_ == 0);
    // read verbose parameter 
    int verb = 0;
    readParameter(paramfile,"verbose",verb, verboseOutput );
    verbose_ = (verb == 1) ? true : false;

    if( ! readParameter(paramfile,"OutputPath",path_, verboseOutput) )
    {
      path_ = ".";
    }
    else
    {
      // create path if not already exists 
      IOInterface :: createGlobalPath ( grid_.comm(), path_ );
    }

    // copy parameter file 
    {
      std::string cmd("cp ");
      cmd += paramfile;
      cmd += " ";
      cmd += path_;
      system( cmd.c_str() );
    }

    {
      std::string dummyfile;
      readParameter(paramfile,"OutputPrefix",dummyfile, verbose_ );
      datapref_ += dummyfile;
    }

    int outputFormat = 0; 
    readParameter(paramfile,"OutputFormat",outputFormat, verboseOutput);
    switch( outputFormat ) 
    {
      case 0: outputFormat_ = grape; break;
      case 1: outputFormat_ = vtk; break;
      case 2: outputFormat_ = vtkvtx; break;
      default:
        DUNE_THROW(NotImplemented,"DataWriter::init: wrong output format");
    }

    int gpdisp = 0;
    readParameter(paramfile,"GrapeDisplay",gpdisp, verbose_ );
    grapeDisplay_ = (gpdisp == 1) ? true : false;

    // get parameters for data writing
    readParameter(paramfile,"SaveStep",saveStep_, verbose_ );

    // read start time 
    readParameter(paramfile,"StartTime",saveTime_, verbose_ );
    
    // read end time 
    readParameter(paramfile,"EndTime",endTime_, verbose_ );
  }

public:
  /** \copydoc IOInterface::write */
  virtual void write(double time, int timestep) const 
  {
    // only write data time > saveTime  
    if( (time >= saveTime_ ) || (time >= endTime_) )
    {
      // check online display 
      display(); 

      if( outputFormat_ == grape )
      {
        // create new path for time step output 
        std::string timeStepPath = IOInterface::createPath ( grid_.comm(),
            path_, datapref_ , writeStep_ );

        // for structured grids copy grid 
        IOInterface::copyMacroGrid(grid_,path_,timeStepPath,datapref_);

        GrapeDataIO<GridType> dataio;
        // call output of IOTuple 
        IOTuple<OutPutDataType>::output(dataio, 
          grid_ ,time, writeStep_, timeStepPath , datapref_, data_ , verbose_ );
      }
#ifdef USE_VTKWRITER
      else if ( outputFormat_ == vtk || outputFormat_ == vtkvtx )
      {
        // write data in vtk output format 
        writeVTKOutput( Element<0>::get(data_), time );
      }
#endif
      else 
      {
        DUNE_THROW(NotImplemented,"DataWriter::write: wrong output format");
      }

      // only write info for proc 0, otherwise on large number of procs
      // this is to much output 
      if(myRank_ <= 0)
      {
        std::cout << "DataWriter["<<myRank_<<"]::write:  time = "<< time << "  writeData: step number " << writeStep_ << std::endl;
      }
      saveTime_ += saveStep_;
      ++writeStep_;
    }
    return;
  }

protected:
#ifdef USE_VTKWRITER
  template <class DFType> 
  void writeVTKOutput(const DFType* func, double time) const 
  {
    // check whether to use vertex data of discontinuous data 
    const bool vertexData = (outputFormat_ == vtkvtx);

    // generate filename 
    std::string name = genFilename( path_, datapref_, writeStep_ );

#ifdef YASPGRID 
    if( vertexData )
    {
      static bool called = false; 
      if( ! called ) 
      {
        std::cerr << "WARNING: vertexData output no tworking with YaspGrid! \n"; 
        called = true;
      }
    }

#else 
    if( vertexData ) 
    {
      // generate leaf grid part 
      typedef AdaptiveLeafGridPart< GridType > GridPartType; 
      GridPartType gridPart( const_cast<GridType&> (grid_) );

      // create vtk output handler 
      typedef VTKIO < GridPartType > VTKIOType; 
      VTKIOType vtkio ( gridPart, VTKOptions::conforming );

      // add all functions 
      ForEachValue<OutPutDataType> forEach(data_); 
      VTKOutputerLagrange< VTKIOType > io( vtkio );
      forEach.apply( io );

      vtkio.write( name.c_str(), Dune::VTKOptions::ascii );
    }
    else
#endif 
    {
      // get grid part 
      typedef typename DFType :: DiscreteFunctionSpaceType :: GridPartType GridPartType; 
      const GridPartType& gridPart = func->space().gridPart();

      // create vtk output handler 
      typedef VTKIO < GridPartType > VTKIOType; 
      VTKIOType vtkio ( gridPart, VTKOptions::nonconforming );

      // add all functions 
      ForEachValue<OutPutDataType> forEach(data_); 
      VTKOutputerDG< VTKIOType > io( vtkio );
      forEach.apply( io );

      // write all data 
      vtkio.write( name.c_str(), Dune::VTKOptions::ascii );
    }
  }
#endif

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

public:  
  /** \brief save structured macro grid file 
    \param macroFileName filename which contains the macro file 
  */
  virtual void saveMacroGrid(const std::string macroFileName) const 
  {
    IOInterface :: writeMacroGrid( grid_, macroFileName, path_, datapref_);
  }
  
}; // end class DataWriter 
  

/** @ingroup Checkpointing 
   \brief Implementation of the IOInterface. 
   This class is managing our checkpointing. 
   The data will be stored in GRAPE output format, meaning that every
   checkpoint also is a visualizable data set. 
*/    
template <class GridImp, 
          class DataImp> 
class CheckPointer : public DataWriter<GridImp,DataImp> 
{
  typedef DataWriter<GridImp,DataImp> BaseType;
  
  typedef GridImp GridType;
  typedef DataImp OutPutDataType; 
  

  int checkPointStep_;
  mutable int checkPointNumber_;

  std::string runfile_;
  std::string checkPointFile_;

  // pointer to adaptation manager 
  const LoadBalancerInterface* lb_;
  
  int balanceRecover_;

  static const char * checkPointPrefix()
  {
    return "checkpoint";
  }
  
public: 
  /** \brief Constructor generating a cechkpointer 
    \param grid corresponding grid 
    \param gridName name of macro grid file (for structured grids)
    \param data Tuple containing discrete functions to write 
    \param paramfile parameter file containing parameters for data writer
    \param checkFile filename for restoring state of program from
           previous runs (default is zero, which means start from new)
    \param lb LoadBalancer instance 
      (default is zero, which means start  from new)

    \note In Addition to the parameters read by DataWriter this class 
          reads the following parameters: 

    # write checkpoint every `CheckPointStep' time step
    CheckPointStep: 500 

    # store checkpoint information to file `CheckPointFile'
    CheckPointFile: checkpoint
  */
  CheckPointer(const GridType & grid, 
               const std::string& gridName, 
               OutPutDataType& data, 
               const std::string paramfile,
               const char * checkFile = 0,
               const LoadBalancerInterface* lb = 0) 
    : BaseType(grid,gridName,data,paramfile) 
    , checkPointStep_(500)
    , checkPointNumber_(1)
    , lb_(lb)
    , balanceRecover_(0)
  {
    this->datapref_ = checkPointPrefix();
    readParameter(paramfile,"CheckPointStep",checkPointStep_, this->verbose_ );

    if( checkFile )
    {
      checkPointFile_ = checkFile;

      // read last counter 
      bool ok = readCheckPoint();

      // read name of check point file 
      std::string dummyfile;
      readParameter(paramfile,"CheckPointFile",dummyfile, this->verbose_ );
      checkPointFile_ = this->path_; checkPointFile_ += "/"; checkPointFile_ += dummyfile;

      // if check point couldn't be opened, try again   
      if(!ok)
      {
        ok = readCheckPoint();
        if( ! ok )
        {
          std::cerr <<"ERROR: unable to open checkpoint file! \n";
          exit(1);
        }
      }

      // copy parameter file 
      std::string cmd("cp ");
      cmd += paramfile; cmd += " "; cmd += this->path_;

      // copy parameter file to outpath 
      system (cmd.c_str());
    }
    else
    {
      // read name of check point file 
      std::string dummyfile;
      readParameter(paramfile,"CheckPointFile",dummyfile, this->verbose_ );
      checkPointFile_ = this->path_; checkPointFile_ += "/"; 
      checkPointFile_ += dummyfile;

      // read last counter 
      readCheckPoint();
    }
  }

  /** \brief Constructor generating a cechkpointer 
    \param grid corresponding grid 
    \param data Tuple containing discrete functions to write 
    \param paramfile parameter file containing parameters for data writer
    \param checkFile filename for restoring state of program from
           previous runs (default is zero, which means start from new)
    \param lb LoadBalancer instance 
      (default is zero, which means start  from new)

    \note In Addition to the parameters read by DataWriter this class 
          reads the following parameters: 

    # write checkpoint every `CheckPointStep' time step
    CheckPointStep: 500 

    # store checkpoint information to file `CheckPointFile'
    CheckPointFile: checkpoint
  */
  CheckPointer(const GridType & grid, OutPutDataType& data, 
               const std::string paramfile,
               const char * checkFile = 0,
               const LoadBalancerInterface* lb = 0) DUNE_DEPRECATED
    : BaseType(grid,data,paramfile) 
    , checkPointStep_(500)
    , checkPointNumber_(1)
    , lb_(lb)
    , balanceRecover_(0)
  {
    this->datapref_ = checkPointPrefix();
    readParameter(paramfile,"CheckPointStep",checkPointStep_, this->verbose_ );

    if( checkFile )
    {
      checkPointFile_ = checkFile;

      // read last counter 
      bool ok = readCheckPoint();

      // read name of check point file 
      std::string dummyfile;
      readParameter(paramfile,"CheckPointFile",dummyfile, this->verbose_ );
      checkPointFile_ = this->path_; checkPointFile_ += "/"; checkPointFile_ += dummyfile;

      // if check point couldn't be opened, try again   
      if(!ok)
      {
        ok = readCheckPoint();
        if( ! ok )
        {
          std::cerr <<"ERROR: unable to open checkpoint file! \n";
          exit(1);
        }
      }

      // copy parameter file 
      std::string cmd("cp ");
      cmd += paramfile; cmd += " "; cmd += this->path_;

      // copy parameter file to outpath 
      system (cmd.c_str());
    }
    else
    {
      // read name of check point file 
      std::string dummyfile;
      readParameter(paramfile,"CheckPointFile",dummyfile, this->verbose_ );
      checkPointFile_ = this->path_; checkPointFile_ += "/"; 
      checkPointFile_ += dummyfile;

      // read last counter 
      readCheckPoint();
    }
  }

  /** \brief restore grid from previous runs 
    \param[in] paramfile parameter filename  
    \param[in] checkFile checkPoint filename 
    \param[in] rank number of my process 
    \param tp TimeProvider to set time and timestep to

    \return Pointer to restored grid 
  */
  template <class TimeProviderImp> 
  static GridType* restoreGrid(const std::string paramfile, 
                               const std::string checkFile,
                               const int rank, 
                               TimeProviderImp& tp)
  {
    GrapeDataIO<GridType> dataio;

    std::string datapref( checkPointPrefix() );
    std::string path;
    std::string checkPointFile;

    int checkPointNumber = 0;
    // if given checkpointfile is not valid use default checkpoint file 
    if( ! readParameter(checkFile,"LastCheckPoint",checkPointNumber) )
    {
      // read default path
      path = IOInterface::readPath(paramfile);
      // set checkpointfile 
      checkPointFile = path;

      // read name of check point file 
      std::string dummyfile;
      readParameter(paramfile,"CheckPointFile",dummyfile);
      checkPointFile += "/"; checkPointFile += dummyfile;
      readParameter(checkPointFile,"LastCheckPoint",checkPointNumber);
    }
    else
    {
      if( ! readParameter(checkFile,"RecoverPath",path) )
      {
        // read default path
        path = IOInterface::readPath(paramfile);
      }
    }

    int timeStep = 0;
    // try to read time step, if fails default value is taken
    readParameter(checkPointFile,"TimeStep",timeStep);

    // now add timestamp and rank 
    path = IOInterface::createRecoverPath(
        path, rank, datapref, checkPointNumber );

    // time is set during grid restore  
    double time = 0.0; 

    GridType* grid = 
      IOTuple<OutPutDataType>::restoreGrid(dataio, time, 
          checkPointNumber , 
          path, datapref);

    tp.setTime( time, timeStep );
    
    return grid;
  }

  /** \brief restores data, assumes that all objects have been created before
  */
  void restoreData()
  {
    // now add timestamp and rank 
    std::string path = IOInterface::createRecoverPath(
        this->path_, this->myRank_ , this->datapref_, checkPointNumber_ );

    GrapeDataIO<GridType> dataio;
    IOTuple<OutPutDataType>::restoreData(this->data_, dataio, 
      this->grid_, checkPointNumber_, path , this->datapref_);

    typedef DofManager<GridType> DMType;
    typedef DofManagerFactory<DMType> DMFactoryType;

    DMType & dm = DMFactoryType::getDofManager(this->grid_);

    // do compress of data 
    dm.dofCompress();
  }

  /** \copydoc IOInterface::write */
  virtual void write(double time, int timestep) const 
  {
    // only write data time > saveTime  
    if( ((timestep % checkPointStep_) == 0) && timestep > 0 )
    {
      // toggle between 0 and 1 
      checkPointNumber_ = (checkPointNumber_ == 0) ? 1 : 0;

      // create new timestep path 
      std::string path = IOInterface::createPath ( this->grid_.comm(),
          this->path_, this->datapref_ , checkPointNumber_ );

      GrapeDataIO<GridType> dataio;
      // call output of IOTuple 
      IOTuple<OutPutDataType>::output(dataio, this->grid_ ,
          time, checkPointNumber_ , path, this->datapref_ , this->data_ , false );
      
      writeCheckPoint(path,time,timestep,
                      checkPointNumber_, 
                      (lb_) ? lb_->balanceCounter() : 0);
    }
    return;
  }

  /** returns actual balance counter, for restoreing LoadBalancer 
    \return Restored load balance counter 
   */
  
  int balanceCounter () const 
  {
    return balanceRecover_;
  }
  
private:
  //! read checkpoint file
  bool readCheckPoint()
  {
    // read Checkpiont file 
    if( readParameter(checkPointFile_,"LastCheckPoint",checkPointNumber_, this->verbose_ ) )
    {
      readParameter(checkPointFile_,"BalanceCounter",balanceRecover_, this->verbose_ );
      return true;
    }
    return false;
  }
  
  // write some info for checkpointing 
  void writeCheckPoint (const std::string& path,
            const double time, int timestep,
            int savestep, int balancecounter) const 
  {   
    // write some needed informantion to current checkpoint file 
    {
      std::string filepref(path);
      filepref += "/";
      filepref += this->datapref_;

      std::string filename = genFilename("",filepref,savestep);

      std::ofstream file (filename.c_str());
      if( file.is_open() )
      {
        file << "Time: "      << std::scientific << time << std::endl;
        file << "SaveCount: " << savestep << std::endl;
      }
      else
      {
        std::cerr << "Couldn't open file `" << filename << "' ! " << std::endl;
      }
    }
    
    // only proc 0 writes global checkpoint file 
    if( this->myRank_ <= 0)
    {
      // write last checkpoint to filename named like the checkpoint files 
      // but with no extentions 
      std::ofstream file (checkPointFile_.c_str());
      if( file.is_open() )
      {
        file << "LastCheckPoint: " << savestep << std::endl;
        file << "Time: " << std::scientific << time << std::endl;
        file << "TimeStep: "  << timestep << std::endl;
        file << "BalanceCounter: " << balancecounter << std::endl;
        file << "# RecoverPath can be edited by hand if data has been moved!" << std::endl;
        file << "RecoverPath: " << this->path_ << std::endl;
        file.close();
      }
      else
      {
        std::cerr << "Couldn't open file `" << checkPointFile_ << "' ! " << std::endl;
      }

      std::cout << "CheckPointer: time = " << time << " , wrote checkpoint path `" << path << "'\n";
    }
  }

}; // end class CheckPointer 
  
} // end namespace DataIO 
#endif
