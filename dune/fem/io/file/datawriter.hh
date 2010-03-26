#ifndef DUNE_DATAWRITER_HH
#define DUNE_DATAWRITER_HH

#include <string>

#include <dune/fem/io/file/iointerface.hh>
#include <dune/fem/io/file/iotuple.hh>
#include <dune/fem/io/parameter.hh>
#include <dune/fem/space/common/loadbalancer.hh>
#include <dune/fem/gridpart/adaptiveleafgridpart.hh>
#include <dune/fem/space/lagrangespace.hh>
#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/quadrature/cachingquadrature.hh>
#include <dune/fem/io/file/persistencemanager.hh>

#ifndef USE_GRAPE 
// define whether to use grape of not 
#define USE_GRAPE HAVE_GRAPE
#endif

#include <dune/grid/io/visual/grapedatadisplay.hh>
#include <dune/fem/io/file/dataoutput.hh>

namespace Dune {

struct DataWriterParameters : public DataOutputParameters
{
  //! base of file name for data file (fem.io.macroGridFile)
  virtual std::string macroGridName () const
  {
    return Parameter::getValue< std::string >( "fem.io.macroGridFile" );
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
class DataWriter : public DataOutput< GridImp, DataImp > 
{
protected:  
  //! \brief type of grid used 
  typedef GridImp GridType;
  //! \brief type of data tuple 
  typedef DataImp OutPutDataType; 

  //! \brief type of this class 
  typedef DataWriter< GridImp, DataImp > ThisType;

  typedef DataOutput< GridImp, DataImp > BaseType;

  using BaseType :: grid_;
  using BaseType :: data_;

  using BaseType :: path_;
  using BaseType :: datapref_;
  using BaseType :: writeStep_;
  using BaseType :: outputFormat_ ;

public: 

 /** \brief Constructor creating data writer 
    \param grid corresponding grid 
    \param gridName corresponding macro grid name (needed for structured grids)
    \param data Tuple containing discrete functions to write 
    \param startTime starting time for a time dependent simulation
    \param endTime final time of a time dependent simulation

    The DataWriter is tuned through \ref Parameter 
    described under \ref DiscFuncIO. 
  */
  DataWriter(const GridType & grid,
             OutPutDataType& data,
             const DataWriterParameters& parameter = DataWriterParameters() )
    : BaseType( grid, data, parameter )
  {
    // save macro grid for structured grids 
    saveMacroGrid( parameter.macroGridName() );
  }

  DataWriter(const GridType & grid,
             OutPutDataType& data,
             const TimeProviderBase& tp,
             const DataWriterParameters& parameter = DataWriterParameters() )
    : BaseType( grid, data, tp, parameter )
  {
    // save macro grid for structured grids 
    saveMacroGrid( parameter.macroGridName() );
  }

  //! destructor 
  virtual ~DataWriter() {}

public:  
  //! print class name 
  virtual const char* myClassName() const { return "DataWriter"; }
    
  //! write binary data 
  virtual void writeBinaryData(const double sequenceStamp) const 
  {
    writeMyData( sequenceStamp, writeStep_ , data_ );
  }

  std::string writeMyData(const double sequenceStamp, 
                          const int step, 
                          Dune::Nil& data) const 
  {

    // create new path for time step output 
    std::string timeStepPath = IOInterface::createPath ( grid_.comm(),
        path_, datapref_ , step );

    // for structured grids copy grid 
    IOInterface::copyMacroGrid(grid_,path_, timeStepPath, datapref_);

    return timeStepPath;
  }

  template <class OutputTupleType>
  std::string writeMyData(const double sequenceStamp, 
                          const int step, 
                          OutputTupleType& data) const 
  {
    // create new path for time step output 
    std::string timeStepPath = IOInterface::createPath ( grid_.comm(),
        path_, datapref_ , step );

    // for structured grids copy grid 
    IOInterface::copyMacroGrid(grid_,path_, timeStepPath, datapref_);

    // create binary io obj 
    BinaryDataIO<GridType> dataio;

    // call output of IOTuple 
    IOTuple<OutputTupleType>::output(dataio, 
       grid_ , sequenceStamp, step, timeStepPath , datapref_, data );

    return timeStepPath;
  }

  /** \brief save structured macro grid file 
    \param macroFileName filename which contains the macro file 
  */
  virtual void saveMacroGrid(const std::string macroFileName) const 
  {
    IOInterface :: writeMacroGrid( grid_, macroFileName, path_, datapref_);
  }
  
}; // end class DataWriter 

//////////////////////////////////////////////////////////////////
//
//  Checkpointer 
//
//////////////////////////////////////////////////////////////////

struct CheckPointerParameters : public DataWriterParameters 
{
  //! base of file name for data file (fem.io.datafileprefix)
  virtual std::string prefix () const
  {
    return checkPointPrefix();
  }

  //! return default value for check point prefix 
  static const char * checkPointPrefix()
  {
    return "checkpoint";
  }
  
};

/** @ingroup Checkpointing 
   \brief Implementation of the IOInterface. 
   This class manages checkpointing. 

   The data will be stored in GRAPE output format, meaning that every
   checkpoint is also a visualizable data set. 
   Constructor and write are simular to the
   Dune::DataWriter, but in addition a
   static method Dune::CheckPointer::restoreGrid
   and a method Dune::CheckPointer::restoreData
   is provided. The template arguments are
   the type of the grid and a tuple type
   of pointers to the discrete functions types
   to be stored.
*/    
template <class GridImp, class DataImp = Dune::Nil > 
class CheckPointer : public DataWriter<GridImp,DataImp> 
{
protected:
  //! type of base class 
  typedef DataWriter<GridImp,DataImp> BaseType;

  using BaseType :: grid_;
  using BaseType :: data_;

  using BaseType :: path_;
  using BaseType :: datapref_;
  using BaseType :: writeStep_;
  using BaseType :: outputFormat_ ;
  using BaseType :: myRank_;
  using BaseType :: verbose_;
  using BaseType :: grapeDisplay_;

  //! type of this class  
  typedef CheckPointer<GridImp,DataImp> ThisType;
  
  //! used grid type 
  typedef GridImp GridType;
  //! used data tuple 
  typedef DataImp OutPutDataType; 

  int checkPointStep_;
  mutable int checkPointNumber_;

  std::string runPrefix_;
  std::string checkPointFile_;

  const OutPutDataType* dataPtr_;

public: 
  /** \brief Constructor generating a checkpointer 
    \param grid corresponding grid 
    \param data Tuple containing discrete functions to write

    \note In Addition to the parameters read by the DataWriter this class 
          reads the following parameters: 

    # write checkpoint every `CheckPointStep' time step
    fem.io.checkpointstep: 500 
    # store checkpoint information to file `CheckPointFile'
    fem.io.checkpointfile: checkpoint
  */
  CheckPointer(const GridType & grid, 
               OutPutDataType& data,
               const CheckPointerParameters& parameter = CheckPointerParameters() ) 
    : BaseType(grid,data,parameter)  
    , checkPointStep_(500)
    , checkPointNumber_(1)
    , dataPtr_( 0 )
  {
    initialize();
  }

  CheckPointer( const GridType & grid, const CheckPointerParameters& parameter = CheckPointerParameters() )
    : BaseType(grid, *( new OutPutDataType () ), parameter )  
    , checkPointStep_(500)
    , checkPointNumber_(1)
    , dataPtr_( &data_ )
  {
    initialize();
  }

  ~CheckPointer() 
  {
    if( dataPtr_ ) 
    {
      delete dataPtr_;
      dataPtr_ = 0;
    }
  }
protected:  
  void initialize() 
  {
    // output format can oinly be binary
    outputFormat_ = BaseType :: binary; 
    // do not display 
    grapeDisplay_ = false ;

    Parameter::get("fem.io.checkpointstep",checkPointStep_,checkPointStep_);

    /*
    // create runfile  prefix 
    {
      runPrefix_ = "run.";
      std::stringstream rankDummy;
      rankDummy << myRank_;
      runPrefix_ += rankDummy.str();
    }
    */
    
    {
      // read name of check point file 
      std::string dummyfile;
      Parameter::get("fem.io.checkpointfile",dummyfile);
      checkPointFile_ = path_; 
      checkPointFile_ += "/"; 
      checkPointFile_ += dummyfile;
    }
    
    Parameter::write("parameter.log");
  }
protected:  
  friend class CheckPointer< GridType , Dune::Nil > ;
  /** \brief Constructor generating a checkpointer to restore data 
    \param grid corresponding grid 
    \param gridName name of macro grid file (for structured grids)
    \param data Tuple containing discrete functions to write 
    \param startTime start time of simulation (needed for next save step)
    \param endTime end time of simulation
    \param checkFile filename for restoring state of program from
           previous runs 

    \note In Addition to the parameters read by the DataWriter this class 
          reads the following parameters: 

    # write checkpoint every `CheckPointStep' time step
    fem.io.checkpointstep: 500 
    # store checkpoint information to file `CheckPointFile'
    fem.io.checkpointfile: checkpoint
  */
  CheckPointer(const GridType & grid, 
               OutPutDataType& data, 
               const char * checkFile)
    : BaseType(grid, data, CheckPointerParameters() ) 
    , checkPointStep_(500)
    , checkPointNumber_(1)
    , dataPtr_( 0 )
  {
    // output format can oinly be binary
    outputFormat_ = BaseType :: binary; 
    // do not display 
    grapeDisplay_ = false ;

    datapref_ = CheckPointerParameters :: checkPointPrefix();
    Parameter::get("fem.io.checkpointstep",checkPointStep_,checkPointStep_);

    /*
    // create runf prefix 
    {
      runPrefix_ = "run.";
      std::stringstream rankDummy;
      rankDummy << myRank_;
      runPrefix_ += rankDummy.str();
    }
    */
    
    if( checkFile )
    {
      checkPointFile_ = checkFile;
      // read last counter 
      bool ok = readCheckPoint();
      // read name of check point file 
      std::string dummyfile;
      Parameter::get("fem.io.checkpointfile",dummyfile);
      checkPointFile_ = path_; 
      checkPointFile_ += "/"; 
      checkPointFile_ += dummyfile;


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
    }
    else
    {
      DUNE_THROW(InvalidStateException,"No CheckPoint file!");
    }
    
    Parameter::write("parameter.log");
  }

public:
  /** \brief restore grid from previous runs 
    \param[in] checkFile checkPoint filename 
    \param[in] rank number of my process (defaults to MPIManager :: rank() )

    \return Pointer to restored grid instance 
  */
  static GridType* restoreGrid(const std::string checkFile,
                               const int rank = MPIManager :: rank() ) 
  {
    std::string datapref( CheckPointerParameters::checkPointPrefix() );
    std::string path;
    std::string checkPointFile;

    const bool verbose = (rank == 0);

    int checkPointNumber = 0;
    // if given checkpointfile is not valid use default checkpoint file 
    if( ! readParameter(checkFile,"LastCheckPoint",checkPointNumber,verbose ) )
    {
      // read default path
      path = IOInterface::readPath();
      // set checkpointfile 
      checkPointFile = path;
      // read name of check point file 
      std::string dummyfile;
      Parameter::get("fem.io.checkpointfile",dummyfile);
      checkPointFile += "/"; checkPointFile += dummyfile;
      readParameter(checkPointFile,"LastCheckPoint",checkPointNumber, verbose);
    }
    else
    {
      if( ! readParameter(checkFile,"RecoverPath",path,verbose) )
      {
        // read default path
        path = IOInterface::readPath();
      }
    }

    // now add timestamp and rank 
    path = IOInterface::createRecoverPath(
        path, rank, datapref, checkPointNumber );

    // time is set during grid restore (not needed here)
    double time = 0.0; 

    BinaryDataIO<GridType> dataio;
    GridType* grid = IOTupleBase::restoreGrid(dataio, time, checkPointNumber, path, datapref);
    assert( grid );
    return grid;
  }

  /** \brief restores data, assumes that all objects have been created before
             this method is called
  */
  template <class InputTupleType>
  static inline void restoreData(
               const GridType& grid, 
               InputTupleType& data,
               const std::string checkFile)
  {
    // check that check point is not empty 
    if( checkFile == "" ) 
    {
      DUNE_THROW(InvalidStateException,"Checkpoint file empty!");
    }
    
    // create temporary check pointer 
    CheckPointer<GridType, InputTupleType> checker( grid, data, checkFile.c_str() );

    // restore data 
    checker.restoreData();
  }

protected:
  /** \brief restores data, assumes that all objects have been created before
   *         this method is called
  */
  void restoreData()
  {
    // now add timestamp and rank 
    std::string path = IOInterface::createRecoverPath(
        path_, myRank_ , datapref_, checkPointNumber_ );

    // restore all persistent values kept by PersistenceManager 
    PersistenceManager::restore( path );

    // restore user data 
    BinaryDataIO<GridType> dataio;
    IOTuple<OutPutDataType>::restoreData(data_, dataio, 
      grid_, checkPointNumber_, path , datapref_ );
  }

public:
  //! print class name 
  virtual const char* myClassName() const { return "CheckPointer"; }
    
  /** \brief returns true if data will be written on next write call
  */
  bool willWrite(const TimeProviderBase& tp) const
  {
    return willWrite( tp.time(), tp.timeStep() );
  }

  /** \copydoc DataWriter::willWrite */
  virtual bool willWrite(double time, int timestep) const 
  {
    // only write data time > saveTime  
    return ( (((timestep % checkPointStep_) == 0) && timestep > 0) );
  }

  virtual void writeBinaryData(const double time) const 
  {
    // toggle between 0 and 1 
    checkPointNumber_ = (checkPointNumber_ == 0) ? 1 : 0;

    // write data 
    std::string path = this->writeMyData( time, checkPointNumber_, data_ );

    // backup all persistent values kept by PersistenceManager 
    PersistenceManager::backup( path );

    // write checkpoint info 
    writeCheckPoint(path, time, 
                    checkPointNumber_ );

    return;
  }

  //! return file name for run file (if needed)
  std::string runFile() const
  {
    return path_ + "/" + runPrefix_;
  }
  
protected:
  //! read checkpoint file
  bool readCheckPoint(const bool warn = true)
  {
    // read Checkpiont file 
    if( readParameter(checkPointFile_,"LastCheckPoint",checkPointNumber_, verbose_, warn ) )
    {
      std::string recoverPath;
      // try to read recover path 
      if( ! readParameter(checkPointFile_,"RecoverPath", recoverPath, verbose_) )
      {
        // default value is output path
        recoverPath = path_;
      }

      // create processor path 
      std::string path = IOInterface :: createPath( grid_.comm(),
          recoverPath, datapref_ , checkPointNumber_);

      // copy store checkpoint run file to run file to 
      // resume from same point 
      path += "/";
      path += runPrefix_;

      // get run file name 
      std::string runFileName( runFile() );
      {
        // clear run file in case we have a restart with 
        // more processors than the checkpoint was written with
        std::ofstream clearfile ( runFileName.c_str() );
        clearfile << "# restarted from checkpoint " << std::endl;
        // close file 
        clearfile.close();
      }
       
      // create test file pointer  
      std::ifstream testfile ( path.c_str() );

      // if file exists then copy it 
      if( testfile.is_open() ) 
      {
        // close test file 
        testfile.close(); 

        std::string cmd("cp ");
        cmd += path;
        cmd += " ";
        cmd += runFileName;
        system ( cmd.c_str() );
      }
      
      // overwrite path with recover path 
      path_ = recoverPath;
      
      return true;
    }
    return false;
  }

  // write some info for checkpointing 
  void writeCheckPoint (const std::string& path,
                        const double time,
                        const int savestep ) const 
  {   
    // write some needed informantion to current checkpoint file 
    {
      std::string filepref(path);
      filepref += "/";
      filepref += datapref_;

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
    if( myRank_ <= 0)
    {
      // write last checkpoint to filename named like the checkpoint files 
      // but with no extentions 
      std::ofstream file (checkPointFile_.c_str());
      if( file.is_open() )
      {
        file << "LastCheckPoint: " << savestep << std::endl;
        file << "Time: " << std::scientific << time << std::endl;
        file << "# RecoverPath can be edited by hand if data has been moved!" << std::endl;
        file << "RecoverPath: " << path_ << std::endl;
        file.close();

        // copy checkpoint file to checkpoint path 
        std::string cmd("cp "); 
        cmd += checkPointFile_;
        cmd += " "; 
        cmd += path; 

        system ( cmd.c_str() );
      }
      else
      {
        std::cerr << "Couldn't open file `" << checkPointFile_ << "' ! " << std::endl;
      }
    }

    /*
    {
      std::string runfilename( runFile() );

      std::ifstream testfile ( runfilename.c_str() );

      // check whether file exists 
      if( testfile.is_open() ) 
      {
        // close first 
        testfile.close();

        // copy run file for recovery  
        std::string cmd;
        cmd += "cp ";
        cmd += runfilename  ;
        cmd += " ";
        cmd += path;
        system ( cmd.c_str() );
      }
    }
    */
  }

}; // end class CheckPointer 
  
} // end namespace DataIO 
#endif
