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

  friend class DataOutput< GridImp, DataImp >;

public: 

  /** \brief Constructor creating data writer 
    \param grid corresponding grid 
    \param data Tuple containing discrete functions to write 
    \param parameter structure for tuning the behavior of the Dune::DataWriter 
                     defaults to Dune::DataWriterParameters
  */
  DataWriter(const GridType & grid,
             OutPutDataType& data,
             const DataWriterParameters& parameter = DataWriterParameters() )
    : BaseType( grid, data, parameter )
  {
    // save macro grid for structured grids 
    saveMacroGrid( parameter.macroGridName() );
  }

  /** \brief Constructor creating data writer 
    \param grid corresponding grid 
    \param data Tuple containing discrete functions to write 
    \param tp   a time provider to set time (e.g. for restart)
    \param parameter structure for tuning the behavior of the Dune::DataWriter
                     defaults to Dune::DataWriterParameters
  */
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

protected:  
  //! print class name 
  virtual const char* myClassName() const { return "DataWriter"; }
    
  //! write binary data 
  virtual void writeBinaryData(const double sequenceStamp) const 
  {
    writeMyBinaryData( sequenceStamp, writeStep_ , data_ );
  }

  std::string writeMyBinaryData(const double sequenceStamp, 
                                const int step, 
                                Dune::Nil& data) const 
  {

    // create new path for time step output 
    std::string timeStepPath = IOInterface::createPath ( grid_.comm(),
        path_, datapref_ , step );

    // for structured grids copy grid 
    IOInterface::copyMacroGrid(grid_,path_, timeStepPath, datapref_);

    // create binary io obj 
    BinaryDataIO<GridType> dataio;

    // call writeGrid of IOTupleBase
    IOTupleBase::writeGrid(dataio, grid_ , sequenceStamp, step, timeStepPath, datapref_ );

    return timeStepPath;
  }

  template <class OutputTupleType>
  std::string writeMyBinaryData(const double sequenceStamp, 
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

  //! return number of timestep to be passed until next checkpoint in written 
  virtual int checkPointStep() const 
  {
    return Parameter::getValue< int > ("fem.io.checkpointstep", 500 );
  }

  //! maximal number of checkpoint stages written (default = 2)
  virtual int maxNumberOfCheckPoints() const 
  {
    return Parameter :: getValue< int > ("fem.io.checkpointmax", 2 );
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

  const int checkPointStep_;
  const int maxCheckPointNumber_;

  std::string checkPointFile_;

  const OutPutDataType* dataPtr_;

public: 
  /** \brief Constructor generating a checkpointer 
    \param grid corresponding grid 
    \param data Tuple containing discrete functions to write 
    \param tp   a time provider to set time (e.g. for restart)
    \param parameter structure for tuning the behavior of the Dune::CheckPointer
                     defaults to Dune::CheckPointerParameters
  */
  CheckPointer(const GridType & grid, 
               OutPutDataType& data,
               const TimeProviderBase& tp,
               const CheckPointerParameters& parameter = CheckPointerParameters() ) 
    : BaseType(grid,data,tp,parameter)  
    , checkPointStep_( parameter.checkPointStep() )
    , maxCheckPointNumber_( parameter.maxNumberOfCheckPoints() )
    , dataPtr_( 0 )
  {
    initialize( parameter );
  }

  /** \brief Constructor generating a checkpointer 
    \param grid corresponding grid 
    \param tp   a time provider to set time (e.g. for restart)
    \param parameter structure for tuning the behavior of the Dune::CheckPointer
                     defaults to Dune::CheckPointerParameters
  */
  CheckPointer( const GridType & grid, 
                const TimeProviderBase& tp,
                const CheckPointerParameters& parameter = CheckPointerParameters() )
    : BaseType(grid, *( new OutPutDataType () ), tp, parameter )  
    , checkPointStep_( parameter.checkPointStep() )
    , maxCheckPointNumber_( parameter.maxNumberOfCheckPoints() )
    , dataPtr_( &data_ )
  {
    initialize( parameter );
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
  void initialize( const CheckPointerParameters& parameter ) 
  {
    // output format can only be binary
    outputFormat_ = BaseType :: binary; 
    // do not display 
    grapeDisplay_ = false ;

    checkPointFile_ = path_; 
    checkPointFile_ += "/"; 
    checkPointFile_ += parameter.prefix();
  }
protected:  
  friend class CheckPointer< GridType , Dune::Nil > ;
  /** \brief Constructor generating a checkpointer to restore data 
    \param grid corresponding grid 
    \param data Tuple containing discrete functions to write 
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
    , checkPointStep_( 0 )
    , maxCheckPointNumber_( 0 )
    , dataPtr_( 0 )
  {
    // output format can oinly be binary
    outputFormat_ = BaseType :: binary; 
    // do not display 
    grapeDisplay_ = false ;

    datapref_ = CheckPointerParameters :: checkPointPrefix();

    if( checkFile )
    {
      // try to read given check point file 
      checkPointFile_ = checkFile;
      // read last counter 
      bool ok = readCheckPoint();

      // if check point couldn't be opened, try again with default  
      if(!ok)
      {
        // read name of check point file 
        checkPointFile_ = path_; 
        checkPointFile_ += "/"; 
        checkPointFile_ += CheckPointerParameters :: checkPointPrefix();

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

    const bool verbose = (rank == 0);

    int checkPointNumber = 0;
    // if given checkpointfile is not valid use default checkpoint file 
    if( ! readParameter(checkFile,"LastCheckPoint",checkPointNumber,verbose ) )
    {
      // read default path
      path = IOInterface::readPath();
      // set checkpointfile 
      std::string checkPointFile = path;
      // try out default checkpoint file 
      checkPointFile += "/"; 
      checkPointFile += CheckPointerParameters :: checkPointPrefix(); 
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

  /** \brief restores data, assumes that all objects have been created and inserted to
   *         PersistenceManager before this method is called
   *
   *  \param grid Grid the data belong to 
   *  \param checkFile check point file 
  */
  static inline 
  void restoreData(const GridType& grid, const std::string checkFile)
  {
    Dune::Nil fakeData;
    restoreData( grid, fakeData, checkFile );
  }

  /** \brief restores data, assumes that all objects have been created and inserted to
   *         PersistenceManager before this method is called
   *
   *  \param grid Grid the data belong to 
   *  \param data tuple of discrete functions to be additionally read during restore  
   *  \param checkFile check point file 
  */
  template <class InputTupleType>
  static inline 
  void restoreData(const GridType& grid, 
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
  std::string restorePersistentData()
  {
    // now add timestamp and rank 
    std::string path = IOInterface::createRecoverPath(
        path_, myRank_ , datapref_, writeStep_ );

    // restore all persistent values kept by PersistenceManager 
    PersistenceManager::restore( path );

    return path;
  }

  void restoreUserData( Dune::Nil& ) 
  {
    // restore persistent data 
    restorePersistentData();
  }

  template <class InputTupleType>
  void restoreUserData(InputTupleType & data) 
  {
    // restore persistent data 
    std::string path = restorePersistentData();

    // restore user data 
    BinaryDataIO<GridType> dataio;
    IOTuple<OutPutDataType>::restoreData(data, dataio, 
      grid_, writeStep_, path , datapref_ );
  }

  void restoreData() 
  {
    restoreUserData( data_ );
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
    // reset writeStep_ when maxCheckPointNumber_ is reached 
    if( writeStep_ >= maxCheckPointNumber_ ) writeStep_ = 0;

    // write data 
    std::string path = this->writeMyBinaryData( time, writeStep_, data_ );

    // backup all persistent values kept by PersistenceManager 
    PersistenceManager::backup( path );

    // write checkpoint info 
    writeCheckPoint(path, time, 
                    writeStep_ );

    return;
  }

protected:
  //! read checkpoint file
  bool readCheckPoint(const bool warn = true)
  {
    // read Checkpiont file 
    if( readParameter(checkPointFile_,"LastCheckPoint",writeStep_, verbose_, warn ) )
    {
      std::string recoverPath;
      // try to read recover path 
      if( ! readParameter(checkPointFile_,"RecoverPath", recoverPath, verbose_) )
      {
        // default value is output path
        recoverPath = path_;
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
  }

}; // end class CheckPointer 
  
} // end namespace DataIO 
#endif
