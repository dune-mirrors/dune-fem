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
#include <dune/fem/io/file/dataoutput.hh>
#include <dune/fem/misc/gridname.hh>
#include <dune/grid/common/backuprestore.hh>

#if USE_GRAPE
#include <dune/grid/io/visual/grapedatadisplay.hh>
#endif

namespace Dune {

struct DataWriterParameters : public DataOutputParameters
{
  //! base of file name for data file (fem.io.macroGridFile)
  virtual std::string macroGridName (const int dim) const
  {
    return Parameter::getValue< std::string >( IOInterface::defaultGridKey( dim ) );
  }

  //! return true if all data should be written to a spearate path per rank 
  virtual bool separateRankPath () const 
  { 
    return false; 
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
  mutable std::stringstream macroGrid_;
  const bool separateRankPath_ ;

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
    : BaseType( grid, data, parameter ),
      separateRankPath_( parameter.separateRankPath() )
  {
    if( parameter.writeMode() ) 
    {
      // save macro grid for structured grids 
      saveMacroGrid( parameter.macroGridName( GridType :: dimension ) );
    }
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
    : BaseType( grid, data, tp, parameter ),
      separateRankPath_( parameter.separateRankPath() )
  {
    if( parameter.writeMode() ) 
    {
      // save macro grid for structured grids 
      saveMacroGrid( parameter.macroGridName( GridType :: dimension ) );
    }
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

  template< class OutputTuple >
  std::string writeMyBinaryData ( const double sequenceStamp, const int step,
                                  OutputTuple &data ) const
  {
    // create new path for time step output
    std::string timeStepPath = IOInterface::createPath( grid_.comm(), path_, datapref_, step, separateRankPath_ );

    typedef IOTuple< OutputTuple > IOTupleType ;

    // if data is given, use BinaryDataIO to write it 
    if( IOTupleType :: length > 0 ) 
    {
      // for structured grids copy grid
      // IOInterface::copyMacroGrid( grid_, macroGrid_.str(), path_, timeStepPath, datapref_ );

      // call output of IOTuple
      IOTupleType :: output( grid_, sequenceStamp, timeStepPath, datapref_, data );
    }

    return timeStepPath;
  }

  /** \brief save structured macro grid file 
    \param macroFileName filename which contains the macro file 
  */
  virtual void saveMacroGrid(const std::string macroFileName) const 
  {
    IOInterface :: writeMacroGrid( grid_, macroGrid_, 
                                   macroFileName, path_, datapref_);
  }
  
}; // end class DataWriter 

//////////////////////////////////////////////////////////////////
//
//  Checkpointer 
//
//////////////////////////////////////////////////////////////////

struct CheckPointerParameters : public DataWriterParameters 
{
protected:
  bool writeMode_;

public:  
  explicit CheckPointerParameters( const bool writeMode = true ) :
    writeMode_( writeMode ) 
  {}

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

  //! writeMode, true when checkpointer is in backup mode 
  virtual bool writeMode() const 
  {
    return writeMode_;
  }

  //! return true if all data should be written to a spearate path per rank 
  virtual bool separateRankPath () const 
  { 
    return false; 
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
template< class GridImp, class DataImp = tuple<> > 
class CheckPointer
: public DataWriter< GridImp, DataImp >
{
protected:
  //! used grid type 
  typedef GridImp GridType;

  typedef DofManager< GridType > DofManagerType ;

  //! call appropriate backup and restore methods on the grid class 
  struct GridPersistentObject : public PersistentObject 
  {
    const GridType& grid_ ;
    const std::string name_;
    
    //! constructor storing grid 
    GridPersistentObject( const GridType& grid ) 
      : grid_( grid ),
        name_( Fem :: gridName( grid_ ) )
    {
      // we need to push the grid at the 
      // front position of the PersistenceManager list 
      // since we have to read it first on restore
      const bool pushFront = true ;
      // add grid at first position 
      PersistenceManager::insert( *this, pushFront );
    }

    //! destructor removing grid object 
    ~GridPersistentObject() 
    { 
      // remove myself
      PersistenceManager::remove( *this ); 
    }

    //! backup grid 
    virtual void backup() const 
    {
      // try backup using stream method first 
      try 
      { 
        std::ostream& stream = PersistenceManager :: backupStream().stream();
        Dune::BackupRestoreFacility< GridType > :: backup( grid_, stream );
      } 
      catch ( Dune :: NotImplemented ) 
      {
#ifndef NDEBUG
        if( Parameter :: verbose () )
          std::cerr << "GridPersistentObject::backup: cannot use stream backup." << std::endl;
#endif

        // try method given a filename 
        try {
          std::string filename( PersistenceManager :: uniqueFileName( name_ ) );
          Dune::BackupRestoreFacility< GridType > :: backup( grid_, filename ); 
        }
        catch ( Dune :: NotImplemented )
        {
          std::cerr << "ERROR: GridPersistentObject::backup: not possible!" << std::endl;
        }
      }

      // backup dof manager 
      DofManagerType :: instance( grid_ ).backup();
    }

    //! restore grid 
    virtual void restore () 
    {
      // restore of the grid is done below in method restoreGrid 
      // here we only need to restore the DofManager 
      DofManagerType :: instance( grid_ ).restore();
    } 
  };

  //! type of base class 
  typedef DataWriter<GridImp,DataImp> BaseType;

  using BaseType :: grid_;
  using BaseType :: data_;

  using BaseType :: path_;
  using BaseType :: datapref_;
  using BaseType :: writeStep_;
  using BaseType :: outputFormat_ ;
  using BaseType :: grapeDisplay_;
  using BaseType :: separateRankPath_;

  // friendship for restoreData calls 
  friend class CheckPointer< GridImp >;

  //! type of this class  
  typedef CheckPointer<GridImp,DataImp> ThisType;
  
  //! used data tuple 
  typedef DataImp OutPutDataType; 

  typedef GridPersistentObject PersistentGridObjectType;
  PersistentGridObjectType* persistentGridObject_ ;

  const int checkPointStep_;
  const int maxCheckPointNumber_;
  int myRank_;

  std::string checkPointFile_;

  bool takeCareOfPersistenceManager_; 

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
    , persistentGridObject_( new PersistentGridObjectType( grid_ ) ) 
    , checkPointStep_( parameter.checkPointStep() )
    , maxCheckPointNumber_( parameter.maxNumberOfCheckPoints() )
    , myRank_( grid.comm().rank() )  
    , takeCareOfPersistenceManager_( true )
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
    , persistentGridObject_( new PersistentGridObjectType( grid_ ) ) 
    , checkPointStep_( parameter.checkPointStep() )
    , maxCheckPointNumber_( parameter.maxNumberOfCheckPoints() )
    , myRank_( grid.comm().rank() )  
    , takeCareOfPersistenceManager_( true )
  {
    initialize( parameter );
  }

  ~CheckPointer() 
  {
    // remove persistent grid object from PersistenceManager
    delete persistentGridObject_;
    persistentGridObject_ = 0;
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
  /** \brief Constructor generating a checkpointer to restore data 
    \param grid corresponding grid 
    \param myRank rank of process 
    \param data Tuple containing discrete functions to write 
    \param checkFile filename for restoring state of program from
           previous runs 
    \param takeCareOfPersistenceManager flag whether to keep persistent 
           values by PersistenceManager (default true)
    \note In Addition to the parameters read by the DataWriter this class 
          reads the following parameters: 

    # write checkpoint every `CheckPointStep' time step
    fem.io.checkpointstep: 500 
    # store checkpoint information to file `CheckPointFile'
    fem.io.checkpointfile: checkpoint
  */
  CheckPointer(const GridType & grid, 
               const int myRank,
               OutPutDataType& data, 
               const char * checkFile,
               const bool takeCareOfPersistenceManager = true,
               const int writeStep = 0 )
    : BaseType(grid, data, CheckPointerParameters( checkFile == 0 ) ) // checkFile != 0 means read mode 
    , persistentGridObject_( 0 ) // do not create a persistent object here, since we are in read mode
    , checkPointStep_( 0 )
    , maxCheckPointNumber_( writeStep + 1 )
    , myRank_( myRank )  
    , takeCareOfPersistenceManager_( takeCareOfPersistenceManager )
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
      initialize( CheckPointerParameters( checkFile == 0 ) );
    }

    // set write step counter to value given in constructor
    writeStep_ = writeStep ;
  }

public:
  /** \brief restore grid from previous runs 
    \param[in] checkFile checkPoint filename 
    \param[in] givenRank number of my process (defaults to MPIManager :: rank() )

    \return Pointer to restored grid instance 
  */
  static GridType* restoreGrid(const std::string checkFile,
                               const int givenRank = -1,
                               const CheckPointerParameters& parameter = CheckPointerParameters() )
  {
    const int rank = ( givenRank < 0 ) ? MPIManager :: rank() : givenRank ;
    std::string datapref( parameter.checkPointPrefix() );
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
      checkPointFile += parameter.checkPointPrefix(); 
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

    // now add timestamp (and rank) 
    path = IOInterface::createRecoverPath(
        path, rank, datapref, checkPointNumber, parameter.separateRankPath() );

    // initialize PersistenceManager 
    PersistenceManager :: startRestore ( path ); 

    GridType* grid = 0;
    try 
    {
      std::istream& stream = PersistenceManager :: restoreStream().stream();
      grid = Dune::BackupRestoreFacility< GridType > :: restore( stream );
    }
    catch ( Dune :: NotImplemented ) 
    {
#ifndef NDEBUG
      if( Parameter :: verbose () )
        std::cerr << "GridPersistentObject::restore: cannot use stream restore." << std::endl;
#endif
      try {
        std::string name ( Fem :: gridName( *grid ) );
        std::string filename( PersistenceManager :: uniqueFileName( name ) );
        grid = Dune::BackupRestoreFacility< GridType > :: restore( filename ); 
      }
      catch ( Dune :: NotImplemented )
      {
        std::cerr << "ERROR: GridPersistentObject::restore: not possible!" << std::endl;
      }
    }

    if( grid == 0 ) 
    {
      DUNE_THROW(InvalidStateException,"Could not recover grid");
    }

    return grid;
  }

  /** \brief restores data, assumes that all objects have been created and inserted to
   *         PersistenceManager before this method is called
   *
   *  \param grid Grid the data belong to 
   *  \param checkFile check point file 
  */
  static inline 
  void restoreData ( const GridType &grid, const std::string checkFile )
  {
    tuple<> fakeData;
    restoreData( grid, fakeData, checkFile );
  }

  /** \brief restores data, assumes that all objects have been created and inserted to
   *         PersistenceManager before this method is called
   *
   *  \param grid Grid the data belong to 
   *  \param data tuple of discrete functions to be additionally read during restore  
   *  \param checkFile check point file 
   *  \param rank rank of process (defaults to grid.comm().rank())
  */
  template <class InputTupleType>
  static inline 
  void restoreData(const GridType& grid, 
                   InputTupleType& data,
                   const std::string checkFile,
                   const int rank = -1 )
  {
    // make rank exchangable 
    const int myRank = ( rank < 0 ) ? grid.comm().rank() : rank ;

    // check that check point is not empty 
    if( checkFile == "" ) 
    {
      DUNE_THROW(InvalidStateException,"Checkpoint file empty!");
    }
    
    // create temporary check pointer 
    CheckPointer<GridType, InputTupleType> checker( grid, myRank, data, checkFile.c_str() );

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
        path_, myRank_ , datapref_, writeStep_, separateRankPath_ );

    // if true also restore PersistenceManager 
    if( takeCareOfPersistenceManager_ ) 
    {
      // restore all persistent values kept by PersistenceManager 
      PersistenceManager::restore( path );
    }

    return path;
  }

  template< class InputTuple >
  void restoreUserData ( InputTuple &data )
  {
    // restore persistent data 
    std::string path = restorePersistentData( );

    typedef IOTuple< InputTuple > IOTupleType ;

    // restore user data 
    if( IOTupleType :: length > 0 ) 
    {
      IOTupleType :: restoreData( data, grid_, path , datapref_ );
    }

    // make data consecutive at the end of the restore process 
    DofManagerType :: instance( grid_ ).compress();
  }

  void restoreData( ) 
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
    const int timestep = tp.timeStep();
    // only write data time > saveTime  
    return ( (checkPointStep_ > 0) && (((timestep % checkPointStep_) == 0) && timestep > 0) );
  }

  template <class OutputTuple> 
  static void writeSingleCheckPoint(const GridType& grid, 
                                    OutputTuple& data,
                                    const double time,
                                    const bool storePersistenceManager, 
                                    const int writeStep = 0 )
  {
    CheckPointer< GridType, OutputTuple > checkPointer( grid, grid.comm().rank(),
                                                        data, 0, storePersistenceManager, writeStep );
    checkPointer.writeBinaryData( time );
  }

  virtual void writeBinaryData(const double time) const 
  {
    // reset writeStep_ when maxCheckPointNumber_ is reached 
    if( writeStep_ >= maxCheckPointNumber_ ) writeStep_ = 0;

    // write data 
    std::string path = writeMyBinaryData( time, writeStep_, data_ );

    // if true also backup PersistenceManager 
    if( takeCareOfPersistenceManager_ ) 
    {
      // backup all persistent values kept by PersistenceManager 
      PersistenceManager::backup( path );
    }

    // write checkpoint info 
    writeCheckPoint(path, time, 
                    writeStep_ );

    return;
  }

protected:
  //! read checkpoint file
  bool readCheckPoint(const bool warn = true)
  {
    const bool verbose = Parameter::verbose(); 

    // read Checkpiont file 
    if( readParameter(checkPointFile_,"LastCheckPoint",writeStep_, verbose, warn ) )
    {
      std::string recoverPath;
      // try to read recover path 
      if( ! readParameter(checkPointFile_,"RecoverPath", recoverPath, verbose) )
      {
        // default value is output path
        recoverPath = path_;
      }

      int storedPersistentManager = 0;
      if( readParameter(checkPointFile_,"PersistenceManager", storedPersistentManager, verbose) )
      {
        takeCareOfPersistenceManager_ = ( storedPersistentManager > 0 );
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
    // only proc 0 writes the global checkpoint file 
    if( myRank_ <= 0 )
    {
      std::string checkpointstr ; 
      {
        std::stringstream checkpoint; 
        checkpoint << "LastCheckPoint: " << savestep << std::endl;
        checkpoint.precision( 16 );
        checkpoint << "Time: " << std::scientific << time << std::endl;
        checkpoint << "SaveCount: " << savestep << std::endl;
        checkpoint << "PersistenceManager: " << takeCareOfPersistenceManager_ << std::endl;
        checkpoint << "NumberProcessors: " << grid_.comm().size() << std::endl;
        checkpoint << "# RecoverPath can be edited by hand if data has been moved!" << std::endl;
        checkpoint << "RecoverPath: " << path_ << std::endl;
        checkpointstr = checkpoint.str();
      }

      // overwrite the last checkpoint file 
      {
        std::ofstream file (checkPointFile_.c_str());
        if( file.is_open() )
        {
          file << checkpointstr;
        }
      }

      // write check point file for this checkpoint
      {
        std::string checkPointStepFile( path );
        checkPointStepFile += "/" + datapref_;

        std::ofstream file ( checkPointStepFile.c_str() );
        if( file.is_open() )
        {
          file << checkpointstr;
        }
      }
    }
  }

}; // end class CheckPointer 
  
} // end namespace DataIO 
#endif
