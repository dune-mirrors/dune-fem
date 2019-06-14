#ifndef DUNE_FEM_DATAWRITER_HH
#define DUNE_FEM_DATAWRITER_HH

#include <string>
#include <tuple>
#include <limits>
#include <memory>

#include <dune/fem/io/file/asciiparser.hh>
#include <dune/fem/io/file/iointerface.hh>
#include <dune/fem/io/file/iotuple.hh>
#include <dune/fem/io/parameter.hh>
#include <dune/fem/space/common/loadbalancer.hh>
#include <dune/fem/gridpart/adaptiveleafgridpart.hh>
#include <dune/fem/space/lagrange.hh>
#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/quadrature/cachingquadrature.hh>
#include <dune/fem/io/file/persistencemanager.hh>
#include <dune/fem/io/file/dataoutput.hh>
#include <dune/fem/misc/gridname.hh>
#include <dune/grid/common/backuprestore.hh>

#if USE_GRAPE
#include <dune/grid/io/visual/grapedatadisplay.hh>
#endif

namespace Dune
{

  namespace Fem
  {

    struct DataWriterParameters : public DataOutputParameters
    {
      explicit DataWriterParameters ( std::string keyPrefix, const ParameterReader &parameter = Parameter::container() )
        : DataOutputParameters( keyPrefix, parameter )
      {}

      explicit DataWriterParameters ( const ParameterReader &parameter = Parameter::container() )
        : DataOutputParameters( parameter )
      {}

      //! base of file name for data file (fem.io.macroGridFile)
      virtual std::string macroGridName (const int dim) const
      {
        return parameter().getValue< std::string >( IOInterface::defaultGridKey( dim, parameter() ) );
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
          // call output of IOTuple
          IOTupleType :: output( grid_, sequenceStamp, timeStepPath, datapref_, data );
        }

        return timeStepPath;
      }
    }; // end class DataWriter

    //////////////////////////////////////////////////////////////////
    //
    //  Checkpointer
    //
    //////////////////////////////////////////////////////////////////

    /**
     * \class CheckPointerParameters
     * \brief local parameter collection for CheckPointer
     *
     * \note For now, the CheckPointer only works with the singleton parameter
     *       class.
     **/
    struct CheckPointerParameters : public DataWriterParameters
    {
    protected:
      bool writeMode_;

    public:
      CheckPointerParameters( const bool writeMode, const std::string keyPrefix = "fem.io." ) :
        DataWriterParameters( keyPrefix, Parameter::container() ),
        writeMode_( writeMode )
      {}

      explicit CheckPointerParameters( const std::string keyPrefix = "fem.io." ) :
        DataWriterParameters( keyPrefix, Parameter::container() ),
        writeMode_( true )
      {}

      //! base of file name for data file (fem.io.datafileprefix)
      virtual std::string prefix () const
      {
        return checkPointPrefix();
      }

      //! return number of timestep to be passed until next checkpoint in written
      virtual int checkPointStep() const
      {
        return parameter().getValue< int > ( keyPrefix_ + "checkpointstep", 500 );
      }

      //! maximal number of checkpoint stages written (default = 2)
      virtual int maxNumberOfCheckPoints() const
      {
        return parameter().getValue< int > ( keyPrefix_ + "checkpointmax", 2 );
      }

      //! return default value for check point prefix
      virtual std::string checkPointPrefix() const
      {
        return parameter().getValue< std::string > ( keyPrefix_ +  "checkpointfile", "checkpoint" );
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

      virtual int outputformat() const
      {
        return 5; // i.e. none
      }
    };

    /** @ingroup Checkpointing
       \brief Implementation of the IOInterface.
       This class manages checkpointing.

       All data that was registered to PersistenceManager will be stored in binary output format.
       The derivation from DataWriter is simply to use the writeStep method. The
       binary output of DataWriter is not used anymore and does not work for
       checkpointing.
    */
    template< class GridImp >
    class CheckPointer
    : public DataWriter< GridImp, std::tuple<> >
    {
      typedef std::tuple<> DataImp;
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
          // this feature is available in dune-grid 2.3.x and later

          // try backup using stream method first
          try
          {
            std::string gridBackup ;
            // get backup stream from grid facility
            {
              std::stringstream gridBackupStream;
              Dune::BackupRestoreFacility< GridType > :: backup( grid_, gridBackupStream );
              gridBackup = gridBackupStream.str();
            }
            PersistenceManager :: backupStream() << gridBackup;
          }
          catch ( const Dune :: NotImplemented& )
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
            catch ( const Dune :: NotImplemented& )
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

      //! type of this class
      typedef CheckPointer< GridImp > ThisType;

      //! used data tuple
      typedef DataImp OutPutDataType;
      OutPutDataType fakeData_; // empty tuple

      typedef GridPersistentObject PersistentGridObjectType;
      std::unique_ptr< PersistentGridObjectType > persistentGridObject_ ;

      OutPutDataType* dataPtr_ ;

      const int checkPointStep_;
      const int maxCheckPointNumber_;
      int myRank_;

      std::string checkPointFile_;

      bool takeCareOfPersistenceManager_;

    public:
      /** \brief Constructor generating a checkpointer
        \param grid corresponding grid
        \param parameter structure for tuning the behavior of the Dune::CheckPointer
                         defaults to Dune::CheckPointerParameters
      */
      CheckPointer(const GridType & grid,
                   const CheckPointerParameters& parameter = CheckPointerParameters() )
        : BaseType(grid, fakeData_, parameter)
        , persistentGridObject_( new PersistentGridObjectType( grid_ ) )
        , dataPtr_ ( 0 )
        , checkPointStep_( parameter.checkPointStep() )
        , maxCheckPointNumber_( parameter.maxNumberOfCheckPoints() )
        , myRank_( grid.comm().rank() )
        , takeCareOfPersistenceManager_( true )
      {
        initialize( parameter );
      }

      /** \brief Constructor generating a checkpointer
        \param grid corresponding grid
        \param data Tuple containing discrete functions to write
        \param tp   a time provider to set time (e.g. for restart)
        \param parameter structure for tuning the behavior of the Dune::CheckPointer
                         defaults to Dune::CheckPointerParameters
      */
      CheckPointer(const GridType & grid,
                   const TimeProviderBase& tp,
                   const CheckPointerParameters& parameter = CheckPointerParameters() )
        : BaseType(grid, fakeData_,tp,parameter)
        , persistentGridObject_( new PersistentGridObjectType( grid_ ) )
        , checkPointStep_( parameter.checkPointStep() )
        , maxCheckPointNumber_( parameter.maxNumberOfCheckPoints() )
        , myRank_( grid.comm().rank() )
        , takeCareOfPersistenceManager_( true )
      {
        initialize( parameter );
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

        // write parameter file
        Parameter::write("parameter.log");
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
        \param writeStep initial counter value, default 0
        \note In Addition to the parameters read by the DataWriter this class
              reads the following parameters:

        # write checkpoint every `CheckPointStep' time step
        fem.io.checkpointstep: 500
        # store checkpoint information to file `CheckPointFile'
        fem.io.checkpointfile: checkpoint
      */
      CheckPointer(const GridType & grid,
                   const int myRank,
                   const char * checkFile,
                   const bool takeCareOfPersistenceManager = true,
                   const int writeStep = 0 )
        : BaseType(grid, fakeData_, CheckPointerParameters( checkFile == 0 ) ) // checkFile != 0 means read mode
        , persistentGridObject_( ) // do not create a persistent object here, since we are in read mode
        , checkPointStep_( 0 )
        , maxCheckPointNumber_( writeStep + 1 )
        , myRank_( myRank )
        , takeCareOfPersistenceManager_( takeCareOfPersistenceManager )
      {
        // output format can oinly be binary
        outputFormat_ = BaseType :: binary;
        // do not display
        grapeDisplay_ = false ;

        CheckPointerParameters parameter( checkFile == 0 );
        datapref_ = parameter.checkPointPrefix();

        if( checkFile )
        {
          // try to read given check point file
          checkPointFile_ = checkFile;

          // read last counter, don't issue warning
          bool ok = readCheckPoint( false );

          // if check point couldn't be opened, try again with default
          if(!ok)
          {

            // read name of check point file
            checkPointFile_ = path_;
            checkPointFile_ += "/";
            checkPointFile_ += parameter.checkPointPrefix();

            const bool warn = (myRank == 0);
            if( warn )
            {
              std::cerr << "WARNING: Coudn't open file `" << checkFile << "' trying file `" << checkPointFile_ << "' instead!" << std::endl;
            }

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
        \param[in] parameter Parameterclass which provides informations about the
        checkpoint

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
        if( ! readParameter(checkFile,"LastCheckPoint",checkPointNumber, verbose, false ) )
        {
          // read default path
          path = IOInterface::readPath();
          // set checkpointfile
          std::string checkPointFile = path;
          // try out default checkpoint file
          checkPointFile += "/";
          checkPointFile += parameter.checkPointPrefix();
          if ( verbose )
          {
            std::cerr << "WARNING: Coudn't open file `" << checkFile << "' trying file `" << checkPointFile << "' instead!" << std::endl;
          }
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
        // this is only available in dune-grid 2.3.x and later
        try
        {
          std::string gridData;
          PersistenceManager :: restoreStream() >> gridData;

          // copy data to stream
          std::stringstream gridStream( gridData );
          // clear grid data
          gridData = std::string();

          // perform restore using grid stream only
          grid = Dune::BackupRestoreFacility< GridType > :: restore( gridStream );
        }
        catch ( const Dune :: NotImplemented& )
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
          catch ( const Dune :: NotImplemented& )
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
      void restoreData ( const GridType &grid, const std::string checkFile, const int rank=-1  )
      {
        // make rank exchangable
        const int myRank = ( rank < 0 ) ? grid.comm().rank() : rank ;

        // check that check point is not empty
        if( checkFile == "" )
        {
          DUNE_THROW(InvalidStateException,"Checkpoint file empty!");
        }

        // create temporary check pointer
        CheckPointer<GridType> checker( grid, myRank, checkFile.c_str() );

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
        // restore of the grid is done in method restoreGrid
        // here we only need to restore the DofManager since now
        // all spaces and index sets have been created
        DofManagerType :: instance( grid_ ).restore();

        // restore persistent data
        std::string path = restorePersistentData( );

        // make data consecutive at the end of the restore process
        // and communicate data
        DofManagerType :: instance( grid_ ).compress();
      }

      void restoreData( )
      {
        restoreUserData( data_ );
      }

    public:
      //! print class name
      virtual const char* myClassName() const { return "CheckPointer"; }

      using BaseType :: willWrite ;
      /** \brief returns true if data will be written on next write call
      */
      bool willWrite(const TimeProviderBase& tp) const
      {
        const int timestep = tp.timeStep();
        // only write data time > saveTime
        return ( (checkPointStep_ > 0) && (((timestep % checkPointStep_) == 0) && timestep > 0) );
      }

      static void writeSingleCheckPoint(const GridType& grid,
                                        const double time,
                                        const bool storePersistenceManager,
                                        const int writeStep = 0 )
      {
        CheckPointer< GridType > checkPointer( grid );
        checkPointer.writeBinaryData( time );
      }

      virtual void writeBinaryData(const double time) const
      {
        // reset writeStep_ when maxCheckPointNumber_ is reached
        if( writeStep_ >= maxCheckPointNumber_ ) writeStep_ = 0;

        // write data
        std::string path = this->writeMyBinaryData( time, writeStep_, data_ );

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

  } // end namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_DATAWRITER_HH
