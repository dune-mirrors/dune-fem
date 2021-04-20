#ifndef DUNE_FEM_PERSISTENCEMANAGER_HH
#define DUNE_FEM_PERSISTENCEMANAGER_HH

#include <fstream>
#include <iostream>
#include <list>
#include <sstream>
#include <string>
#include <type_traits>
#include <utility>

#include <dune/fem/io/file/iointerface.hh>
#include <dune/fem/io/parameter.hh>
#include <dune/fem/io/streams/binarystreams.hh>
#include <dune/fem/io/streams/tuples.hh>
#include <dune/fem/storage/singleton.hh>

namespace Dune
{

  namespace Fem
  {

  /** @addtogroup Checkpointing

      The Dune::Fem::PersistenceManager manages a list of persistent
      object the state should be recovered during a restart process.
      The singleton instance is available using persistenceManager.

      Variables of any type can be added to the manager
      using operator<< and removed through operator>>.
      For general types the state is written to an ascii file
      \c persistentobjects through the call to PersistenceManager::backup(path).

      Finer control is available for user defined classes by
      deriving these from the Dune::PersistentObject class and
      implementing the virtual functions \c backup and \c restore.
      Values can then either be written to the ascii file
      \c persistentobjects using the method
      Dune::Fem::PersistenceManager::backupValue giving a token and
      the value to store (a structure simular to the Dune::Parameter
      files is used).
      If data is to be written to a different stream a unique filename
      can be obtained from the Dune::PersistenceManager using the
      \c uniqueFileName method.

      If a class is derived from Dune::AutoPersistenceObject
      then instances of this class are automatically added to the
      PersistenceManager on construction and removed on destruction;
      again the virtual methods backup/restore must be provided.

      The restart process is activated through the method \c restore.
      Note that care must be taken to add the objects in exactly
      the same order when using the PersistenceManager for creating
      the backup file as when used to restore the saved states.

      The Dune::Fem::Parameter class is automatically backuped and restored.
      For the restore process the saved parameter file is read first and
      existing values for parameters are replaced by these values.
      Therefore it is possible to reread the state of a variable
      orignally defined through a run time parameter using the Parameter
      class:
      \code
      class A : public Dune::AutoPersistentObject {
        int var_,b_;
        ComplexType U_;
        A() : var_(Dune::Parameter::getValue("token",0,var_) {}
        virtual backup() const {
          // note var_ is saved in the parameter file

          // write b to ascii file
          PersistenceManager::backupValue("b",b);

          // write U to stream
          ofstream out((PersistenceManager::uniqueFileName()+"classA").c_str());
          U.write(out);
        }
        virtual restore() {
          // set var_ through parameter file
          var_= Dune::Parameter::getValue<int>("token",var_);

          // read b from ascii file
          PersistenceManager::restoreValue("b",b);

          // read U from stream
          ifstream in((PersistenceManager::uniqueFileName()+"classA").c_str());
          U.read(in);
        }
      };
      \endcode

    */

    class PersistenceManager;

    /** \class   PersistentObject
        \ingroup Checkpointing
        \brief   base class for persistent objects
     */
    class PersistentObject
    {
      typedef PersistentObject ThisType;

      friend class PersistenceManager;

    public:
      virtual ~PersistentObject() {}
      /** \brief backup persistent object */
      virtual void backup () const = 0;
      /** \brief restore persistent object */
      virtual void restore () = 0;

    protected:
      /** \brief insert possible sub data of object */
      virtual void insertSubData() {}
      /** \brief remove possible sub data of object */
      virtual void removeSubData() {}

      virtual void *pointer ()
      {
        return this;
      }
    };



    template< class ObjectType >
    struct IsPersistent
    {
      static const bool value = std::is_convertible< ObjectType *, PersistentObject * >::value;
    };



    /** \class   PersistenceManager
     *  \ingroup Checkpointing
     *  \brief   class with singleton instance managing all
     *           persistent objects
     */
    class PersistenceManager
    {
      typedef PersistenceManager ThisType;
      template <class ObjectType,bool isPersistent>
      struct WrapObject;

    public:
      // make backup and restore streams exchangeable
#ifdef FEM_PERSISTENCEMANAGERSTREAMTRAITS
      typedef FEM_PERSISTENCEMANAGERSTREAMTRAITS :: BackupStreamType  BackupStreamType;
      typedef FEM_PERSISTENCEMANAGERSTREAMTRAITS :: RestoreStreamType RestoreStreamType;
      static const bool singleBackupRestoreFile = FEM_PERSISTENCEMANAGERSTREAMTRAITS ::
        singleBackupRestoreFile ;
#else
      typedef Fem :: BinaryFileOutStream  BackupStreamType ;
      typedef Fem :: BinaryFileInStream   RestoreStreamType ;
      static const bool singleBackupRestoreFile = false ;
#endif

    private:
      typedef std::list< std::pair< PersistentObject *, unsigned int > > PersistentType;
      typedef PersistentType::iterator IteratorType;

      PersistentType objects_;
      int fileCounter_,lineNo_;
      std::string path_;
      std::ifstream inAsciStream_;
      std::ofstream outAsciStream_;
      bool closed_,invalid_;

      BackupStreamType*  backupStream_;
      RestoreStreamType* restoreStream_;

      PersistenceManager ()
      : fileCounter_( 0 ),
        lineNo_(),
        path_(),
        closed_( false ),
        invalid_( false ),
        backupStream_( 0 ),
        restoreStream_( 0 )
      {}

      PersistenceManager ( const ThisType & );
      ThisType &operator= ( const ThisType & );

      BackupStreamType& backupStreamObj ()
      {
        assert( backupStream_ );
        return *backupStream_ ;
      }

      RestoreStreamType& restoreStreamObj ()
      {
        assert( restoreStream_ );
        return *restoreStream_ ;
      }

    public:
      template< class ObjectType >
      void insertObject( ObjectType &object, const bool pushFront = false  )
      {
        IteratorType end = objects_.end();
        for( IteratorType it = objects_.begin(); it != end; ++it )
        {
          if( it->first->pointer() != &object )
            continue;
          ++it->second;
          return;
        }

        // for objects like the grid
        // we need to allow to add this later
        if ( closed_ && ! pushFront )
        {
#ifndef NDEBUG
          std::cerr << "WARNING: new object added to PersistenceManager "
                    << "although backup/restore has been called - "
                    << "Object will be ignored!" << std::endl;
#endif
          return;
        }

        PersistentObject *obj =
          WrapObject< ObjectType, IsPersistent< ObjectType > :: value >
          :: apply( object );

        // insert possible sub data
        obj->insertSubData();

        if( pushFront )
          objects_.push_front( std :: make_pair( obj, 1 ) );
        else
          objects_.push_back( std :: make_pair( obj, 1 ) );
      }

      template< class ObjectType >
      void removeObject ( ObjectType &object )
      {
        IteratorType end = objects_.end();
        for( IteratorType it = objects_.begin(); it != end; ++it )
        {
          if( it->first->pointer() != &object )
            continue;

          --it->second;
          if( it->second == 0 )
          {
            if (closed_) invalid_=true;
            PersistentObject *obj = it->first;
            // remove possible sub data
            obj->removeSubData();
            objects_.erase( it );
            if( !IsPersistent< ObjectType > :: value )
              delete obj;
          }
          return;
        }
      }

      void backupObjects ( const std::string& path )
      {
        if( invalid_ )
        {
#ifndef NDEBUG
          std::cerr << "WARNING: backup called although objects "
                    << "have been removed from the PersistenceManager! "
                    << "Backup ignored!" << std::endl;
#endif
          return;
        }
        closed_ = true;
        startBackup( path );
        typedef PersistentType::iterator IteratorType;

        for( IteratorType it = objects_.begin(); it != objects_.end(); ++it )
          it->first->backup();

        closeStreams();
      }

      void restoreObjects ( const std::string &path )
      {
        if (invalid_) {
#ifndef NDEBUG
          std::cerr << "WARNING: restore called although objects "
                    << "have been removed from the PersistenceManager! "
                    << "Restore ignored!" << std::endl;
#endif
          return;
        }
        closed_ = true;
        startRestoreImpl( path );
        typedef PersistentType :: iterator IteratorType;

        for( IteratorType it = objects_.begin(); it != objects_.end(); ++it )
          it->first->restore( );

        closeStreams();
      }

      std::string getUniqueFileName ( const std::string &tag )
      {
        return generateFilename( path_ + "/" + tag, ++fileCounter_ );
      }

      std::string getUniqueTag ( const std::string &tag )
      {
        return generateFilename( tag, ++fileCounter_ );
      }

      template< class T >
      void backup ( const std::string &token, const T &value )
      {
        backupStreamObj() << token;
        backupStreamObj() << value;
      }

      template< class T >
      void restore ( const std::string &token, T &value )
      {
        std::string readToken ;
        restoreStreamObj() >> readToken;
        restoreStreamObj() >> value;
        if( token != readToken )
        {
          DUNE_THROW(InvalidStateException,"wrong object restored in PersistenceManager" << token << " " << readToken );
        }
      }

      //! clear all objects registered to PersistenceManager
      void reset()
      {
        closeStreams();   // flush and close streams
        objects_.clear(); // clear all objects
        path_.clear();
        lineNo_ = 0;
        closed_ =  false;
        invalid_ = false;
      }

    public:
      friend class Dune::Fem::Singleton< PersistenceManager >;

      static PersistenceManager &instance ()
      {
        return Singleton< PersistenceManager >::instance();
      }

      static BackupStreamType& backupStream()
      {
        return instance ().backupStreamObj();
      }

      static RestoreStreamType& restoreStream()
      {
        return instance ().restoreStreamObj();
      }

      static void insert ( PersistentObject &object, const bool pushFront = false )
      {
        instance().insertObject( object, pushFront );
      }

      static void remove ( PersistentObject &object )
      {
        instance().removeObject( object );
      }

      static void backup ( const std::string& path )
      {
        instance().backupObjects( path );
      }

      static void restore ( const std::string& path )
      {
        instance().restoreObjects( path );
      }

      static void startRestore ( const std::string& path )
      {
        instance().startRestoreImpl( path );
      }

      static std::string uniqueFileName(const std::string& tag = "" )
      {
        return instance().getUniqueFileName( tag );
      }

      static std::string uniqueTag(const std::string& tag = "" )
      {
        return instance().getUniqueTag( tag );
      }

      template< class T >
      static void backupValue ( const std::string &token, const T &value )
      {
        instance().backup( token, value );
      }

      template< class T >
      static void restoreValue ( const std::string &token, T &value )
      {
        instance().restore( token, value );
      }

    private:
      const char* myTag() const { return "persistentobjects"; }

      // create filename for persistent objects
      std::string createFilename( const std::string& path,
                                  const int rank,
                                  const int size ) const
      {
        std::stringstream s;
        const int number = ( singleBackupRestoreFile ) ? size : rank ;
        s << path << myTag() << "." << number ;
        return s.str();
      }

      void startBackup ( const std::string &path )
      {
        path_ = path + "/";

        if( createDirectory( path_ ) )
        {
          const int rank = MPIManager :: rank() ;
          const int size = MPIManager :: size() ;
          std::string filename( createFilename( path_, rank, size ) );

          assert( backupStream_ == 0 );
          backupStream_ = Fem :: StreamFactory<BackupStreamType> :: create( filename );

          if( rank == 0 )
          {
            std::ofstream paramfile( (path_ + "parameter").c_str() );
            if( paramfile )
            {
              // write parameters on rank 0
              Parameter::write( paramfile, true );
            }
          }
        }
        else
          std::cerr << "Error: Unable to create '" << path_ << "'" << std::endl;
      }

      void startRestoreImpl ( const std::string &path )
      {
        if( restoreStream_ == 0 )
        {
          path_ = path + "/";
          const int rank = MPIManager :: rank();
          const int size = MPIManager :: size();
          std::string filename( createFilename( path_, rank, size ) );
          // create strema with stream factory
          restoreStream_ = Fem :: StreamFactory<RestoreStreamType> :: create( filename );

          if( Parameter :: verbose () )
            std::cout << "Restore from " << filename << std::endl;

          if( ! restoreStream_ )
          {
            std::cout << "Error opening global stream: " << path_+myTag()
                      << std::endl;
            abort();
          }

          // restore parameter
          Parameter::container().clear();
          Parameter::append(path_ + "parameter");
        }
      }

      void closeStreams ()
      {
        if( backupStream_ )
        {
          backupStream_->flush();
          delete backupStream_;
          backupStream_ = 0;
        }

        if( restoreStream_ )
        {
          delete restoreStream_;
          restoreStream_ = 0;
        }
      }
    };


    // !!!! not accessable outside namespace Dune::Fem ?!?!?!
    namespace
    {
      PersistenceManager &persistenceManager __attribute__((unused)) = PersistenceManager::instance();
    }


    template< class ObjectType >
    inline PersistenceManager &
    operator<< ( PersistenceManager &pm, ObjectType &object )
    {
      static_assert( !std::is_pointer< ObjectType >::value, "Do not add pointers to PersistenceManager." );
      pm.insertObject( object );
      return pm;
    }


    template< class ObjectType >
    inline PersistenceManager &
    operator>> ( PersistenceManager &pm, ObjectType &object )
    {
      pm.removeObject( object );
      return pm;
    }


    template< class ObjectType >
    struct PersistenceManager::WrapObject< ObjectType, true >
    {
      static PersistentObject *apply( ObjectType &obj )
      {
        return &obj;
      }
    };


    template< class ObjectType >
    struct PersistenceManager::WrapObject< ObjectType, false >
    : public PersistentObject
    {
      typedef WrapObject< ObjectType, false > ThisType;
      typedef PersistentObject BaseType;

    protected:
      ObjectType& obj_;
      std::string token_;

      WrapObject( ObjectType &obj )
      : obj_( obj ),
        // store unique token of this object
        token_( "_token"+PersistenceManager::uniqueTag() )
      {}

    public:
      virtual ~WrapObject ()
      {}

      virtual void backup () const
      {
        PersistenceManager::backupValue( token_, obj_ );
      }

      virtual void restore ()
      {
        PersistenceManager::restoreValue( token_, obj_ );
      }

    protected:
      virtual void *pointer ()
      {
        return &obj_;
      }

    public:
      static PersistentObject *apply ( ObjectType &obj )
      {
        return new ThisType( obj );
      }
    };



    /** \class   AutoPersistentObject
     *  \ingroup Checkpointing
     *  \brief   base class for auto persistent objects
     */
    class AutoPersistentObject
    : public PersistentObject
    {
      typedef AutoPersistentObject ThisType;
      typedef PersistentObject BaseType;

    protected:
      AutoPersistentObject ()
      {
        PersistenceManager::insert( *this );
      }

      AutoPersistentObject ( const ThisType & )
      {
        PersistenceManager::insert( *this );
      }

      virtual ~AutoPersistentObject ()
      {
        PersistenceManager::remove( *this );
      }
    };

  } // end namespace Fem

} // end namespace Dune

#endif // #ifndef DUNE_FEM_PERSISTENCEMANAGER_HH
