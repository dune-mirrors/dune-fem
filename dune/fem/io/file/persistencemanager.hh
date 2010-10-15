#ifndef DUNE_FEM_PERSISTENCEMANAGER_HH
#define DUNE_FEM_PERSISTENCEMANAGER_HH

#include <fstream>
#include <list>

#include <dune/common/static_assert.hh>
#include <dune/common/typetraits.hh>

#include <dune/fem/io/streams/virtualstreams.hh>
#include <dune/fem/io/file/asciiparser.hh>
#include <dune/fem/io/file/iointerface.hh>
#include <dune/fem/io/parameter.hh>

namespace Dune
{
  
/** @addtogroup Checkpointing
 *  
 *  The Dune::PersistenceManager manages a list of persistent
 *  object the state should be recovered during a restart process.
 *  The singleton instance is available using persistenceManager.
 *
 *  Variables of any type can be added to the manager
 *  using operator<< and removed through operator>>.
 *  For general types the state is written to an ascii file
 *  \c persistentobjects through the call to PersistenceManager::backup(path).
 *  
 *  Finer control is available for user defined classes by 
 *  deriving these from the Dune::PersistentObject class and
 *  implementing the virtual functions \c backup and \c restore.
 *  Values can then either be written to the ascii file
 *  \c persistentobjects using the method
 *  Dune::PersistenceManager::backupValue giving a token and
 *  the value to store (a structure simular to the Dune::Parameter
 *  files is used).
 *  If data is to be written to a different stream a unique filename
 *  can be obtained from the Dune::PersistenceManager using the
 *  \c uniqueFileName method.
 *
 *  If a class is derived from Dune::AutoPersistenceObject 
 *  then instances of this class are automatically added to the
 *  PersistenceManager on construction and removed on destruction;
 *  again the virtual methods backup/restore must be provided.
 *
 *  The restart process is activated through the method \c restore.
 *  Note that care must be taken to add the objects in exactly
 *  the same order when using the PersistenceManager for creating
 *  the backup file as when used to restore the saved states.
 *
 *  The Dune::Parameter class is automatically backuped and restored.
 *  For the restore process the saved parameter file is read first and
 *  existing values for parameters are replaced by these values.
 *  Therefore it is possible to reread the state of a variable 
 *  orignally defined through a run time parameter using the Parameter
 *  class:
 *  \code 
 *  class A : public Dune::AutoPersistentObject {
 *    int var_,b_;
 *    ComplexType U_;
 *    A() : var_(Dune::Parameter::getValue("token",0,var_) {}
 *    virtual backup() const {
 *      // note var_ is saved in the parameter file
 *
 *      // write b to ascii file 
 *      PersistenceManager::backupValue("b",b);
 *      
 *      // write U to stream
 *      ofstream out((PersistenceManager::uniqueFileName()+"classA").c_str());
 *      U.write(out);
 *    }
 *    virtual restore() {
 *      // set var_ through parameter file
 *      var_= Dune::Parameter::getValue<int>("token",var_);
 *      
 *      // read b from ascii file 
 *      PersistenceManager::restoreValue("b",b);
 *      
 *      // read U from stream
 *      ifstream in((PersistenceManager::uniqueFileName()+"classA").c_str());
 *      U.read(in);
 *    }
 *  };
 *  \endcode
 *  
 **/

  class PersistenceManager;
  
  /** \class   PersistentObject
   *  \ingroup Checkpointing
   *  \brief   base class for persistent objects
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
    static const bool value = Dune::Conversion< ObjectType *, PersistentObject * >::exists;
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

  private:
    typedef std::list< std::pair< PersistentObject *, unsigned int > > PersistentType;
    typedef PersistentType::iterator IteratorType;

    PersistentType objects_;
    int fileCounter_,lineNo_;
    std::string path_;
    std::ifstream inAsciStream_;
    std::ofstream outAsciStream_;
    bool closed_,invalid_;
    
    PersistenceManager ()
    : fileCounter_( 0 ),
      lineNo_(),
      path_(), 
      closed_( false ),
      invalid_( false )
    {}

    PersistenceManager ( const ThisType & );
    ThisType &operator= ( const ThisType & );

  public:
    template< class ObjectType >
    void insertObject( ObjectType &object )
    {
      IteratorType end = objects_.end();
      for( IteratorType it = objects_.begin(); it != end; ++it )
      {
        if( it->first->pointer() != &object )
          continue;
        ++it->second;
        return;
      }
      if (closed_) {
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

      closeAscii();
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
      closed_=true;
      startRestore(path);
      typedef PersistentType :: iterator IteratorType;

      for( IteratorType it = objects_.begin(); it != objects_.end(); ++it )
        it->first->restore( );

      closeAscii();
    }

    std::string getUniqueFileName (const std::string& tag )
    { 
      return genFilename( path_, tag, ++fileCounter_ );
    }

    template< class T >
    void backup ( const std::string &token, const T &value )
    {
      outAsciStream_ << token << ": " << value << std::endl;
    }

    template< class T >
    void restore ( const std::string &token, T &value )
    {
      std::string linestring;
      int startLine=lineNo_;
      bool found;
      do { 
        std::getline(inAsciStream_,linestring);
        std::stringstream line(linestring);
        /*
        std::cout << "  " << lineNo_ 
                  << " : " << linestring 
                  << " , " << line << std::endl;
        */
        if (!(inAsciStream_)) {
          std::cout << "Error in restore process! " 
                    << "The token " << token << " was not found "
                    << "in file " << path_+myTag() << " "
                    << "after line number " << startLine
                    << std::endl;
          abort();
        }
        lineNo_++;
        found = Dune::readParameter(line,token,value,false,false);
      } while (!found);
    }

  public:
    static PersistenceManager &instance ()
    {
      static PersistenceManager theInstance;
      return theInstance;
    }

    static void insert ( PersistentObject &object )
    {
      instance().insertObject( object );
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

    static std::string uniqueFileName(const std::string& tag = "" )
    {
      return instance().getUniqueFileName( tag );
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

    void startBackup ( const std::string &path )
    {
      path_ = path + "/";
      fileCounter_ = 0;
      lineNo_ = 0;

      if( createDirectory( path_ ) )
      {
        outAsciStream_.open((path_+ myTag()).c_str());  
        outAsciStream_ << std::scientific;
        outAsciStream_.precision(16);
        outAsciStream_ << "Persistent Objects" << std::endl;
        // write parameters 
        Parameter::write( path_, "parameter", true );
      }
      else
        std::cerr << "Error: Unable to create '" << path_ << "'" << std::endl;
    }

    void startRestore ( const std::string &path )
    {
      path_=path + "/";
      fileCounter_ = 0;
      lineNo_ = 0;
      inAsciStream_.open((path_+myTag()).c_str());  
      if (!inAsciStream_) {
        std::cout << "Error opening global stream: " << path_+myTag()
                  << std::endl;
        abort();
      }
      std::string tmp;
      inAsciStream_ >> tmp;
      std::cout << tmp << " ";
      inAsciStream_ >> tmp;
      std::cout << tmp << " ";
      std::cout << std::endl;

      Parameter::clear();
      Parameter::append(path_+"parameter");
    }

    void closeAscii ()
    {
      if( outAsciStream_.is_open() )
        outAsciStream_.close();
      if( inAsciStream_.is_open() )
        inAsciStream_.close();
    }
  };
  
  namespace
  {
    PersistenceManager &persistenceManager = PersistenceManager::instance();
  }


  template< class ObjectType >
  inline PersistenceManager &
  operator<< ( PersistenceManager &pm, ObjectType &object )
  {
    dune_static_assert( !TypeTraits< ObjectType >::isPointer, "Do not add pointers to PersistenceManager." );
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

  private:
    ObjectType& obj_;
    
    WrapObject( ObjectType &obj )
    : obj_( obj )
    {}
    
  public:
    virtual ~WrapObject ()
    {}
    
    virtual void backup () const
    {
      PersistenceManager::backupValue( "_token"+PersistenceManager::uniqueFileName(), obj_ );
    }
    
    virtual void restore ()
    {
      PersistenceManager::restoreValue( "_token"+PersistenceManager::uniqueFileName(), obj_ );
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

}

#endif // #ifndef DUNE_FEM_PERSISTENCEMANAGER_HH
