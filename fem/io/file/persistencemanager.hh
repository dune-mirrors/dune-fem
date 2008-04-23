#ifndef DUNE_FEM_PERSISTENCE_HH
#define DUNE_FEM_PERSISTENCE_HH

#include <map>
#include <fstream>

#include <dune/fem/io/parameter.hh>
#include <dune/fem/io/streams/virtualstreams.hh>
#include <dune/fem/io/file/asciiparser.hh>

namespace Dune
{

  class PersistenceManager;
  
  class PersistentObject 
  {
    typedef PersistentObject ThisType;

  public:
    virtual void backup ( ) const {}
    virtual void restore ( ) {}
  };
  template <class ObjectType,bool isPersistent>
  struct WrappObject;  
  template <class ObjectType> 
  struct WrappObject<ObjectType,true>
  {
    static PersistentObject* apply(ObjectType& obj) {
      return &obj;
    }
  };
  template <class ObjectType> 
  struct WrappObject<ObjectType,false> 
  {
    struct ObjectWrapper : public PersistentObject {
      ObjectWrapper(ObjectType& obj) :
        obj_(obj) {}
      virtual void backup ( ) const;
      virtual void restore ( );
      ObjectType& obj_;
    };
    static PersistentObject* apply(ObjectType& obj) {
      return new ObjectWrapper(obj);
    }
  };
  class PersistenceManager
  {
    typedef PersistenceManager ThisType;

  private:
    typedef std :: vector< std::pair<PersistentObject *, unsigned int > > PersistentType;
    typedef PersistentType :: iterator IteratorType;

  private:
    PersistentType objects_;
    
  private:
    inline PersistenceManager () :
      fileCounter_(0), lineNo_(), path_()
    {}

    PersistenceManager ( const ThisType & );
    ThisType &operator= ( const ThisType & );

  public:
    template <class ObjectType>
    inline void insertObject( ObjectType& object) {
      PersistentObject* obj = 
        WrappObject<ObjectType,
                    Dune::Conversion<ObjectType,PersistentObject>::exists>
        ::apply(object);
      for (int i=0;i<objects_.size();i++) {
        if (objects_[i].first==obj) {
          objects_[i].second++;
          return;
        }
      }
      objects_.push_back(std::make_pair(obj,1));
    }

    inline void removeObject ( PersistentObject &object )
    {
      /*
      typedef PersistentMapType :: iterator IteratorType;
      
      IteratorType it = objects_.find( &object );
      if( it != objects_.end() )
      {
        if( (--it->second) == 0 )
          erase( it );
      }
      */
    }

    inline void backupObjects ( const std::string& path ) 
    {
      startBackup(path);
      typedef PersistentType :: iterator IteratorType;

      for( IteratorType it = objects_.begin(); it != objects_.end(); ++it )
        it->first->backup( );

      closeAsci();
    }
    
    inline void restoreObjects ( const std::string& path) 
    {
      startRestore(path);
      typedef PersistentType :: iterator IteratorType;

      for( IteratorType it = objects_.begin(); it != objects_.end(); ++it )
        it->first->restore( );

      closeAsci();
    }

    inline std::string getUniqueFileName() { 
      fileCounter_++;
      return genFilename(path_,"",fileCounter_);
    }
    template <class T>
    inline void backup(const std::string& token,const T& value) {
      outAsciStream_ << token << ": " << value << std::endl;
    }
    template <class T>
    inline void restore(const std::string& token,T& value) {
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
                    << "in file " << path_+"checkpoint" << " "
                    << "after line number " << startLine
                    << std::endl;
          abort();
        }
        lineNo_++;
        found = Dune::readParameter(line,token,value,false,false);
      } while (!found);
    }

  public:
    inline static PersistenceManager &instance ()
    {
      static PersistenceManager theInstance;
      return theInstance;
    }

    inline static void insert ( PersistentObject &object )
    {
      instance().insertObject( object );
    }

    inline static void remove ( PersistentObject &object )
    {
      instance().removeObject( object );
    }

    inline static void backup ( const std::string& path )
    {
      instance().backupObjects( path );
    }

    inline static void restore ( const std::string& path )
    {
      instance().restoreObjects( path );
    }

    inline static std::string uniqueFileName() { 
      return instance().getUniqueFileName();
    }
    template <class T>
    inline static void backupValue(const std::string& token,const T& value) {
      instance().backup(token,value);
    }
    template <class T>
    inline static void restoreValue(const std::string& token,T& value) {
      instance().restore(token,value);
    }
  private:
    void startBackup(const std::string& path) {
      path_=path+"/";
      fileCounter_=0;
      lineNo_=0;
      outAsciStream_.open((path_+"checkpoint").c_str());  
      outAsciStream_ << "Persistent Objects" << std::endl;
    }
    void startRestore(const std::string& path) {
      path_=path+"/";
      fileCounter_=0;
      lineNo_=0;
      inAsciStream_.open((path_+"checkpoint").c_str());  
      if (!inAsciStream_) {
        std::cout << "Error opening global checkpoint stream!"
                  << std::endl;
        abort();
      }
      std::string tmp;
      inAsciStream_ >> tmp;
      std::cout << tmp << " ";
      inAsciStream_ >> tmp;
      std::cout << tmp << " ";
      std::cout << std::endl;
    }
    void closeAsci() {
      if (outAsciStream_.is_open())
        outAsciStream_.close();
      if (inAsciStream_.is_open())
        inAsciStream_.close();
    }
    int fileCounter_,lineNo_;
    std::string path_;
    std::ifstream inAsciStream_;
    std::ofstream outAsciStream_;
  };
  
  namespace {
  PersistenceManager& persistenceManager = 
    PersistenceManager::instance();
  }

  template <class ObjectType>
  inline PersistenceManager &operator<< ( PersistenceManager &pm,
                                          ObjectType &object )
  {
    pm.insertObject( object );
    return pm;
  }

  inline PersistenceManager &operator>> ( PersistenceManager &pm,
                                          PersistentObject &object )
  {
    pm.removeObject( object );
    return pm;
  }


  class AutoPersistentObject
  : public PersistentObject
  {
    typedef AutoPersistentObject ThisType;
    typedef PersistentObject BaseType;

  protected:
    inline AutoPersistentObject ()
    {
      PersistenceManager :: insert( *this );
    }

    inline AutoPersistentObject ( const ThisType & )
    {
      PersistenceManager :: insert( *this );
    }

    inline ~AutoPersistentObject ()
    {
      PersistenceManager :: remove( *this );
    }
  };

  template <class ObjectType> 
  void WrappObject<ObjectType,false> ::
       ObjectWrapper :: backup ( ) const {
    PersistenceManager::backupValue(
        "_token"+PersistenceManager::uniqueFileName(),obj_
    );
  }
  template <class ObjectType> 
  void WrappObject<ObjectType,false> ::
       ObjectWrapper :: restore ( ) {
    PersistenceManager::restoreValue(
        "_token"+PersistenceManager::uniqueFileName(),obj_
    );
  }

}

#endif
