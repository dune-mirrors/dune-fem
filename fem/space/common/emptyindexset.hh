#ifndef DUNE_EMPTYINDEXSET_HH
#define DUNE_EMPTYINDEXSET_HH

//- system includes 
#include <iostream>
#include <string> 
#include <rpc/xdr.h>
#include <cassert>

#include <dune/grid/common/defaultindexsets.hh>

/** @file
 @brief Provides default empty index set class for persistent index sets. 
*/

namespace Dune {

/*!
  The EmptyIndexSet implements all additional method of a DUNE fem index set with 
  an empty default implementation. 
*/
class EmptyIndexSet
{
  // dummy value 
  enum { myType = -1 };
public:  
  //! return false mean the no memory has to be allocated 
  //! and no compress of date has to be done 
  bool compress () { 
    return false; 
  }

  //! return true if the index set is consecutive 
  inline bool consecutive() const { return false; }

  //! return true if the index set is persistent 
  inline bool persistent() const { return false; }

  //! do nothing here, because fathers index should already exist 
  template <class EntityType> 
  inline void insertEntity( const EntityType & en ) {}

  //! do nothing here, because fathers index should already exist 
  template <class EntityType> 
  inline void removeEntity( const EntityType & en ) {}

  //! nothing to do here 
  inline void resize () {}

  //! no extra memory for restriction is needed
  inline int additionalSizeEstimate () const { return 0; }

  static int type() { return myType; }

  //! we have no old size 
  inline int numberOfHoles ( const int codim ) const { return 0; }
  
  //! return old index, for dof manager only 
  inline int oldIndex (const int hole, const int codim ) const { return 0; }
  
  //! return new index, for dof manager only 
  inline int newIndex (const int hole, const int codim ) const { return 0; }
  
  //! write index set to xdr file 
  inline bool write_xdr(const std::string filename , int timestep) 
  {
    FILE  *file;
    XDR   xdrs;
    const char *path = "";

    std::string fnstr  = genFilename(path,filename, timestep);
    const char * fn = fnstr.c_str();
    file = fopen(fn, "wb");
    if (!file)
    {
        std::cerr << "\aERROR in DefaultGridIndexSet::write_xdr(..): could not open <"
                  << filename << ">!" << std::endl;
      return false;
    }

    xdrstdio_create(&xdrs, file, XDR_ENCODE);
    this->processXdr(&xdrs);

    xdr_destroy(&xdrs);
    fclose(file);

    return true;
  }
  
  //! read index set to xdr file 
  inline bool read_xdr(const std::string filename, int timestep) 
  {
    FILE   *file;
    XDR     xdrs;
    const char *path = "";

    std::string fnstr = genFilename(path,filename, timestep);
    const char * fn  = fnstr.c_str();
    std::cout << "Reading <" << fn << "> \n";
    file = fopen(fn, "rb");
    if(!file)
    {
      std::cerr << "\aERROR in DefaultGridIndexSet::read_xdr(..): could not open <" 
                << filename << ">!" << std::endl;
      return(false);
    }

    // read xdr 
    xdrstdio_create(&xdrs, file, XDR_DECODE);
    this->processXdr(&xdrs);

    xdr_destroy(&xdrs);
    fclose(file);
    return true;
  }

protected: 
  // read/write from/to xdr stream 
  inline bool processXdr(XDR *xdrs)
  {
    int type = myType; 
    xdr_int ( xdrs, &type);
    if(type != myType)
    {
      std::cerr << "\nERROR: DefaultGridIndex: wrong type choosen! \n\n";
      assert(type == myType);
    }
    return true;
  }
};

} // end namespace Dune 
#endif
