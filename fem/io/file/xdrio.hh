#ifndef DUNE_XDRIO_HH
#define DUNE_XDRIO_HH

//- system headers  
#include <iostream>
#include <rpc/xdr.h>
#include <cassert>

//- Dune includes 
#include <dune/common/fvector.hh>

namespace Dune { 

//! wrapper class for xdr calls specified by data type
template <class T> struct XdrIO;
  
// xdr method for double  
template <> 
struct XdrIO<double>
{
  // read/write data to xdr stream 
  static inline int io(XDR * xdrs,double& value)
  {
    assert( xdrs );
    return xdr_double(xdrs, &value);
  }
};

// xdr method for float  
template <> 
struct XdrIO<float>
{
  // read/write data to xdr stream 
  static inline int io(XDR * xdrs,float& value)
  {
    assert( xdrs );
    return xdr_float(xdrs, &value);
  }
};

// xdr method for int   
template <> 
struct XdrIO<int>
{
  // read/write data to xdr stream 
  static inline int io(XDR * xdrs,int& value)
  {
    assert( xdrs );
    return xdr_int(xdrs, &value);
  }
};

// xdr method for int   
template <> 
struct XdrIO<bool>
{
  // read/write data to xdr stream 
  static inline int io(XDR * xdrs,bool& value)
  {
    assert( xdrs );
    // convert to char  
    char val = (value == true) ? 1 : 0;
    int ret = xdr_char(xdrs, &val);
    // convert back 
    value = (val == 1) ? true : false;
    return ret;
  }
};

// xdr method for int   
template <> 
struct XdrIO<long>
{
  // read/write data to xdr stream 
  static inline int io(XDR * xdrs,long& value)
  {
    assert( xdrs );
    return xdr_long(xdrs, &value);
  }
};

// xdr method for char   
template <> 
struct XdrIO<char>
{
  // read/write data to xdr stream 
  static inline int io(XDR * xdrs,char& value)
  {
    assert( xdrs );
    return xdr_char(xdrs, &value);
  }
};

// xdr method for FieldVector   
template <class T, int n> 
struct XdrIO<FieldVector<T,n> >
{
  typedef FieldVector<T,n> VectorType;
  // read/write data to xdr stream 
  static inline int io(XDR * xdrs,VectorType& vec)
  {
    assert( xdrs );
    int result = 1;
    for(int i=0; i<n; ++i)
    {
      result |= XdrIO<T>::io(xdrs,vec[i]);
    }
    return result;
  }
};

//! xdr stream class 
class XDRStream 
{
protected:  
  // file handler 
  FILE*   file_;
  // xdr struct 
  XDR     xdrs_; 

  const bool in_;
private:
  // do not copy
  XDRStream(const  XDRStream&);
public:
  //! constructor 
  XDRStream(bool in) : file_(0) , in_(in) {}
    
  //! destructor 
  ~XDRStream() 
  {
    xdr_destroy(xdrs());
    fclose(file_);
  }

  //! return true if stream is valid
  bool operator ! () const 
  {
    return file_ != 0;
  }
  
  //! write value to stream 
  template <class T> 
  XDRStream& operator << (const T& value)
  {
    assert( !in_ );
    XdrIO<T>::io(xdrs(),value);
    return *this;
  }
  
  //! read value from stream 
  template <class T> 
  XDRStream& operator >> (T& value)  
  {
    assert( in_ );
    XdrIO<T>::io(xdrs(),value);
    return *this;
  }

  //! read/write value for stream 
  template <class T> 
  int inout(T& value)
  {
    return XdrIO<T>::io(xdrs(),value);
  }

protected:  
  //! return pointer to xdr 
  XDR* xdrs() { return &xdrs_; }
};

//! xdr stream for writing  
class XDRWriteStream : public XDRStream  
{
  typedef XDRStream BaseType;
  // do not copy
  XDRWriteStream(const  XDRWriteStream&);
public:  
  XDRWriteStream(const std::string& filename) :
    BaseType(false)
  {
    this->file_ = fopen(filename.c_str(),"wb");
    if(! this->file_)
    { 
      std::cerr <<"\aXDRWriteStream:: couldn't open file `"<<filename<<"' !\n";
      std::cerr.flush();
      return ;
    } 
    xdrstdio_create(this->xdrs(), this->file_, XDR_ENCODE);
  }
};

//! xdr stream for reading 
class XDRReadStream :  public XDRStream 
{
  typedef XDRStream BaseType;
  
  // do not copy
  XDRReadStream(const  XDRReadStream&);
public:  
  XDRReadStream(const std::string& filename) 
    : BaseType(true)
  {
    this->file_ = fopen(filename.c_str(),"rb");
    if(! this->file_)
    { 
      std::cerr <<"\aXDRReadStream:: couldn't open file `"<<filename<<"' !\n";
      std::cerr.flush();
      return ;
    } 
    xdrstdio_create(this->xdrs(), this->file_, XDR_DECODE);
  }
};

} // end namespace Dune 
#endif
