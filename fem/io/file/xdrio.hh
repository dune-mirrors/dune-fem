#ifndef DUNE_XDRIO_HH
#define DUNE_XDRIO_HH

//- system headers  
#include <rpc/xdr.h>
#include <cassert>

namespace Dune { 

//
template <class T> struct XdrIO;
  
// xdr method for double  
template <> 
struct XdrIO<double>
{
  // read/write data to xdr stream 
  static int io(XDR * xdrs,double& value)
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
  static int io(XDR * xdrs,float& value)
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
  static int io(XDR * xdrs,int& value)
  {
    assert( xdrs );
    return xdr_int(xdrs, &value);
  }
};

} // end namespace Dune 
#endif
