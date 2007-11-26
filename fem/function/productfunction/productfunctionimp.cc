#ifndef DUNE_PRODUCTFUNCTION_CC
#define DUNE_PRODUCTFUNCTION_CC



#include <algorithm>
// BUG: return values are initialized with *= 0. Problematic for nan
// values!
namespace Dune {

// Constructor making discrete function  
template<class DiscreteFunctionSpaceType, class DiscreteFunctionSpace2Type>
  inline ProductDiscreteFunction< DiscreteFunctionSpaceType, DiscreteFunctionSpace2Type>::
  ProductDiscreteFunction(const DiscreteFunctionSpaceType& f, const DiscreteFunctionSpace2Type& f2) :
    // DiscreteFunctionDefaultType ( f ), 
    name_ ("no name"),
    dm_(DofManagerFactoryType::getDofManager(f.grid())),
    memPair_(dm_.addDofSet(&dofVec_, f.mapper(), name_)),
    dofVec_ ( *memPair_.second ),
    functionSpace_(f),
    functionSpace2_ (f2)
{
    std::cout << "Size old = " << dofVec_.size() << "  Size  f2 = " << functionSpace2_.size();
    std::cout << " New size will be " << dofVec_.size() * functionSpace2_.size() << "\n";
    dofVec_.resize ( dofVec_.size() * functionSpace2_.size() );
}
// Desctructor 
template<class DiscreteFunctionSpaceType, class DiscreteFunctionSpace2Type>
inline ProductDiscreteFunction< DiscreteFunctionSpaceType, DiscreteFunctionSpace2Type>::
~ProductDiscreteFunction()
{
  dm_.removeDofSet(*memPair_.first);
}

template<class DiscreteFunctionSpaceType, class DiscreteFunctionSpace2Type>
inline void ProductDiscreteFunction< DiscreteFunctionSpaceType, DiscreteFunctionSpace2Type>::clear ()
{
  const int size = dofVec_.size();
  DofStorageType &vec = dofVec_;
  for(int i=0; i<size; ++i) vec[i] = 0.0; 
  //set ( 0.0 ); 
}

template<class DiscreteFunctionSpaceType, class DiscreteFunctionSpace2Type>
inline void ProductDiscreteFunction< DiscreteFunctionSpaceType, DiscreteFunctionSpace2Type>::
print(std::ostream &s ) const
{
  RangeFieldType sum = 0.;
  ConstDofIteratorType enddof = dend ( );
  for(ConstDofIteratorType itdof = dbegin ( ); itdof != enddof; ++itdof) 
  {
    s << (*itdof) << " DofValue \n";
    sum += std::abs(*itdof);
  } 
  s << "sum = " << sum << "\n";
}
//*************************************************************************
//  Interface Methods 
//*************************************************************************
//return local function for global dof index dofIndex2 of space two
template<class DiscreteFunctionSpaceType, class DiscreteFunctionSpace2Type> 
inline AdaptiveDiscreteFunction< DiscreteFunctionSpaceType>
ProductDiscreteFunction<DiscreteFunctionSpaceType, DiscreteFunctionSpace2Type>::localFunction(int dofIndex2) const
{
 return AdaptiveDiscreteFunction< DiscreteFunctionSpaceType> ( name(), this->functionSpace_, &(dofVec_.leakPointer()[this->functionSpace_.size()*dofIndex2]));
 
}

//return local function for local dof index dofIndex2 of space two and given entity2
template<class DiscreteFunctionSpaceType, class DiscreteFunctionSpace2Type> 
template < class Entity2Type>
inline AdaptiveDiscreteFunction< DiscreteFunctionSpaceType>
ProductDiscreteFunction<DiscreteFunctionSpaceType, DiscreteFunctionSpace2Type>::localFunction(const Entity2Type &en2,int dofIndex2) const
{
 return AdaptiveDiscreteFunction< DiscreteFunctionSpaceType> ( name(), this->functionSpace_, &(dofVec_.leakPointer()[this->functionSpace_.size()*space2().mapToGlobal( en2 ,dofIndex2)]));
 
}
#if 0
//local function for given entity2 and local coordinate of second space and discrete function of first space
template<class DiscreteFunctionSpaceType, class DiscreteFunctionSpace2Type>
template < class Entity2Type, class LocalCoord2Type>
inline void ProductDiscreteFunction<DiscreteFunctionSpaceType, DiscreteFunctionSpace2Type>::
localFunction(const Entity2Type &en2, const LocalCoord2Type &loc2,  DiscreteFunction1Type &discFunc) const
{
 typedef typename DiscreteFunctionSpace2Type::BaseFunctionSetType BaseFunctionSetType;
 typename DiscreteFunctionSpace2Type::RangeType tmp_;
 discFunc.clear();
 const BaseFunctionSetType bSet2 = space2().baseFunctionSet(en2);
 const int numOfDofs = bSet2.numBaseFunctions();
 
 for(int i = 0; i< numOfDofs; i++) 
 {
  const int map = space2().mapToGlobal( en2 , i );  
  bSet2.evaluate(i,loc2,tmp_);
  DiscreteFunction1Type df = this->localFunction(map);
	    discFunc.addScaled(df,tmp_);
 }
	
}
#endif

//local function for given entity2, quadrature type and qudrature point number of second space and discrete function of first space
template<class DiscreteFunctionSpaceType, class DiscreteFunctionSpace2Type>
template < class Entity2Type, class QuadratureType>
inline void ProductDiscreteFunction<DiscreteFunctionSpaceType, DiscreteFunctionSpace2Type>::localFunction(const Entity2Type &en2, const QuadratureType &quad2,  int pointNr, DiscreteFunction1Type &discFunc) const
{
 typedef typename DiscreteFunctionSpace2Type::BaseFunctionSetType BaseFunctionSetType;
 typename DiscreteFunctionSpace2Type::RangeType tmp_;
 discFunc.clear();
 const BaseFunctionSetType bSet2 = space2().baseFunctionSet(en2);
 const int numOfDofs = bSet2.numBaseFunctions();
 
 for(int i = 0; i< numOfDofs; i++) 
 {
  const int map = space2().mapToGlobal( en2 , i );  
  bSet2.evaluate(i,quad2[pointNr],tmp_);
  DiscreteFunction1Type df = this->localFunction(map);
	    discFunc.addScaled(df,tmp_);
 }
	
}

//local function for given entity2, quadrature type and qudrature point number of second space and discrete function of first space
template<class DiscreteFunctionSpaceType, class DiscreteFunctionSpace2Type>
template < class Entity2Type, class PointType>
inline void ProductDiscreteFunction<DiscreteFunctionSpaceType, DiscreteFunctionSpace2Type>::
localFunction(const Entity2Type &en2, const PointType &pt, DiscreteFunction1Type &discFunc) const
{
 typedef typename DiscreteFunctionSpace2Type::BaseFunctionSetType BaseFunctionSetType;
 typename DiscreteFunctionSpace2Type::RangeType tmp_;
 discFunc.clear();
 const BaseFunctionSetType bSet2 = space2().baseFunctionSet(en2);
 const int numOfDofs = bSet2.numBaseFunctions();
 
 for(int i = 0; i< numOfDofs; i++) 
 {
  const int map = space2().mapToGlobal( en2 , i );  
  bSet2.evaluate(pt,tmp_);
  DiscreteFunction1Type df = this->localFunction(map);
	    discFunc.addScaled(df,tmp_);
 }
}

template<class DiscreteFunctionSpaceType, class DiscreteFunctionSpace2Type> 
inline typename ProductDiscreteFunction< DiscreteFunctionSpaceType, DiscreteFunctionSpace2Type>::DofIteratorType 
ProductDiscreteFunction< DiscreteFunctionSpaceType, DiscreteFunctionSpace2Type>::dbegin ( )
{
  return dofVec_.begin();
}

template<class DiscreteFunctionSpaceType, class DiscreteFunctionSpace2Type> 
inline typename ProductDiscreteFunction< DiscreteFunctionSpaceType, DiscreteFunctionSpace2Type>::DofIteratorType 
ProductDiscreteFunction< DiscreteFunctionSpaceType, DiscreteFunctionSpace2Type>::dend ()
{
  return dofVec_.end();
}

template<class DiscreteFunctionSpaceType, class DiscreteFunctionSpace2Type> 
inline typename ProductDiscreteFunction< DiscreteFunctionSpaceType, DiscreteFunctionSpace2Type>::ConstDofIteratorType 
ProductDiscreteFunction< DiscreteFunctionSpaceType, DiscreteFunctionSpace2Type>::dbegin ( ) const
{
  return dofVec_.begin();
}

template<class DiscreteFunctionSpaceType, class DiscreteFunctionSpace2Type> 
inline typename ProductDiscreteFunction< DiscreteFunctionSpaceType, DiscreteFunctionSpace2Type>::ConstDofIteratorType 
ProductDiscreteFunction< DiscreteFunctionSpaceType, DiscreteFunctionSpace2Type>::dend () const 
{
  return dofVec_.end();
}
//**************************************************************************
//  Read and Write Methods 
//**************************************************************************
template<class DiscreteFunctionSpaceType, class DiscreteFunctionSpace2Type>
inline bool ProductDiscreteFunction< DiscreteFunctionSpaceType, DiscreteFunctionSpace2Type>::
write_xdr(const std::string fn) const
{
  // create write stream 
  XDRWriteStream xdr(fn);
  // make sure data is only written in compressed state. 
  /*
  if( dofVec_.size() != spc_.size() )
  {
    DUNE_THROW(InvalidStateException,"DofVector not in compressed state while writing!");
  }
  */

  // write data 
  return dofVec_.processXdr(xdr);
}

template<class DiscreteFunctionSpaceType, class DiscreteFunctionSpace2Type>
inline bool ProductDiscreteFunction< DiscreteFunctionSpaceType, DiscreteFunctionSpace2Type>::
read_xdr(const std::string fn)
{
  XDRReadStream xdr(fn);
  // make sure data is only written in compressed state. 
  /* 
  if( dofVec_.size() != spc_.size() )
  {
    DUNE_THROW(InvalidStateException,"DofVector size does not match while reading!");
  }
  */

  // read data 
  return dofVec_.processXdr(xdr);
}

template<class DiscreteFunctionSpaceType, class DiscreteFunctionSpace2Type>
inline bool ProductDiscreteFunction< DiscreteFunctionSpaceType, DiscreteFunctionSpace2Type>::
write_ascii(const std::string fn) const
{
  std::fstream outfile( fn.c_str() , std::ios::out );
  if (!outfile)
  { 
    printf( "\aERROR in ProductDiscreteFunction::write_ascii(..): couldnot open <%s>!\n", fn.c_str());
    fflush(stderr);
    return false;
  }

  {
    int length = this->functionSpace_.size();
    outfile << length << "\n";
    ConstDofIteratorType enddof = dend ( );
    for(ConstDofIteratorType itdof = dbegin ( );itdof != enddof; ++itdof) 
    {
      outfile << (*itdof) << " ";
    }
    outfile << "\n";
  }

  outfile.close();
  return true;
}


template<class DiscreteFunctionSpaceType, class DiscreteFunctionSpace2Type>
inline bool ProductDiscreteFunction< DiscreteFunctionSpaceType, DiscreteFunctionSpace2Type>::
read_ascii(const std::string fn)
{
  FILE *infile=0;
  infile = fopen( fn.c_str(), "r" );
  assert(infile != 0); 
  {
    int length;
    fscanf(infile,"%d \n",&length); 
    assert(length == this->functionSpace_.size( )); 

    DofIteratorType enddof = dend ( );
    for(DofIteratorType itdof = dbegin ( );itdof != enddof; ++itdof) 
    {
      fscanf(infile,"%le \n",& (*itdof)); 
    }
  }
  fclose(infile);
  return true;
}

template<class DiscreteFunctionSpaceType, class DiscreteFunctionSpace2Type>
inline bool ProductDiscreteFunction< DiscreteFunctionSpaceType, DiscreteFunctionSpace2Type>::
write_pgm(const std::string fn) const
{
  std::ofstream out( fn.c_str() );

  enum { dim = GridType::dimension };
 
  if(out)
  {
    int danz = 129; 
  
    out << "P2\n " << danz << " " << danz <<"\n255\n";
    ConstDofIteratorType enddof = dend ();
    for(ConstDofIteratorType itdof = dbegin (); itdof != enddof; ++itdof) {
      out << (int)((*itdof)*255.) << "\n";
    }
    out.close();
  }
  else 
  {
    std::cerr << "Couldn't open file '"<<fn<<"' \n";
  }
  return true;
}

template<class DiscreteFunctionSpaceType, class DiscreteFunctionSpace2Type>
inline bool ProductDiscreteFunction< DiscreteFunctionSpaceType, DiscreteFunctionSpace2Type>::
read_pgm(const std::string fn)
{
  FILE *in;
  int v;

  in = fopen( fn.c_str(), "r" );
  assert(in);

  fscanf( in, "P2\n%d %d\n%d\n", &v, &v, &v );
  DofIteratorType enddof = dend ();
  for(DofIteratorType itdof = dbegin (); itdof != enddof; ++itdof) {
    fscanf( in, "%d", &v );
    (*itdof) = ((double)v)/255.;
  } 
  fclose( in );
  return true;
}

template<class DiscreteFunctionSpaceType, class DiscreteFunctionSpace2Type>
inline void ProductDiscreteFunction< DiscreteFunctionSpaceType, DiscreteFunctionSpace2Type>::
addScaled( const ThisType &org, const RangeFieldType &c )
{
  dofVec_.axpy(org.dofVec_ , c);
}


} // end namespace 
#endif
