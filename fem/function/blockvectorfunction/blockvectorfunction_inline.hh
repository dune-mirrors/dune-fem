#ifndef DUNE_BLOCKVECTORFUNCTION_INLINE_HH
#define DUNE_BLOCKVECTORFUNCTION_INLINE_HH

#include "blockvectorfunction.hh"

namespace Dune 
{

  template< class DiscreteFunctionSpaceImp >
  inline BlockVectorDiscreteFunction< DiscreteFunctionSpaceImp >
    :: BlockVectorDiscreteFunction ( const DiscreteFunctionSpaceType &f ) 
  : BaseType( "no name ", f, lfFactory_ ), 
    lfFactory_( *this ),
    mapper_( f.blockMapper() ) ,
    dm_(DofManagerFactoryType::getDofManager(f.grid())),
    memPair_(dm_.addDofSet(&dofVec_, mapper_, BaseType :: name())),
    dofVec_( *memPair_.second ),
    leakPtr_(dofVec_),
    localFunc_( *this )
  {}


  template< class DiscreteFunctionSpaceImp >
  inline BlockVectorDiscreteFunction< DiscreteFunctionSpaceImp >
    :: BlockVectorDiscreteFunction ( const std :: string &name,
                                     const DiscreteFunctionSpaceType &f )
  : BaseType( name, f, lfFactory_ ),
    lfFactory_( *this ),
    mapper_( f.blockMapper() ),
    dm_(DofManagerFactoryType::getDofManager(f.grid())),
    memPair_(dm_.addDofSet(&dofVec_, mapper_, name)),
    dofVec_( *memPair_.second ),
    leakPtr_(dofVec_),
    localFunc_( *this )
  {}

  
  template< class DiscreteFunctionSpaceImp >
  inline BlockVectorDiscreteFunction< DiscreteFunctionSpaceImp >
    :: BlockVectorDiscreteFunction ( const std :: string &name,
                                     const DiscreteFunctionSpaceType &f,
                                     const DofStorageType &data )
  : BaseType( name, f , lfFactory_ ),
    lfFactory_( *this ),
    mapper_( f.blockMapper() ),
    dm_(DofManagerFactoryType::getDofManager(f.grid())),
    memPair_(dm_.addDummyDofSet(&dofVec_, mapper_, name, &data)),
    dofVec_( *memPair_.second ),
    leakPtr_(dofVec_),
    localFunc_( *this )
  {}

  
  template< class DiscreteFunctionSpaceImp >
  inline BlockVectorDiscreteFunction< DiscreteFunctionSpaceImp >
    :: BlockVectorDiscreteFunction ( const ThisType &other ) 
  : BaseType( other.name(), other.functionSpace_, lfFactory_ ),
    lfFactory_( *this ),
    mapper_( other.functionSpace_.blockMapper() ),
    dm_( other.dm_ ),
    memPair_( dm_.addDofSet( &dofVec_, mapper_, BaseType :: name() ) ),
    dofVec_( *memPair_.second ),
    leakPtr_( dofVec_ ),
    localFunc_( *this )
  {
    built_ = other.built_; 
    dofVec_ = other.dofVec_;
  } 

  
  template< class DiscreteFunctionSpaceImp >
  inline BlockVectorDiscreteFunction< DiscreteFunctionSpaceImp >
    :: ~BlockVectorDiscreteFunction ()
  {
    dm_.removeDofSet( *memPair_.first );
  }


template<class DiscreteFunctionSpaceType>
inline void BlockVectorDiscreteFunction<DiscreteFunctionSpaceType>::clear ()
{
  const int size = dofVec_.size();
  for(int i=0; i<size; ++i) dofVec_[i] = 0.0; 
}

template< class DiscreteFunctionSpaceType >
inline void BlockVectorDiscreteFunction< DiscreteFunctionSpaceType >
  :: print ( std::ostream &out ) const
{
  out << "BlockVectorDiscreteFunction '" << name() << "'" << std :: endl;
  
  const ConstDofIteratorType end = dend();
  for( ConstDofIteratorType dit = dbegin(); dit != end; ++dit )
    out << (*dit) << std :: endl;
}

//*************************************************************************
//  Interface Methods 
//*************************************************************************

#if 0
template<class DiscreteFunctionSpaceType> template <class EntityType>
inline typename BlockVectorDiscreteFunction<DiscreteFunctionSpaceType>:: LocalFunctionType 
BlockVectorDiscreteFunction<DiscreteFunctionSpaceType>:: localFunction ( const EntityType &en ) const
{
  return LocalFunctionType (en,*this);
}
#endif

template<class DiscreteFunctionSpaceType> 
inline typename BlockVectorDiscreteFunction<DiscreteFunctionSpaceType>::DofIteratorType 
BlockVectorDiscreteFunction<DiscreteFunctionSpaceType>::dbegin ()
{
  return DofIteratorType( dofVec_ , 0 );     
}


template<class DiscreteFunctionSpaceType> 
inline typename BlockVectorDiscreteFunction<DiscreteFunctionSpaceType>::DofIteratorType 
BlockVectorDiscreteFunction<DiscreteFunctionSpaceType>::dend ( )
{
  return DofIteratorType( dofVec_  , dofVec_.size() );     
}

template<class DiscreteFunctionSpaceType> 
inline typename BlockVectorDiscreteFunction<DiscreteFunctionSpaceType>::ConstDofIteratorType 
BlockVectorDiscreteFunction<DiscreteFunctionSpaceType>::dbegin ( ) const
{
  DofIteratorType tmp ( dofVec_ , 0 );     
  return ConstDofIteratorType(tmp);
}

template<class DiscreteFunctionSpaceType> 
inline typename BlockVectorDiscreteFunction<DiscreteFunctionSpaceType>::ConstDofIteratorType 
BlockVectorDiscreteFunction<DiscreteFunctionSpaceType>::dend ( ) const 
{
  DofIteratorType tmp ( dofVec_ , dofVec_.size() );     
  return ConstDofIteratorType(tmp);
}
//**************************************************************************
//  Read and Write Methods 
//**************************************************************************
template<class DiscreteFunctionSpaceType>
inline bool BlockVectorDiscreteFunction<DiscreteFunctionSpaceType>::
processXdrs(XDRStream& xdr) const  
{
  int len = dofVec_.size();
  xdr.inout( len );
 
  // make sure data is only read in compressed state. 
  if( (dofVec_.size() != mapper_.size()) || (len != dofVec_.size()) )
  {
    DUNE_THROW(InvalidStateException,"BlockVectorDiscreteFunction::processXdrs: sizes do not match!");
  }
  
  // write/read vector 
  const int vecSize = dofVec_.size();
  for(int i=0; i<vecSize; ++i) 
  {
#if OLD_XDR_METHOD
    xdr.bytes( dofVec_[i] , DofBlockType :: dimension );
#else 
    xdr.inout( dofVec_[i] );
#endif
  }
  return true;
}

#if DUNE_FEM_COMPATIBILITY
template<class DiscreteFunctionSpaceType>
inline bool BlockVectorDiscreteFunction<DiscreteFunctionSpaceType>::
write_xdr( const std::string filename ) const
{
  // create write stream 
  XDRWriteStream xdr(filename);
  return processXdrs(xdr);
}

template<class DiscreteFunctionSpaceType>
inline bool BlockVectorDiscreteFunction<DiscreteFunctionSpaceType>::
read_xdr( const std::string filename)
{
  // create read stream 
  XDRReadStream xdr(filename);
  return processXdrs(xdr);
}

template<class DiscreteFunctionSpaceType>
inline bool BlockVectorDiscreteFunction<DiscreteFunctionSpaceType>::
write_ascii( const std::string filename ) const
{
  const char * fn = filename.c_str();
  std::fstream outfile( fn , std::ios::out );
  if(outfile)
  {
    int length = this->functionSpace_.size();
    outfile << length << std::endl;
    ConstDofIteratorType enddof = this->dend ();
    for(ConstDofIteratorType itdof = this->dbegin ();itdof != enddof; ++itdof) 
    {
      outfile << (*itdof) << std::endl;
    }
    outfile.close();
  }
  else 
  { 
    fprintf(stderr,"\aERROR in BlockVectorDiscreteFunction::read_xdr(..): couldnot open <%s>!\n", fn);
    fflush(stderr);
    return(false);
  }
  return true;
}


template<class DiscreteFunctionSpaceType>
inline bool BlockVectorDiscreteFunction<DiscreteFunctionSpaceType>::
read_ascii( const std::string filename )
{
  std::ifstream infile ( filename.c_str() );
  if( infile )
  {
    int length;
    infile >> length; 
    if(length != this->functionSpace_.size()) 
    {
      DUNE_THROW(InvalidStateException,"ERROR: wrong number of dofs stored in file!"); 
    }

    DofIteratorType enddof = this->dend ();
    for(DofIteratorType itdof = this->dbegin ();itdof != enddof; ++itdof) 
    {
      infile >> (*itdof); 
    }
    return true;
  }
  return false;
}
#endif

template<class DiscreteFunctionSpaceType>
inline bool BlockVectorDiscreteFunction<DiscreteFunctionSpaceType>::
write_pgm( const std::string filename ) const
{
  const char * fn = filename.c_str();
  std::ofstream out( fn );

  enum { dim = GridType::dimension };
  
  int danz = 129; 
  /*
  int danz = this->functionSpace_.getGrid().size(level_, dim );
  danz = (int) pow (( double ) danz, (double) (1.0/dim) );
  std::cout << danz << " Danz!\n";
  */
  
  out << "P2\n " << danz << " " << danz <<"\n255\n";
  ConstDofIteratorType enddof = this->dend ( );
  for(ConstDofIteratorType itdof = this->dbegin ( ); itdof != enddof; ++itdof) {
    out << (int)((*itdof)*255.) << "\n";
  }
  out.close();
  return true;
}

template<class DiscreteFunctionSpaceType>
inline bool BlockVectorDiscreteFunction<DiscreteFunctionSpaceType>::
read_pgm( const std::string filename )
{
  const char * fn = filename.c_str();
  FILE *in;
  int v;

  in = fopen( fn, "r" );
  assert(in);
  fscanf( in, "P2\n%d %d\n%d\n", &v, &v, &v );
  DofIteratorType enddof = this->dend ( );
  for(DofIteratorType itdof = this->dbegin ( ); itdof != enddof; ++itdof) {
    fscanf( in, "%d", &v );
    (*itdof) = ((double)v)/255.;
  } 
  fclose( in );
  return true;
}


template< class DiscreteFunctionSpaceType >
inline void BlockVectorDiscreteFunction< DiscreteFunctionSpaceType >
  :: addScaled ( const DiscreteFunctionType &g,
                 const RangeFieldType &s )
{
  int length = dofVec_.size();
  const DofStorageType &gvec = g.dofVec_;
  assert(length == gvec.size());
  dofVec_.axpy( s, gvec );
}

template<class DiscreteFunctionSpaceType>
template<class GridIteratorType>
inline void BlockVectorDiscreteFunction<DiscreteFunctionSpaceType>::
substractLocal( GridIteratorType &it , 
 const BlockVectorDiscreteFunction<DiscreteFunctionSpaceType> &g)
{
  localFunction( *it , localFunc_ );
  
  BlockVectorDiscreteFunction<DiscreteFunctionSpaceType> &G = 
      const_cast<BlockVectorDiscreteFunction<DiscreteFunctionSpaceType> &> (g);
  G.localFunction(*it,G.localFunc_);

  int length = localFunc_.numberOfDofs();
  for(int i=0; i<length; ++i)
    localFunc_[i] -= G.localFunc_[i];
}

template<class DiscreteFunctionSpaceType>
template<class GridIteratorType>
inline void BlockVectorDiscreteFunction<DiscreteFunctionSpaceType>::
setLocal( GridIteratorType &it , const RangeFieldType & scalar )
{
  localFunction( *it , localFunc_ );
  int length = localFunc_.numberOfDofs();
  for(int i=0; i<length; ++i)
    localFunc_[i] = scalar;
}

template<class DiscreteFunctionSpaceType>
inline void BlockVectorDiscreteFunction<DiscreteFunctionSpaceType>::
enableDofCompression()
{
  assert( memPair_.first );
  memPair_.first->enableDofCompression();
}

} // end namespace Dune 
#endif
