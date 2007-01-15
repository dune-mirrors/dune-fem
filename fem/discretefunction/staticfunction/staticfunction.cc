#ifndef DUNE_DISCFUNCARRAY_CC
#define DUNE_DISCFUNCARRAY_CC

namespace Dune 
{

// Constructor making discrete function  
template<class DiscreteFunctionSpaceType, class DofStorageImp >
inline StaticDiscreteFunction< DiscreteFunctionSpaceType,DofStorageImp >::
StaticDiscreteFunction(const DiscreteFunctionSpaceType & f) 
: DiscreteFunctionDefaultType ( f )  
  , name_ ( "no name" )
  , mapper_(f.indexSet())
  , dm_(DofManagerFactoryType::getDofManager(f.grid()))
  , memPair_(dm_.addDofSet(&dofVec_, mapper_, name_)) 
  , dofVec_( *memPair_.second ) 
  , localFunc_ ( f , mapper_ , dofVec_ ) 
  , leakPointer_(dofVec_)                             
{
}

// Constructor making discrete function  
template<class DiscreteFunctionSpaceType, class DofStorageImp >
inline StaticDiscreteFunction< DiscreteFunctionSpaceType,DofStorageImp >::
StaticDiscreteFunction(const std::string name, const DiscreteFunctionSpaceType & f ) 
: DiscreteFunctionDefaultType ( f )  
  , name_ ( name )
  , mapper_(f.indexSet(),1)
  , dm_(DofManagerFactoryType::getDofManager(f.grid()))
  , memPair_(dm_.addDofSet(&dofVec_, mapper_, name_)) 
  , dofVec_( *memPair_.second ) 
  , localFunc_ ( f , mapper_, dofVec_ ) 
  , leakPointer_(dofVec_)                             
{
}

// Constructor making discrete function  
template<class DiscreteFunctionSpaceType, class DofStorageImp >
inline StaticDiscreteFunction< DiscreteFunctionSpaceType,DofStorageImp >::
StaticDiscreteFunction(const std::string name, 
    const DiscreteFunctionSpaceType & f, const DofStorageType & data ) 
: DiscreteFunctionDefaultType ( f )  
  , name_ ( name )
  , mapper_(f.indexSet(),1)
  , dm_(DofManagerFactoryType::getDofManager(f.grid()))
  , memPair_(dm_.addDummyDofSet(&dofVec_, mapper_, name_, &data)) 
  , dofVec_( *memPair_.second ) 
  , localFunc_ ( f , mapper_ , dofVec_ ) 
  , leakPointer_(dofVec_)                             
{
}

template<class DiscreteFunctionSpaceType, class DofStorageImp >
inline StaticDiscreteFunction< DiscreteFunctionSpaceType,DofStorageImp >::
StaticDiscreteFunction(const StaticDiscreteFunction< DiscreteFunctionSpaceType,DofStorageImp > & df ) :
  DiscreteFunctionDefaultType ( df.functionSpace_ ) 
  , mapper_(df.functionSpace_.indexSet(),1)
  , dm_(df.dm_)
  , memPair_(dm_.addDofSet(&dofVec_, mapper_, name_)) 
  , dofVec_( *memPair_.second ) 
  , localFunc_ ( this->functionSpace_ , mapper_, dofVec_ )
  , leakPointer_(dofVec_)        
{
  name_ = df.name_;
  built_ = df.built_; 

  dofVec_ = df.dofVec_;
} 

// Desctructor 
template<class DiscreteFunctionSpaceType, class DofStorageImp >
inline StaticDiscreteFunction< DiscreteFunctionSpaceType,DofStorageImp >::
~StaticDiscreteFunction() 
{
  dm_.removeDofSet(*memPair_.first);
}

template<class DiscreteFunctionSpaceType, class DofStorageImp >
inline void StaticDiscreteFunction< DiscreteFunctionSpaceType,DofStorageImp >::clear ()
{
  const int size = dofVec_.size();
  for(int i=0; i<size; ++i) dofVec_[i] = 0.0; 
}

template<class DiscreteFunctionSpaceType, class DofStorageImp >
inline void StaticDiscreteFunction< DiscreteFunctionSpaceType,DofStorageImp >::print(std::ostream &s ) const
{
  s << "StaticDiscreteFunction '" << name_ << "'\n";
  ConstDofIteratorType enddof = this->dend ();
  for(ConstDofIteratorType itdof = this->dbegin (); itdof != enddof; ++itdof) 
  {
    s << (*itdof) << " \n";
  } 
}
//*************************************************************************
//  Interface Methods 
//*************************************************************************
template<class DiscreteFunctionSpaceType, class DofStorageImp > template <class EntityType>
inline void
StaticDiscreteFunction< DiscreteFunctionSpaceType,DofStorageImp >::
localFunction ( const EntityType &en , LocalFunctionType &lf )
{
  lf.init ( en );
}

template<class DiscreteFunctionSpaceType, class DofStorageImp >
inline typename StaticDiscreteFunction<
DiscreteFunctionSpaceType,DofStorageImp >:: LocalFunctionImp *
StaticDiscreteFunction< DiscreteFunctionSpaceType,DofStorageImp >::
newLocalFunctionObject () const

{
    return new LocalFunctionImp ( this->functionSpace_ , mapper_, dofVec_ );
}

template<class DiscreteFunctionSpaceType, class DofStorageImp > template <class EntityType>
inline typename StaticDiscreteFunction< DiscreteFunctionSpaceType,DofStorageImp >:: LocalFunctionType 
StaticDiscreteFunction< DiscreteFunctionSpaceType,DofStorageImp >:: localFunction ( const EntityType &en ) const
{
  return LocalFunctionType (en,*this);
}

template<class DiscreteFunctionSpaceType, class DofStorageImp > 
inline typename StaticDiscreteFunction< DiscreteFunctionSpaceType,DofStorageImp >::DofIteratorType 
StaticDiscreteFunction< DiscreteFunctionSpaceType,DofStorageImp >::dbegin ()
{
  DofIteratorType tmp ( dofVec_ , 0 );     
  return tmp;
}


template<class DiscreteFunctionSpaceType, class DofStorageImp > 
inline typename StaticDiscreteFunction< DiscreteFunctionSpaceType,DofStorageImp >::DofIteratorType 
StaticDiscreteFunction< DiscreteFunctionSpaceType,DofStorageImp >::dend ( )
{
  DofIteratorType tmp ( dofVec_  , dofVec_.size() );     
  return tmp;
}

template<class DiscreteFunctionSpaceType, class DofStorageImp > 
inline typename StaticDiscreteFunction< DiscreteFunctionSpaceType,DofStorageImp >::ConstDofIteratorType 
StaticDiscreteFunction< DiscreteFunctionSpaceType,DofStorageImp >::dbegin ( ) const
{
  DofIteratorType tmp ( dofVec_ , 0 );     
  ConstDofIteratorType tmp2(tmp);
  return tmp2;
}

template<class DiscreteFunctionSpaceType, class DofStorageImp > 
inline typename StaticDiscreteFunction< DiscreteFunctionSpaceType,DofStorageImp >::ConstDofIteratorType 
StaticDiscreteFunction< DiscreteFunctionSpaceType,DofStorageImp >::dend ( ) const 
{
  DofIteratorType tmp ( dofVec_ , dofVec_.size() );     
  ConstDofIteratorType tmp2(tmp);
  return tmp2;
}
//**************************************************************************
//  Read and Write Methods 
//**************************************************************************
template<class DiscreteFunctionSpaceType, class DofStorageImp >
inline bool StaticDiscreteFunction< DiscreteFunctionSpaceType,DofStorageImp >::
write_xdr( std::string filename ) const
{
  const char * fn = filename.c_str();
  FILE  *file;
  XDR   xdrs;

  file = fopen(fn, "wb");
  if (!file)
  { 
    printf( "\aERROR in StaticDiscreteFunction::write_xdr(..): couldnot open <%s>!\n", fn);
    fflush(stderr);
    return false;
  }

  xdrstdio_create(&xdrs, file, XDR_ENCODE);   
  
  writeXdrs(&xdrs);

  xdr_destroy(&xdrs);
  fclose(file);
  
  return true;
}

template<class DiscreteFunctionSpaceType, class DofStorageImp >
inline bool StaticDiscreteFunction< DiscreteFunctionSpaceType,DofStorageImp >::
writeXdrs(XDR * xdrs) const  
{
  int len = dofVec_.size();
  xdr_int( xdrs, &len );
 
  enum { blockSize = DofBlockType :: dimension };
  for(int i=0; i<len; ++i) 
  {
    DofBlockType &dof = dofVec_[i];
    xdr_vector(xdrs,(char *) &dof[0], blockSize , sizeof(RangeFieldType) ,(xdrproc_t)xdr_double);
  }
  return true;
}

template<class DiscreteFunctionSpaceType, class DofStorageImp >
inline bool StaticDiscreteFunction< DiscreteFunctionSpaceType,DofStorageImp >::
readXdrs(XDR * xdrs)   
{
  int len = 0 ;
  xdr_int( xdrs, &len );
  assert( len == (int) dofVec_.size() );
 
  enum { blockSize = DofBlockType :: dimension };
  for(int i=0; i<len; ++i) 
  {
    DofBlockType &dof = dofVec_[i];
    xdr_vector(xdrs,(char *) &dof[0], blockSize , sizeof(RangeFieldType) ,(xdrproc_t)xdr_double);
  }
  return true;
}

template<class DiscreteFunctionSpaceType, class DofStorageImp >
inline bool StaticDiscreteFunction< DiscreteFunctionSpaceType,DofStorageImp >::
read_xdr( std::string filename)
{
  const char * fn = filename.c_str();
  FILE   *file;
  XDR     xdrs;

  std::cout << "Reading <" << fn << "> \n";
  file = fopen(fn, "rb");
  if(!file)
  { 
    printf( "\aERROR in StaticDiscreteFunction::read_xdr(..): couldnot open <%s>!\n", fn);
    fflush(stderr);
    return(false);
  }

  // read xdr 
  xdrstdio_create(&xdrs, file, XDR_DECODE);     

  readXdrs(&xdrs);
  
  xdr_destroy(&xdrs);
  fclose(file);
  return true;
}

template<class DiscreteFunctionSpaceType, class DofStorageImp >
inline bool StaticDiscreteFunction< DiscreteFunctionSpaceType,DofStorageImp >::
write_ascii( std::string filename ) const
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
    fprintf(stderr,"\aERROR in StaticDiscreteFunction::read_xdr(..): couldnot open <%s>!\n", fn);
    fflush(stderr);
    return(false);
  }
  return true;
}


template<class DiscreteFunctionSpaceType, class DofStorageImp >
inline bool StaticDiscreteFunction< DiscreteFunctionSpaceType,DofStorageImp >::
read_ascii( std::string filename )
{
  const char * fn = filename.c_str();
  FILE *infile=NULL;
  infile = fopen( fn, "r" );
  if(!infile)
  {
    std::cerr << "Couldnt open file! "<< fn << "\n";
    abort();
  }
  {
    int length;
    fscanf(infile,"%d \n",&length); 
    if(length != this->functionSpace_.size()) 
    {
      std::cerr << "ERROR: wrong number of dofs stored in file!\n"; 
      abort();
    }
    DofIteratorType enddof = this->dend ();
    for(DofIteratorType itdof = this->dbegin ();itdof != enddof; ++itdof) 
    {
      fscanf(infile,"%le \n",& (*itdof)); 
    }
  }
  fclose(infile);
  return true;
}

template<class DiscreteFunctionSpaceType, class DofStorageImp >
inline bool StaticDiscreteFunction< DiscreteFunctionSpaceType,DofStorageImp >::
write_pgm( std::string filename ) const
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

template<class DiscreteFunctionSpaceType, class DofStorageImp >
inline bool StaticDiscreteFunction< DiscreteFunctionSpaceType,DofStorageImp >::
read_pgm( std::string filename )
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


template<class DiscreteFunctionSpaceType, class DofStorageImp >
inline void StaticDiscreteFunction< DiscreteFunctionSpaceType,DofStorageImp >::
addScaled( const StaticDiscreteFunction< DiscreteFunctionSpaceType,DofStorageImp > &g, 
           const RangeFieldType &scalar )
{
  int length = dofVec_.size();
  const DofStorageImp &gvec = g.dofVec_;
  assert(length == gvec.size());
  
  for(int i=0; i<length; i++)
    dofVec_[i] += scalar*gvec[i];
}


template<class DiscreteFunctionSpaceType, class DofStorageImp >
template<class GridIteratorType>
inline void StaticDiscreteFunction< DiscreteFunctionSpaceType,DofStorageImp >::
substractLocal( GridIteratorType &it , 
 const StaticDiscreteFunction< DiscreteFunctionSpaceType,DofStorageImp > &g)
{
  localFunction( *it , localFunc_ );
  
  StaticDiscreteFunction< DiscreteFunctionSpaceType,DofStorageImp > &G = 
      const_cast<StaticDiscreteFunction< DiscreteFunctionSpaceType,DofStorageImp > &> (g);
  G.localFunction(*it,G.localFunc_);

  int length = localFunc_.numberOfDofs();
  for(int i=0; i<length; i++)
    localFunc_[i] -= G.localFunc_[i];
}

template<class DiscreteFunctionSpaceType, class DofStorageImp >
template<class GridIteratorType>
inline void StaticDiscreteFunction< DiscreteFunctionSpaceType,DofStorageImp >::
setLocal( GridIteratorType &it , const RangeFieldType & scalar )
{
  localFunction( *it , localFunc_ );
  int length = localFunc_.numberOfDofs();
  for(int i=0; i<length; i++)
    localFunc_[i] = scalar;
}
//**********************************************************************
//  --StaticDiscreteLocalFunction 
//**********************************************************************
template<class DiscreteFunctionType >
inline StaticDiscreteLocalFunction < DiscreteFunctionType >::
StaticDiscreteLocalFunction( 
    const DiscreteFunctionSpaceType &f , 
    const MapperType& mapper,
    DofStorageType & dofVec )
 : fSpace_ ( f )
 , mapper_(mapper)
 , en_(0), dofVec_ ( dofVec ) 
 , uniform_(! (fSpace_.multipleGeometryTypes()))
 , init_(false)
 , baseSet_(0)
{
  // only works for discontinuous spaces at the moment 
  assert( ! f.continuous() );
}
      
template<class DiscreteFunctionType >
inline StaticDiscreteLocalFunction < DiscreteFunctionType >::~StaticDiscreteLocalFunction() 
{
}

template<class DiscreteFunctionType >
inline typename StaticDiscreteLocalFunction < DiscreteFunctionType >::RangeFieldType & 
StaticDiscreteLocalFunction < DiscreteFunctionType >::operator [] (int num) 
{
  return (* (values_[num]));
}

template<class DiscreteFunctionType >
inline const typename StaticDiscreteLocalFunction < DiscreteFunctionType >::RangeFieldType & 
StaticDiscreteLocalFunction < DiscreteFunctionType >::operator [] (int num) const
{
  return (* (values_[num]));
}

template<class DiscreteFunctionType >
inline int StaticDiscreteLocalFunction < DiscreteFunctionType >::
numDofs () const 
{
  return numOfDof_;
}

template<class DiscreteFunctionType >
inline void StaticDiscreteLocalFunction < DiscreteFunctionType >::
evaluate (EntityType &en, const DomainType & x, RangeType & ret) const 
{
  ret = 0.0;
  xtmp_ = en.geometry().local(x);
  evaluateLocal(en, xtmp_, ret);
}

template<class DiscreteFunctionType >
inline void StaticDiscreteLocalFunction < DiscreteFunctionType >::
evaluate (const DomainType & local, RangeType & ret) const 
{
  enum { dimRange = DiscreteFunctionSpaceType::DimRange };
  assert(init_);
  assert(en().geometry().checkInside(local));
  ret = 0.0;
  const BaseFunctionSetType& bSet = baseFunctionSet();

  for (int i = 0; i < bSet.numBaseFunctions(); ++i)
  {
    bSet.eval(i, local , tmp_);
    for (int l = 0; l < dimRange; ++l) {
      ret[l] += (*values_[i]) * tmp_[l];
    }
  }
}

template<class DiscreteFunctionType >
inline void StaticDiscreteLocalFunction < DiscreteFunctionType >::
evaluateLocal (EntityType &en, const DomainType & local, RangeType & ret) const
{
  evaluate(local,ret);
}

template<class DiscreteFunctionType >
inline void StaticDiscreteLocalFunction < DiscreteFunctionType >::
evaluateLocal(const DomainType & x, RangeType & ret) const
{
  enum { dimRange = DiscreteFunctionSpaceType::DimRange };
  assert(init_);
  assert(en().geometry().checkInside(x));
  ret = 0.0;
  const BaseFunctionSetType& bSet = baseFunctionSet();

  for (int i = 0; i < bSet.numBaseFunctions(); ++i)
  {
    bSet.eval(i, x, tmp_);
    for (int l = 0; l < dimRange; ++l) 
    {
      ret[l] += (*values_[i]) * tmp_[l];
    }
  }
}

template<class DiscreteFunctionType > 
template <class QuadratureType> 
inline void StaticDiscreteLocalFunction < DiscreteFunctionType >::
evaluate (EntityType &en, QuadratureType &quad, int quadPoint, RangeType & ret) const 
{
  evaluateLocal(en, quad.point(quadPoint), ret);
}

template<class DiscreteFunctionType > 
template <class QuadratureType> 
inline void StaticDiscreteLocalFunction < DiscreteFunctionType >::
evaluate (QuadratureType &quad, int quadPoint, RangeType & ret) const 
{
  evaluateLocal(quad.point(quadPoint), ret);
}

template<class DiscreteFunctionType >
template <class QuadratureType>
inline void StaticDiscreteLocalFunction < DiscreteFunctionType >::
jacobian (EntityType &en, QuadratureType &quad, int quadPoint, JacobianRangeType & ret) const
{
  jacobianLocal(en, quad.point(quadPoint), ret);
}

template<class DiscreteFunctionType >
inline void StaticDiscreteLocalFunction < DiscreteFunctionType >::
jacobian(EntityType& en, const DomainType& x,
         JacobianRangeType& ret) const
{
  ret *= 0.0;
  xtmp_ = en.geometry().local(x);
  jacobianLocal(en, xtmp_, ret);
}

template<class DiscreteFunctionType >
inline void StaticDiscreteLocalFunction < DiscreteFunctionType >::
jacobianLocal(EntityType& en, const DomainType& x,
              JacobianRangeType& ret) const
{
  assert(init_);
  enum { dim = EntityType::dimension };
  enum { dimRange = DiscreteFunctionSpaceType::DimRange };

  ret = 0.0;
  const BaseFunctionSetType& bSet = baseFunctionSet();

  for (int i = 0; i < bSet.numBaseFunctions(); ++i) {
    tmpGrad_ *= 0.0;
    bSet.jacobian(i, x, tmpGrad_);

    for (int l = 0; l < dimRange; ++l) {
      tmpGrad_[l] *= *(values_[i]);
        // * umtv or umv?
      en.geometry().jacobianInverseTransposed(x).umv(tmpGrad_[l], ret[l]);
    }
  }
}

template<class DiscreteFunctionType >
inline
const typename
StaticDiscreteLocalFunction < DiscreteFunctionType >:: BaseFunctionSetType&
StaticDiscreteLocalFunction < DiscreteFunctionType >::
baseFunctionSet() const 
{
  assert(init_ && baseSet_);
  return *baseSet_;
}

template<class DiscreteFunctionType >
template <class EntityImp>
inline void StaticDiscreteLocalFunction < DiscreteFunctionType >::
init (const EntityImp &en ) const
{
  if(!uniform_ || !init_)
  {
    baseSet_  = &fSpace_.baseFunctionSet(en);
    numOfDof_ = baseSet_->numBaseFunctions();

    if(numOfDof_ > values_.size())
      values_.resize( numOfDof_ );

    init_ = true;
  }

  // cache entity 
  en_ = &en;

  // cache local dofs 
  DofBlockType& dofs = dofVec_[mapper_.mapToGlobal(en,0)] ;
  assert( numOfDof_ == DofBlockType :: dimension );
  for(int i=0; i<numOfDof_; i++) values_ [i] = &dofs[i]; 
  return ;
} 

} // end namespace Dune 
#endif
