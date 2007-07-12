#ifndef DUNE_BLOCKVECTORFUNCTION_INLINE_HH
#define DUNE_BLOCKVECTORFUNCTION_INLINE_HH

namespace Dune 
{

// Constructor making discrete function  
template<class DiscreteFunctionSpaceType>
inline BlockVectorDiscreteFunction<DiscreteFunctionSpaceType>::
BlockVectorDiscreteFunction(const DiscreteFunctionSpaceType & f) 
: DiscreteFunctionDefaultType ( f )  
  , name_ ( "no name" )
  , mapper_(f.indexSet())
  , dm_(DofManagerFactoryType::getDofManager(f.grid()))
  , memPair_(dm_.addDofSet(&dofVec_, mapper_, name_)) 
  , dofVec_( *memPair_.second ) 
  , localFunc_ ( f , mapper_ , dofVec_ ) 
{
}

// Constructor making discrete function  
template<class DiscreteFunctionSpaceType>
inline BlockVectorDiscreteFunction<DiscreteFunctionSpaceType>::
BlockVectorDiscreteFunction(const std::string name, const DiscreteFunctionSpaceType & f ) 
: DiscreteFunctionDefaultType ( f )  
  , name_ ( name )
  , mapper_(f.indexSet(),1)
  , dm_(DofManagerFactoryType::getDofManager(f.grid()))
  , memPair_(dm_.addDofSet(&dofVec_, mapper_, name_)) 
  , dofVec_( *memPair_.second ) 
  , localFunc_ ( f , mapper_, dofVec_ ) 
{
}

// Constructor making discrete function  
template<class DiscreteFunctionSpaceType>
inline BlockVectorDiscreteFunction<DiscreteFunctionSpaceType>::
BlockVectorDiscreteFunction(const std::string name, 
    const DiscreteFunctionSpaceType & f, const DofStorageType & data ) 
: DiscreteFunctionDefaultType ( f )  
  , name_ ( name )
  , mapper_(f.indexSet(),1)
  , dm_(DofManagerFactoryType::getDofManager(f.grid()))
  , memPair_(dm_.addDummyDofSet(&dofVec_, mapper_, name_, &data)) 
  , dofVec_( *memPair_.second ) 
  , localFunc_ ( f , mapper_ , dofVec_ ) 
{
}

template<class DiscreteFunctionSpaceType>
inline BlockVectorDiscreteFunction<DiscreteFunctionSpaceType>::
BlockVectorDiscreteFunction(const BlockVectorDiscreteFunction<DiscreteFunctionSpaceType> & df ) :
  DiscreteFunctionDefaultType ( df.functionSpace_ ) 
  , mapper_(df.functionSpace_.indexSet(),1)
  , dm_(df.dm_)
  , memPair_(dm_.addDofSet(&dofVec_, mapper_, name_)) 
  , dofVec_( *memPair_.second ) 
  , localFunc_ ( this->functionSpace_ , mapper_, dofVec_ )
{
  name_ = df.name_;
  built_ = df.built_; 

  dofVec_ = df.dofVec_;
} 

// Desctructor 
template<class DiscreteFunctionSpaceType>
inline BlockVectorDiscreteFunction<DiscreteFunctionSpaceType>::
~BlockVectorDiscreteFunction() 
{
  dm_.removeDofSet(*memPair_.first);
}

template<class DiscreteFunctionSpaceType>
inline void BlockVectorDiscreteFunction<DiscreteFunctionSpaceType>::clear ()
{
  const int size = dofVec_.size();
  for(int i=0; i<size; ++i) dofVec_[i] = 0.0; 
}

template<class DiscreteFunctionSpaceType>
inline void BlockVectorDiscreteFunction<DiscreteFunctionSpaceType>::print(std::ostream &s ) const
{
  s << "BlockVectorDiscreteFunction '" << name_ << "'\n";
  ConstDofIteratorType enddof = this->dend ();
  for(ConstDofIteratorType itdof = this->dbegin (); itdof != enddof; ++itdof) 
  {
    s << (*itdof) << " \n";
  } 
}
//*************************************************************************
//  Interface Methods 
//*************************************************************************
template<class DiscreteFunctionSpaceType> template <class EntityType>
inline void
BlockVectorDiscreteFunction<DiscreteFunctionSpaceType>::
localFunction ( const EntityType &en , LocalFunctionType &lf )
{
  lf.init ( en );
}

template<class DiscreteFunctionSpaceType>
inline typename BlockVectorDiscreteFunction<
DiscreteFunctionSpaceType>:: LocalFunctionImp *
BlockVectorDiscreteFunction<DiscreteFunctionSpaceType>::
newObject () const

{
  return new LocalFunctionImp ( this->functionSpace_ , mapper_, dofVec_ );
}

template<class DiscreteFunctionSpaceType> template <class EntityType>
inline typename BlockVectorDiscreteFunction<DiscreteFunctionSpaceType>:: LocalFunctionType 
BlockVectorDiscreteFunction<DiscreteFunctionSpaceType>:: localFunction ( const EntityType &en ) const
{
  return LocalFunctionType (en,*this);
}

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
write_xdr( std::string filename ) const
{
  const char * fn = filename.c_str();
  XDR   xdrs;
  FILE* file = fopen(fn, "wb");

  if (!file)
  { 
    printf( "\aERROR in BlockVectorDiscreteFunction::write_xdr(..): couldnot open <%s>!\n", fn);
    fflush(stderr);
    return false;
  }

  xdrstdio_create(&xdrs, file, XDR_ENCODE);   
  
  writeXdrs(&xdrs);

  xdr_destroy(&xdrs);
  fclose(file);
  
  return true;
}

template<class DiscreteFunctionSpaceType>
inline bool BlockVectorDiscreteFunction<DiscreteFunctionSpaceType>::
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

template<class DiscreteFunctionSpaceType>
inline bool BlockVectorDiscreteFunction<DiscreteFunctionSpaceType>::
readXdrs(XDR * xdrs)   
{
  int len = 0 ;
  xdr_int( xdrs, &len );

  // for parallel runs len might be bigger than 
  // size we read we use Grape 
  assert( len >= dofVec_.size() );
 
  enum { blockSize = DofBlockType :: dimension };
  const int vecSize = dofVec_.size();
  for(int i=0; i<vecSize; ++i) 
  {
    DofBlockType &dof = dofVec_[i];
    xdr_vector(xdrs,(char *) &dof[0], blockSize , sizeof(RangeFieldType) ,(xdrproc_t)xdr_double);
  }
  return true;
}

template<class DiscreteFunctionSpaceType>
inline bool BlockVectorDiscreteFunction<DiscreteFunctionSpaceType>::
read_xdr( std::string filename)
{
  const char * fn = filename.c_str();
  XDR xdrs;

  std::cout << "Reading <" << fn << "> \n";
  FILE* file = fopen(fn, "rb");
  if(!file)
  { 
    printf( "\aERROR in BlockVectorDiscreteFunction::read_xdr(..): couldnot open <%s>!\n", fn);
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

template<class DiscreteFunctionSpaceType>
inline bool BlockVectorDiscreteFunction<DiscreteFunctionSpaceType>::
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
    fprintf(stderr,"\aERROR in BlockVectorDiscreteFunction::read_xdr(..): couldnot open <%s>!\n", fn);
    fflush(stderr);
    return(false);
  }
  return true;
}


template<class DiscreteFunctionSpaceType>
inline bool BlockVectorDiscreteFunction<DiscreteFunctionSpaceType>::
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

template<class DiscreteFunctionSpaceType>
inline bool BlockVectorDiscreteFunction<DiscreteFunctionSpaceType>::
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

template<class DiscreteFunctionSpaceType>
inline bool BlockVectorDiscreteFunction<DiscreteFunctionSpaceType>::
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


template<class DiscreteFunctionSpaceType>
inline void BlockVectorDiscreteFunction<DiscreteFunctionSpaceType>::
addScaled( const BlockVectorDiscreteFunction<DiscreteFunctionSpaceType> &g, 
           const RangeFieldType &scalar )
{
  int length = dofVec_.size();
  const DofStorageType &gvec = g.dofVec_;
  assert(length == gvec.size());
  dofVec_.axpy(scalar,gvec);
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
 , en_(0)
 , dofVec_ ( dofVec ) 
 , multipleBaseFunctionSets_((fSpace_.multipleBaseFunctionSets()))
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
evaluate (const DomainType & local, RangeType & ret) const 
{
  enum { dimRange = DiscreteFunctionSpaceType::DimRange };
  assert(init_);
  assert(en().geometry().checkInside(local));
  ret = 0.0;

  const int numBase = baseFunctionSet().numBaseFunctions();
  for (int i = 0; i < numBase; ++i)
  {
    baseFunctionSet().evaluate(i, local , tmp_);
    for (int l = 0; l < dimRange; ++l) {
      ret[l] += (*values_[i]) * tmp_[l];
    }
  }
}

template<class DiscreteFunctionType > 
template <class QuadratureType> 
inline void StaticDiscreteLocalFunction < DiscreteFunctionType >::
evaluate (const QuadratureType &quad, const int quadPoint, RangeType & ret) const 
{
  enum { dimRange = DiscreteFunctionSpaceType::DimRange };
  assert(init_);

  ret = 0.0;
  const int numBaseFunctions = baseFunctionSet().numBaseFunctions();
  for (int i = 0; i < numBaseFunctions; ++i)
  {
    baseFunctionSet().evaluate(i, quad, quadPoint , tmp_);
    tmp_ *= (*values_[i]);
    ret += tmp_;
  }
}

template<class DiscreteFunctionType >
inline void StaticDiscreteLocalFunction < DiscreteFunctionType >::
jacobian(const DomainType& x,
         JacobianRangeType& ret) const
{
  assert(init_);
  enum { dim = EntityType::dimension };
  enum { dimRange = DiscreteFunctionSpaceType::DimRange };
  typedef typename DiscreteFunctionSpaceType::GridType::ctype ctype;

  ret = 0.0;
  const FieldMatrix<ctype,dim,dim>& inv =
    en().geometry().jacobianInverseTransposed(x);

  const int numBaseFct = baseFunctionSet().numBaseFunctions();
  for (int i = 0; i < numBaseFct; ++i) 
  {
    baseFunctionSet().jacobian(i, x , tmpGrad_);
    for (int l = 0; l < dimRange; ++l) 
    {
      tmpGrad_[l] *= *(values_[i]);
      inv.umv(tmpGrad_[l], ret[l]);
    }
  }
}

template<class DiscreteFunctionType >
template <class QuadratureType>
inline void StaticDiscreteLocalFunction < DiscreteFunctionType >::
jacobian (const QuadratureType &quad, const int quadPoint, JacobianRangeType & ret) const
{
  assert(init_);
  enum { dim = EntityType::dimension };
  enum { dimRange = DiscreteFunctionSpaceType::DimRange };
  typedef typename DiscreteFunctionSpaceType::GridType::ctype ctype;

  ret = 0.0;
  const FieldMatrix<ctype,dim,dim>& inv =
    en().geometry().jacobianInverseTransposed(quad.point(quadPoint));

  const int numBaseFct = baseFunctionSet().numBaseFunctions();
  for (int i = 0; i < numBaseFct; ++i) 
  {
    baseFunctionSet().jacobian(i, quad, quadPoint , tmpGrad_);
    for (int l = 0; l < dimRange; ++l) 
    {
      tmpGrad_[l] *= *(values_[i]);
      inv.umv(tmpGrad_[l], ret[l]);
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
  assert(init_);
  return baseSet_;
}

//! axpy operation for factor 
template<class DiscreteFunctionType >
template <class QuadratureType>
inline void 
StaticDiscreteLocalFunction < DiscreteFunctionType >::
axpy(const QuadratureType& quad, const int quadPoint, const RangeType& factor)
{
  const int numDof = this->numDofs();
  for(int i=0; i<numDof; ++i)
  {
    this->baseFunctionSet().evaluate( i , quad, quadPoint, tmp_ );
    (*values_[i]) += tmp_ * factor;
  }
}
  
//! axpy operation for factor 
template<class DiscreteFunctionType >
template <class QuadratureType>
inline void 
StaticDiscreteLocalFunction < DiscreteFunctionType >::
axpy(const QuadratureType& quad, const int quadPoint, const JacobianRangeType& factor)
{
  const int numDof = this->numDofs();

  const JacobianInverseType& jti =
    en().geometry().jacobianInverseTransposed(quad.point(quadPoint));
  rightMultiply( factor, jti, factorInv_ );

  for(int i=0; i<numDof; ++i)
  {
    // evaluate gradient on reference element
    this->baseFunctionSet().jacobian(i, quad, quadPoint , tmpGrad_);
    for (int l = 0; l < dimRange; ++l)
    {
      (*values_[i]) += tmpGrad_[l] * factorInv_[l];
    }
  }
}
  
//! axpy operation for factor 
template<class DiscreteFunctionType >
template <class QuadratureType>
inline void 
StaticDiscreteLocalFunction < DiscreteFunctionType >::
axpy(const QuadratureType& quad, 
     const int quadPoint, 
     const RangeType& factor1, 
     const JacobianRangeType& factor2)
{
  const int numDof = this->numDofs();

  const JacobianInverseType& jti =
    en().geometry().jacobianInverseTransposed(quad.point(quadPoint));
  rightMultiply( factor2, jti, factorInv_ );

  for(int i=0; i<numDof; ++i)
  {
    // evaluate gradient on reference element
    this->baseFunctionSet().evaluate(i, quad, quadPoint, tmp_ );
    (*values_[i]) += tmp_ * factor1;
    this->baseFunctionSet().jacobian(i, quad, quadPoint, tmpGrad_);
    for (int l = 0; l < dimRange; ++l)
    {
      (*values_[i]) += tmpGrad_[l] * factorInv_[l];
    }
  }
}

template<class DiscreteFunctionType >
inline void 
StaticDiscreteLocalFunction < DiscreteFunctionType >::
rightMultiply(const JacobianRangeType& factor,
              const JacobianInverseType& jInv,
              JacobianRangeType& result) const
{
  enum { rows = JacobianRangeType :: rows };
  enum { cols = JacobianInverseType :: rows };
  for (int i=0; i<rows; ++i)
  {
    for (int j=0; j<cols; ++j)
    {
      result[i][j] = 0;
      for (int k=0; k<cols; ++k)
      {
        result[i][j] += factor[i][k] * jInv[k][j];
      }
    }
  }
}


// --init
template<class DiscreteFunctionType >
template <class EntityImp>
inline void StaticDiscreteLocalFunction < DiscreteFunctionType >::
init (const EntityImp &en ) const
{
  // in not initilized or we have more then one base function set
  if(multipleBaseFunctionSets_ || !init_)
  {
    baseSet_  = fSpace_.baseFunctionSet(en);
    init_ = true;
    numOfDof_ = baseFunctionSet().numBaseFunctions();

    if(numOfDof_ > values_.size())
      values_.resize( numOfDof_ );
  }

  // cache entity 
  en_ = &en;

  // cache local dofs 
  DofBlockType& dofs = dofVec_[mapper_.mapToGlobal(en,0)] ;
  assert( numOfDof_ == DofBlockType :: dimension );
  for(int i=0; i<numOfDof_; ++i) 
  {
    // assert that mapper matches with space mapping 
    assert( (mapper_.mapToGlobal(en,0) * DofBlockType::dimension + i) ==
             fSpace_.mapToGlobal(en,i) );
    values_ [i] = &dofs[i]; 
  }
  return ;
} 

} // end namespace Dune 
#endif
