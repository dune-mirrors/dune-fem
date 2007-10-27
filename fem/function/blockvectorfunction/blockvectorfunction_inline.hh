#ifndef DUNE_BLOCKVECTORFUNCTION_INLINE_HH
#define DUNE_BLOCKVECTORFUNCTION_INLINE_HH

namespace Dune 
{

// Constructor making discrete function  
template<class DiscreteFunctionSpaceType>
inline BlockVectorDiscreteFunction<DiscreteFunctionSpaceType>::
BlockVectorDiscreteFunction(const DiscreteFunctionSpaceType & f) 
  : DiscreteFunctionDefaultType ( f , lfFactory_ ) 
  , lfFactory_( *this )
  , name_ ( "no name" )
  , mapper_( f.blockMapper() ) 
  , dm_(DofManagerFactoryType::getDofManager(f.grid()))
  , memPair_(dm_.addDofSet(&dofVec_, mapper_, name_)) 
  , dofVec_( *memPair_.second ) 
  , leakPtr_(dofVec_)
#ifdef NEW_LOCALFUNCTION
  , localFunc_( *this )
#else
  , localFunc_( f, mapper_, dofVec_, leakPtr_ ) 
#endif
{
}

// Constructor making discrete function  
template<class DiscreteFunctionSpaceType>
inline BlockVectorDiscreteFunction<DiscreteFunctionSpaceType>::
BlockVectorDiscreteFunction(const std::string name, const DiscreteFunctionSpaceType & f ) 
  : DiscreteFunctionDefaultType ( f , lfFactory_ ) 
  , lfFactory_( *this )
  , name_ ( name )
  , mapper_( f.blockMapper() ) 
  , dm_(DofManagerFactoryType::getDofManager(f.grid()))
  , memPair_(dm_.addDofSet(&dofVec_, mapper_, name_)) 
  , dofVec_( *memPair_.second ) 
  , leakPtr_(dofVec_)
#ifdef NEW_LOCALFUNCTION
  , localFunc_( *this )
#else
  , localFunc_( f, mapper_, dofVec_, leakPtr_ ) 
#endif
{
}

// Constructor making discrete function  
template<class DiscreteFunctionSpaceType>
inline BlockVectorDiscreteFunction<DiscreteFunctionSpaceType>::
BlockVectorDiscreteFunction(const std::string name, 
    const DiscreteFunctionSpaceType & f, const DofStorageType & data ) 
  : DiscreteFunctionDefaultType ( f , lfFactory_ ) 
  , lfFactory_( *this )
  , name_ ( name )
  , mapper_( f.blockMapper() ) 
  , dm_(DofManagerFactoryType::getDofManager(f.grid()))
  , memPair_(dm_.addDummyDofSet(&dofVec_, mapper_, name_, &data)) 
  , dofVec_( *memPair_.second ) 
  , leakPtr_(dofVec_)
#ifdef NEW_LOCALFUNCTION
  , localFunc_( *this )
#else
  , localFunc_( f, mapper_, dofVec_, leakPtr_ ) 
#endif
{
}

template< class DiscreteFunctionSpaceType >
inline BlockVectorDiscreteFunction< DiscreteFunctionSpaceType >
  :: BlockVectorDiscreteFunction ( const ThisType &df ) 
: DiscreteFunctionDefaultType( df.functionSpace_, lfFactory_ )
  , lfFactory_( *this )
  , mapper_( df.functionSpace_.blockMapper() ) 
  , dm_(df.dm_)
  , memPair_(dm_.addDofSet(&dofVec_, mapper_, name_)) 
  , dofVec_( *memPair_.second ) 
  , leakPtr_(dofVec_)
#ifdef NEW_LOCALFUNCTION
  , localFunc_( *this )
#else
  , localFunc_( this->functionSpace_, mapper_, dofVec_, leakPtr_ ) 
#endif
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

template< class DiscreteFunctionSpaceType >
inline void BlockVectorDiscreteFunction< DiscreteFunctionSpaceType >
  :: print ( std::ostream &out ) const
{
  out << "BlockVectorDiscreteFunction '" << name_ << "'" << std :: endl;
  
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
//**********************************************************************
//  --BlockVectorLocalFunction 
//**********************************************************************
template< class TraitsImp >
inline BlockVectorLocalFunction < TraitsImp >::
BlockVectorLocalFunction( 
    const DiscreteFunctionSpaceType &f , 
    const MapperType& mapper,
    DofStorageType & dofVec , 
    LeakPointerType& leakPtr )
 : fSpace_ ( f ),
   mapper_(mapper),
   en_(0),
   dofVec_ ( dofVec ),
   leakPtr_ ( leakPtr ),
   needCheckGeometry_( true ),
   baseSet_()
{
}
      
template< class TraitsImp >
inline BlockVectorLocalFunction < TraitsImp >::~BlockVectorLocalFunction() 
{
}

#if 0
template< class TraitsImp >
inline typename BlockVectorLocalFunction< TraitsImp > :: RangeFieldType &
BlockVectorLocalFunction < TraitsImp >::operator [] (int num) 
{
  return (* (values_[num]));
}

template< class TraitsImp >
inline const typename BlockVectorLocalFunction < TraitsImp >::RangeFieldType & 
BlockVectorLocalFunction < TraitsImp >::operator [] (int num) const
{
  return (* (values_[num]));
}

template< class TraitsImp >
inline int BlockVectorLocalFunction < TraitsImp >::
numDofs () const 
{
  return numDofs_;
}
#endif

template< class TraitsImp >
inline void BlockVectorLocalFunction < TraitsImp >::
evaluate (const DomainType & local, RangeType & ret) const 
{
  enum { dimRange = DiscreteFunctionSpaceType::DimRange };
  assert( en_ != 0 );
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

template< class TraitsImp > 
template <class QuadratureType> 
inline void BlockVectorLocalFunction < TraitsImp >::
evaluate (const QuadratureType &quad, const int quadPoint, RangeType & ret) const 
{
  enum { dimRange = DiscreteFunctionSpaceType::DimRange };
  assert( en_ != 0 );

  ret = 0.0;
  const int numBaseFunctions = baseFunctionSet().numBaseFunctions();
  for (int i = 0; i < numBaseFunctions; ++i)
  {
    baseFunctionSet().evaluate(i, quad, quadPoint , tmp_);
    tmp_ *= (*values_[i]);
    ret += tmp_;
  }
}

template< class TraitsImp >
inline void BlockVectorLocalFunction < TraitsImp >::
jacobian(const DomainType& x,
         JacobianRangeType& ret) const
{
  assert( en_ != 0 );
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

template< class TraitsImp >
template <class QuadratureType>
inline void BlockVectorLocalFunction < TraitsImp >::
jacobian (const QuadratureType &quad, const int quadPoint, JacobianRangeType & ret) const
{
  assert( en_ != 0 );
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

template< class TraitsImp >
inline
const typename
BlockVectorLocalFunction < TraitsImp >:: BaseFunctionSetType&
BlockVectorLocalFunction < TraitsImp >::
baseFunctionSet() const 
{
  assert( en_ != 0 );
  return baseSet_;
}

//! axpy operation for factor 
template< class TraitsImp >
template <class QuadratureType>
inline void 
BlockVectorLocalFunction < TraitsImp >::
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
template< class TraitsImp >
template <class QuadratureType>
inline void 
BlockVectorLocalFunction < TraitsImp >::
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
template< class TraitsImp >
template <class QuadratureType>
inline void 
BlockVectorLocalFunction < TraitsImp >::
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

template< class TraitsImp >
inline void 
BlockVectorLocalFunction < TraitsImp >::
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
template< class DiscreteFunctionType >
inline void BlockVectorLocalFunction< DiscreteFunctionType >
  :: init ( const EntityType &entity ) const
{
  const bool multipleBaseSets = fSpace_.multipleBaseFunctionSets();

  if( multipleBaseSets || needCheckGeometry_ )
  {
    // if multiple base sets skip geometry call
    bool updateBaseSet = true;
    if( !multipleBaseSets && (en_ != 0) )
      updateBaseSet = (baseSet_.geometryType() != entity.geometry().type());

    if( multipleBaseSets || updateBaseSet )
    {
      baseSet_  = fSpace_.baseFunctionSet( entity );
      // do not use baseFunctionSet() here, entity might not have been set 
      numDofs_ = baseSet_.numBaseFunctions();

      if( numDofs_ > values_.size() )
        values_.resize( numDofs_ );
    }
  }
  
  // cache entity 
  en_ = &entity;

  // cache local dofs 
  //DofBlockType& dofs = dofVec_[mapper_.mapToGlobal(entity,0)] ;
  //assert( numDofs_ == DofBlockType :: dimension );
  for(int i=0; i<numDofs_; ++i) 
  {
    // assert that mapper matches with space mapping 
    /*
    assert( (mapper_.mapToGlobal(entity,0) * DofBlockType::dimension + i) ==
             fSpace_.mapToGlobal(entity,i) );
    values_ [i] = &dofs[i]; 
    */
    values_ [i] = &leakPtr_[ fSpace_.mapToGlobal(entity,i) ]; 
  }
} 

} // end namespace Dune 
#endif
