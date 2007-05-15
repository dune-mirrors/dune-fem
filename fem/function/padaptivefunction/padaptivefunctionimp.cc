#ifndef DUNE_PADAPTIVEFUNCTION_CC
#define DUNE_PADAPTIVEFUNCTION_CC



#include <algorithm>
// BUG: return values are initialized with *= 0. Problematic for nan
// values!
namespace Dune {

// Constructor making discrete function  
template<class DiscreteFunctionSpaceType, class DiscreteFunctionSpace2Type>
  inline PAdaptiveDiscreteFunction< DiscreteFunctionSpaceType, DiscreteFunctionSpace2Type>::
  PAdaptiveDiscreteFunction(const DiscreteFunctionSpaceType& f, const DiscreteFunctionSpace2Type& f2) :
    DiscreteFunctionDefaultType ( f ), 
    name_ ("no name"),
    dm_(DofManagerFactoryType::getDofManager(f.grid())),
    memPair_(dm_.addDofSet(&dofVec_, f.mapper(), name_)),
    dofVec_ ( *memPair_.second ),
    functionSpace2_ (f2)
{
    std::cout << "Size old = " << dofVec_.size() << "  Size  f2 = " << functionSpace2_.size();
    dofVec_.resize ( dofVec_.size() * functionSpace2_.size() );
    std::cout << "  Size new = " << dofVec_.size() << std::endl;
   
}
// Desctructor 
template<class DiscreteFunctionSpaceType, class DiscreteFunctionSpace2Type>
inline PAdaptiveDiscreteFunction< DiscreteFunctionSpaceType, DiscreteFunctionSpace2Type>::
~PAdaptiveDiscreteFunction()
{
  dm_.removeDofSet(*memPair_.first);
}


template<class DiscreteFunctionSpaceType, class DiscreteFunctionSpace2Type>
inline void PAdaptiveDiscreteFunction< DiscreteFunctionSpaceType, DiscreteFunctionSpace2Type>::set ( RangeFieldType x )
{
  int size = dofVec_.size();
  DofStorageType &vec = dofVec_;
  for(int i=0; i<size; ++i) vec[i] = x; 
}  

template<class DiscreteFunctionSpaceType, class DiscreteFunctionSpace2Type>
inline void PAdaptiveDiscreteFunction< DiscreteFunctionSpaceType, DiscreteFunctionSpace2Type>::clear ()
{
  const int size = dofVec_.size();
  DofStorageType &vec = dofVec_;
  for(int i=0; i<size; ++i) vec[i] = 0.0; 
  //set ( 0.0 ); 
}

template<class DiscreteFunctionSpaceType, class DiscreteFunctionSpace2Type>
inline void PAdaptiveDiscreteFunction< DiscreteFunctionSpaceType, DiscreteFunctionSpace2Type>::
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
template<class DiscreteFunctionSpaceType, class DiscreteFunctionSpace2Type> 
inline ParaLocalFunctionAdapt<PAdaptiveDiscreteFunction< DiscreteFunctionSpaceType, DiscreteFunctionSpace2Type> > *
PAdaptiveDiscreteFunction< DiscreteFunctionSpaceType, DiscreteFunctionSpace2Type>::
newLocalFunctionObject ( ) const
{
  return new ParaLocalFunctionAdapt<ThisType> ( *this );
}

template<class DiscreteFunctionSpaceType, class DiscreteFunctionSpace2Type> 
inline AdaptiveDiscreteFunction< DiscreteFunctionSpaceType>
PAdaptiveDiscreteFunction<DiscreteFunctionSpaceType, DiscreteFunctionSpace2Type>::localFunction(int dofIndex2) const
{
 return AdaptiveDiscreteFunction< DiscreteFunctionSpaceType> ( name(), this->functionSpace_, &(dofVec_.leakPointer()[this->functionSpace_.size()*dofIndex2]));
 
}

template<class DiscreteFunctionSpaceType, class DiscreteFunctionSpace2Type> 
template < class Entity2Type>
inline AdaptiveDiscreteFunction< DiscreteFunctionSpaceType>
PAdaptiveDiscreteFunction<DiscreteFunctionSpaceType, DiscreteFunctionSpace2Type>::localFunction(const Entity2Type &en2,int dofIndex2) const
{
 return AdaptiveDiscreteFunction< DiscreteFunctionSpaceType> ( name(), this->functionSpace_, &(dofVec_.leakPointer()[this->functionSpace_.size()*space2().mapToGlobal( en2 ,dofIndex2)]));
 
}

template<class DiscreteFunctionSpaceType, class DiscreteFunctionSpace2Type>
template < class Entity2Type, class LocalCoord2Type>
inline void PAdaptiveDiscreteFunction<DiscreteFunctionSpaceType, DiscreteFunctionSpace2Type>::localFunction(const Entity2Type &en2, const LocalCoord2Type &loc2,  DiscreteFunction1Type &discFunc) const
{
 typedef typename DiscreteFunctionSpace2Type::BaseFunctionSetType BaseFunctionSetType;
 typename DiscreteFunctionSpace2Type::RangeType tmp_;
 discFunc.clear();
 const BaseFunctionSetType& bSet2 = space2().baseFunctionSet(en2);
 const int numOfDofs = bSet2.numBaseFunctions();
 
 for(int i = 0; i< numOfDofs; i++) 
 {
  const int map = space2().mapToGlobal( en2 , i );  
  bSet2.evaluate(i,loc2,tmp_);
  DiscreteFunction1Type df = this->localFunction(map);
	    discFunc.addScaled(df,tmp_);
 }
	
}

template<class DiscreteFunctionSpaceType, class DiscreteFunctionSpace2Type> 
inline typename PAdaptiveDiscreteFunction< DiscreteFunctionSpaceType, DiscreteFunctionSpace2Type>::DofIteratorType 
PAdaptiveDiscreteFunction< DiscreteFunctionSpaceType, DiscreteFunctionSpace2Type>::dbegin ( )
{
  return dofVec_.begin();
}

template<class DiscreteFunctionSpaceType, class DiscreteFunctionSpace2Type> 
inline typename PAdaptiveDiscreteFunction< DiscreteFunctionSpaceType, DiscreteFunctionSpace2Type>::DofIteratorType 
PAdaptiveDiscreteFunction< DiscreteFunctionSpaceType, DiscreteFunctionSpace2Type>::dend ()
{
  return dofVec_.end();
}

template<class DiscreteFunctionSpaceType, class DiscreteFunctionSpace2Type> 
inline typename PAdaptiveDiscreteFunction< DiscreteFunctionSpaceType, DiscreteFunctionSpace2Type>::ConstDofIteratorType 
PAdaptiveDiscreteFunction< DiscreteFunctionSpaceType, DiscreteFunctionSpace2Type>::dbegin ( ) const
{
  return dofVec_.begin();
}

template<class DiscreteFunctionSpaceType, class DiscreteFunctionSpace2Type> 
inline typename PAdaptiveDiscreteFunction< DiscreteFunctionSpaceType, DiscreteFunctionSpace2Type>::ConstDofIteratorType 
PAdaptiveDiscreteFunction< DiscreteFunctionSpaceType, DiscreteFunctionSpace2Type>::dend () const 
{
  return dofVec_.end();
}
//**************************************************************************
//  Read and Write Methods 
//**************************************************************************
template<class DiscreteFunctionSpaceType, class DiscreteFunctionSpace2Type>
inline bool PAdaptiveDiscreteFunction< DiscreteFunctionSpaceType, DiscreteFunctionSpace2Type>::
write_xdr(std::string fn) const
{
  FILE  *file;
  XDR   xdrs;
  file = fopen(fn.c_str(), "wb");
  if (!file)
  { 
    printf( "\aERROR in PAdaptiveDiscreteFunction::write_xdr(..): couldnot open <%s>!\n", fn.c_str());
    fflush(stderr);
    return false;
  }

  xdrstdio_create(&xdrs, file, XDR_ENCODE);   
  dofVec_.processXdr(&xdrs);

  xdr_destroy(&xdrs);
  fclose(file);
  
  return true;
}

template<class DiscreteFunctionSpaceType, class DiscreteFunctionSpace2Type>
inline bool PAdaptiveDiscreteFunction< DiscreteFunctionSpaceType, DiscreteFunctionSpace2Type>::
read_xdr(std::string fn)
{
  FILE   *file;
  XDR     xdrs;
  std::cout << "Reading <" << fn << "> \n";
  file = fopen(fn.c_str() , "rb");
  if(!file)
  { 
    printf( "\aERROR in PAdaptiveDiscreteFunction::read_xdr(..): couldnot open <%s>!\n", fn.c_str());
    fflush(stderr);
    return(false);
  }

  // read xdr 
  xdrstdio_create(&xdrs, file, XDR_DECODE);     
  dofVec_.processXdr(&xdrs);
  
  xdr_destroy(&xdrs);
  fclose(file);
  return true;
}

template<class DiscreteFunctionSpaceType, class DiscreteFunctionSpace2Type>
inline bool PAdaptiveDiscreteFunction< DiscreteFunctionSpaceType, DiscreteFunctionSpace2Type>::
write_ascii(std::string fn) const
{
  std::fstream outfile( fn.c_str() , std::ios::out );
  if (!outfile)
  { 
    printf( "\aERROR in PAdaptiveDiscreteFunction::write_ascii(..): couldnot open <%s>!\n", fn.c_str());
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
inline bool PAdaptiveDiscreteFunction< DiscreteFunctionSpaceType, DiscreteFunctionSpace2Type>::
read_ascii(std::string fn)
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
inline bool PAdaptiveDiscreteFunction< DiscreteFunctionSpaceType, DiscreteFunctionSpace2Type>::
write_pgm(std::string fn) const
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
inline bool PAdaptiveDiscreteFunction< DiscreteFunctionSpaceType, DiscreteFunctionSpace2Type>::
read_pgm(std::string fn)
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
inline void PAdaptiveDiscreteFunction< DiscreteFunctionSpaceType, DiscreteFunctionSpace2Type>::
addScaled( const PAdaptiveDiscreteFunction< DiscreteFunctionSpaceType, DiscreteFunctionSpace2Type> &g, 
           const RangeFieldType &scalar )
{
  int length = dofVec_.size();

  DofStorageType &v = dofVec_;
  const DofStorageType &gvec = g.dofVec_;
  
  for(int i=0; i<length; i++)
    v[i] += scalar*gvec[i];
}

//**********************************************************************
//  --ParaLocalFunctionAdapt 
//**********************************************************************
template<class DiscreteFunctionType>
inline ParaLocalFunctionAdapt<DiscreteFunctionType>::
ParaLocalFunctionAdapt(const DiscreteFunctionType &df) :
  tmp_(0.0),
  xtmp_(0.0),
  tmpGrad_(0.0),
  numOfDof_(-1),
  df_(df),
  fSpace_ ( df.space() ),
  values_ (),
  dofVec_ ( df.dofVec_ ),
  baseSet_(0),
  init_(false),
  geoType_(0) // init as Vertex
{}
      

template<class DiscreteFunctionType>
inline ParaLocalFunctionAdapt<DiscreteFunctionType>::~ParaLocalFunctionAdapt() 
{
}

template<class DiscreteFunctionType>
inline std::string  
ParaLocalFunctionAdapt<DiscreteFunctionType>::name() const 
{
  std::ostringstream nombre; 
  nombre << df_.name() << "_local"; 
  return nombre.str();
}

template<class DiscreteFunctionType>
inline typename ParaLocalFunctionAdapt<DiscreteFunctionType>::RangeFieldType & 
ParaLocalFunctionAdapt<DiscreteFunctionType>::operator [] (int num) 
{
  // check that storage (dofVec_) and mapper are in sync:
  assert(dofVec_.size() >= fSpace_.size());
  return (* (values_[num]));
}

template<class DiscreteFunctionType>
inline const typename ParaLocalFunctionAdapt<DiscreteFunctionType>::RangeFieldType & 
ParaLocalFunctionAdapt<DiscreteFunctionType>::operator [] (int num) const
{ 
  // check that storage (dofVec_) and mapper are in sync:
  assert(dofVec_.size() >= fSpace_.size());
  return (* (values_[num]));
}


template<class DiscreteFunctionType>
inline int ParaLocalFunctionAdapt<DiscreteFunctionType>::
numDofs () const 
{
  return numOfDof_;
}

template<class DiscreteFunctionType> 
inline 
const typename 
ParaLocalFunctionAdapt<DiscreteFunctionType>::BaseFunctionSetType& 
ParaLocalFunctionAdapt<DiscreteFunctionType>::baseFunctionSet() const {
  assert(init_ && baseSet_);
  return *baseSet_;
}


template<class DiscreteFunctionType> 
template <class EntityType> 
inline void ParaLocalFunctionAdapt<DiscreteFunctionType>::
init (const EntityType &en) const
{
  // NOTE: if init is false, then LocalFunction has been create before. 
  // if fSpace_.multipleGeometryTypes() is true, then grid has elements 
  // of different geometry type (hybrid grid) and we have to check geometry
  // type again, if not we skip this part, because calling the entity's
  // geometry method is not a cheep call 
  
  if( !init_ || ( ! fSpace_.multipleGeometryTypes() ) )
  {
    if(geoType_ != en.geometry().type()) 
    {
      baseSet_  = &fSpace_.baseFunctionSet(en);
      numOfDof_ = baseSet_->numBaseFunctions();
    
      if(numOfDof_ > this->values_.size())
      {
        this->values_.resize( numOfDof_ );
      }

      init_ = true;
      geoType_ = en.geometry().type();
    }
  }

  // make sure that the base functions set we have is for right geom type
  assert( geoType_ == en.geometry().type() );

  for(int i=0; i<numOfDof_; ++i)
  {
    values_ [i] = &(this->dofVec_[ fSpace_.mapToGlobal ( en , i) ]);
  }
  return;
} 

template<class DiscreteFunctionType> 
inline void ParaLocalFunctionAdapt<DiscreteFunctionType>::
assign(int numDof, const RangeType& dofs) 
{
  assert(false); // untested and most probably wrong
  for (size_t i = 0; i < dimrange; ++i) {
    *(values_[numDof + dimrange*i]) = dofs[i]; 
  }
}

} // end namespace 
#endif
