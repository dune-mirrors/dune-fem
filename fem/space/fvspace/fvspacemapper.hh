#ifndef DUNE_FVSPACEMAPPER_HH
#define DUNE_FVSPACEMAPPER_HH

//- Dune includes 
#include <dune/fem/space/common/dofmapperinterface.hh>

namespace Dune {
 
/** \brief 
 The FVSpaceMapper maps local to global dofs for the FVSpace. 
  */
template <class IndexSetImp, int polOrd, int dimrange>
class FiniteVolumeMapper; 

//*****************************************************************
//
// specialisation for polynom order 0 and arbitray dimrange 
//
//*****************************************************************
template <class IndexSetImp, int dimrange>
class FiniteVolumeMapper<IndexSetImp,0,dimrange>
: public DofMapperDefault < FiniteVolumeMapper <IndexSetImp,0,dimrange> > 
{
  enum { numCodims = IndexSetImp::ncodim };

  int numberOfDofs_;
  const IndexSetImp & indexSet_;
    
public:
  typedef IndexSetImp IndexSetType;
  
  FiniteVolumeMapper ( const IndexSetType  & is , int numDofs ) 
    : numberOfDofs_ (numDofs) , indexSet_ (is) 
  {
    assert(numberOfDofs_ == dimrange);  
  }

  //! return size of function space, here number of elements 
  int size () const
  {
    return dimrange * indexSet_.size(0);
  }
  
  //! map Entity an local Dof number to global Dof number 
  //! for polOrd = 0
  template <class EntityType>
  int mapToGlobal (EntityType &en, int localNum ) const
  {
    return (dimrange * indexSet_.template index<0> (en,localNum)) + localNum;
  }

  //! return old index for hole 
  int oldIndex (const int hole, int ) const
  {   
    // corresponding number of set is newn 
    const int newn  = static_cast<int> (hole/dimrange);
    // local number of dof is local 
    const int local = (hole % dimrange); 
    return (dimrange * indexSet_.oldIndex(newn,0)) + local;
  }

  //! return new index for hole 
  int newIndex (const int hole, int ) const
  {
    // corresponding number of set is newn 
    const int newn = static_cast<int> (hole / dimrange);
    // local number of dof is local 
    const int local = (hole % dimrange); 
    return (dimrange * indexSet_.newIndex(newn,0)) + local;
  }

  //! return numbber of exsiting hole 
  int numberOfHoles ( int ) const
  {   
    // this index set works only for codim = 0 at the moment
    return dimrange * indexSet_.numberOfHoles(0);
  }
  
  // is called once and calcs the insertion points too
  int newSize() const 
  {
    return this->size();
  }

  //! return number of dof per element 
  int numDofs () const 
  {
    assert( numberOfDofs_ == dimrange );
    return numberOfDofs_;
  }

  //! return the sets needsCompress 
  bool needsCompress () const { return indexSet_.needsCompress(); }
};

template <class IndexSetImp>
class FiniteVolumeMapper<IndexSetImp,0,1>
: public DofMapperDefault < FiniteVolumeMapper <IndexSetImp,0,1> > 
{
  // corresp. index set 
  const IndexSetImp & indexSet_;
public:
  typedef IndexSetImp IndexSetType;
  
  FiniteVolumeMapper ( const IndexSetType  & is , int numDofs ) : indexSet_ (is) {}

  // we have virtual function ==> virtual destructor 
  virtual ~FiniteVolumeMapper () {}

  //! return size of function space, here number of elements 
  int size () const
  {
    return indexSet_.size(0);
  }
  
  //! map Entity an local Dof number to global Dof number 
  //! for polOrd = 0
  template <class EntityType>
  int mapToGlobal (EntityType &en, int localNum ) const
  {
    return indexSet_.template index<0> (en,localNum);
  }

  //! return old index of hole 
  int oldIndex (const int hole, int ) const
  {   
    return indexSet_.oldIndex(hole,0);
  }

  //! return new index of hole 
  int newIndex (const int hole, int ) const
  {
    return indexSet_.newIndex(hole,0);
  }

  //! return number of holes 
  int numberOfHoles ( int ) const
  {   
    // this index set works only for codim = 0 at the moment
    return indexSet_.numberOfHoles(0);
  }
  
  // is called once and calcs the insertion points too
  int newSize() const 
  {
    return this->size();
  }

  //! return number of dof per entity, here this method returns 1
  int numDofs () const 
  {
    return 1;
  }

  //! return the sets needsCompress 
  bool needsCompress () const { return indexSet_.needsCompress(); }
};

} // end namespace Dune 
#endif
