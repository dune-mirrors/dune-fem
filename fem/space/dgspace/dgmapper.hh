#ifndef DUNE_DGMAPPER_HH
#define DUNE_DGMAPPER_HH

#include "../common/dofmapperinterface.hh"

namespace Dune {

//***************************************************************************
//
//!  DG Mapper for mapping of local dof numbers to global dof numbers, 
//!  i.e. the entry in the vector of unknowns
//
//***************************************************************************
template <class IndexSetType, int polOrd, int dimRange>
class DGMapper
: public DofMapperDefault < DGMapper <IndexSetType,polOrd, dimRange > >
{
  // index set of grid, i.e. the element indices 
  const IndexSetType &indexSet_;

  // number of dofs on element 
  const int numberOfDofs_;
public:
  //! Constructor 
  DGMapper(const IndexSetType& iset , int numDof) :
    indexSet_ (iset), numberOfDofs_ (numDof)  {}

  //! return size of function space 
  //! see dofmanager.hh for definition of IndexSet, which 
  //! is a wrapper for en.index 
  int size () const
  {
    // return number of dofs * number of elements 
    return (numberOfDofs_ * indexSet_.size( 0 ));
  }

  //! map Entity an local Dof number to global Dof number 
  //! see dofmanager.hh for definition of IndexSet, which 
  //! is a wrapper for en.index 
  template <class EntityType>
  int mapToGlobal (EntityType &en, int localNum ) const
  {
    // indexSet_.index<0>(en,0) return 0 entity of codim the codim 0
    // entitys of en which is the element index 
    return ((indexSet_.template index<0>(en,0) * numberOfDofs_) + localNum);
  };

  //! default implementation if not overlaoded 
  int numDofs () const
  {
    return numberOfDofs_;
  }

  //! only called once, if grid was adapted 
  int newSize() const
  {
    return this->size();
  }

  //! return old index, for dof manager only 
  int oldIndex (const int hole, int ) const
  {
    // corresponding number of set is newn 
    const int newn  = static_cast<int> (hole / numberOfDofs_);
    // local number of dof is local 
    const int local = (hole % numberOfDofs_);
    return (numberOfDofs_ * indexSet_.oldIndex(newn,0)) + local;
  }

  //! return new index, for dof manager only 
  int newIndex (const int hole, int ) const
  {
    // corresponding number of set is newn 
    const int newn = static_cast<int> (hole / numberOfDofs_);
    // local number of dof is local 
    const int local = (hole % numberOfDofs_);
    return (numberOfDofs_ * indexSet_.newIndex(newn,0)) + local;
  }

  //! return size of grid entities per level and codim 
  //! for dof mapper 
  int numberOfHoles ( int ) const
  {
    // this index set works only for codim = 0 at the moment
    return numberOfDofs_ * indexSet_.numberOfHoles(0);
  }

  //! return whether the data has to be compressed or not
  bool needsCompress() const 
  {
    return indexSet_.needsCompress();
  }
};

} // end namespace Dune 
#endif
