#ifndef DUNE_FVSPACE_INLINE_HH
#define DUNE_FVSPACE_INLINE_HH

//- system includes 
#include <algorithm>

//- Dune includes 
#include <dune/fem/space/common/allgeomtypes.hh>

namespace Dune {

template <class FunctionSpaceImp, class GridPartImp, int polOrd, template <class> class BaseFunctionStorageImp >
inline FiniteVolumeSpace<FunctionSpaceImp, GridPartImp, polOrd, BaseFunctionStorageImp>::
FiniteVolumeSpace (GridPartType & gridPart) :
    DefaultType(gridPart),
    baseFuncSet_(),
    gridPart_(gridPart), 
    mapper_(0),
    dm_(DofManagerFactoryType::getDofManager(gridPart.grid()))
{
  makeFunctionSpace(gridPart);
}

template <class FunctionSpaceImp, class GridPartImp, int polOrd, template <class> class BaseFunctionStorageImp >
inline void FiniteVolumeSpace<FunctionSpaceImp, GridPartImp, polOrd, BaseFunctionStorageImp>::
makeFunctionSpace (GridPartType& gridPart) 
{
  // add index set to list of indexset of dofmanager 
  typedef DofManager<GridType> DofManagerType;
  typedef DofManagerFactory<DofManagerType> DofManagerFactoryType;
  DofManagerType & dm = DofManagerFactoryType::getDofManager(gridPart.grid());
  dm.addIndexSet(gridPart.grid(), 
      const_cast<typename GridPartType::IndexSetType&>(gridPart.indexSet()));
  
  // search the macro grid for diffrent element types 
  AllGeomTypes< typename GridPartType::IndexSetType ,
                typename GridPartType::GridType > 
            allGeomTypes(gridPart.indexSet());
  
  const std::vector<GeometryType>& geomTypes = allGeomTypes.geomTypes(0);

  for(size_t i=0; i<geomTypes.size(); ++i)
  {
    const GeometryType & geoType = geomTypes[i];
    if(baseFuncSet_.find( geoType ) == baseFuncSet_.end() ) 
    {
      // get pointer to base function set 
      const BaseFunctionSetImp* baseSet = 
        & SingletonProviderType::getObject(geoType);
      
      baseFuncSet_[ geoType ] = baseSet; 
      const int numDofs = baseSet->numBaseFunctions();
      
      MapperSingletonKeyType key(gridPart.indexSet(),numDofs);
      mapper_ = & MapperProviderType::getObject(key);

      // make sure we got the right mapper 
      assert( mapper_->numDofs() == numDofs );
    }
  }
}
  
template <class FunctionSpaceImp, class GridPartImp, int polOrd, template <class> class BaseFunctionStorageImp >
inline FiniteVolumeSpace<FunctionSpaceImp, GridPartImp, polOrd, BaseFunctionStorageImp>::
~FiniteVolumeSpace () 
{
  typedef typename BaseFunctionMapType :: iterator iterator;
  iterator end = baseFuncSet_.end();
  for(iterator it = baseFuncSet_.begin(); it != end; ++it)
  {
    const BaseFunctionSetImp * set = (*it).second; 
    SingletonProviderType::removeObject(*set);
  }

  baseFuncSet_.clear();
  MapperProviderType::removeObject(*mapper_);
}  

template <class FunctionSpaceImp, class GridPartImp, int polOrd, template <class> class BaseFunctionStorageImp >
inline DFSpaceIdentifier FiniteVolumeSpace<FunctionSpaceImp, GridPartImp, polOrd, BaseFunctionStorageImp>::
type () const 
{
  return FiniteVolumeSpace_id;
}

template <class FunctionSpaceImp, class GridPartImp, int polOrd, template <class> class BaseFunctionStorageImp >
template <class EntityType> 
inline const 
typename FiniteVolumeSpace<FunctionSpaceImp, GridPartImp, polOrd, BaseFunctionStorageImp>::BaseFunctionSetType  
FiniteVolumeSpace<FunctionSpaceImp, GridPartImp, polOrd, BaseFunctionStorageImp>::
baseFunctionSet (const EntityType &en) const 
{
  return this->baseFunctionSet(en.geometry().type());
}

template <class FunctionSpaceImp, class GridPartImp, int polOrd, template <class> class BaseFunctionStorageImp >
inline const 
typename FiniteVolumeSpace<FunctionSpaceImp, GridPartImp, polOrd, BaseFunctionStorageImp>::BaseFunctionSetType  
FiniteVolumeSpace<FunctionSpaceImp, GridPartImp, polOrd, BaseFunctionStorageImp>::
baseFunctionSet (const GeometryType& geomType) const 
{
  assert( baseFuncSet_.find(geomType) != baseFuncSet_.end());
  return BaseFunctionSetType(baseFuncSet_[geomType]);
}

template <class FunctionSpaceImp, class GridPartImp, int polOrd, template <class> class BaseFunctionStorageImp >
inline int FiniteVolumeSpace<FunctionSpaceImp, GridPartImp, polOrd, BaseFunctionStorageImp>::
size () const
{
  assert( mapper_ );
  return mapper_->size ();
}

template <class FunctionSpaceImp, class GridPartImp, int polOrd, template <class> class BaseFunctionStorageImp >
template< class EntityType> 
inline int FiniteVolumeSpace<FunctionSpaceImp, GridPartImp, polOrd, BaseFunctionStorageImp>::
mapToGlobal ( EntityType &en, int localNum ) const
{
  assert( mapper_ );
  return mapper_->mapToGlobal ( en , localNum );
}

template <class FunctionSpaceImp, class GridPartImp, int polOrd, template <class> class BaseFunctionStorageImp >
typename FiniteVolumeSpace<FunctionSpaceImp, GridPartImp, polOrd, BaseFunctionStorageImp>::MapperType&
FiniteVolumeSpace<FunctionSpaceImp, GridPartImp, polOrd, BaseFunctionStorageImp>::mapper() const {
  assert(mapper_);
  return *mapper_;
}
   
} // end namespace Dune 
#endif
