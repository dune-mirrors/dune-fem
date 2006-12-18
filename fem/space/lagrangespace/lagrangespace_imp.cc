#ifndef DUNE_LAGRANGESPACE_CC
#define DUNE_LAGRANGESPACE_CC

//- system includes 
#include <algorithm>

//- Dune includes 
#include <dune/fem/space/common/allgeomtypes.hh>

namespace Dune {

template <class FunctionSpaceImp, class GridPartImp, int polOrd, template <class> class BaseFunctionStorageImp >
inline LagrangeDiscreteFunctionSpace<FunctionSpaceImp, GridPartImp, polOrd, BaseFunctionStorageImp>::
LagrangeDiscreteFunctionSpace (GridPartType & gridPart) :
    DefaultType(gridPart),
    baseFuncSet_(GeometryIdentifier::numTypes,0),
    grid_(gridPart), 
    mapper_(0),
    dm_(DofManagerFactoryType::getDofManager(gridPart.grid()))
{
  makeFunctionSpace(gridPart);
}

template <class FunctionSpaceImp, class GridPartImp, int polOrd, template <class> class BaseFunctionStorageImp >
inline void LagrangeDiscreteFunctionSpace<FunctionSpaceImp, GridPartImp, polOrd, BaseFunctionStorageImp>::
makeFunctionSpace (GridPartType& gridPart) 
{
  // add index set to list of indexset of dofmanager 
  typedef DofManager<GridType> DofManagerType;
  typedef DofManagerFactory<DofManagerType> DofManagerFactoryType;
  DofManagerType & dm = DofManagerFactoryType::getDofManager(gridPart.grid());
  dm.addIndexSet(gridPart.grid(), 
      const_cast<typename GridPartType::IndexSetType&>(gridPart.indexSet()));
  
  //std::cout << "Constructor of LagrangeDiscreteFunctionSpace! \n";
  // search the macro grid for diffrent element types 
  AllGeomTypes< typename GridPartType::IndexSetType ,
                typename GridPartType::GridType > 
            allGeomTypes(gridPart.indexSet());
  
  const std::vector<GeometryType>& geomTypes = allGeomTypes.geomTypes(0);

  for(size_t i=0; i<geomTypes.size(); ++i)
  {
    GeometryIdentifier::IdentifierType id =
                GeometryIdentifier::fromGeo(geomTypes[i]);

    if(baseFuncSet_[id] == 0) 
    {
      baseFuncSet_[id] = & setBaseFuncSetPointer(geomTypes[i]);
      int numDofs = baseFuncSet_[id]->numBaseFunctions();
      
      MapperSingletonKeyType key(gridPart.indexSet(),numDofs);
      mapper_ = & MapperProviderType::getObject(key);

      // make sure we got the right mapper 
      assert( mapper_->numDofs() == numDofs );
    }
  }
}
  
template <class FunctionSpaceImp, class GridPartImp, int polOrd, template <class> class BaseFunctionStorageImp >
inline LagrangeDiscreteFunctionSpace<FunctionSpaceImp, GridPartImp, polOrd, BaseFunctionStorageImp>::
~LagrangeDiscreteFunctionSpace () 
{
  for(unsigned int i=0; i<baseFuncSet_.size(); i++)
  {
    if (baseFuncSet_[i] != 0)
    {
      BaseFunctionSetType * set = (BaseFunctionSetType *) baseFuncSet_[i]; 
      SingletonProviderType::removeObject(*set);
    }
  }

  MapperProviderType::removeObject(*mapper_);
}  

template <class FunctionSpaceImp, class GridPartImp, int polOrd, template <class> class BaseFunctionStorageImp >
inline DFSpaceIdentifier LagrangeDiscreteFunctionSpace<FunctionSpaceImp, GridPartImp, polOrd, BaseFunctionStorageImp>::
type () const 
{
  return LagrangeSpace_id;
}

template <class FunctionSpaceImp, class GridPartImp, int polOrd, template <class> class BaseFunctionStorageImp >
template <class EntityType> 
inline const 
typename LagrangeDiscreteFunctionSpace<FunctionSpaceImp, GridPartImp, polOrd, BaseFunctionStorageImp>::BaseFunctionSetType &  
LagrangeDiscreteFunctionSpace<FunctionSpaceImp, GridPartImp, polOrd, BaseFunctionStorageImp>::
baseFunctionSet (const EntityType &en) const 
{
  GeometryIdentifier::IdentifierType id =
    GeometryIdentifier::fromGeometry(en.geometry());
  return this->baseFunctionSet(id);
}

template <class FunctionSpaceImp, class GridPartImp, int polOrd, template <class> class BaseFunctionStorageImp >
inline const 
typename LagrangeDiscreteFunctionSpace<FunctionSpaceImp, GridPartImp, polOrd, BaseFunctionStorageImp>::BaseFunctionSetType &  
LagrangeDiscreteFunctionSpace<FunctionSpaceImp, GridPartImp, polOrd, BaseFunctionStorageImp>::
baseFunctionSet (const GeometryIdentifier::IdentifierType id) const 
{
  assert(id < (int) baseFuncSet_.size());
  assert(id >= 0);
  assert(baseFuncSet_[id]);
  return *baseFuncSet_[id];
}

template <class FunctionSpaceImp, class GridPartImp, int polOrd, template <class> class BaseFunctionStorageImp >
inline int LagrangeDiscreteFunctionSpace<FunctionSpaceImp, GridPartImp, polOrd, BaseFunctionStorageImp>::
size () const
{
  return mapper_->size ();
}

template <class FunctionSpaceImp, class GridPartImp, int polOrd, template <class> class BaseFunctionStorageImp >
template< class EntityType> 
inline int LagrangeDiscreteFunctionSpace<FunctionSpaceImp, GridPartImp, polOrd, BaseFunctionStorageImp>::
mapToGlobal ( EntityType &en, int localNum ) const
{
  return mapper_->mapToGlobal ( en , localNum );
}

template <class FunctionSpaceImp, class GridPartImp, int polOrd, template <class> class BaseFunctionStorageImp >
const typename LagrangeDiscreteFunctionSpace<FunctionSpaceImp, GridPartImp, polOrd, BaseFunctionStorageImp>::MapperType&
LagrangeDiscreteFunctionSpace<FunctionSpaceImp, GridPartImp, polOrd, BaseFunctionStorageImp>::mapper() const {
  assert(mapper_);
  return *mapper_;
}
   
template <class FunctionSpaceImp, class GridPartImp, int polOrd, template <class> class BaseFunctionStorageImp >
typename LagrangeDiscreteFunctionSpace<FunctionSpaceImp, GridPartImp, polOrd, BaseFunctionStorageImp>::
BaseFunctionSetType& 
LagrangeDiscreteFunctionSpace<FunctionSpaceImp, GridPartImp, polOrd, BaseFunctionStorageImp>::
setBaseFuncSetPointer(GeometryType type) 
{
  return SingletonProviderType::getObject(type);
}

} // end namespace Dune 
#endif
