#ifndef DUNE_DGADAPTMANAGERIMP_HH
#define DUNE_DGADAPTMANAGERIMP_HH

//- local includes  
#include <dune/fem/space/common/adaptmanager.hh>
#include <dune/fem/quadrature/cachequad.hh>

//- local includes 
#include "dgspace.hh"

/************************************************************
1) Gewichte zwischen Vater/Sohn default Implementieren auf Gitter
   (father.weight(son) oder so
2) Caching der Basisfunktionen fuer Vater-BasisFunktionen fuer
   Kinderquadraturen
*************************************************************/

namespace Dune{
/** @ingroup RestrictProlongImpl
    @{
**/

//***********************************************************************
/** \brief This is a restriction/prolongation operator for DG data. 
 */
template <class DiscreteFunctionImp, int polOrd> 
class RestrictProlongDiscontinuousSpace
: public RestrictProlongInterfaceDefault<RestrictProlongTraits< 
  RestrictProlongDiscontinuousSpace<DiscreteFunctionImp,polOrd> > >
{
  typedef RestrictProlongInterfaceDefault<RestrictProlongTraits< 
  RestrictProlongDiscontinuousSpace<DiscreteFunctionImp,polOrd> > > BaseType;

public:
  typedef DiscreteFunctionImp DiscreteFunctionType;
  typedef typename DiscreteFunctionType::FunctionSpaceType FunctionSpaceType;
  typedef typename FunctionSpaceType :: GridPartType GridPartType;
  typedef typename FunctionSpaceType :: GridType GridType;
  typedef typename DiscreteFunctionType::LocalFunctionType LocalFunctionType;

  typedef typename DiscreteFunctionType::RangeFieldType RangeFieldType;
  typedef typename DiscreteFunctionType::DomainType DomainType;
  typedef CachingQuadrature<GridPartType,0> QuadratureType;
  typedef typename GridType::template Codim<0>::Entity::LocalGeometry LocalGeometry;

protected:
  using BaseType :: calcWeight;
  using BaseType :: entitiesAreCopies;
  
public:  
  //! Constructor
  explicit RestrictProlongDiscontinuousSpace( DiscreteFunctionType &df )
  : df_( df ),
    quadord_( 2 * df.space().order() ),
    weight_( -1.0 )
  {}

  /** \brief explicit set volume ratio of son and father
   *
   *  \param[in]  weight  volume of son / volume of father
   *
   *  \note If this ratio is set, it is assume to be constant.
   */
  void setFatherChildWeight ( const RangeFieldType &weight ) const
  {
    weight_ = weight;
  }

  //! restrict data to father 
  template< class EntityType >
  void restrictLocal ( const EntityType &father, const EntityType &son, bool initialize ) const
  {
    // if father and son are copies, do nothing
    if( entitiesAreCopies( df_.space().indexSet(), father, son ) )
      return;
    
    typename FunctionSpaceType::RangeType ret (0.0);
    typename FunctionSpaceType::RangeType phi (0.0);
    assert( !father.isLeaf() );
    const RangeFieldType weight = (weight_ < 0.0) ? calcWeight( father, son ) : weight_;

    LocalFunctionType vati_ = df_.localFunction( father);
    LocalFunctionType sohn_ = df_.localFunction( son   );

    const typename FunctionSpaceType::BaseFunctionSetType & baseset =
      vati_.baseFunctionSet();
    const LocalGeometry& geometryInFather = son.geometryInFather();

    const int vati_numDofs = vati_.numDofs();
    if( initialize )
    {
      for(int i=0; i<vati_numDofs; ++i) 
        vati_[i] = 0.0;
    }
    
    QuadratureType quad(son,quadord_);
    const int nop = quad.nop();
    for( int qP = 0; qP < nop; ++qP )
    {
      sohn_.evaluate(quad[qP],ret);
      for(int i=0; i<vati_numDofs; ++i) 
      {
        baseset.evaluate(i,geometryInFather.global(quad.point(qP)),phi);
        vati_[i] += quad.weight(qP) * weight * (ret * phi) ;
      }
    }
  }

  //! prolong data to children 
  template< class EntityType >
  void prolongLocal ( const EntityType &father, const EntityType &son, bool initialize ) const
  {
    // if father and son are copies, do nothing
    if( entitiesAreCopies( df_.space().indexSet(), father, son ) )
      return;
    
    typename FunctionSpaceType::RangeType ret (0.0);
    typename FunctionSpaceType::RangeType phi (0.0);

    LocalFunctionType vati_ = df_.localFunction( father);
    LocalFunctionType sohn_ = df_.localFunction( son   );

    const int sohn_numDofs = sohn_.numDofs();
    for(int i=0; i<sohn_numDofs; ++i) sohn_[i] = 0.;

    const typename FunctionSpaceType::BaseFunctionSetType &baseset
      = sohn_.baseFunctionSet();
    const LocalGeometry& geometryInFather = son.geometryInFather();

    QuadratureType quad(son,quadord_);
    const int nop = quad.nop();
    for( int qP = 0; qP < nop; ++qP )
    {
      vati_.evaluate(geometryInFather.global(quad.point(qP)),ret);
      
      for( int i = 0; i < sohn_numDofs; ++i )
      {
        baseset.evaluate(i,quad[qP],phi);
        sohn_[i] += quad.weight(qP) * (ret * phi) ;
      }
    }
  }

  //! add discrete function to communicator 
  template <class CommunicatorImp>
  void addToList(CommunicatorImp& comm)
  {
    comm.addToList(df_);
  }

private:
  mutable DiscreteFunctionType & df_;
  const int quadord_;
  mutable RangeFieldType weight_;
};

/** \brief This is a restriction/prolongation operator for DG data of order zero. 
 */
  template <class DiscreteFunctionImp> 
class RestrictProlongDiscontinuousSpace<DiscreteFunctionImp,0> :
public RestrictProlongPieceWiseConstantData<DiscreteFunctionImp> 
{
public:
  typedef DiscreteFunctionImp DiscreteFunctionType;
  typedef RestrictProlongPieceWiseConstantData<DiscreteFunctionImp> BaseType; 
public:  
  //! Constructor
  RestrictProlongDiscontinuousSpace( DiscreteFunctionType & df ) : 
    BaseType ( df ) 
  {
  }
};

/** \brief specialization of RestrictProlongDefault for
    DiscontinuousGalerkinSpace.
*/
template <class DiscFunc,
          class FunctionSpaceImp, 
          class GridPartImp, 
          int polOrd, 
          template <class> class StorageImp> 
class RestrictProlongDefaultImplementation< 
  DiscFunc,
  DiscontinuousGalerkinSpace<FunctionSpaceImp, GridPartImp, polOrd,StorageImp> 
  > 
: public RestrictProlongDiscontinuousSpace<
  DiscFunc,polOrd >
{
public:
  //! type of discrete function 
  typedef DiscFunc DiscreteFunctionType;
  //! type of base class  
  typedef RestrictProlongDiscontinuousSpace<DiscreteFunctionType,polOrd> 
          BaseType;
public:  
  //! Constructor
  RestrictProlongDefaultImplementation ( DiscreteFunctionType & df ) : 
    BaseType(df) 
  {
  }
};

/** \brief specialization of RestrictProlongDefault for
    LegendreDiscontinuousGalerkinSpace.
*/
template <class DiscFunc,
          class FunctionSpaceImp, 
          class GridPartImp, 
          int polOrd, 
          template <class> class StorageImp> 
class RestrictProlongDefaultImplementation<DiscFunc,
 LegendreDiscontinuousGalerkinSpace<FunctionSpaceImp, GridPartImp, polOrd,StorageImp> > 
: public RestrictProlongDiscontinuousSpace<DiscFunc,polOrd >
{
public:
  //! type of discrete function 
  typedef DiscFunc DiscreteFunctionType;
  //! type of base class  
  typedef RestrictProlongDiscontinuousSpace<DiscreteFunctionType,polOrd> BaseType;
public:  
  //! Constructor
  RestrictProlongDefaultImplementation( DiscreteFunctionType & df ) : 
    BaseType(df) 
  {
  }
};

///@}
} // end namespace Dune 
#endif
