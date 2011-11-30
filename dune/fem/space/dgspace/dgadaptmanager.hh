#ifndef DUNE_DGADAPTMANAGERIMP_HH
#define DUNE_DGADAPTMANAGERIMP_HH

//- local includes  
#include <dune/fem/space/common/adaptmanager.hh>
#include <dune/fem/quadrature/cachingquadrature.hh>
#include <dune/fem/function/localfunction/temporarylocalfunction.hh>

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
  RestrictProlongDiscontinuousSpace<DiscreteFunctionImp,polOrd>,
  typename DiscreteFunctionImp::DiscreteFunctionSpaceType::GridType::ctype > >
{
  typedef RestrictProlongInterfaceDefault<RestrictProlongTraits< 
    RestrictProlongDiscontinuousSpace<DiscreteFunctionImp,polOrd>,
    typename DiscreteFunctionImp::DiscreteFunctionSpaceType::GridType::ctype > > BaseType;

  typedef typename BaseType::DomainFieldType DomainFieldType;

public:
  typedef DiscreteFunctionImp DiscreteFunctionType;
  typedef typename DiscreteFunctionType::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
  typedef typename DiscreteFunctionSpaceType :: GridPartType GridPartType;
  typedef typename DiscreteFunctionSpaceType :: GridType GridType;

  typedef typename DiscreteFunctionType::LocalFunctionType LocalFunctionType;
  typedef ConstLocalFunction< DiscreteFunctionType > ConstLocalFunctionType;

  typedef typename DiscreteFunctionSpaceType :: RangeType RangeType;
  typedef typename DiscreteFunctionSpaceType :: DomainType DomainType;
  typedef CachingQuadrature<GridPartType,0> QuadratureType;

  // grid part entity 
  typedef typename GridPartType::template Codim<0>::EntityType  EntityType ;
  // grid entity 
  typedef typename GridType::template Codim<0>::Entity  GridEntityType ;

  typedef typename EntityType::LocalGeometry LocalGeometry;
protected:
  using BaseType :: calcWeight;
  using BaseType :: entitiesAreCopies;
  
public:  
  //! Constructor
  explicit RestrictProlongDiscontinuousSpace( DiscreteFunctionType &df )
  : df_( df ),
    constFct_( df_ ),
    quadord_( 2 * df.space().order() ),
    weight_( -1.0 )
  {}

  /** \brief explicit set volume ratio of son and father
   *
   *  \param[in]  weight  volume of son / volume of father
   *
   *  \note If this ratio is set, it is assume to be constant.
   */
  void setFatherChildWeight ( const DomainFieldType &weight ) const
  {
    weight_ = weight;
  }

  //! restrict data to father 
  void restrictLocal ( const GridEntityType &dad, const GridEntityType &filius, bool initialize ) const
  {
    // convert from grid entities to grid part entities 
    const EntityType& father = df_.space().gridPart().convert( dad ); 
    const EntityType& son    = df_.space().gridPart().convert( filius ); 
      
    // if father and son are copies, do nothing
    if( entitiesAreCopies( df_.space().indexSet(), father, son ) )
      return;

    //std::cout << "Restrict " << df_.space().indexSet().index( father ) << "  "
    //          << df_.space().indexSet().index( son ) << std::endl;
    
    assert( !father.isLeaf() );

    // get local function of son as a constant function
    // in case father and son are not children (GeoGridPart)
    constFct_.init( son );

    // get father function 
    LocalFunctionType fatherLf = df_.localFunction( father);

    restrictLocal( son.geometryInFather(), constFct_, fatherLf, initialize );
  }

  //! restrict data to father 
  template < class LocalGeom, class LFSon, class LFFather>
  void restrictLocal ( const LocalGeom& geometryInFather,
                       const LFSon& sonLf, LFFather& fatherLf, bool initialize ) const 
  {
    const DomainFieldType weight = (weight_ < 0.0) ? calcWeight( fatherLf.entity(), sonLf.entity() ) : weight_;

    // clear father on initialize 
    if( initialize )
    {
      fatherLf.clear();
    }
    
    RangeType value ;
    QuadratureType quad( sonLf.entity(), quadord_ );
    const int nop = quad.nop();
    for( int qP = 0; qP < nop; ++qP )
    {
      // evaluate son function 
      sonLf.evaluate( quad[ qP ], value );

      // apply weight 
      value *= quad.weight( qP ) * weight ;

      // add to father 
      fatherLf.axpy( geometryInFather.global(quad.point(qP) ), value );
    }
  }

  //! prolong data to children 
  void prolongLocal ( const GridEntityType &dad, const GridEntityType &filius, bool initialize ) const
  {
    // convert from grid entities to grid part entities 
    const EntityType& father = df_.space().gridPart().convert( dad ); 
    const EntityType& son    = df_.space().gridPart().convert( filius ); 
      
    // if father and son are copies, do nothing
    if( entitiesAreCopies( df_.space().indexSet(), father, son ) )
      return;
    
    //std::cout << "Prolong  " << df_.space().indexSet().index( father ) << "  "
    //          << df_.space().indexSet().index( son ) << std::endl;
    
    constFct_.init( father );
    LocalFunctionType sonLf = df_.localFunction( son );

    // call prolongLocal method 
    prolongLocal( son.geometryInFather(), constFct_, sonLf, initialize );
  }

  //! prolong data to children 
  template < class LocalGeom, class LFFather, class LFSon>
  void prolongLocal ( const LocalGeom& geometryInFather,
                      const LFFather& fatherLf, LFSon& sonLf, bool initialize ) const 
  {
    // clear son 
    sonLf.clear();

    RangeType value ;
    QuadratureType quad( sonLf.entity(), quadord_ );
    const int nop = quad.nop();
    for( int qP = 0; qP < nop; ++qP )
    {
      // evaluate father 
      fatherLf.evaluate(geometryInFather.global(quad.point(qP)), value );
      
      // apply weight 
      value *= quad.weight( qP );

      // add to son 
      sonLf.axpy( quad[ qP ], value );
    }
  }

  //! add discrete function to communicator 
  template< class CommunicatorImp >
  void addToList ( CommunicatorImp &comm )
  {
    comm.addToList( df_ );
  }

protected:
  DiscreteFunctionType &df_;
  mutable ConstLocalFunctionType constFct_;
  const int quadord_;
  mutable DomainFieldType weight_;
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

/** \brief specialization of RestrictProlongDefault for
    LagrangeDiscontinuousGalerkinSpace.
*/
template <class DiscFunc,
          class FunctionSpaceImp, 
          class GridPartImp, 
          int polOrd, 
          template <class> class StorageImp> 
class RestrictProlongDefaultImplementation<DiscFunc,
 LagrangeDiscontinuousGalerkinSpace<FunctionSpaceImp, GridPartImp, polOrd,StorageImp> > 
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
