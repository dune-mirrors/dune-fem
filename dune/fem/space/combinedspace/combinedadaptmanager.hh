#ifndef DUNE_COMBINEDADAPTMANAGERIMP_HH
#define DUNE_COMBINEDADAPTMANAGERIMP_HH

//- local includes  
#include <dune/fem/space/common/adaptmanager.hh>
#include <dune/fem/quadrature/cachingquadrature.hh>
#include <dune/fem/space/common/restrictprolonginterface.hh>

//- local includes 
#include "combinedspace.hh"

//----------------------------------------------------------------------
//-   1) Gewichte zwischen Vater/Sohn default Implementieren auf Gitter
//-      (father.weight(son) oder so
//-   2) Caching der Basisfunktionen fuer Vater-BasisFunktionen fuer
//-      Kinderquadraturen
//----------------------------------------------------------------------

namespace Dune
{

/** @ingroup RestrictProlongImpl
    @{
**/

//***********************************************************************
/** \brief This is a restriction/prolongation operator for combined DG data. 
 */
template< class DiscreteFunctionImp, int polOrd >
class RestrictProlongCombinedSpace
: public RestrictProlongInterfaceDefault
  < RestrictProlongTraits< RestrictProlongCombinedSpace< DiscreteFunctionImp,polOrd > > >
{
  typedef RestrictProlongInterfaceDefault
    < RestrictProlongTraits< RestrictProlongCombinedSpace< DiscreteFunctionImp, polOrd > > >
    BaseType;

public:
  typedef DiscreteFunctionImp DiscreteFunctionType;

  typedef typename DiscreteFunctionType::DiscreteFunctionSpaceType
    DiscreteFunctionSpaceType;
  typedef typename DiscreteFunctionType::LocalFunctionType LocalFunctionType;

  typedef typename DiscreteFunctionSpaceType::GridPartType GridPartType;
  typedef typename DiscreteFunctionSpaceType::GridType GridType;

  typedef typename DiscreteFunctionSpaceType::RangeFieldType RangeFieldType;
  typedef typename DiscreteFunctionSpaceType::RangeType  RangeType;
  typedef typename DiscreteFunctionSpaceType::DomainType DomainType;

  typedef CachingQuadrature<GridPartType,0> QuadratureType;
  typedef typename GridType::template Codim<0>::Entity::LocalGeometry LocalGeometry;

  static const int dimRange = DiscreteFunctionSpaceType::dimRange;

protected:
  using BaseType :: calcWeight;
  using BaseType :: entitiesAreCopies;

public:  
  //! Constructor
  explicit RestrictProlongCombinedSpace( DiscreteFunctionType &df )
  : df_( df ),
    quadord_( 2*df.space().order() ),
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
    
    typedef typename DiscreteFunctionSpaceType::DomainFieldType DomainFieldType;

    assert( !father.isLeaf() );
    const DomainFieldType weight = (weight_ < 0.0) ? calcWeight( father, son) : weight_;

    LocalFunctionType vati = df_.localFunction( father);
    LocalFunctionType sohn = df_.localFunction( son   );

    if(initialize) 
    {
      vati.clear();
    }
    
    const LocalGeometry& geometryInFather = son.geometryInFather();

    RangeType value ;

    QuadratureType quad( son, quadord_);

    const int nop = quad.nop();
    for(int qP = 0; qP < nop; ++qP) 
    {
      sohn.evaluate( quad[ qP ], value );
      // calculate factor 
      const DomainFieldType intel = quad.weight(qP) * weight;

      // apply weight 
      value *= intel;

      // apply axpy on father 
      vati.axpy( geometryInFather.global(quad.point(qP)), value );
    }
  }

  //! prolong data to children 
  template <class EntityType>
  void prolongLocal ( EntityType &father, EntityType &son, bool initialize ) const
  {
    // if father and son are copies, do nothing 
    if( this->entitiesAreCopies( df_.space().indexSet(), father, son ) ) return ; 
    
    // get local functions 
    LocalFunctionType vati = df_.localFunction( father);
    LocalFunctionType sohn = df_.localFunction( son   );

    // set sohn to zero
    sohn.clear();

    // get quadrature 
    QuadratureType quad( son, quadord_ );

    // get geometry 
    const LocalGeometry& geometryInFather = son.geometryInFather();
    
    RangeType value ;

    // get number of points 
    const int nop = quad.nop();
    for(int qP = 0; qP < nop; ++qP) 
    {
      // evaluate father 
      vati.evaluate( geometryInFather.global( quad.point(qP) ), value );

      // apply weight 
      value *= quad.weight(qP);

      // add to son 
      sohn.axpy( quad[ qP ], value );
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
  int quadord_;
  mutable RangeFieldType weight_;
};

/** \brief This is a restriction/prolongation operator for 
    combined DG data with polynomial order 0. 
*/
template <class DiscreteFunctionImp> 
class RestrictProlongCombinedSpace<DiscreteFunctionImp,0>
: public RestrictProlongPieceWiseConstantData<DiscreteFunctionImp> 
{
public:
  typedef DiscreteFunctionImp DiscreteFunctionType;
  typedef RestrictProlongPieceWiseConstantData<DiscreteFunctionImp>  BaseType; 
public:  
  //! Constructor
  RestrictProlongCombinedSpace( DiscreteFunctionType & df ) :
    BaseType(df)
  {}
};

/** \brief specialization of RestrictProlongDefault for
    CombinedSpace.
*/
template <class DiscFunc,
          class DiscreteFunctionSpaceImp, 
          int N, 
          DofStoragePolicy policy> 
class RestrictProlongDefaultImplementation< 
  DiscFunc,CombinedSpace<DiscreteFunctionSpaceImp,N,policy> > 
: public RestrictProlongCombinedSpace<
  DiscFunc,DiscreteFunctionSpaceImp :: polynomialOrder >
{
public:
  //! type of discrete function 
  typedef DiscFunc DiscreteFunctionType;
  //! type of base class  
  typedef RestrictProlongCombinedSpace<
      DiscFunc, DiscreteFunctionSpaceImp :: polynomialOrder > BaseType;
public:  
  //! Constructor
  RestrictProlongDefaultImplementation ( DiscreteFunctionType & df ) : 
    BaseType(df) 
  {
  }
};

///@}
} // end namespace Dune 
#endif
