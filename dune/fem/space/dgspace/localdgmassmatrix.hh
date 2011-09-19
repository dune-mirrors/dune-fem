#ifndef DUNE_FEM_LOCALDGMASSMATRIX_HH
#define DUNE_FEM_LOCALDGMASSMATRIX_HH

#include <dune/fem/space/dgspace.hh>
#include <dune/fem/operator/1order/localmassmatrix.hh>

namespace Dune {

/** \brief Local Mass Matrix for DG space */
template <class FunctionSpaceImp, class GridPartImp, int polOrd,
          template<class> class BaseFunctionStorageImp, 
          class VolumeQuadratureImp>
class LocalMassMatrix< 
  DiscontinuousGalerkinSpace< FunctionSpaceImp, GridPartImp, polOrd, BaseFunctionStorageImp >,
  VolumeQuadratureImp > 
  : public LocalMassMatrixImplementationDgOrthoNormal< 
      DiscontinuousGalerkinSpace< FunctionSpaceImp, GridPartImp, polOrd, BaseFunctionStorageImp >, VolumeQuadratureImp >
{
  typedef DiscontinuousGalerkinSpace< FunctionSpaceImp, GridPartImp, polOrd, BaseFunctionStorageImp > DiscreteFunctionSpaceImp;
  typedef LocalMassMatrixImplementationDgOrthoNormal< DiscreteFunctionSpaceImp, VolumeQuadratureImp > BaseType;
public:
  LocalMassMatrix( const DiscreteFunctionSpaceImp& spc, const int volQuadOrd = -1 )
    : BaseType( spc, volQuadOrd )
  {}
};

/** \brief Local Mass Matrix for Legendre space */
template <class FunctionSpaceImp, class GridPartImp, int polOrd,
          template<class> class BaseFunctionStorageImp, 
          class VolumeQuadratureImp>
class LocalMassMatrix< 
  LegendreDiscontinuousGalerkinSpace< FunctionSpaceImp, GridPartImp, polOrd, BaseFunctionStorageImp >,
  VolumeQuadratureImp > 
  : public LocalMassMatrixImplementationDgOrthoNormal< 
      LegendreDiscontinuousGalerkinSpace< FunctionSpaceImp, GridPartImp, polOrd, BaseFunctionStorageImp >, VolumeQuadratureImp >
{
  typedef LegendreDiscontinuousGalerkinSpace< FunctionSpaceImp, GridPartImp, polOrd, BaseFunctionStorageImp > DiscreteFunctionSpaceImp;
  typedef LocalMassMatrixImplementationDgOrthoNormal< DiscreteFunctionSpaceImp, VolumeQuadratureImp > BaseType;
public:
  LocalMassMatrix( const DiscreteFunctionSpaceImp& spc, const int volQuadOrd = -1 )
    : BaseType( spc, volQuadOrd )
  {}
};

} // end namespace Dune 
#endif // end DUNE_FEM_LOCALDGMASSMATRIX_HH
