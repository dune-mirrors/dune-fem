#ifndef DUNE_FEM_SPACE_DISCONTINUOUSGALERKIN_LOCALDGMASSMATRIX_HH
#define DUNE_FEM_SPACE_DISCONTINUOUSGALERKIN_LOCALDGMASSMATRIX_HH

// dune-fem includes
#include <dune/fem/operator/1order/localmassmatrix.hh>

// local includes
#include "declaration.hh"


namespace Dune
{

  namespace Fem 
  {

    // LocalMassMatrix
    // ---------------

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

  } // namespace Fem 

} // namespace Dune 

#endif // #ifndef DUNE_FEM_SPACE_DISCONTINUOUSGALERKIN_LOCALDGMASSMATRIX_HH
