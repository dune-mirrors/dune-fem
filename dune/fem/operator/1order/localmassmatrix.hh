#ifndef DUNE_FEM_LOCALMASSMATRIX_HH
#define DUNE_FEM_LOCALMASSMATRIX_HH

//- Dune includes 
#include <dune/common/fmatrix.hh>
#include <dune/common/dynmatrix.hh>
#include <dune/common/geometrytypeindex.hh>

#include <dune/fem/version.hh>

#include <dune/fem/quadrature/cachingquadrature.hh>
#include <dune/fem/space/common/allgeomtypes.hh>
#include <dune/fem/misc/checkgeomaffinity.hh>

namespace Dune {

/*! @addtogroup PassHyp
 *  ** @{
*/

/** \brief Local Mass Matrix inversion implementation, select the correct method in your
    implementation */ 
template <class DiscreteFunctionSpaceImp, class VolumeQuadratureImp> 
class LocalMassMatrixImplementation  
{
public:  
  typedef DiscreteFunctionSpaceImp DiscreteFunctionSpaceType;
  typedef typename DiscreteFunctionSpaceType :: RangeFieldType ctype;
  typedef typename DiscreteFunctionSpaceType :: RangeFieldType RangeFieldType;
  typedef typename DiscreteFunctionSpaceType :: RangeType RangeType;

  enum { dgNumDofs = DiscreteFunctionSpaceType :: localBlockSize };

  typedef FieldMatrix<ctype, dgNumDofs, dgNumDofs > DGMatrixType;
  typedef FieldVector<ctype, dgNumDofs >            DGVectorType;

  typedef typename DiscreteFunctionSpaceType :: GridPartType GridPartType;
  typedef typename DiscreteFunctionSpaceType :: BaseFunctionSetType BaseFunctionSetType; 
  typedef typename GridPartType :: GridType GridType;
  typedef typename GridType :: template Codim<0> :: Entity EntityType;
  typedef typename GridType :: template Codim<0> :: Geometry Geometry;

  typedef VolumeQuadratureImp VolumeQuadratureType;

  //! is true if grid is structured grid 
  enum { StructuredGrid = Capabilities::isCartesian< GridType >::v };

  typedef AllGeomTypes< typename GridPartType :: IndexSetType,GridType> GeometryInformationType;
  typedef typename GeometryInformationType :: DomainType DomainType;
  typedef DynamicMatrix< RangeFieldType > MatrixType; 

protected:  
  const DiscreteFunctionSpaceType& spc_;

  GeometryInformationType geoInfo_;
  const int volumeQuadOrd_;
  const bool affine_;

  mutable DGMatrixType dgMatrix_;
  mutable DGVectorType dgX_, dgRhs_;

  mutable DynamicVector< RangeFieldType > rhs_;
  mutable DynamicVector< RangeFieldType > x_;
  mutable MatrixType matrix_;

  mutable std::vector< RangeType > phi_;
  mutable std::vector< RangeType > phiMass_;

  typedef std::map< const int, MatrixType* > MassMatrixStorageType ;
  typedef std::vector< MassMatrixStorageType > LocalInverseMassMatrixStorageType;

  mutable LocalInverseMassMatrixStorageType localInverseMassMatrix_;

  //! dummy caller 
  struct NoMassDummyCaller
  {
    enum { dimRange = DiscreteFunctionSpaceType::dimRange };
    typedef FieldMatrix<ctype, dimRange, dimRange> MassFactorType;
    // return false since we don;t have a mass term
    bool hasMass() const { return false; }
    void mass(const EntityType&,
              const VolumeQuadratureType&,
              const int, 
              MassFactorType&) const
    {
    }
  };

  template <class BaseFunctionSetType> 
  MatrixType& getLocalMassMatrix(const EntityType& en, 
                                 const Geometry& geo,
                                 const BaseFunctionSetType& baseSet,
                                 const int numBaseFct ) const 
  {
    const GeometryType geomType = geo.type();
    typedef typename MassMatrixStorageType :: iterator iterator ;
    MassMatrixStorageType& massMap = localInverseMassMatrix_[ GlobalGeometryTypeIndex :: index( geomType ) ];

    iterator it = massMap.find( numBaseFct ); 
    if( it == massMap.end() ) 
    {
      MatrixType* matrix = new MatrixType( numBaseFct, numBaseFct, 0.0 );
      massMap[ numBaseFct ] = matrix ;

      // create quadrature 
      VolumeQuadratureType volQuad(en, volumeQuadOrd_ );

      // setup matrix 
      buildMatrixNoMassFactor(en, geo, baseSet, volQuad, numBaseFct, *matrix, false );

      return *matrix;
    }

    return *((*it).second );
  }

public:
  //! constructor taking space and volume quadrature order 
  LocalMassMatrixImplementation(const DiscreteFunctionSpaceType& spc, const int volQuadOrd = -1 ) 
    : spc_(spc) 
    , geoInfo_( spc.indexSet() ) 
    , volumeQuadOrd_ ( ( volQuadOrd < 0 ) ? ( spc.order() * 2 ) : volQuadOrd )
    , affine_ ( setup() )
    , rhs_(), x_(), matrix_() 
    , phi_( spc_.mapper().maxNumDofs() )
    , phiMass_( spc_.mapper().maxNumDofs() )
    , localInverseMassMatrix_( GlobalGeometryTypeIndex :: size( GridType::dimension ) )
  {
  }

  //! copy constructor 
  LocalMassMatrixImplementation(const LocalMassMatrixImplementation& org) 
    : spc_(org.spc_),
      geoInfo_( org.geoInfo_),
      volumeQuadOrd_( org.volumeQuadOrd_ ),
      affine_( org.affine_ ),
      rhs_( org.rhs_ ), x_( org.x_ ), matrix_( org.matrix_ ),
      phi_( org.phi_ ),
      phiMass_( org.phiMass_ ),
      localInverseMassMatrix_( GlobalGeometryTypeIndex :: size( GridType::dimension ) )
  {
  }

  ~LocalMassMatrixImplementation()
  {
    typedef typename MassMatrixStorageType :: iterator iterator ;
    for (int i=0;i< localInverseMassMatrix_.size();++i)
    {
      const iterator end = localInverseMassMatrix_[i].end();
      for ( iterator it =  localInverseMassMatrix_[i].begin();it != end; ++it )
        delete (*it).second;
    }
  }

public:  
  //! returns true if geometry mapping is affine 
  bool affine () const { return affine_; }

  //! return mass factor for diagonal mass matrix 
  double getAffineMassFactor(const Geometry& geo) const 
  {
    return geoInfo_.referenceVolume( geo.type() ) / geo.volume(); 
  }

  //! apply local dg mass matrix to local function lf
  //! using the massFactor method of the caller 
  template <class MassCallerType, class LocalFunctionType> 
  void applyInverse(const MassCallerType& caller, 
                    const EntityType& en, 
                    LocalFunctionType& lf) const 
  {
    // get geometry 
    const Geometry& geo = en.geometry();

    if( geo.affine() && ! caller.hasMass() ) 
    {
      applyInverseLocally( caller, en, geo, lf );
    }
    else 
    {
      applyInverseDefault( caller, en, geo, lf );
    }
  }

  //! apply local dg mass matrix to local function lf without mass factor 
  template <class LocalFunctionType> 
  void applyInverse(const EntityType& en, 
                    LocalFunctionType& lf) const 
  {
    NoMassDummyCaller caller;
    applyInverse(caller, en, lf );
  }

  //! apply local dg mass matrix to local function lf without mass factor 
  template <class LocalFunctionType> 
  void applyInverse(LocalFunctionType& lf) const 
  {
    applyInverse( lf.entity(), lf );
  }

  /////////////////////////////////////////////
  // end of public methods 
  /////////////////////////////////////////////
  
  
protected:  
  ///////////////////////////////////////////////////////////
  //  applyInverse for DG spaces 
  ///////////////////////////////////////////////////////////
  template <class MassCallerType, class LocalFunctionType> 
  void applyInverseDgOrthoNormalBasis(const MassCallerType& caller, 
                                      const EntityType& en, 
                                      LocalFunctionType& lf) const 
  {
    // get geometry 
    const Geometry& geo = en.geometry();
    assert( dgNumDofs == lf.numDofs() );

    // affine_ can be a static information 
    const bool isAffine = ( affine_ ) ? true : geo.affine() ;

    // in case of affine mappings we only have to multiply with a factor 
    if( isAffine && ! caller.hasMass() )
    {
      const double massVolInv = getAffineMassFactor( geo );

      // apply inverse mass matrix  
      for(int l=0; l<dgNumDofs; ++l) 
      {
        lf[l] *= massVolInv;
      }

      return; 
    }
    else 
    {
      // copy local function to right hand side 
      for(int l=0; l<dgNumDofs; ++l) 
      {
        // copy 
        dgRhs_[ l ] = lf[ l ];
      }

      // setup local mass matrix 
      buildMatrix( caller, en, geo, lf.baseFunctionSet(), dgNumDofs, dgMatrix_ );

      // solve linear system  
      dgMatrix_.solve( dgX_, dgRhs_ );

      // copy back to local function 
      for(int l=0; l<dgNumDofs; ++l)
      {
        lf[ l ] = dgX_[ l ];
      }

      return; 
    }
  }

  ///////////////////////////////////////////////////////////
  //  standard applyInverse method 
  ///////////////////////////////////////////////////////////
  //! apply local mass matrix to local function lf
  //! using the massFactor method of the caller 
  template <class MassCallerType, class LocalFunctionType> 
  void applyInverseDefault(MassCallerType& caller, 
                           const EntityType& en, 
                           const Geometry& geo,
                           LocalFunctionType& lf) const 
  {
    const int numDofs = lf.numDofs();

    rhs_.resize( numDofs );
    x_.resize( numDofs );

    // copy local function to right hand side 
    for(int l=0; l<numDofs; ++l) 
    {
      // copy 
      rhs_[ l ] = lf[ l ];
    }

    // resize matrix 
    matrix_.resize( numDofs, numDofs );

    // setup local mass matrix 
    buildMatrix( caller, en, geo, lf.baseFunctionSet(), numDofs, matrix_ );

    // solve linear system  
    matrix_.solve( x_, rhs_ );

    for(int l=0; l<numDofs; ++l) 
    {
      // copy back
      lf[ l ] = x_[ l ];
    }
    return; 
  }

  ///////////////////////////////////////////////////////////
  //  local applyInverse method for affine geometries 
  ///////////////////////////////////////////////////////////
  //! apply local mass matrix to local function lf
  //! using the massFactor method of the caller 
  template <class MassCallerType, class LocalFunctionType> 
  void applyInverseLocally(MassCallerType& caller, 
                           const EntityType& en, 
                           const Geometry& geo,
                           LocalFunctionType& lf) const 
  {
    const int numDofs = lf.numDofs();

    // get local inverted mass matrix 
    MatrixType& invMassMatrix = 
      getLocalMassMatrix( en, geo, lf.baseFunctionSet(), numDofs );

    // resize vectors 
    rhs_.resize( numDofs );
    x_.resize( numDofs );

    const double massVolInv = getAffineMassFactor( geo );

    // copy local function to right hand side 
    for(int l=0; l<numDofs; ++l) 
    {
      // copy 
      rhs_[ l ] = lf[ l ];
    }

    // apply local inverse mass matrix
    invMassMatrix.mv( rhs_, x_ );

    for(int l=0; l<numDofs; ++l) 
    {
      // copy back
      lf[ l ] = x_[ l ] * massVolInv ;
    }

    return; 
  }

  //! setup and return affinity 
  bool setup() const 
  {
    // for structured grids this is always true 
    if( StructuredGrid ) return true;

    // get types for codim 0 
    const std::vector<Dune::GeometryType>& geomTypes = geoInfo_.geomTypes(0);

    // for simplices we also have affine mappings 
    if( (geomTypes.size() == 1) && geomTypes[0].isSimplex() ) 
    {
      return true;
    }

    bool affine = checkAffinity();
    //std::cout << "found affine grid? = " << affine << "\n";
    return affine;
  }

  //! check whether all geometry mappings are affine 
  bool checkAffinity() const 
  {
    return GeometryAffinityCheck<VolumeQuadratureType> :: 
      checkAffinity( spc_.begin(), spc_.end(), volumeQuadOrd_);
  }

  //! build local mass matrix 
  template <class MassCallerType, class Matrix> 
  void buildMatrix(const MassCallerType& caller,
                   const EntityType& en,
                   const Geometry& geo, 
                   const BaseFunctionSetType& set,
                   const int numDofs,
                   Matrix& matrix) const 
  {
    assert( numDofs == set.numBaseFunctions() );

    // clear matrix 
    matrix = 0;

    // create quadrature 
    VolumeQuadratureType volQuad(en, volumeQuadOrd_ );

    if( caller.hasMass() )
    {
      // build matix with calling 
      buildMatrixWithMassFactor(caller, en, geo, set, volQuad, numDofs, matrix);
    }
    else 
    {
      buildMatrixNoMassFactor(en, geo, set, volQuad, numDofs, matrix);
    }
  }

  //! build local mass matrix with mass factor 
  template <class Matrix>
  void buildMatrixNoMassFactor(
                   const EntityType& en,
                   const Geometry& geo, 
                   const BaseFunctionSetType& set,
                   const VolumeQuadratureType& volQuad,
                   const int numDofs,
                   Matrix& matrix,
                   const bool applyIntegrationElement = true ) const 
  {
    const int volNop = volQuad.nop();
    for(int qp=0; qp<volNop; ++qp) 
    {
      // calculate integration weight 
      const double intel = ( applyIntegrationElement ) ? 
          ( volQuad.weight(qp) * geo.integrationElement(volQuad.point(qp)) ) : volQuad.weight(qp) ;

      // eval base functions 
      set.evaluateAll(volQuad[qp], phi_);

      for(int m=0; m<numDofs; ++m)
      {
        const RangeType& phi_m = phi_[m];
        const ctype val = intel * (phi_m * phi_m);
        matrix[m][m] += val;
       
        for(int k=m+1; k<numDofs; ++k) 
        {
          const ctype val = intel * (phi_m * phi_[k]);
          matrix[m][k] += val;
          matrix[k][m] += val;
        }
      }
    }
  }

  //! build local mass matrix with mass factor 
  template <class MassCallerType, class Matrix> 
  void buildMatrixWithMassFactor(
                   const MassCallerType& caller,
                   const EntityType& en,
                   const Geometry& geo, 
                   const BaseFunctionSetType& set,
                   const VolumeQuadratureType& volQuad,
                   const int numDofs,
                   Matrix& matrix) const 
  {
    typedef typename MassCallerType :: MassFactorType MassFactorType;
    MassFactorType mass;

    const int volNop = volQuad.nop();
    for(int qp=0; qp<volNop; ++qp) 
    {
      // calculate integration weight 
      const double intel = volQuad.weight(qp)
         * geo.integrationElement(volQuad.point(qp));

      // eval base functions 
      set.evaluateAll( volQuad[qp], phi_);

      // call mass factor 
      caller.mass( en, volQuad, qp, mass);

      // apply mass matrix to all base functions 
      for(int m=0; m<numDofs; ++m)
      {
        mass.mv( phi_[m], phiMass_[m] );
      }

      // add values to matrix 
      for(int m=0; m<numDofs; ++m)
      {
        for(int k=0; k<numDofs; ++k) 
        {
          matrix[m][k] += intel * (phiMass_[m] * phi_[k]);
        }
      }
    }
  }
}; // end class LocalMassMatrixImplementation

///////////////////////////////////////////////////////////////////
//
//  default LocalMassMatrix Implementation 
//
///////////////////////////////////////////////////////////////////

/** \brief Local Mass Matrix for arbitrary spaces */ 
template <class DiscreteFunctionSpaceImp, class VolumeQuadratureImp> 
class LocalMassMatrix 
  : public LocalMassMatrixImplementation< DiscreteFunctionSpaceImp, VolumeQuadratureImp >
{
  typedef LocalMassMatrixImplementation< DiscreteFunctionSpaceImp, VolumeQuadratureImp > BaseType;
public:
  LocalMassMatrix( const DiscreteFunctionSpaceImp& spc, const int volQuadOrd = -1 )
    : BaseType( spc, volQuadOrd )
  {}
};

/** \brief Local Mass Matrix for DG Operators (deprecated, use Fem::LocalMassMatrix)*/ 
template <class DiscreteFunctionSpaceImp, class VolumeQuadratureImp> 
class LocalDGMassMatrix 
  : public LocalMassMatrixImplementation< DiscreteFunctionSpaceImp, VolumeQuadratureImp >
{
  typedef LocalMassMatrixImplementation< DiscreteFunctionSpaceImp, VolumeQuadratureImp > BaseType;
public:
  DUNE_VERSION_DEPRECATED(1,3,remove) 
  LocalDGMassMatrix( const DiscreteFunctionSpaceImp& spc, const int volQuadOrd = -1 )
    : BaseType( spc, volQuadOrd )
  {}
};

///////////////////////////////////////////////////////////////////
//
//  DG LocalMassMatrix Implementation 
//
///////////////////////////////////////////////////////////////////

/** \brief DG Local Mass Matrix for arbitrary spaces */ 
template <class DiscreteFunctionSpaceImp, class VolumeQuadratureImp> 
class LocalMassMatrixImplementationDgOrthoNormal
  : public LocalMassMatrixImplementation< DiscreteFunctionSpaceImp, VolumeQuadratureImp >
{
  typedef LocalMassMatrixImplementation< DiscreteFunctionSpaceImp, VolumeQuadratureImp > BaseType;
public:
  typedef typename BaseType :: EntityType  EntityType;

  LocalMassMatrixImplementationDgOrthoNormal( const DiscreteFunctionSpaceImp& spc, const int volQuadOrd = -1 )
    : BaseType( spc, volQuadOrd )
  {}

  //! apply local dg mass matrix to local function lf
  //! using the massFactor method of the caller 
  template <class MassCallerType, class LocalFunctionType> 
  void applyInverse(const MassCallerType& caller, 
                    const EntityType& en, 
                    LocalFunctionType& lf) const 
  {
    BaseType :: applyInverseDgOrthoNormalBasis( caller, en, lf );
  }

  //! apply local dg mass matrix to local function lf without mass factor 
  template <class LocalFunctionType> 
  void applyInverse(const EntityType& en, 
                    LocalFunctionType& lf) const 
  {
    typename BaseType :: NoMassDummyCaller caller;
    applyInverse(caller, en, lf );
  }

  //! apply local dg mass matrix to local function lf without mass factor 
  template <class LocalFunctionType> 
  void applyInverse(LocalFunctionType& lf) const 
  {
    applyInverse( lf.entity(), lf );
  }

};

//! @}  
} // end namespace Dune 
#endif
