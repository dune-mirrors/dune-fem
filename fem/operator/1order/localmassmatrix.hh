#ifndef DUNE_MASSMATRIXINVERTER_HH
#define DUNE_MASSMATRIXINVERTER_HH

//- Dune includes 
#include <dune/common/fmatrix.hh>
#include <dune/fem/quadrature/cachequad.hh>
#include <dune/fem/space/common/allgeomtypes.hh>
#include <dune/fem/misc/checkgeomaffinity.hh>

namespace Dune {

/*! @addtogroup PassHyp
 *  ** @{
*/

/** \brief Local Mass Matrix for DG Operators */ 
template <class DiscreteFunctionSpaceImp, class VolumeQuadratureImp> 
class LocalDGMassMatrix 
{
public:  
  typedef DiscreteFunctionSpaceImp DiscreteFunctionSpaceType;
  enum { numDofs_ = DiscreteFunctionSpaceType :: localBlockSize };
  typedef typename DiscreteFunctionSpaceType :: RangeFieldType ctype;
  typedef typename DiscreteFunctionSpaceType :: RangeType RangeType;
  typedef FieldMatrix<ctype, numDofs_, numDofs_ > MatrixType;
  typedef FieldVector<ctype, numDofs_ > VectorType;

  typedef typename DiscreteFunctionSpaceType :: GridPartType GridPartType;
  typedef typename DiscreteFunctionSpaceType :: BaseFunctionSetType BaseFunctionSetType; 
  typedef typename GridPartType :: GridType GridType;
  typedef typename GridType :: template Codim<0> :: Entity EntityType;
  typedef typename GridType :: template Codim<0> :: Geometry Geometry;

  typedef VolumeQuadratureImp VolumeQuadratureType;

  //! is true if grid is structured grid 
  enum { StructuredGrid = ! Capabilities::IsUnstructured<GridType>::v };

  typedef AllGeomTypes< typename GridPartType :: IndexSetType,GridType> GeometryInformationType;
  typedef typename GeometryInformationType :: DomainType DomainType;
protected:  
  const DiscreteFunctionSpaceType& spc_;
  const int volumeQuadOrd_;

  GeometryInformationType geoInfo_;

  mutable MatrixType matrix_;
  mutable VectorType x_, rhs_;

  mutable RangeType phi_[numDofs_];
  mutable RangeType phiMass_[numDofs_];

  const bool affine_;

  //! dummy caller 
  struct NoMassDummyCaller
  {
    enum { dimRange = DiscreteFunctionSpaceType::DimRange };
    typedef FieldMatrix<ctype, dimRange, dimRange> MassFactorType;
    // return false since we don;t have a mass term
    bool hasMass() const { return false; }
    void mass(const EntityType&,
              const VolumeQuadratureType&,
              const int, 
              const MassFactorType&) const
    {
    }
  };

public:
  //! constructor taking space and volume quadrature order 
  LocalDGMassMatrix(const DiscreteFunctionSpaceType& spc, const int volQuadOrd ) 
    : spc_(spc) 
    , volumeQuadOrd_ ( volQuadOrd )
    , geoInfo_( spc.indexSet() ) 
    , affine_ ( setup() )
  {
    // only for DG spaces at the moment 
    assert( spc_.continuous() == false );
  }

  //! copy constructor 
  LocalDGMassMatrix(const LocalDGMassMatrix& org) 
    : spc_(org.spc_),
      volumeQuadOrd_( org.volumeQuadOrd_ ),
      geoInfo_( org.geoInfo_),
      affine_( org.affine_ )
  {
  }

public:  
  //! apply local dg mass matrix to local function lf
  //! using the massFactor method of the caller 
  template <class MassCallerType, class LocalFunctionType> 
  void applyInverse(MassCallerType& caller, 
                    const EntityType& en, 
                    LocalFunctionType& lf) const 
  {
    // get geometry 
    const Geometry& geo = en.geometry();

    assert( numDofs_ == lf.numDofs() );
    // in case of affine mappings we only have to multiply with a factor 
    if( affine_ && ! caller.hasMass() )
    {
      const double massVolInv = geoInfo_.referenceVolume( geo.type() ) / geo.volume(); 

      // apply inverse mass matrix  
      for(int l=0; l<numDofs_; ++l) 
      {
        lf[l] *= massVolInv;
      }

      return; 
    }
    else 
    {
      // copy local function to right hand side 
      for(int l=0; l<numDofs_; ++l) 
      {
        // copy 
        rhs_[l] = lf[l];
        x_[l] = 0;
      }

      // setup local mass matrix 
      buildMatrix( caller, en, geo, lf.baseFunctionSet(), matrix_ );

      // invert matrix 
      matrix_.invert();
      // apply matrix 
      matrix_.umv( rhs_ , x_ );

      // copy back to local function 
      for(int l=0; l<numDofs_; ++l) 
      {
        lf[l] = x_[l];
      }

      return; 
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

public:  
  //! returns true if geometry mapping is affine 
  bool affine () const { return affine_; }

protected:
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
  template <class MassCallerType> 
  void buildMatrix(MassCallerType& caller,
                   const EntityType& en,
                   const Geometry& geo, 
                   const BaseFunctionSetType& set,
                   MatrixType& matrix) const 
  {
    assert( numDofs_ == set.numBaseFunctions() );

    matrix = 0;
    VolumeQuadratureType volQuad(en, volumeQuadOrd_ );

    if( caller.hasMass() )
    {
      // build matix with calling 
      buildMatrixWithMassFactor(caller, en, geo, set, volQuad, matrix);
    }
    else 
    {
      buildMatrixNoMassFactor(en, geo, set, volQuad, matrix);
    }
  }

  //! build local mass matrix with mass factor 
  void buildMatrixNoMassFactor(
                   const EntityType& en,
                   const Geometry& geo, 
                   const BaseFunctionSetType& set,
                   const VolumeQuadratureType& volQuad,
                   MatrixType& matrix) const 
  {
    const int volNop = volQuad.nop();
    for(int qp=0; qp<volNop; ++qp) 
    {
      // calculate integration weight 
      const double intel = volQuad.weight(qp)
         * geo.integrationElement(volQuad.point(qp));

      for(int m=0; m<numDofs_; ++m)
      {  
        // eval base functions 
        set.evaluate(m, volQuad[qp], phi_[m] );
      }

      for(int m=0; m<numDofs_; ++m)
      {
        {
          const ctype val = intel * (phi_[m] * phi_[m]);
          matrix[m][m] += val;
        }
        
        for(int k=m+1; k<numDofs_; ++k) 
        {
          const ctype val = intel * (phi_[m] * phi_[k]);
          matrix[m][k] += val;
          matrix[k][m] += val;
        }
      }
    }
  }

  //! build local mass matrix with mass factor 
  template <class MassCallerType> 
  void buildMatrixWithMassFactor(
                   MassCallerType& caller,
                   const EntityType& en,
                   const Geometry& geo, 
                   const BaseFunctionSetType& set,
                   const VolumeQuadratureType& volQuad,
                   MatrixType& matrix) const 
  {
    typedef typename MassCallerType :: MassFactorType MassFactorType;
    MassFactorType mass;

    const int volNop = volQuad.nop();
    for(int qp=0; qp<volNop; ++qp) 
    {
      // calculate integration weight 
      const double intel = volQuad.weight(qp)
         * geo.integrationElement(volQuad.point(qp));

      for(int m=0; m<numDofs_; ++m)
      {  
        // eval base functions 
        set.evaluate(m, volQuad[qp], phi_[m] );
      }

      // call mass factor 
      caller.mass( en, volQuad, qp, mass);

      // apply mass matrix to all base functions 
      for(int m=0; m<numDofs_; ++m)
      {
        mass.mv( phi_[m], phiMass_[m] );
      }

      // add values to matrix 
      for(int m=0; m<numDofs_; ++m)
      {
        for(int k=0; k<numDofs_; ++k) 
        {
          matrix[m][k] += intel * (phiMass_[m] * phi_[k]);
        }
      }
    }
  }
}; // end class LocalDGMassMatrix

//! @}  
} // end namespace Dune 
#endif
