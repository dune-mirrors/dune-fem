#ifndef DUNE_FEM_CHECKGEOMETRYAFFINITY_HH
#define DUNE_FEM_CHECKGEOMETRYAFFINITY_HH

#include <dune/common/fvector.hh>
#include <dune/grid/common/gridenums.hh>

#include <dune/fem/space/common/allgeomtypes.hh>

namespace Dune
{

  namespace Fem
  {

    /*! @addtogroup HelperClasses
     ** @{
     */

    //! Helper class to check affinity of the grid's geometries
    template <class QuadratureType>
    struct GeometryAffinityCheck
    {
      // return true if geometry is affine
      template< class EntityType, class ElementGeometryType >
      static bool checkGeometry( const EntityType& entity,
                                 const ElementGeometryType& geo,
                                 const int quadOrd )
      {
        // if method tells that geometry is not affine
        // then check it carefully
        if( ! geo.affine() )
        {
          // get quadrature of desired order
          QuadratureType volQuad( entity, quadOrd );
          const int nop = volQuad.nop();

          // check all integration elements against the first
          const double oldIntel = geo.integrationElement( volQuad.point(0) );
          for(int l=1; l<nop; ++l)
          {
            const double intel = geo.integrationElement( volQuad.point(l) );
            if( std::abs( oldIntel - intel ) > 1e-12 )
              return false;
          }
        }
        return true ;
      }

      // return true if geometry is affine
      template< class EntityType >
      static bool checkGeometry( const EntityType& entity, const int quadOrd )
      {
        return checkGeometry( entity, entity.geometry(), quadOrd );
      }

      //! check whether all geometry mappings are affine
      template <class IteratorType>
      static inline bool checkAffinity(const IteratorType& begin,
                                       const IteratorType& endit,
                                       const int quadOrd)
      {
        for(IteratorType it = begin; it != endit; ++it)
        {
          if( ! checkGeometry( *it, quadOrd ) ) return false ;
        }
        return true;
      }

      //! check whether all geometry mappings are affine
      template <class GridPartType, class Vector >
      static inline void checkElementAffinity(const GridPartType& gridPart,
                                              const int quadOrd,
                                              Vector& affineGeomtryVec )
      {
        typedef typename GridPartType :: template Codim< 0 > :: IteratorType  IteratorType;
        typedef typename GridPartType :: template Codim< 0 > :: EntityType    EntityType;
        const IteratorType endit = gridPart.template end<0> ();
        affineGeomtryVec.resize( gridPart.indexSet().size( 0 ) );
        for(IteratorType it = gridPart.template begin<0>(); it != endit; ++it)
        {
          const EntityType& entity = *it ;
          const int index = gridPart.indexSet().index( entity );
          affineGeomtryVec[ index ] = checkGeometry( entity, quadOrd );
        }

        //for( size_t i=0; i<affineGeomtryVec.size(); ++ i)
        //  std::cout << "geo is " << affineGeomtryVec[ i ] << std::endl;
      }
    };

    //! Helper class to check whether grid is structured or not
    template <class GridPartType>
    struct CheckCartesian
    {
      typedef typename GridPartType :: GridType GridType ;
    protected:
      //! check whether this is a cartesian grid or not
      template <class IndexSetType>
      static inline bool doCheck(const GridType& grid, const IndexSetType& indexSet)
      {
        typedef typename GridType :: template Codim<0> :: LevelIterator MacroIteratorType;
        typedef typename GridType :: template Codim<0> :: Entity  EntityType;
        typedef typename GridType :: template Codim<0> :: Geometry Geometry;

        typedef typename GridType :: LevelGridView  MacroGridView ;

        // get macro grid view
        MacroGridView macroView = grid.levelGridView( 0 );

        const MacroIteratorType endit = macroView.template end<0> ();
        MacroIteratorType it = macroView.template begin<0> ();

        // empty grids are considerd as cartesian
        if( it == endit ) return true;

        typedef AllGeomTypes< IndexSetType, GridType> GeometryInformationType;
        GeometryInformationType geoInfo( indexSet );

        // if we have more than one geometry Type return false
        if( geoInfo.geomTypes(0).size() > 1 ) return false;

        // if this type is not cube return false
        if( ! geoInfo.geomTypes(0)[0].isCube() ) return false;

        typedef typename GridType :: ctype ctype;
        enum { dimension = GridType :: dimension };
        enum { dimworld  = GridType :: dimensionworld };

        // grid width
        FieldVector<ctype, dimension> h(0);
        FieldVector<ctype, dimension-1> mid(0.5);

        const int map[3] = {1, 2, 4};
        {
          const Geometry geo = it->geometry();
          if ( ! geo.type().isCube() ) return false;

          // calculate grid with
          for(int i=0; i<dimension; ++i)
          {
            h[i] = (geo.corner( 0 ) - geo.corner( map[i] )).two_norm();
          }
        }

        // loop over all macro elements
        for(MacroIteratorType it = macroView.template begin<0> ();
            it != endit; ++it)
        {
          const EntityType& en = *it;
          const Geometry geo = en.geometry();

          const FieldVector<ctype, dimworld> enBary =
            geo.global( geoInfo.localCenter( geo.type() ));

          typedef typename MacroGridView :: IntersectionIterator IntersectionIteratorType;

          // if geometry is not cube, it's not a cartesian grid
          if ( ! geo.type().isCube() ) return false;

          for(int i=0; i<dimension; ++i)
          {
            const ctype w = (geo.corner(0) - geo.corner( map[i] )).two_norm();
            if( std::abs( h[i] - w ) > 1e-15 ) return false;
          }

          IntersectionIteratorType endnit = macroView.iend( en );
          for(IntersectionIteratorType nit = macroView.ibegin( en );
              nit != endnit; ++nit)
          {
            const typename IntersectionIteratorType::Intersection& inter=*nit;
            if( ! checkIntersection(inter) ) return false;

            if( inter.neighbor() )
            {
              EntityType nb = inter.outside();
              Geometry nbGeo = nb.geometry();

              FieldVector<ctype, dimworld> diff = nbGeo.global( geoInfo.localCenter( nbGeo.type() ));
              diff -= enBary;
              assert( diff.two_norm() > 1e-15 );
              diff /= diff.two_norm();

              // scalar product should be 1 or -1
              if( std::abs(std::abs((diff * inter.unitOuterNormal(mid))) - 1.0) > 1e-12 ) return false;
            }
          }
        }
        return true;
      }

      template <class ctype, int dim>
      class ReferenceNormals
      {
        const FieldVector<ctype,dim-1> mid_;
        enum { numberOfNormals = 2 * dim };
        FieldVector<ctype,dim> refNormal_[numberOfNormals];
      public:
        ReferenceNormals () : mid_(0.5)
        {
          for(int i=0; i<numberOfNormals; ++i)
          {
            // get i-th reference normal
            FieldVector<ctype,dim>& refNormal = refNormal_[i];
            // reset normal
            refNormal = 0;
            // set one component
            int comp = ((int) i/2);
            refNormal[comp] = ((i % 2) == 0) ? -1 : 1;
          }
        }

        const FieldVector<ctype,dim>& referenceNormal(const int i) const
        {
          assert( i >= 0 && i< numberOfNormals );
          return refNormal_[i];
        }

        const FieldVector<ctype,dim-1>& faceMidPoint() const
        {
          return mid_;
        }
      };

    public:
      // check that element is provided following the DUNE reference cube
      template <class IntersectionType>
      static bool checkIntersection(const IntersectionType& nit)
      {
        if ( ! nit.type().isCube() ) return false;

        typedef typename IntersectionType :: Entity EntityType;
        typedef typename EntityType :: Geometry :: ctype ctype;
        enum { dimworld = EntityType :: Geometry :: coorddimension };

        // get reference normals
        static const ReferenceNormals<ctype,dimworld> normals;

        // get current normal
        FieldVector<ctype,dimworld> unitNormal = nit.unitOuterNormal(normals.faceMidPoint());
        unitNormal -= normals.referenceNormal( nit.indexInInside() );

        // if normals are not equal grid is not cartesian
        if( unitNormal.infinity_norm() > 1e-10 ) return false;

        return true;
      }

      //! check whether all the is grid is a cartesian grid
      static inline bool check(const GridPartType& gridPart)
      {
        bool cartesian = doCheck( gridPart.grid() , gridPart.indexSet() );
        int val = (cartesian) ? 1 : 0;
        // take global minimum
        val = gridPart.comm().min( val );
        return (val == 1) ? true : false;
      }
    };

    //! @}

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_CHECKGEOMETRYAFFINITY_HH
