#ifndef DUNE_FEM_SPACE_LAGRANGE_RESTRICTPROLONG_HH
#define DUNE_FEM_SPACE_LAGRANGE_RESTRICTPROLONG_HH

// C++ includes
#include <map>

// dune-geometry includes
#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/type.hh>

// dune-fem includes
#include <dune/fem/gridpart/leafgridpart.hh>

// local includes
#include "lagrangepoints.hh"


namespace Dune
{

  namespace Fem
  {

    template< class G, int ord >
    struct LagrangeLocalRestrictProlong
    {
      typedef G GridType;

      typedef typename GridType::ctype ctype;
      static const int dimension = GridType::dimension;

      typedef FieldVector< ctype, dimension > DomainVector;

      typedef LagrangePointSet< LeafGridPart< GridType >, ord > LagrangePointSetType;

    private:
      typedef typename LagrangePointSetType::template Codim< 0 >::SubEntityIteratorType
        EntityDofIterator;

      typedef std::map< const GeometryType, const LagrangePointSetType * > LagrangePointSetMapType;

    public:
      ~LagrangeLocalRestrictProlong ()
      {
        typedef typename LagrangePointSetMapType::iterator Iterator;
        const Iterator end = lagrangePointSet_.end();
        for( Iterator it = lagrangePointSet_.begin(); it != end; ++it )
          delete it->second;
      }

      template< class DomainField >
      void setFatherChildWeight ( const DomainField &weight ) {}

      template< class LFFather, class LFSon, class LocalGeometry >
      void restrictLocal ( LFFather &lfFather,
                           const LFSon &lfSon,
                           const LocalGeometry &geometryInFather,
                           const bool initialize ) const
      {
        static const int dimRange = LFSon::dimRange;

        const auto &refSon = Dune::ReferenceElements< ctype, dimension >::general( lfSon.entity().type() );

        const LagrangePointSetType &pointSet = lagrangePointSet( lfFather.entity(), lfFather.order() );

        const EntityDofIterator send = pointSet.template endSubEntity< 0 >( 0 );
        for( EntityDofIterator sit = pointSet.template beginSubEntity< 0 >( 0 ); sit != send; ++sit )
        {
          const unsigned int dof = *sit;
          const DomainVector &pointInFather = pointSet.point( dof );
          const DomainVector pointInSon = geometryInFather.local( pointInFather );
          if( refSon.checkInside( pointInSon ) )
          {
            typename LFSon::RangeType phi;
            lfSon.evaluate( pointInSon, phi );
            for( int coordinate = 0; coordinate < dimRange; ++coordinate )
              lfFather[ dimRange * dof + coordinate ] = phi[ coordinate ];
          }
        }
      }
      template< class LFFather >
      void restrictFinalize ( LFFather &lfFather ) const
      {}

      template< class LFFather, class LFSon, class LocalGeometry >
      void prolongLocal ( const LFFather &lfFather, LFSon &lfSon,
                          const LocalGeometry &geometryInFather,
                          bool initialize ) const
      {
        static const int dimRange = LFFather::dimRange;

        const LagrangePointSetType &pointSet = lagrangePointSet( lfSon.entity(), lfSon.order() );

        const EntityDofIterator send = pointSet.template endSubEntity< 0 >( 0 );
        for( EntityDofIterator sit = pointSet.template beginSubEntity< 0 >( 0 ); sit != send; ++sit )
        {
          const unsigned int dof = *sit;
          const DomainVector &pointInSon = pointSet.point( dof );
          const DomainVector pointInFather = geometryInFather.global( pointInSon );

          typename LFFather::RangeType phi;
          lfFather.evaluate( pointInFather, phi );
          for( int coordinate = 0; coordinate < dimRange; ++coordinate )
            lfSon[ dimRange * dof + coordinate ] = phi[ coordinate ];
        }
      }

      bool needCommunication () const { return true; }

    protected:
      template< class Entity >
      const LagrangePointSetType &lagrangePointSet ( const Entity &entity, const int order ) const
      {
        return lagrangePointSet( entity.type(), order );
      }

      const LagrangePointSetType &lagrangePointSet ( const GeometryType &type, const int order ) const
      {
        typedef typename LagrangePointSetMapType::iterator Iterator;
        Iterator it = lagrangePointSet_.find( type );
        if( it == lagrangePointSet_.end() )
          it = lagrangePointSet_.insert( it, std::make_pair( type, new LagrangePointSetType( type, order ) ) );
        assert( it->second != 0 );
        return *(it->second);
      }

    private:
      mutable LagrangePointSetMapType lagrangePointSet_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SPACE_LAGRANGE_RESTRICTPROLONG_HH
