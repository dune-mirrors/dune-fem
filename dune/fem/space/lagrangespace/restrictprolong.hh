#ifndef DUNE_LAGRANGESPACE_RESTRICTPROLONG_HH
#define DUNE_LAGRANGESPACE_RESTRICTPROLONG_HH

#include <map>

#include <dune/common/geometrytype.hh>

#include <dune/grid/common/genericreferenceelements.hh>

#include <dune/fem/gridpart/gridpart.hh>
#include <dune/fem/function/localfunction/localfunction.hh>
#include <dune/fem/space/lagrangespace/lagrangepoints.hh>

namespace Dune
{

  template< class G, int ord >
  struct LagrangeLocalRestrictProlong
  {
    typedef G Grid;

    typedef typename Grid::ctype ctype;
    static const int dimensionworld = Grid::dimensionworld;
    typedef FieldVector< ctype, dimensionworld > DomainVector;

    typedef typename Grid::template Codim< 0 >::Entity Entity;

    typedef Dune::LagrangePointSet< LeafGridPart< Grid >, ord > LagrangePointSet;

  private:
    typedef typename Entity::LocalGeometry LocalGeometry;

    typedef typename LagrangePointSet::template Codim< 0 >::SubEntityIteratorType
      EntityDofIterator;

    typedef std::map< const GeometryType, const LagrangePointSet * > LagrangePointSetMap;

  public:
    ~LagrangeLocalRestrictProlong ()
    {
      typedef typename LagrangePointSetMap::iterator Iterator;
      const Iterator end = lagrangePointSet_.end();
      for( Iterator it = lagrangePointSet_.begin(); it != end; ++it )
        delete it->second;
    }

    template< class FT, class ST >
    void restrictLocal ( LocalFunction< FT > &fatherFunction,
                         const LocalFunction< ST > &sonFunction,
                         const bool initialize ) const
    {
      static const int dimRange = LocalFunction< ST >::dimRange;

      const Entity &father = fatherFunction.entity();
      const Entity &son = sonFunction.entity();

      const GenericReferenceElement< ctype, dimensionworld > &refSon
        = GenericReferenceElements< ctype, dimensionworld >::general( son.type() );

      const LagrangePointSet &pointSet = lagrangePointSet( father );
      const LocalGeometry &geometryInFather = son.geometryInFather();

      const EntityDofIterator send = pointSet.template endSubEntity< 0 >( 0 );
      for( EntityDofIterator sit = pointSet.template beginSubEntity< 0 >( 0 ); sit != send; ++sit )
      {
        const unsigned int dof = *sit;
        const DomainVector &pointInFather = pointSet.point( dof );
        const DomainVector pointInSon = geometryInFather.local( pointInFather );
        if( refSon.checkInside( pointInSon ) )
        {
          typename LocalFunction< ST >::RangeType phi;
          sonFunction.evaluate( pointInSon, phi );
          for( int coordinate = 0; coordinate < dimRange; ++coordinate )
            fatherFunction[ dimRange * dof + coordinate ] = phi[ coordinate ];
        }
      }
    }

    template< class FT, class ST >
    void prolongLocal ( const LocalFunction< FT > &fatherFunction,
                        LocalFunction< ST > &sonFunction ) const
    {
      static const int dimRange = LocalFunction< FT >::dimRange;

      const Entity &son = sonFunction.entity();

      const LagrangePointSet &pointSet = lagrangePointSet( son );
      const LocalGeometry &geometryInFather = son.geometryInFather();

      const EntityDofIterator send = pointSet.template endSubEntity< 0 >( 0 );
      for( EntityDofIterator sit = pointSet.template beginSubEntity< 0 >( 0 ); sit != send; ++sit )
      {
        const unsigned int dof = *sit;
        const DomainVector &pointInSon = pointSet.point( dof );
        const DomainVector pointInFather = geometryInFather.global( pointInSon );
        
        typename LocalFunction< FT >::RangeType phi;
        fatherFunction.evaluate( pointInFather, phi );
        for( int coordinate = 0; coordinate < dimRange; ++coordinate )
          sonFunction[ dimRange * dof + coordinate ] = phi[ coordinate ];
      }
    }

    const LagrangePointSet &lagrangePointSet ( const Entity &entity ) const
    {
      return lagrangePointSet( entity.type() );
    }

    const LagrangePointSet &lagrangePointSet ( const GeometryType &type ) const
    {
      typedef typename LagrangePointSetMap::iterator Iterator;
      Iterator it = lagrangePointSet_.find( type );
      if( it == lagrangePointSet_.end() )
        it = lagrangePointSet_.insert( it, std::make_pair( type, new LagrangePointSet( type ) ) );
      assert( it->second != 0 );
      return *(it->second);
    }

  private:
    mutable LagrangePointSetMap lagrangePointSet_;
  };

}

#endif // #ifndef DUNE_LAGRANGESPACE_RESTRICTPROLONG_HH
