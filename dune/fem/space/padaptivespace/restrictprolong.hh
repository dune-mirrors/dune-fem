#ifndef DUNE_FEM_PLAGRANGESPACE_RESTRICTPROLONG_HH
#define DUNE_EFM_PLAGRANGESPACE_RESTRICTPROLONG_HH

#include <map>

#include <dune/geometry/type.hh>
#include <dune/geometry/referenceelements.hh>

#include <dune/fem/function/localfunction/localfunction.hh>
#include <dune/fem/space/lagrangespace/lagrangepoints.hh>

namespace Dune
{

  namespace Fem
  {

    template< class G, class LagrangePointSetProvider >
    struct PLagrangeLocalRestrictProlong
    {
      typedef G Grid;

      typedef typename Grid::ctype ctype;
      static const int dimension = Grid::dimension;
      typedef FieldVector< ctype, dimension > DomainVector;

      typedef typename Grid::template Codim< 0 >::Entity Entity;

      typedef typename LagrangePointSetProvider :: LagrangePointSetType  LagrangePointSet;

    private:
      typedef typename Entity::LocalGeometry LocalGeometry;

      typedef typename LagrangePointSet::template Codim< 0 >::SubEntityIteratorType
        EntityDofIterator;

    public:
      PLagrangeLocalRestrictProlong ( const LagrangePointSetProvider &lpsProvider )
      : lpsProvider_( lpsProvider )
      {}

      template< class DomainField >
      void setFatherChildWeight ( const DomainField &weight ) {}

      template< class FT, class ST, class LocalGeometry >
      void restrictLocal ( LocalFunction< FT > &fatherFunction,
                           const LocalFunction< ST > &sonFunction,
                          const LocalGeometry &geometryInFather,
                           const bool initialize ) const
      {
        static const int dimRange = LocalFunction< ST >::dimRange;

        const Entity &father = fatherFunction.entity();
        const Entity &son = sonFunction.entity();

        const GenericReferenceElement< ctype, dimension > &refSon
          = GenericReferenceElements< ctype, dimension >::general( son.type() );

        const LagrangePointSet &pointSet = lagrangePointSet( father );

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

      template< class FT, class ST, class LocalGeometry >
      void prolongLocal ( const LocalFunction< FT > &fatherFunction,
                          LocalFunction< ST > &sonFunction,
                          const LocalGeometry &geometryInFather,
                          bool initialize ) const
      {
        static const int dimRange = LocalFunction< FT >::dimRange;

        const Entity &son = sonFunction.entity();

        const LagrangePointSet &pointSet = lagrangePointSet( son );

        const EntityDofIterator send = pointSet.template endSubEntity< 0 >( 0 );
        for( EntityDofIterator sit = pointSet.template beginSubEntity< 0 >( 0 ); sit != send; ++sit )
        {
          const unsigned int dof = *sit;
          const DomainVector &pointInSon = pointSet.point( dof );
          const DomainVector pointInFather = geometryInFather.global( pointInSon );
          
          typename LocalFunction< FT >::RangeType phi;
          fatherFunction.evaluate( pointInFather, phi );
          for( int coordinate = 0; coordinate < dimRange; ++coordinate )
          {
            const int idx = dimRange * dof + coordinate  ;
            sonFunction[ idx ] = phi[ coordinate ];
          }
        }
      }

      template< class FT, class ST >
      void localInterpolation ( const LocalFunction< FT > &oldFunction,
                                LocalFunction< ST > &newFunction ) const
      {
        static const int dimRange = LocalFunction< ST >::dimRange;

        const Entity &entity = newFunction.entity();

        const LagrangePointSet &pointSet = lagrangePointSet( entity );

        const EntityDofIterator send = pointSet.template endSubEntity< 0 >( 0 );
        for( EntityDofIterator sit = pointSet.template beginSubEntity< 0 >( 0 ); sit != send; ++sit )
        {
          const unsigned int dof = *sit;
          const DomainVector &localPoint = pointSet.point( dof );
          
          typename LocalFunction< FT >::RangeType phi;
          oldFunction.evaluate( localPoint, phi );
          for( int coordinate = 0; coordinate < dimRange; ++coordinate )
          {
            const int idx = dimRange * dof + coordinate  ;
            newFunction[ idx ] = phi[ coordinate ];
          }
        }
      }

      bool needCommunication () const { return false; }

      const LagrangePointSet &lagrangePointSet ( const Entity &entity ) const
      {
        return lpsProvider_.lagrangePointSet( entity );
      }

    protected:
      const LagrangePointSetProvider& lpsProvider_;
    };

  } // namespace Fem

} // namespace Dune 

#endif // #ifndef DUNE_FEM_LAGRANGESPACE_RESTRICTPROLONG_HH
