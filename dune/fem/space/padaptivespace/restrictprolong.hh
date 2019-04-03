#ifndef DUNE_FEM_SPACE_PADAPTIVE_RESTRICTPROLONG_HH
#define DUNE_FEM_SPACE_PADAPTIVE_RESTRICTPROLONG_HH

#include <dune/geometry/referenceelements.hh>

#include <dune/fem/function/localfunction/localfunction.hh>
#include <dune/fem/space/lagrange/lagrangepoints.hh>


namespace Dune
{

  namespace Fem
  {

    // PLagrangeLocalRestrictProlong
    // -----------------------------

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

      template< class LFFather, class LFSon, class LocalGeometry >
      void restrictLocal ( LFFather &lfFather, const LFSon &lfSon,
                           const LocalGeometry &geometryInFather, bool initialize ) const
      {
        static const int dimRange = LFSon::dimRange;

        const Entity &father = lfFather.entity();
        const Entity &son = lfSon.entity();

        auto refSon = referenceElement< ctype, dimension >( son.type() );

        const LagrangePointSet &pointSet = lagrangePointSet( father );

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
                          const LocalGeometry &geometryInFather, bool initialize ) const
      {
        static const int dimRange = LFFather::dimRange;

        const Entity &son = lfSon.entity();

        const LagrangePointSet &pointSet = lagrangePointSet( son );

        const EntityDofIterator send = pointSet.template endSubEntity< 0 >( 0 );
        for( EntityDofIterator sit = pointSet.template beginSubEntity< 0 >( 0 ); sit != send; ++sit )
        {
          const unsigned int dof = *sit;
          const DomainVector &pointInSon = pointSet.point( dof );
          const DomainVector pointInFather = geometryInFather.global( pointInSon );

          typename LFFather::RangeType phi;
          lfFather.evaluate( pointInFather, phi );
          for( int coordinate = 0; coordinate < dimRange; ++coordinate )
          {
            const int idx = dimRange * dof + coordinate  ;
            lfSon[ idx ] = phi[ coordinate ];
          }
        }
      }

      template< class ArgLocal, class DestLocal >
      void localInterpolation ( const ArgLocal &argLocal,
                                DestLocal &destLocal ) const
      {
        static const int dimRange = DestLocal::dimRange;

        const Entity &entity = destLocal.entity();

        const LagrangePointSet &pointSet = lagrangePointSet( entity );

        const EntityDofIterator send = pointSet.template endSubEntity< 0 >( 0 );
        for( EntityDofIterator sit = pointSet.template beginSubEntity< 0 >( 0 ); sit != send; ++sit )
        {
          const unsigned int dof = *sit;
          const DomainVector &localPoint = pointSet.point( dof );

          typename ArgLocal::RangeType phi;
          argLocal.evaluate( localPoint, phi );
          for( int coordinate = 0; coordinate < dimRange; ++coordinate )
          {
            const int idx = dimRange * dof + coordinate  ;
            destLocal[ idx ] = phi[ coordinate ];
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

#endif // #ifndef DUNE_FEM_SPACE_PADAPTIVE_RESTRICTPROLONG_HH
