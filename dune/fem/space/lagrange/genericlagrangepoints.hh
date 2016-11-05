#ifndef DUNE_FEM_SPACE_LAGRANGE_GENERICLAGRANGEPOINTS_HH
#define DUNE_FEM_SPACE_LAGRANGE_GENERICLAGRANGEPOINTS_HH

// C++ includes
#include <iostream>

// dune-common includes
#include <dune/common/fvector.hh>

// local includes
#include "genericgeometry.hh"


namespace Dune
{

  namespace Fem
  {

    template< class ImplType, unsigned int order, bool bottom = true >
    class GenericLagrangePoint;


    template< unsigned int order, bool bottom >
    class GenericLagrangePoint< PointGeometry, order, bottom >
    {
    public:
      typedef PointGeometry GeometryType;
      static const unsigned int dimension = GeometryType::dimension;

      typedef LocalCoordinate< GeometryType, unsigned int > DofCoordinateType;

      static const unsigned int polynomialOrder = order;

      template< class, unsigned int, bool >
      friend class GenericLagrangePoint;

      template< class, class, unsigned int >
      friend class GenericLagrangeBaseFunction;

    private:
      typedef GenericLagrangePoint< GeometryType, polynomialOrder > ThisType;

    public:
      static const unsigned int numLagrangePoints = 1;

    protected:
      DofCoordinateType dofCoordinate_;

    public:
      template< unsigned int codim >
      struct Codim;

    public:
      GenericLagrangePoint ( unsigned int index )
      {
        dofCoordinate( index, dofCoordinate_ );
      }

      GenericLagrangePoint ( const ThisType &point )
      : dofCoordinate_( point.dofCoordinate_ )
      {}

      template< class LocalCoordinateType >
      static void dofSubEntity ( LocalCoordinateType &coordinate,
                                 unsigned int &codim,
                                 unsigned int &subEntity )
      {
        codim = 0;
        subEntity = 0;
      }

      template< class LocalCoordinateType >
      static void dofSubEntity ( LocalCoordinateType &coordinate,
                                 unsigned int &codim,
                                 unsigned int &subEntity,
                                 unsigned int &dofNumber )
      {
        codim = 0;
        subEntity = 0;
        dofNumber = 0;
      }

      void dofSubEntity ( unsigned int &codim,
                          unsigned int &subEntity )
      {
        dofSubEntity( dofCoordinate_, codim, subEntity );
      }

      void dofSubEntity ( unsigned int &codim,
                          unsigned int &subEntity,
                          unsigned int &dofNumber )
      {
        dofSubEntity( dofCoordinate_, codim, subEntity, dofNumber );
      }

      static unsigned int entityDofNumber ( unsigned int codim,
                                            unsigned int subEntity,
                                            unsigned int dof )
      {
        //assert( (codim == 0) && (subEntity == 0) && (dof == 0) );
        return 0;
      }

      template< class LocalCoordinateType >
      static unsigned int height ( LocalCoordinateType &coordinate )
      {
        return polynomialOrder;
      }

      unsigned int height ()
      {
        return height( dofCoordinate_ );
      }

      template< class FieldType >
      void local ( FieldVector< FieldType, dimension > &coordinate ) const
      {
        const FieldType factor = FieldType( 1 ) / FieldType( polynomialOrder );
        for( unsigned int i = 0; i < dimension; ++i )
          coordinate[ i ] = factor * dofCoordinate_[ i ];
      }

      /** \brief obtain the maximal number of DoFs in one entity of a codimension
       *
       *  \param[in]  codim  codimension, the information is desired for
       *
       *  \returns maximal number of DoFs for one entity in the codimension
       */
      static unsigned int maxDofs ( unsigned int codim )
      {
        return ((codim == 0) ? 1 : 0);
      }

      /** \brief obtain the number of DoFs on one entity
       *
       *  \param[in]  codim      codimension of the entity
       *  \param[in]  subEntity  number of the subentity (of the given codimension)
       *
       *  \returns the number of DoFs associated with the specified entity
       */
      static unsigned int numDofs ( unsigned int codim, unsigned int subEntity )
      {
        return ((codim == 0) ? 1 : 0);
      }

      /** \brief obtain the total number of DoFs in a codimension
       *
       *  \param[in]  codim      codimension the information is desired for
       *
       *  \returns the number of DoFs associated with the codimension
       */
      static unsigned int numDofs ( unsigned int codim )
      {
        return ((codim == 0) ? 1 : 0);
      }

    protected:
      template< class LocalCoordinateType >
      static inline void dofCoordinate ( unsigned int index,
                                         LocalCoordinateType &coordinate )
      {
        assert( index <= numLagrangePoints );
        coordinate = 0;
      }
    };


    /** \cond */
    template< unsigned int order, bool bottom >
    template< unsigned int codim >
    struct GenericLagrangePoint< PointGeometry, order, bottom >::Codim
    {
      static unsigned int maxDofs ()
      {
        return ((codim == 0) ? 1 : 0);
      }
    };
    /** \endcond */


    template< class BaseGeometry, bool bottom >
    class GenericLagrangePoint< PyramidGeometry< BaseGeometry >, 0, bottom >
    {
    public:
      typedef BaseGeometry BaseGeometryType;
      typedef PyramidGeometry< BaseGeometryType > GeometryType;

      static const unsigned int dimension = GeometryType::dimension;

      typedef LocalCoordinate< GeometryType, unsigned int > DofCoordinateType;

      static const unsigned int polynomialOrder = 0;

      template< class, unsigned int, bool >
      friend class GenericLagrangePoint;

      template< class, class, unsigned int >
      friend class GenericLagrangeBaseFunction;

    private:
      typedef GenericLagrangePoint< GeometryType, polynomialOrder > ThisType;

    public:
      static const unsigned int numLagrangePoints = 1;

    protected:
      DofCoordinateType dofCoordinate_;

    public:
      template< unsigned int codim >
      struct Codim;

    public:
      GenericLagrangePoint ( unsigned int index )
      {
        dofCoordinate( index, dofCoordinate_ );
      }

      GenericLagrangePoint ( const ThisType &point )
      : dofCoordinate_( point.dofCoordinate_ )
      {}

      template< class LocalCoordinateType >
      static void dofSubEntity ( LocalCoordinateType &coordinate,
                                 unsigned int &codim,
                                 unsigned int &subEntity )
      {
        codim = (bottom ? 0 : dimension);
        subEntity = 0;
      }

      template< class LocalCoordinateType >
      static void dofSubEntity ( LocalCoordinateType &coordinate,
                                 unsigned int &codim,
                                 unsigned int &subEntity,
                                 unsigned int &dofNumber )
      {
        codim = (bottom ? 0 : dimension);
        subEntity = 0;
        dofNumber = 0;
      }

      void dofSubEntity ( unsigned int &codim,
                          unsigned int &subEntity )
      {
        dofSubEntity( dofCoordinate_, codim, subEntity );
      }

      void dofSubEntity ( unsigned int &codim,
                          unsigned int &subEntity,
                          unsigned int &dofNumber )
      {
        dofSubEntity( dofCoordinate_, codim, subEntity, dofNumber );
      }

      static unsigned int entityDofNumber ( unsigned int codim,
                                            unsigned int subEntity,
                                            unsigned int dof )
      {
        return 0;
      }

      template< class LocalCoordinateType >
      static unsigned int height ( LocalCoordinateType &coordinate )
      {
        return polynomialOrder;
      }

      unsigned int height ()
      {
        return height( dofCoordinate_ );
      }

      template< class FieldType >
      void local ( FieldVector< FieldType, dimension > &coordinate ) const
      {
        const FieldType factor = FieldType( 1 ) / FieldType( polynomialOrder );
        for( unsigned int i = 0; i < dimension; ++i )
          coordinate[ i ] = factor * dofCoordinate_[ i ];
      }

      /** \brief obtain the maximal number of DoFs in one entity of a codimension
       *
       *  \param[in]  codim  codimension, the information is desired for
       *
       *  \returns maximal number of DoFs for one entity in the codimension
       */
      static unsigned int maxDofs ( unsigned int codim )
      {
        if( bottom )
          return ((codim == 0) ? 1 : 0);
        else
          return ((codim == dimension) ? 1 : 0);
      }

      /** \brief obtain the number of DoFs on one entity
       *
       *  \param[in]  codim      codimension of the entity
       *  \param[in]  subEntity  number of the subentity (of the given codimension)
       *
       *  \returns the number of DoFs associated with the specified entity
       */
      static unsigned int numDofs ( unsigned int codim, unsigned int subEntity )
      {
        if( bottom )
          return ((codim == 0) ? 1 : 0);
        else
          return ((codim == dimension) ? 1 : 0);
      }

      /** \brief obtain the total number of DoFs in a codimension
       *
       *  \param[in]  codim      codimension the information is desired for
       *
       *  \returns the number of DoFs associated with the codimension
       */
      static unsigned int numDofs ( unsigned int codim )
      {
        if( bottom )
          return ((codim == 0) ? 1 : 0);
        else
          return ((codim == dimension) ? 1 : 0);
      }

    protected:
      template< class LocalCoordinateType >
      static void dofCoordinate ( unsigned int index,
                                  LocalCoordinateType &coordinate )
      {
        assert( index <= numLagrangePoints );
        coordinate = 0;
      }
    };


    /** \cond */
    template< class BaseGeometry, bool bottom >
    template< unsigned int codim >
    struct GenericLagrangePoint< PyramidGeometry< BaseGeometry >, 0, bottom >::Codim
    {
      static inline unsigned int maxDofs ()
      {
        if( bottom )
          return ((codim == 0) ? 1 : 0);
        else
          return ((codim == dimension) ? 1 : 0);
      }
    };
    /** \endcond */


    template< class BaseGeometry, unsigned int order, bool bottom >
    class GenericLagrangePoint< PyramidGeometry< BaseGeometry >, order, bottom >
    {
    public:
      typedef BaseGeometry BaseGeometryType;
      typedef PyramidGeometry< BaseGeometryType > GeometryType;
      static const unsigned int dimension = GeometryType :: dimension;

      typedef LocalCoordinate< GeometryType, unsigned int > DofCoordinateType;

      static const unsigned int polynomialOrder = order;

      template< class, unsigned int, bool >
      friend class GenericLagrangePoint;

      template< class, class, unsigned int >
      friend class GenericLagrangeBaseFunction;

    private:
      typedef GenericLagrangePoint< GeometryType, polynomialOrder > ThisType;

      typedef GenericLagrangePoint< GeometryType, polynomialOrder - 1, false >
        OrderReductionType;
      typedef GenericLagrangePoint< BaseGeometryType, polynomialOrder >
        DimensionReductionType;

    public:
      static const unsigned int numLagrangePoints
        = DimensionReductionType::numLagrangePoints + OrderReductionType::numLagrangePoints;

    protected:
      DofCoordinateType dofCoordinate_;

    public:
      template< unsigned int codim >
      struct Codim;

    public:
      GenericLagrangePoint ( unsigned int index )
      {
        dofCoordinate( index, dofCoordinate_ );
      }

      GenericLagrangePoint ( const ThisType &point )
      : dofCoordinate_( point.dofCoordinate_ )
      {}

      template< class LocalCoordinateType >
      static void dofSubEntity ( LocalCoordinateType &coordinate,
                                 unsigned int &codim,
                                 unsigned int &subEntity )
      {
        if( !useDimReduction( coordinate ) )
        {
          --(*coordinate);
          OrderReductionType::dofSubEntity( coordinate, codim, subEntity );
          ++(*coordinate);

          if( bottom && (codim > 0) )
            subEntity += BaseGeometryType::numSubEntities( codim - 1 );
        }
        else
        {
          DimensionReductionType::dofSubEntity( coordinate.base(), codim, subEntity );
          if( bottom )
            ++codim;
        }
      }

      template< class LocalCoordinateType >
      static void dofSubEntity ( LocalCoordinateType &coordinate,
                                 unsigned int &codim,
                                 unsigned int &subEntity,
                                 unsigned int &dofNumber )
      {
        if( !useDimReduction( coordinate ) )
        {
          --(*coordinate);
          OrderReductionType::template dofSubEntity( coordinate, codim, subEntity, dofNumber );
          ++(*coordinate);

          if( bottom )
            subEntity += (codim > 0 ? BaseGeometryType :: numSubEntities( codim - 1 ) : 0);
          else
            dofNumber += DimensionReductionType::numDofs( codim, subEntity );
        }
        else
        {
          DimensionReductionType::dofSubEntity( coordinate.base(), codim, subEntity, dofNumber );
          if( bottom )
            ++codim;
        }
      }


      void dofSubEntity ( unsigned int &codim,
                          unsigned int &subEntity )
      {
        dofSubEntity( dofCoordinate_, codim, subEntity );
      }

      void dofSubEntity ( unsigned int &codim,
                          unsigned int &subEntity,
                          unsigned int &dofNumber )
      {
        dofSubEntity( dofCoordinate_, codim, subEntity, dofNumber );
      }

      static unsigned int entityDofNumber ( unsigned int codim,
                                            unsigned int subEntity,
                                            unsigned int dof )
      {
        if( bottom )
        {
          if( codim == 0 )
            return OrderReductionType::entityDofNumber( codim, subEntity, dof )
                   + DimensionReductionType::numLagrangePoints;

          const unsigned int numBaseSubEntities
            = BaseGeometryType::numSubEntities( codim - 1 );
          if( subEntity >= numBaseSubEntities )
            return OrderReductionType::entityDofNumber( codim, subEntity - numBaseSubEntities, dof )
                   + DimensionReductionType::numLagrangePoints;
          else
            return DimensionReductionType::entityDofNumber( codim - 1, subEntity, dof );
        }
        else
        {
          const unsigned int numBaseEntityDofs
            = DimensionReductionType::numDofs( codim, subEntity );
          if( dof >=  numBaseEntityDofs )
            return OrderReductionType::entityDofNumber( codim, subEntity, dof - numBaseEntityDofs )
                   + DimensionReductionType::numLagrangePoints;
          else
            return DimensionReductionType::entityDofNumber( codim, subEntity, dof );
        }
      }

      template< class LocalCoordinateType >
      static unsigned int height ( LocalCoordinateType &coordinate )
      {
        if( !useDimReduction( coordinate ) ) {
          --(*coordinate);
          unsigned int h = OrderReductionType :: height( coordinate );
          ++(*coordinate);
          return h;
        } else
          return DimensionReductionType :: height( coordinate.base() );
      }

      unsigned int height ()
      {
        return height( dofCoordinate_ );
      }

      template< class FieldType >
      void local ( FieldVector< FieldType, dimension > &coordinate ) const
      {
        const FieldType factor = FieldType( 1 ) / FieldType( polynomialOrder );
        for( unsigned int i = 0; i < dimension; ++i )
          coordinate[ i ] = factor * dofCoordinate_[ i ];
      }

      /** \brief obtain the maximal number of DoFs in one entity of a codimension
       *
       *  \param[in]  codim  codimension, the information is desired for
       *
       *  \returns maximal number of DoFs for one entity in the codimension
       */
      static unsigned int maxDofs ( unsigned int codim )
      {
        const unsigned int maxOrderDofs
          = OrderReductionType::maxDofs( codim );

        if( bottom && (codim == 0) )
          return maxOrderDofs;

        const unsigned int maxDimDofs
          = DimensionReductionType::maxDofs( bottom ? codim - 1 : codim );

        return (bottom ? std::max( maxDimDofs, maxOrderDofs )
                       : maxDimDofs + maxOrderDofs);
      }

      /** \brief obtain the number of DoFs on one entity
       *
       *  \param[in]  codim      codimension of the entity
       *  \param[in]  subEntity  number of the subentity (of the given codimension)
       *
       *  \returns the number of DoFs associated with the specified entity
       */
      static unsigned int numDofs ( unsigned int codim, unsigned int subEntity )
      {
        if( bottom )
        {
          if( bottom && (codim == 0) )
            return OrderReductionType::numDofs( codim, subEntity );

          const unsigned int numBaseSubEntities
            = BaseGeometryType::numSubEntities( codim - 1 );
          if( subEntity < numBaseSubEntities )
            return DimensionReductionType::numDofs( codim - 1, subEntity );
          else
            return OrderReductionType::numDofs( codim, subEntity - numBaseSubEntities );
        }
        else
        {
          return DimensionReductionType::numDofs( codim, subEntity )
                 + OrderReductionType::numDofs( codim, subEntity );
        }
      }

      /** \brief obtain the total number of DoFs in a codimension
       *
       *  \param[in]  codim      codimension the information is desired for
       *
       *  \returns the number of DoFs associated with the codimension
       */
      static unsigned int numDofs ( unsigned int codim )
      {
        if( bottom )
        {
          const unsigned int orderDofs = OrderReductionType::numDofs( codim );
          if( codim > 0 )
            return orderDofs + DimensionReductionType::numDofs( codim - 1 );
          else
            return orderDofs;
        }
        else
        {
          return DimensionReductionType::numDofs( codim )
                 + OrderReductionType::numDofs( codim );
        }
      }

      template< class LocalCoordinateType >
      static bool useDimReduction ( const LocalCoordinateType &coordinate )
      {
        return (*coordinate == 0);
      }

    protected:
      template< class LocalCoordinateType >
      static void dofCoordinate ( unsigned int index,
                                  LocalCoordinateType &coordinate )
      {
        assert( index <= numLagrangePoints );

        if( index < DimensionReductionType::numLagrangePoints )
        {
          (*coordinate) = 0;
          DimensionReductionType :: dofCoordinate( index, coordinate.base() );
        }
        else
        {
          const int orderIndex = index - DimensionReductionType::numLagrangePoints;
          OrderReductionType::dofCoordinate( orderIndex, coordinate );
          ++(*coordinate);
        }
      }
    };


    /** \cond */
    template< class BaseGeometry, unsigned int order, bool bottom >
    template< unsigned int codim >
    struct GenericLagrangePoint< PyramidGeometry< BaseGeometry >, order, bottom >::Codim
    {
      static unsigned int maxDofs ()
      {
        const unsigned int maxOrderDofs
          = OrderReductionType::template Codim< codim >::maxDofs();

        const unsigned int maxDimDofs
          = DimensionReductionType::template Codim< (bottom ? codim - 1 : codim) >::maxDofs();

        if( bottom )
          return std::max( maxDimDofs, maxOrderDofs );
        else
          return maxDimDofs + maxOrderDofs;
      }
    };
    /** \endcond */


    template< class FirstGeometry, class SecondGeometry, unsigned int order, bool bottom >
    class GenericLagrangePoint< ProductGeometry< FirstGeometry, SecondGeometry >, order, bottom >
    {
    public:
      typedef ProductGeometry< FirstGeometry, SecondGeometry > GeometryType;
      static const unsigned int dimension = GeometryType::dimension;

      typedef LocalCoordinate< GeometryType, unsigned int > DofCoordinateType;

      static const unsigned int polynomialOrder = order;

      template< class, unsigned int, bool >
      friend class GenericLagrangePoint;

      template< class, class, unsigned int >
      friend class GenericLagrangeBaseFunction;

    private:
      typedef GenericLagrangePoint< GeometryType, polynomialOrder > ThisType;

      typedef GenericLagrangePoint< FirstGeometry, polynomialOrder > FirstReductionType;
      typedef GenericLagrangePoint< SecondGeometry, polynomialOrder > SecondReductionType;

    public:
      static const unsigned int numLagrangePoints
        = FirstReductionType::numLagrangePoints * SecondReductionType::numLagrangePoints;

    protected:
      DofCoordinateType dofCoordinate_;

    private:
      template< unsigned int codim, unsigned int i >
      struct CodimIterator;

    public:
      template< unsigned int codim >
      struct Codim;

    public:
      GenericLagrangePoint ( unsigned int index )
      {
        dofCoordinate( index, dofCoordinate_ );
      }

      GenericLagrangePoint ( const ThisType &point )
      : dofCoordinate_( point.dofCoordinate_ )
      {}

      template< class LocalCoordinateType >
      static void dofSubEntity( LocalCoordinateType &coordinate,
                                unsigned int &codim,
                                unsigned int &subEntity )
      {
        unsigned int firstCodim, secondCodim;
        unsigned int firstSubEntity, secondSubEntity;

        FirstReductionType::dofSubEntity( coordinate.first(), firstCodim, firstSubEntity );
        SecondReductionType::dofSubEntity( coordinate.second(), secondCodim, secondSubEntity );

        codim = firstCodim + secondCodim;

        subEntity = 0;
        for( unsigned int i = 0; i < secondCodim; ++i )
          subEntity += FirstGeometry::numSubEntities( codim-i  ) * SecondGeometry::numSubEntities( i );
        subEntity += firstSubEntity + secondSubEntity * FirstGeometry::numSubEntities( firstCodim );
      }

      template< class LocalCoordinateType >
      static void dofSubEntity( LocalCoordinateType &coordinate,
                                unsigned int &codim,
                                unsigned int &subEntity,
                                unsigned int &dofNumber )
      {
        unsigned int firstCodim, secondCodim;
        unsigned int firstSubEntity, secondSubEntity;
        unsigned int firstDofNumber, secondDofNumber;

        FirstReductionType::dofSubEntity( coordinate.first(), firstCodim, firstSubEntity, firstDofNumber );
        SecondReductionType::dofSubEntity( coordinate.second(), secondCodim, secondSubEntity, secondDofNumber );

        codim = firstCodim + secondCodim;

        subEntity = 0;
        for( unsigned int i = 0; i < secondCodim; ++i )
          subEntity += FirstGeometry::numSubEntities( codim-i  ) * SecondGeometry::numSubEntities( i );
        subEntity += firstSubEntity + secondSubEntity * FirstGeometry::numSubEntities( firstCodim );

        dofNumber = firstDofNumber + secondDofNumber * FirstReductionType::numDofs( firstCodim, firstSubEntity );
      }

      void dofSubEntity ( unsigned int &codim, unsigned int &subEntity )
      {
        dofSubEntity( dofCoordinate_, codim, subEntity );
      }

      void dofSubEntity ( unsigned int &codim,
                          unsigned int &subEntity,
                          unsigned int &dofNumber )
      {
        dofSubEntity( dofCoordinate_, codim, subEntity, dofNumber );
      }

      static unsigned int entityDofNumber ( unsigned int codim,
                                            unsigned int subEntity,
                                            unsigned int dofNumber )
      {
        unsigned int firstCodim = codim;
        unsigned int secondCodim = 0;
        for( ; secondCodim < codim; --firstCodim, ++secondCodim )
        {
          const unsigned int num
            = FirstGeometry::numSubEntities( firstCodim ) * SecondGeometry::numSubEntities( secondCodim );

          if( subEntity < num )
            break;
          subEntity -= num;
        }

        const unsigned int n = FirstGeometry::numSubEntities( firstCodim );
        const unsigned int firstSubEntity = subEntity % n;
        const unsigned int secondSubEntity = subEntity / n;

        const unsigned int m = FirstReductionType::numDofs( firstCodim, firstSubEntity );
        const unsigned int firstDofNumber = dofNumber % m;
        const unsigned int secondDofNumber = dofNumber / m;

        const unsigned int firstEntityDofNumber
          = FirstReductionType::entityDofNumber( firstCodim, firstSubEntity, firstDofNumber );
        const unsigned int secondEntityDofNumber
          = SecondReductionType::entityDofNumber( secondCodim, secondSubEntity, secondDofNumber );

        return firstEntityDofNumber + secondEntityDofNumber * FirstReductionType::numLagrangePoints;
      }

      template< class LocalCoordinateType >
      static unsigned int height ( LocalCoordinateType &coordinate )
      {
        const unsigned int firstHeight = FirstReductionType::height( coordinate.first() );
        const unsigned int secondHeight = SecondReductionType::height( coordinate.second() );

        return ((firstHeight < secondHeight) ? firstHeight : secondHeight);
      }

      unsigned int height ()
      {
        return height( dofCoordinate_ );
      }

      template< class FieldType >
      void local ( FieldVector< FieldType, dimension > &coordinate ) const
      {
        const FieldType factor = FieldType( 1 ) / FieldType( polynomialOrder );
        for( unsigned int i = 0; i < dimension; ++i )
          coordinate[ i ] = factor * dofCoordinate_[ i ];
      }

      /** \brief obtain the maximal number of DoFs in one entity of a codimension
       *
       *  \param[in]  codim  codimension, the information is desired for
       *
       *  \returns maximal number of DoFs for one entity in the codimension
       */
      static unsigned int maxDofs ( unsigned int codim )
      {
        unsigned int max = 0;
        for( unsigned int i = 0; i <= codim; ++i ) {
          const unsigned int n
            = FirstReductionType :: maxDofs( codim - i )
            * SecondReductionType :: maxDofs( i );
          max = (max >= n) ? max : n;
        }
        return max;
      }

      /** \brief obtain the number of DoFs on one entity
       *
       *  \param[in]  codim      codimension of the entity
       *  \param[in]  subEntity  number of the subentity (of the given codimension)
       *
       *  \returns the number of DoFs associated with the specified entity
       */
      static unsigned int numDofs ( unsigned int codim, unsigned int subEntity )
      {
        unsigned int firstCodim = codim;
        unsigned int secondCodim = 0;
        for( ; secondCodim <= codim; --firstCodim, ++secondCodim )
        {
          const unsigned int numSubEntities
            = FirstGeometry::numSubEntities( firstCodim ) * SecondGeometry :: numSubEntities( secondCodim );

          if( subEntity < numSubEntities )
            break;
          subEntity -= numSubEntities;
        }

        const unsigned int n = FirstGeometry::numSubEntities( firstCodim );
        const unsigned int firstSubEntity = subEntity % n;
        const unsigned int secondSubEntity = subEntity / n;

        return FirstReductionType :: numDofs( firstCodim, firstSubEntity )
               * SecondReductionType :: numDofs( secondCodim, secondSubEntity );
      }

      /** \brief obtain the total number of DoFs in a codimension
       *
       *  \param[in]  codim      codimension the information is desired for
       *
       *  \returns the number of DoFs associated with the codimension
       */
      static unsigned int numDofs ( unsigned int codim )
      {
        unsigned int count = 0;

        unsigned int firstCodim = codim;
        unsigned int secondCodim = 0;
        for( ; secondCodim <= codim; --firstCodim, ++secondCodim )
          count += FirstReductionType :: numDofs( firstCodim )
                 * SecondReductionType :: numDofs( secondCodim );

        return count;
      }

    protected:
      template< class LocalCoordinateType >
      static void dofCoordinate ( unsigned int index,
                                  LocalCoordinateType &coordinate )
      {
        assert( index <= numLagrangePoints );

        const unsigned int firstIndex
          = index % FirstReductionType :: numLagrangePoints;
        const unsigned int secondIndex
          = index / FirstReductionType :: numLagrangePoints;

        FirstReductionType :: dofCoordinate( firstIndex, coordinate.first() );
        SecondReductionType :: dofCoordinate( secondIndex, coordinate.second() );
      }
    };


    /** \cond */
    template< class FirstGeometry, class SecondGeometry, unsigned int order, bool bottom >
    template< unsigned int codim >
    struct GenericLagrangePoint< ProductGeometry< FirstGeometry, SecondGeometry >, order, bottom >::Codim
    {
      static unsigned int maxDofs ()
      {
        return CodimIterator< codim, codim >::maxDofs();
      }
    };
    /** \endcond */


    /** \cond */
    template< class FirstGeometry, class SecondGeometry, unsigned int order, bool bottom >
    template< unsigned int codim, unsigned int i >
    struct GenericLagrangePoint< ProductGeometry< FirstGeometry, SecondGeometry >, order, bottom >::CodimIterator
    {
      static unsigned int maxDofs ()
      {
        typedef typename FirstReductionType::template Codim< codim-i > FirstCodim;
        typedef typename SecondReductionType::template Codim< i > SecondCodim;
        const unsigned int n = FirstCodim::maxDofs() * SecondCodim::maxDofs();

        const unsigned int m = CodimIterator< codim, i-1 >::maxDofs();
        return std::max( m, n );
      }
    };

    template< class FirstGeometry, class SecondGeometry, unsigned int order, bool bottom >
    template< unsigned int codim>
    struct GenericLagrangePoint< ProductGeometry< FirstGeometry, SecondGeometry >, order, bottom >::CodimIterator<codim, 0>
    {
      static unsigned int maxDofs ()
      {
        typedef typename FirstReductionType::template Codim< codim > FirstCodim;
        typedef typename SecondReductionType::template Codim< 0 > SecondCodim;
        return FirstCodim::maxDofs() * SecondCodim::maxDofs();
      }
    };
    /** \endcond */

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SPACE_LAGRANGE_GENERICLAGRANGEPOINTS_HH
