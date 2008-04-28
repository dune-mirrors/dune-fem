#ifndef DUNE_FEM_DIRICHLETBOUNDARYOPERATOR_HH
#define DUNE_FEM_DIRICHLETBOUNDARYOPERATOR_HH

#include <dune/common/exceptions.hh>

#include <dune/fem/storage/array.hh>

#include <dune/fem/operator/common/operator.hh>
#include <dune/fem/operator/model/boundarymodel.hh>

namespace Dune
{

  template< class WrappedOperatorImp, class BoundaryModelImp >
  class DirichletBoundaryOperator;


  
  template< class WrappedOperatorImp, class BoundaryModelImp >
  class DirichletBoundaryRangeProjection
  : public Operator< typename WrappedOperatorImp :: RangeProjectionType :: DomainFieldType,
                     typename WrappedOperatorImp :: RangeProjectionType :: RangeFieldType,
                     typename WrappedOperatorImp :: RangeProjectionType :: DomainType,
                     typename WrappedOperatorImp :: RangeProjectionType :: RangeType >
  {
  public:
    typedef WrappedOperatorImp WrappedOperatorType;
    typedef BoundaryModelImp BoundaryModelType;

    typedef DirichletBoundaryOperator< WrappedOperatorType, BoundaryModelType > OperatorType;

    typedef typename WrappedOperatorType :: RangeProjectionType WrappedProjectionType;
 
    typedef typename WrappedProjectionType :: DomainFunctionType DomainFunctionType;
    typedef typename WrappedProjectionType :: RangeFunctionType RangeFunctionType;

    typedef typename WrappedProjectionType :: DomainFieldType DomainFieldType;
    typedef typename WrappedProjectionType :: RangeFieldType RangeFieldType;
 
  private:
    typedef DirichletBoundaryRangeProjection< WrappedOperatorType, BoundaryModelType >
      ThisType;
    typedef Operator< DomainFieldType, RangeFieldType, DomainFunctionType, RangeFunctionType >
      BaseType;

    friend class DirichletBoundaryOperator< WrappedOperatorType, BoundaryModelType >;

  public:
    typedef typename WrappedProjectionType :: RangeFunctionSpaceType
      RangeFunctionSpaceType;
  
  protected:
    const OperatorType &operator_;
    const WrappedProjectionType wrappedProjection_;

  protected:
    inline explicit DirichletBoundaryRangeProjection ( const OperatorType &op )
    : operator_( op ),
      wrappedProjection_( operator_.wrappedOperator_.rangeProjection() )
    {
    }

  public:
    inline void operator() ( const DomainFunctionType &u, RangeFunctionType &w ) const
    {
      typedef typename RangeFunctionType :: DofIteratorType DofIteratorType;
      
      wrappedProjection_( u, w );

      const DynamicArray< bool > &isDirichletDof = operator_.isDirichletDof_;
      const DynamicArray< RangeFieldType > &dirichletValues = operator_.dirichletValues_;
      
      const DofIteratorType end = w.dend();
      unsigned int globalDofNumber = 0;
      for( DofIteratorType it = w.dbegin(); it != end; ++it )
      {
        if( isDirichletDof[ globalDofNumber ] )
          *it = dirichletValues[ globalDofNumber ];
        ++globalDofNumber;
      }
    }
  };



  /** \class DirichletBoundaryOperator
   *
   *  \note The wrapped operator must be an endomorphism.
   */
  template< class WrappedOperatorImp, class BoundaryModelImp >
  class DirichletBoundaryOperator
  : public Operator< typename WrappedOperatorImp :: DomainFieldType,
                     typename WrappedOperatorImp :: RangeFieldType,
                     typename WrappedOperatorImp :: DomainType,
                     typename WrappedOperatorImp :: RangeType >
  {
  public:
    typedef WrappedOperatorImp WrappedOperatorType;
    typedef BoundaryModelImp BoundaryModelType;
    
    typedef typename WrappedOperatorType :: DomainFieldType DomainFieldType;
    typedef typename WrappedOperatorType :: RangeFieldType RangeFieldType;

    typedef typename WrappedOperatorType :: DomainFunctionType DomainFunctionType;
    typedef typename WrappedOperatorType :: RangeFunctionType RangeFunctionType;

  private:
    typedef DirichletBoundaryOperator< WrappedOperatorType, BoundaryModelType > ThisType;
    typedef Operator< DomainFieldType, RangeFieldType, DomainFunctionType, RangeFunctionType >
      BaseType;

    friend class DirichletBoundaryRangeProjection< WrappedOperatorType, BoundaryModelType >;

  public:
    typedef typename WrappedOperatorType :: DomainFunctionSpaceType
      DomainFunctionSpaceType;
    typedef typename WrappedOperatorType :: RangeFunctionSpaceType
      RangeFunctionSpaceType;

    typedef typename RangeFunctionSpaceType :: GridPartType GridPartType;
    typedef typename RangeFunctionSpaceType :: LagrangePointSetType LagrangePointSetType;

    typedef typename WrappedOperatorType :: DomainProjectionType DomainProjectionType;
    typedef DirichletBoundaryRangeProjection< WrappedOperatorImp, BoundaryModelType >
      RangeProjectionType;

  protected:
    const WrappedOperatorType &wrappedOperator_;
    const BoundaryModelType &boundaryModel_;

    const unsigned int numDofs_;
    DynamicArray< bool > isDirichletDof_;
    DynamicArray< RangeFieldType > dirichletValues_;

  public:
    inline DirichletBoundaryOperator ( const WrappedOperatorType &wrappedOperator,
                                       const BoundaryModelType &boundaryModel )
    : wrappedOperator_( wrappedOperator ),
      boundaryModel_( boundaryModel ),
      numDofs_( wrappedOperator_.rangeFunctionSpace().size() ),
      isDirichletDof_( numDofs_ ),
      dirichletValues_( numDofs_ )
    {
      rebuild();
    }

    inline void operator() ( const DomainFunctionType &u, RangeFunctionType &w ) const
    {
      typedef typename DomainFunctionType :: ConstDofIteratorType DomainDofIteratorType;
      typedef typename RangeFunctionType :: DofIteratorType RangeDofIteratorType;
      
      wrappedOperator_( u, w );

      DomainDofIteratorType uIt = u.dbegin();
      RangeDofIteratorType wIt = w.dbegin();
      const RangeDofIteratorType wEnd = w.dend();
      
      unsigned int globalDofNumber = 0;
      for( ; wIt != wEnd; ++wIt, ++uIt )
      {
        if( isDirichletDof_[ globalDofNumber ] )
          *wIt = *uIt;
        ++globalDofNumber;
      }
    }

    template< class MatrixType >
    inline void assembleMatrix ( MatrixType &matrix ) const
    {
      wrappedOperator_.assembleMatrix( matrix );

      for( unsigned int i = 0; i < numDofs_; ++i )
        if( isDirichletDof_[ i ] )
          matrix.unitRow( i );
    }

    inline const DomainFunctionSpaceType &domainFunctionSpace () const
    {
      return wrappedOperator_.domainFunctionSpace();
    }

    void rebuild ()
    {
      typedef typename RangeFunctionSpaceType :: IteratorType IteratorType;

      isDirichletDof_.assign( false );
      dirichletValues_.assign( 0 );
      
      const RangeFunctionSpaceType &dfSpace = rangeFunctionSpace();
      
      const IteratorType end = dfSpace.end();
      for( IteratorType it = dfSpace.begin(); it != end; ++it )
        markBoundaryDofs( *it, dfSpace );
    }
    
    inline const RangeFunctionSpaceType &rangeFunctionSpace () const
    {
      return wrappedOperator_.rangeFunctionSpace();
    }

    inline const DomainProjectionType domainProjection ()
    {
      return wrappedOperator_.domainProjection();
    }

    inline const RangeProjectionType rangeProjection () const
    {
      return RangeProjectionType( *this );
    }

  protected:
    template< class EntityType >
    inline void markBoundaryDofs ( const EntityType &entity,
                                   const RangeFunctionSpaceType &dfSpace )
    {
      enum { faceCodim = 1 };
      typedef typename LagrangePointSetType :: template Codim< faceCodim >
                                            :: SubEntityIteratorType
        FaceDofIteratorType;

      typedef typename GridPartType :: IntersectionIteratorType IntersectionIteratorType;

      if( !entity.hasBoundaryIntersections() )
        return;
      
      const LagrangePointSetType &lagrangePointSet = dfSpace.lagrangePointSet( entity );
      
      const GridPartType &gridPart = dfSpace.gridPart();
      const IntersectionIteratorType end = gridPart.iend( entity );
      for( IntersectionIteratorType it = gridPart.ibegin( entity ); it != end; ++it )
      {
        if( !it.boundary() )
          continue;
        if( boundaryModel_.boundaryType( it ) != BoundaryModelType :: Dirichlet )
          continue;

        const int faceNumber = it.numberInSelf();
        FaceDofIteratorType faceDofIt
          = lagrangePointSet.template beginSubEntity< faceCodim >( faceNumber );
        const FaceDofIteratorType faceDofEnd
          = lagrangePointSet.template endSubEntity< faceCodim >( faceNumber );
        for( ; faceDofIt != faceDofEnd; ++faceDofIt )
        {
          const int localDofNumber = *faceDofIt;
          const int globalDofNumber = dfSpace.mapToGlobal( entity, localDofNumber );

          bool &isDirichletDof = isDirichletDof_[ globalDofNumber ];
          if( isDirichletDof )
            continue;
          isDirichletDof = true;

          typename BoundaryModelType :: RangeType g;
          boundaryModel_.dirichletValues( it, lagrangePointSet, localDofNumber, g );
          dirichletValues_[ globalDofNumber ] = g[ 0 ];
        }
      }
    }
  };

}

#endif
