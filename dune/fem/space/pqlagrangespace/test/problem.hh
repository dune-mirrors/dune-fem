#ifndef DUNE_FEM_PQLAGRANGESPACE_PROBLEM_HH
#define DUNE_FEM_PQLAGRANGESPACE_PROBLEM_HH

namespace Dune
{

  template< class EntityType >
  struct BoundaryCheck
  {
    typedef typename EntityType::Geometry::LocalCoordinate DomainType;
    typedef FieldVector< typename EntityType::Geometry::ctype, 1 > RangeType;

    BoundaryCheck( const EntityType &entity )
    : entity_( entity )
    {}

    //! check wether x is on boundary element
    void evaluate ( const DomainType &x, RangeType &out ) const
    {
      DomainType xGlobal = entity_.geometry().global( x );

      const double eps = 1e-8;

      if( (xGlobal[ 0 ] < eps) || (xGlobal[ 1 ] < eps) || (xGlobal[ 0 ] > 1-eps) || (xGlobal[ 1 ] > 1-eps) )
        out[ 0 ] = 1;
      else
        out[ 0 ] = 0;
    }

  private:
    const EntityType &entity_;
  };



  // right hand side of governing problem 
  template< class FunctionSpaceImp >
  class RHSFunction
  : public Fem::Function< FunctionSpaceImp, RHSFunction< FunctionSpaceImp > >
  {
    typedef RHSFunction< FunctionSpaceImp > ThisType;
    typedef Fem::Function< FunctionSpaceImp, ThisType > BaseType;

  public:
    typedef FunctionSpaceImp                                         FunctionSpaceType;

    typedef typename FunctionSpaceType :: DomainType                 DomainType;
    typedef typename FunctionSpaceType :: RangeType                  RangeType;

    typedef typename FunctionSpaceType :: DomainFieldType            DomainFieldType;
    typedef typename FunctionSpaceType :: RangeFieldType             RangeFieldType;

    // f( x, y, z ) = 2 \sum_{i=1}^{dimworld} \prod_{j \neq i} (x_j-x_j^2)
    void evaluate( const DomainType &x , RangeType &phi ) const
    {
      enum { dimension = DomainType :: dimension };
      
      phi = 0;
      for( int i = 0; i < dimension; ++i ) { 
        RangeType tmp = 2;
        for( int j = 0; j < dimension; ++j ) {
          if( i == j )
            continue;
          const DomainFieldType &x_j = x[ j ];
          tmp *= x_j - SQR( x_j );
        }
        phi += tmp;
      }
    }
  }; // end class RHSFunction



  //! the exact solution to the problem for EOC calculation 
  template< class FunctionSpaceImp >
  class ExactSolution
  : public Fem::Function < FunctionSpaceImp, ExactSolution< FunctionSpaceImp > >
  {
    typedef ExactSolution< FunctionSpaceImp > ThisType;
    typedef Fem::Function< FunctionSpaceImp, ThisType > BaseType;

  public:
    typedef FunctionSpaceImp                                         FunctionSpaceType;

    typedef typename FunctionSpaceType :: DomainType                 DomainType;
    typedef typename FunctionSpaceType :: RangeType                  RangeType;
    typedef typename FunctionSpaceType :: JacobianRangeType          JacobianRangeType;

    enum { dimDomain = DomainType :: dimension };

    typedef typename FunctionSpaceType :: DomainFieldType            DomainFieldType;
    typedef typename FunctionSpaceType :: RangeFieldType             RangeFieldType;

    // u( x, y, z ) = \prod_{i=1}^{dimworld} (x_i - x_i^2).
    void evaluate ( const DomainType &x , RangeType &phi ) const
    {
      phi = 1;
      for( int i = 0; i < dimDomain; ++i )
      {
        const DomainFieldType &x_i = x[ i ];
        phi *= x_i - SQR( x_i );
      }
    }

    void jacobian ( const DomainType &x, JacobianRangeType &ret ) const
    {
      for( unsigned int i = 0; i < dimDomain; ++i )
      {
        const DomainFieldType &x_i = x[ i ];
        ret[ 0 ][ i ] = 1.0 - 2.0 * x_i;

        for( unsigned int j = 0; j < dimDomain; ++j )
        {
          const DomainFieldType &x_j = x[ j ];
          if( i != j )
            ret[ 0 ][ i ] *= x_j - SQR( x_j );
        }
      }
    }
   
    void evaluate ( const DomainType &x , RangeFieldType time, RangeType &phi ) const
    {
      evaluate( x, phi );
    }
  }; // end class ExactSolution



  // diffusion coefficient for this problem
  template< class FunctionSpaceImp >
  class Tensor
  : public Fem::Function< FunctionSpaceImp, Tensor< FunctionSpaceImp > >
  {
    typedef Tensor< FunctionSpaceImp > ThisType;
    typedef Fem::Function< FunctionSpaceImp, ThisType > BaseType;

  public:
    typedef FunctionSpaceImp                                         FunctionSpaceType;

    typedef typename FunctionSpaceType :: DomainType                 DomainType;
    typedef typename FunctionSpaceType :: RangeType                  RangeType;

    typedef typename FunctionSpaceType :: DomainFieldType            DomainFieldType;
    typedef typename FunctionSpaceType :: RangeFieldType             RangeFieldType;

    void evaluate ( int i, int j, const DomainType &x, RangeType &phi ) const
    {
      evaluate( x, phi );
    }
    
    void evaluate ( const DomainType &x, RangeType &phi ) const
    {
      phi = 1;
    }

    void evaluate ( const DomainType &x, RangeFieldType time, RangeType &phi ) const
    {
      evaluate( x, phi );
    }
  }; //end class Tensor

}

#endif // #ifndef DUNE_FEM_PQLAGRANGESPACE_PROBLEM_HH
