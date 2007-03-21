namespace Dune
{
    
  template< class ct >
  LagrangeLineQuadrature< ct >
    :: LagrangeLineQuadrature( int order, size_t id )
    : QuadratureImp< ct, 1 >( id ),
      order_( (order < 0) ? 0 : order )
  {
    const ctype volume = 1.0;
      
    CoordinateType refpoints[ 3 ];
    
    switch( order_ ) {
    case 0:
      refpoints[ 0 ][ 0 ] = 0.5;
      
      this->addQuadraturePoint( refpoints[ 0 ], volume );
      break;

    case 1:
      refpoints[ 0 ][ 0 ] = 0;
      refpoints[ 1 ][ 0 ] = 1;

      for( int i = 0; i < 2; ++i )
        this->addQuadraturePoint( refpoints[ i ], volume / 2 );
      break;
        
    case 2:
      refpoints[ 0 ][ 0 ] = 0;
      refpoints[ 1 ][ 0 ] = 0.5;
      refpoints[ 2 ][ 0 ] = 1;
    
      for( int i = 0; i < 3; ++i )
        this->addQuadraturePoint( refpoints[ i ], volume / 3 );
      break;
    }
  }


    
  template< class ct >
  LagrangeTriangleQuadrature< ct >
    :: LagrangeTriangleQuadrature( int order, size_t id )
    : QuadratureImp< ct, 2 >( id ),
      order_( (order < 0) ? 0 : order )
  {
    const ctype volume = 0.5;
      
    CoordinateType refpoints[ 6 ];
    
    switch( order_ ) {
    case 0:
      refpoints[ 0 ][ 0 ] = 1.0 / 3.0;
      refpoints[ 0 ][ 1 ] = 1.0 / 3.0;
      
      this->addQuadraturePoint( refpoints[ 0 ], volume );
      break;

    case 1:
      refpoints[ 0 ][ 0 ] = 0;
      refpoints[ 0 ][ 1 ] = 0;

      refpoints[ 1 ][ 0 ] = 1;
      refpoints[ 1 ][ 1 ] = 0;

      refpoints[ 2 ][ 0 ] = 0;
      refpoints[ 2 ][ 1 ] = 1;
      
      for( int i = 0; i < 3; ++i )
        this->addQuadraturePoint( refpoints[ i ], volume / 3 );
      break;
        
    case 2:
      refpoints[ 0 ][ 0 ] = 0;
      refpoints[ 0 ][ 1 ] = 0;

      refpoints[ 1 ][ 0 ] = 1;
      refpoints[ 1 ][ 1 ] = 0;

      refpoints[ 2 ][ 0 ] = 0;
      refpoints[ 2 ][ 1 ] = 1;
        
      refpoints[ 3 ][ 0 ] = 0.5;
      refpoints[ 3 ][ 1 ] = 0.5;

      refpoints[ 4 ][ 0 ] = 0;
      refpoints[ 4 ][ 1 ] = 0.5;

      refpoints[ 5 ][ 0 ] = 0.5;
      refpoints[ 5 ][ 1 ] = 0;
    
      for( int i = 0; i < 6; ++i )
        this->addQuadraturePoint( refpoints[ i ], volume / 6 );
      break;
    }
  }



  template< class ct >
  LagrangeQuadrilateralQuadrature< ct >
    :: LagrangeQuadrilateralQuadrature( int order, size_t id )
    : QuadratureImp< ct, 2 >( id ),
      order_( (order < 0) ? 0 : order )
  {
    const ct volume = 1.0;
   
    CoordinateType refpoints[ 8 ];
    
    switch( order_ ) {
    case 0:
      refpoints[ 0 ][ 0 ] = 0.5;
      refpoints[ 0 ][ 1 ] = 0.5;
      
      this->addQuadraturePoint( refpoints[ 0 ], volume );
      break;

    case 1:
      refpoints[ 0 ][ 0 ] = 0;
      refpoints[ 0 ][ 1 ] = 0;

      refpoints[ 1 ][ 0 ] = 1;
      refpoints[ 1 ][ 1 ] = 0;

      refpoints[ 2 ][ 0 ] = 0;
      refpoints[ 2 ][ 1 ] = 1;

      refpoints[ 3 ][ 0 ] = 1;
      refpoints[ 3 ][ 1 ] = 1;
     
      for( int i = 0; i < 4; ++i )
        this->addQuadraturePoint( refpoints[ i ], volume / 4 );
      break;
        
    case 2:
      refpoints[ 0 ][ 0 ] = 0;
      refpoints[ 0 ][ 1 ] = 0;

      refpoints[ 1 ][ 0 ] = 1;
      refpoints[ 1 ][ 1 ] = 0;

      refpoints[ 2 ][ 0 ] = 0;
      refpoints[ 2 ][ 1 ] = 1;
      
      refpoints[ 3 ][ 0 ] = 1;
      refpoints[ 3 ][ 1 ] = 1;
        
      refpoints[ 4 ][ 0 ] = 0;
      refpoints[ 4 ][ 1 ] = 0.5;

      refpoints[ 5 ][ 0 ] = 1;
      refpoints[ 5 ][ 1 ] = 0.5;

      refpoints[ 6 ][ 0 ] = 0.5;
      refpoints[ 6 ][ 1 ] = 0;
      
      refpoints[ 7 ][ 0 ] = 0.5;
      refpoints[ 7 ][ 1 ] = 1;
    
      for( int i = 0; i < 8; ++i )
        this->addQuadraturePoint( refpoints[ i ], volume / 8 );
      break;
    }
  }

}
