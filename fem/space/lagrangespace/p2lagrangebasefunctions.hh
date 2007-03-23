#ifndef DUNE_P2_LAGRANGE_BASE_FUNCTIONS_HH
#define DUNE_P2_LAGRANGE_BASE_FUNCTIONS_HH

namespace Dune {

//*****************************************************************
//
//! Lagrange base for lines and polynom order = 2
//! (0) 0-----1 (1)
//
//*****************************************************************
template<class FunctionSpaceType>
class LagrangeBaseFunction < FunctionSpaceType , GeometryIdentifier::Line , 2 >  
: public BaseFunctionInterface<FunctionSpaceType> 
{
  typedef typename FunctionSpaceType::DomainType DomainType;
  typedef typename FunctionSpaceType::RangeType RangeType;
  typedef typename FunctionSpaceType::RangeFieldType RangeFieldType;  
  RangeFieldType factor[3];
  
public:
  LagrangeBaseFunction ( int baseNum ) 
    : BaseFunctionInterface< FunctionSpaceType >() 
  {
      switch (baseNum) {
      case 0:
          // 2x^2 - 3x +1
          factor[0] =  2;
          factor[1] = -3;
          factor[2] =  1;
          break;
      case 1:
          // -4x^2 + 4x
          factor[0] = -4;
          factor[1] =  4;
          factor[2] =  0;
          break;
      case 2:
          // 2x^2 - x
          factor[0] =  2;
          factor[1] = -1;
          factor[2] =  0;
          break;
      }
  }  

  //! evaluate the function
    virtual void evaluate ( const FieldVector<deriType, 0> &diffVariable, 
                          const DomainType & x, RangeType & phi) const 
  {
    phi  = factor[0] * x[0] * x[0];
    phi += factor[1] * x[0];
    phi += factor[2];
  }

  //! evaluate first derivative 
    virtual void evaluate ( const FieldVector<deriType, 1> &diffVariable, 
                          const DomainType & x, RangeType & phi) const 
  {
    phi = 2*factor[0] * x[0] + factor[1];
  }

  //! evaluate second derivative 
    virtual void evaluate ( const FieldVector<deriType, 2> &diffVariable, 
                          const DomainType & x, RangeType & phi) const 
  {
    phi = 2*factor[0];
  }

};



/** \brief Lagrange base functions for triangles and polynom order = 2
 *
 *  \internal The base functions are given as follows:
 *  \f{eqnarray*}
 *  \varphi_0( x, y ) &=& 2 (x + y)^2 - 3 (x + y) + 1 \\
 *  \varphi_1( x, y ) &=& 2 x^2 - x \\
 *  \varphi_2( x, y ) &=& 2 y^2 - y \\
 *  \varphi_3( x, y ) &=& 4 x y \\
 *  \varphi_4( x, y ) &=& 4 y (1 - x - y) \\
 *  \varphi_5( x, y ) &=& 4 x (1 - x - y)
 *  \f}
 */
template< class FunctionSpaceType >
class LagrangeBaseFunction < FunctionSpaceType , GeometryIdentifier :: Triangle, 2 >
: public BaseFunctionInterface< FunctionSpaceType >
{
public:
  typedef typename FunctionSpaceType :: DomainType DomainType;
  typedef typename FunctionSpaceType :: DomainFieldType DomainFieldType;
  typedef typename FunctionSpaceType :: RangeType RangeType;
  typedef typename FunctionSpaceType :: RangeFieldType RangeFieldType;
  
private:
  RangeFieldType factor[ 6 ];
  
public:
  LagrangeBaseFunction ( int dof )
    : BaseFunctionInterface< FunctionSpaceType >()
  {
    switch( dof ) {
      case 0: // 2 (x + y)^2 - 3 (x + y) + 1
        factor[ 0 ] = 1;
        factor[ 1 ] = -3;
        factor[ 2 ] = -3;
        factor[ 3 ] = 2;
        factor[ 4 ] = 4;
        factor[ 5 ] = 2;
        break;

      case 1: // 2 x^2 - x
        factor[ 0 ] = 0;
        factor[ 1 ] = -1;
        factor[ 2 ] = 0;
        factor[ 3 ] = 2;
        factor[ 4 ] = 0;
        factor[ 5 ] = 0;
        break;

      case 2: // 2 y^2 - y
        factor[ 0 ] = 0;
        factor[ 1 ] = 0;
        factor[ 2 ] = -1;
        factor[ 3 ] = 0;
        factor[ 4 ] = 0;
        factor[ 5 ] = 2;
        break;

      case 3: // 4 x y
        factor[ 0 ] = 0;
        factor[ 1 ] = 0;
        factor[ 2 ] = 0;
        factor[ 3 ] = 0;
        factor[ 4 ] = 4;
        factor[ 5 ] = 0;
        break;

      case 4: // -4 y^2 - 4 x y + 4 y
        factor[ 0 ] = 0;
        factor[ 1 ] = 0;
        factor[ 2 ] = 4;
        factor[ 3 ] = 0;
        factor[ 4 ] = -4;
        factor[ 5 ] = -4;
        break;

      case 5: // -4 x^2 - 4 x y + 4 x
        factor[ 0 ] = 0;
        factor[ 1 ] = 4;
        factor[ 2 ] = 0;
        factor[ 3 ] = -4;
        factor[ 4 ] = -4;
        factor[ 5 ] = 0;
        break;

      default:
        DUNE_THROW( RangeError, "Invalid local DoF number." );
    }
  }  

  //! evaluate the function
  virtual void evaluate ( const FieldVector< deriType, 0 > &diffVariable, 
                          const DomainType &x,
                          RangeType &phi
                        ) const 
  {
    const DomainFieldType &x1 = x[ 0 ];
    const DomainFieldType &x2 = x[ 1 ];
      
    phi = factor[ 0 ];
    phi += factor[ 1 ] * x1 + factor[ 2 ] * x2;
    phi += factor[ 3 ] * x1 * x1 + factor[ 4 ] * x1 * x2 + factor[ 5 ] * x2 * x2;
  }

  //! evaluate first derivative 
  virtual void evaluate ( const FieldVector< deriType, 1 > &diffVariable, 
                          const DomainType &x,
                          RangeType &phi
                        ) const 
  {
    const DomainFieldType &x1 = x[ 0 ];
    const DomainFieldType &x2 = x[ 1 ];
 
    switch( diffVariable[ 0 ] ) {
      case 0:
        phi = factor[ 1 ];
        phi += 2 * factor[ 3 ] * x1 + factor[ 4 ] * x2;
        break;

      case 1:
        phi = factor[ 2 ];
        phi += 2 * factor[ 5 ] * x2 + factor[ 4 ] * x1;
        break;

      default:
        DUNE_THROW( RangeError, "Invalid derivative requested." );
    }
  }

  //! evaluate second derivative 
  virtual void evaluate ( const FieldVector< deriType, 2 > &diffVariable, 
                          const DomainType &x,
                          RangeType &phi
                        ) const 
  {
    assert( (diffVariable[ 0 ] >= 0) && (diffVariable[ 0 ] < 2) );
    assert( (diffVariable[ 1 ] >= 0) && (diffVariable[ 1 ] < 2) );

    const deriType d = diffVariable[ 0 ] + diffVariable[ 1 ];
    phi = factor[ 3 + d ];
    if( d != 1 )
        phi *= 2;
  }
};



/** \brief Lagrange base functions for quadrilaterals and polynom order = 2
 *
 *  \internal The base functions are given as follows:
 *  \f{eqnarray*}
 *  \varphi_0( x, y ) &=& (1 - 2 (x +y)) x y + 2 (x + y)^2 - 3 (x + y) + 1 \\
 *  \varphi_1( x, y ) &=& 2 (1 - y) x^2 + x y (2 y - 1) - x \\
 *  \varphi_2( x, y ) &=& 2 (1 - x) y^2 + x y (2 x - 1) - y \\
 *  \varphi_3( x, y ) &=& x y (2 x + 2 y - 3) \\
 *  \varphi_4( x, y ) &=& 4 y (1 - x) (1 - y) \\
 *  \varphi_5( x, y ) &=& 4 x y (1 - y) \\
 *  \varphi_6( x, y ) &=& 4 x (1 - x) (1 - y) \\
 *  \varphi_7( x, y ) &=& 4 x y (1 - x)
 *  \f}
 */
template< class FunctionSpaceType >
class LagrangeBaseFunction < FunctionSpaceType , GeometryIdentifier :: Quadrilateral, 2 >
: public BaseFunctionInterface< FunctionSpaceType >
{
public:
  typedef typename FunctionSpaceType :: DomainType DomainType;
  typedef typename FunctionSpaceType :: DomainFieldType DomainFieldType;
  typedef typename FunctionSpaceType :: RangeType RangeType;
  typedef typename FunctionSpaceType :: RangeFieldType RangeFieldType;  

private:
  RangeFieldType factor[ 8 ];
  
public:
  LagrangeBaseFunction ( int dof )
    : BaseFunctionInterface< FunctionSpaceType >()
  {
    switch( dof )
    {
      case 0: // (1 - 2(x+y)) x y + 2 (x+y)^2 - 3(x+y) + 1
        factor[ 0 ] = 1;   // 1
        factor[ 1 ] = -3;  // x
        factor[ 2 ] = -3;  // y
        factor[ 3 ] = 2;   // x^2
        factor[ 4 ] = 5;   // x y
        factor[ 5 ] = 2;   // y^2
        factor[ 6 ] = -2;  // x^2 y
        factor[ 7 ] = -2;  // x y^2
        break;

      case 1: // 2(1 - y) x^2 + x y (2y - 1) - x
        factor[ 0 ] = 0;   // 1
        factor[ 1 ] = -1;  // x
        factor[ 2 ] = 0;   // y
        factor[ 3 ] = 2;   // x^2
        factor[ 4 ] = -1;  // x y
        factor[ 5 ] = 0;   // y^2
        factor[ 6 ] = -2;  // x^2 y
        factor[ 7 ] = 2;   // x y^2
        break;

      case 2: // 2(1 - x) y^2 + x y (2x - 1) - y
        factor[ 0 ] = 0;   // 1
        factor[ 1 ] = 0;   // x
        factor[ 2 ] = -1;  // y
        factor[ 3 ] = 0;   // x^2
        factor[ 4 ] = -1;  // x y
        factor[ 5 ] = 2;   // y^2
        factor[ 6 ] = 2;   // x^2 y
        factor[ 7 ] = -2;  // x y^2
        break;

      case 3: // x y (2x + 2y - 3)
        factor[ 0 ] = 0;   // 1
        factor[ 1 ] = 0;   // x
        factor[ 2 ] = 0;   // y
        factor[ 3 ] = 0;   // x^2
        factor[ 4 ] = -3;  // x y
        factor[ 5 ] = 0;   // y^2
        factor[ 6 ] = 2;   // x^2 y
        factor[ 7 ] = 2;   // x y^2
        break;

      case 4: // 4 y (1 -x) - 4 y^2 (1 - x)
        factor[ 0 ] = 0;   // 1
        factor[ 1 ] = 0;   // x
        factor[ 2 ] = 4;   // y
        factor[ 3 ] = 0;   // x^2
        factor[ 4 ] = -4;  // x y
        factor[ 5 ] = -4;  // y^2
        factor[ 6 ] = 0;   // x^2 y
        factor[ 7 ] = 4;   // x y^2
        break;

      case 5: // 4 x y (1 - y)
        factor[ 0 ] = 0;   // 1
        factor[ 1 ] = 0;   // x
        factor[ 2 ] = 0;   // y
        factor[ 3 ] = 0;   // x^2
        factor[ 4 ] = 4;   // x y
        factor[ 5 ] = 0;   // y^2
        factor[ 6 ] = 0;   // x^2 y
        factor[ 7 ] = -4;  // x y^2
        break;

      case 6: // 4 x (1 - y) - 4 x^2 (1 - y)
        factor[ 0 ] = 0;   // 1
        factor[ 1 ] = 4;   // x
        factor[ 2 ] = 0;   // y
        factor[ 3 ] = -4;  // x^2
        factor[ 4 ] = -4;  // x y
        factor[ 5 ] = 0;   // y^2
        factor[ 6 ] = 4;   // x^2 y
        factor[ 7 ] = 0;   // x y^2
        break;

      case 7: // 4 x y (1 - x)
        factor[ 0 ] = 0;   // 1
        factor[ 1 ] = 0;   // x
        factor[ 2 ] = 0;   // y
        factor[ 3 ] = 0;   // x^2
        factor[ 4 ] = 4;   // x y
        factor[ 5 ] = 0;   // y^2
        factor[ 6 ] = -4;  // x^2 y
        factor[ 7 ] = 0;   // x y^2
        break;

      default:
        DUNE_THROW( RangeError, "Invalid local DoF number." );
    }
  }  

  //! evaluate the function
  virtual void evaluate ( const FieldVector< deriType, 0 > &diffVariable, 
                          const DomainType &x,
                          RangeType &phi
                        ) const
  {
    const DomainFieldType &x1 = x[ 0 ];
    const DomainFieldType &x2 = x[ 1 ];
      
    phi = factor[ 0 ];
    phi += factor[ 1 ] * x1 + factor[ 2 ] * x2;
    phi += factor[ 3 ] * x1 * x1 + factor[ 5 ] * x2 * x2;
    phi += x1 * x2 * (factor[ 4 ] + factor[ 6 ] * x1 + factor[ 7 ] * x2);
  }

  //! evaluate first derivative 
  virtual void evaluate ( const FieldVector< deriType, 1 > &diffVariable, 
                          const DomainType &x,
                          RangeType &phi
                        ) const
  {
    const DomainFieldType &x1 = x[ 0 ];
    const DomainFieldType &x2 = x[ 1 ];
 
    switch( diffVariable[ 0 ] ) {
      case 0:
        phi = factor[ 1 ];
        phi += 2 * factor[ 3 ] * x1 + factor[ 4 ] * x2;
        phi += 2 * factor[ 6 ] * x1 * x2 + factor[ 7 ] * x2 * x2;
        break;

      case 1:
        phi = factor[ 2 ];
        phi += 2 * factor[ 5 ] * x2 + factor[ 4 ] * x1;
        phi += 2 * factor[ 7 ] * x1 * x2 + factor[ 6 ] * x1 * x1;
        break;

      default:
        DUNE_THROW( RangeError, "Invalid derivative requested." );
    }
  }

  //! evaluate second derivative 
  virtual void evaluate ( const FieldVector< deriType, 2 > &diffVariable, 
                          const DomainType &x,
                          RangeType &phi
                        ) const
  {
    const DomainFieldType &x1 = x[ 0 ];
    const DomainFieldType &x2 = x[ 1 ];
 
    switch( diffVariable[ 0 ] + diffVariable[ 1 ] ) {
      case 0:
        assert( diffVariable[ 1 ] == 2 );
        phi = 2 * factor[ 5 ] + factor[ 7 ] * x1;
        break;

      case 1:
        assert( diffVariable[ 1 ] == 1 );
        phi = factor[ 4 ] + 2 * (factor[ 6 ] * x1 + factor[ 7 ] * x2);
        break;

      case 2:
        assert( diffVariable[ 1 ] == 0 );
        phi = 2 * factor[ 3 ] + factor[ 6 ] * x2;
        break;

      default:
        DUNE_THROW( RangeError, "Invalid derivative requested:" );
    }
  }
};



/** \brief Lagrange base functions for tetrahedra and polynom order = 2
 *
 *  \internal The base functions are given as follows:
 *  \f{eqnarray*}
 *  \varphi_0( x, y ) &=& 2 (x + y + z)^2 - 3 (x + y + z) + 1 \\
 *  \varphi_1( x, y ) &=& 2 x^2 - x \\
 *  \varphi_2( x, y ) &=& 2 y^2 - y \\
 *  \varphi_3( x, y ) &=& 2 z^2 - z \\
 *  \varphi_4( x, y ) &=& 4 x (1 - x - y - z) \\
 *  \varphi_5( x, y ) &=& 4 x y \\
 *  \varphi_6( x, y ) &=& 4 y (1 - x - y - z) \\
 *  \varphi_7( x, y ) &=& 4 z (1 - x - y - z) \\
 *  \varphi_8( x, y ) &=& 4 x z \\
 *  \varphi_9( x, y ) &=& 4 y z
 *  \f}
 */
template< class FunctionSpaceType >
class LagrangeBaseFunction < FunctionSpaceType , GeometryIdentifier :: Tetrahedron, 2 >
: public BaseFunctionInterface< FunctionSpaceType >
{
public:
  typedef typename FunctionSpaceType :: DomainType DomainType;
  typedef typename FunctionSpaceType :: DomainFieldType DomainFieldType;
  typedef typename FunctionSpaceType :: RangeType RangeType;
  typedef typename FunctionSpaceType :: RangeFieldType RangeFieldType;
  
private:
  RangeFieldType factor[ 10 ];
  
public:
  LagrangeBaseFunction ( int dof )
    : BaseFunctionInterface< FunctionSpaceType >()
  {
    switch( dof ) {
      case 0: // 2 (x + y + z)^2 - 3 (x + y + z) + 1
        factor[ 0 ] = 1;  // 1
        factor[ 1 ] = -3; // x
        factor[ 2 ] = -3; // y
        factor[ 3 ] = -3; // z
        factor[ 4 ] = 2;  // x^2
        factor[ 5 ] = 2;  // y^2
        factor[ 6 ] = 2;  // z^2
        factor[ 7 ] = 4;  // x y
        factor[ 8 ] = 4;  // x z
        factor[ 9 ] = 4;  // y z
        break;

      case 1: // 2 x^2 - x
        factor[ 0 ] = 0;  // 1
        factor[ 1 ] = -1; // x
        factor[ 2 ] = 0;  // y
        factor[ 3 ] = 0;  // z
        factor[ 4 ] = 2;  // x^2
        factor[ 5 ] = 0;  // y^2
        factor[ 6 ] = 0;  // z^2
        factor[ 7 ] = 0;  // x y
        factor[ 8 ] = 0;  // x z
        factor[ 9 ] = 0;  // y z
        break;

      case 2: // 2 y^2 - y
        factor[ 0 ] = 0;  // 1
        factor[ 1 ] = 0;  // x
        factor[ 2 ] = -1; // y
        factor[ 3 ] = 0;  // z
        factor[ 4 ] = 0;  // x^2
        factor[ 5 ] = 2;  // y^2
        factor[ 6 ] = 0;  // z^2
        factor[ 7 ] = 0;  // x y
        factor[ 8 ] = 0;  // x z
        factor[ 9 ] = 0;  // y z
        break;
        
      case 3: // 2 z^2 - z
        factor[ 0 ] = 0;  // 1
        factor[ 1 ] = 0;  // x
        factor[ 2 ] = 0;  // y
        factor[ 3 ] = -1; // z
        factor[ 4 ] = 0;  // x^2
        factor[ 5 ] = 0;  // y^2
        factor[ 6 ] = 2;  // z^2
        factor[ 7 ] = 0;  // x y
        factor[ 8 ] = 0;  // x z
        factor[ 9 ] = 0;  // y z
        break;

      case 4: // 4 x (1 - x - y - z)
        factor[ 0 ] = 0;  // 1
        factor[ 1 ] = 4;  // x
        factor[ 2 ] = 0;  // y
        factor[ 3 ] = 0;  // z
        factor[ 4 ] = -4; // x^2
        factor[ 5 ] = 0;  // y^2
        factor[ 6 ] = 0;  // z^2
        factor[ 7 ] = -4; // x y
        factor[ 8 ] = -4; // x z
        factor[ 9 ] = 0;  // y z
        break;
        
      case 5: // 4 x y
        factor[ 0 ] = 0;  // 1
        factor[ 1 ] = 0;  // x
        factor[ 2 ] = 0;  // y
        factor[ 3 ] = 0;  // z
        factor[ 4 ] = 0;  // x^2
        factor[ 5 ] = 0;  // y^2
        factor[ 6 ] = 0;  // z^2
        factor[ 7 ] = 4;  // x y
        factor[ 8 ] = 0;  // x z
        factor[ 9 ] = 0;  // y z
        break;
       
      case 6: // 4 y (1 - x - y - z)
        factor[ 0 ] = 0;  // 1
        factor[ 1 ] = 0;  // x
        factor[ 2 ] = 4;  // y
        factor[ 3 ] = 0;  // z
        factor[ 4 ] = 0;  // x^2
        factor[ 5 ] = -4; // y^2
        factor[ 6 ] = 0;  // z^2
        factor[ 7 ] = -4; // x y
        factor[ 8 ] = 0;  // x z
        factor[ 9 ] = -4; // y z
        break;

      case 7: // 4 z (1 - x - y - z)
        factor[ 0 ] = 0;  // 1
        factor[ 1 ] = 0;  // x
        factor[ 2 ] = 0;  // y
        factor[ 3 ] = 4;  // z
        factor[ 4 ] = 0;  // x^2
        factor[ 5 ] = 0;  // y^2
        factor[ 6 ] = -4; // z^2
        factor[ 7 ] = 0;  // x y
        factor[ 8 ] = -4; // x z
        factor[ 9 ] = -4; // y z
        break;

      case 8: // 4 x z
        factor[ 0 ] = 0;  // 1
        factor[ 1 ] = 0;  // x
        factor[ 2 ] = 0;  // y
        factor[ 3 ] = 0;  // z
        factor[ 4 ] = 0;  // x^2
        factor[ 5 ] = 0;  // y^2
        factor[ 6 ] = 0;  // z^2
        factor[ 7 ] = 0;  // x y
        factor[ 8 ] = 4;  // x z
        factor[ 9 ] = 0;  // y z
        break;
 
      case 9: // 4 y z
        factor[ 0 ] = 0;  // 1
        factor[ 1 ] = 0;  // x
        factor[ 2 ] = 0;  // y
        factor[ 3 ] = 0;  // z
        factor[ 4 ] = 0;  // x^2
        factor[ 5 ] = 0;  // y^2
        factor[ 6 ] = 0;  // z^2
        factor[ 7 ] = 0;  // x y
        factor[ 8 ] = 0;  // x z
        factor[ 9 ] = 4;  // y z
        break;

      default:
        DUNE_THROW( RangeError, "Invalid local DoF number." );
    }
  }  

  //! evaluate the function
  virtual void evaluate ( const FieldVector< deriType, 0 > &diffVariable, 
                          const DomainType &x,
                          RangeType &phi
                        ) const 
  {
    phi = factor[ 0 ];
    for( int i = 0; i < 3; ++i ) {
        const DomainFieldType &xi = x[ i ];
        phi += (factor[ i + 4 ] * xi + factor[ i + 1 ]) * xi;
    }
    phi += factor[ 7 ] * x[ 0 ] * x[ 1 ];
    phi += factor[ 8 ] * x[ 0 ] * x[ 2 ];
    phi += factor[ 9 ] * x[ 1 ] * x[ 2 ];
  }

  //! evaluate first derivative 
  virtual void evaluate ( const FieldVector< deriType, 1 > &diffVariable, 
                          const DomainType &x,
                          RangeType &phi
                        ) const 
  {
    const DomainFieldType &x1 = x[ 0 ];
    const DomainFieldType &x2 = x[ 1 ];
    const DomainFieldType &x3 = x[ 2 ];
 
    switch( diffVariable[ 0 ] ) {
      case 0:
        phi = factor[ 1 ];
        phi += 2 * factor[ 4 ] * x1 + factor[ 7 ] * x2 + factor[ 8 ] * x3;
        break;

      case 1:
        phi = factor[ 2 ];
        phi += 2 * factor[ 5 ] * x2 + factor[ 7 ] * x1 + factor[ 9 ] * x3;
        break;

      case 2:
        phi = factor[ 3 ];
        phi += 2 * factor[ 6 ] * x3 + factor[ 8 ] * x1 + factor[ 9 ] * x2;
        break;

      default:
        DUNE_THROW( RangeError, "Invalid derivative requested." );
    }
  }

  //! evaluate second derivative 
  virtual void evaluate ( const FieldVector< deriType, 2 > &diffVariable, 
                          const DomainType &x,
                          RangeType &phi
                        ) const 
  {
    assert( (diffVariable[ 0 ] >= 0) && (diffVariable[ 0 ] < 3) );
    assert( (diffVariable[ 1 ] >= 0) && (diffVariable[ 1 ] < 3) );

    if( diffVariable[ 0 ] == diffVariable[ 1 ] )
        phi = 2 * factor[ 4 + diffVariable[ 0 ] ];
    else
        phi = factor[ 6 + diffVariable[ 0 ] + diffVariable[ 1 ] ];
  }
};



//*****************************************************************
//
//! Lagrange base for tetrahedra and polynom order = 2
//
//*****************************************************************
template< class FunctionSpaceType >
class LagrangeBaseFunction < FunctionSpaceType, GeometryIdentifier :: Pyramid, 2 >  
: public BaseFunctionInterface< FunctionSpaceType >
{
  typedef typename FunctionSpaceType::DomainType DomainType;
  typedef typename FunctionSpaceType::RangeType RangeType;
  typedef typename FunctionSpaceType::RangeFieldType RangeFieldType;  
  RangeFieldType factor[3];
  
public:
  LagrangeBaseFunction ( int baseNum )
    : BaseFunctionInterface<FunctionSpaceType>()
  {
      DUNE_THROW(NotImplemented, "Second order Lagrange elements for pyramids are not implemented yet!");
  }  

  //! evaluate the function
    virtual void evaluate ( const FieldVector<deriType, 0> &diffVariable, 
                          const DomainType & x, RangeType & phi) const 
    {}

  //! evaluate first derivative 
    virtual void evaluate ( const FieldVector<deriType, 1> &diffVariable, 
                          const DomainType & x, RangeType & phi) const 
    {}

    //! evaluate second derivative 
    virtual void evaluate ( const FieldVector<deriType, 2> &diffVariable, 
                            const DomainType & x, RangeType & phi) const 
    {}

};



//*****************************************************************
//
//! Lagrange base for prisms and polynom order = 2
//
//*****************************************************************
template< class FunctionSpaceType >
class LagrangeBaseFunction < FunctionSpaceType, GeometryIdentifier :: Prism, 2 >  
: public BaseFunctionInterface< FunctionSpaceType >
{
  typedef typename FunctionSpaceType::DomainType DomainType;
  typedef typename FunctionSpaceType::RangeType RangeType;
  typedef typename FunctionSpaceType::RangeFieldType RangeFieldType;  
  RangeFieldType factor[3];
  
public:
  LagrangeBaseFunction ( int baseNum ) 
    : BaseFunctionInterface< FunctionSpaceType >()
  {
      DUNE_THROW(NotImplemented, "Second order Lagrange elements for tetrahedra are not implemented yet!");
  }  

  //! evaluate the function
    virtual void evaluate ( const FieldVector<deriType, 0> &diffVariable, 
                          const DomainType & x, RangeType & phi) const 
    {}

  //! evaluate first derivative 
    virtual void evaluate ( const FieldVector<deriType, 1> &diffVariable, 
                          const DomainType & x, RangeType & phi) const 
    {}

    //! evaluate second derivative 
    virtual void evaluate ( const FieldVector<deriType, 2> &diffVariable, 
                            const DomainType & x, RangeType & phi) const 
    {}

};

//*****************************************************************
//
//! Lagrange base for hexahedra and polynom order = 2
//
//*****************************************************************
template< class FunctionSpaceType >
class LagrangeBaseFunction < FunctionSpaceType, GeometryIdentifier :: Hexahedron, 2 >  
: public BaseFunctionInterface< FunctionSpaceType >
{
  typedef typename FunctionSpaceType::DomainType DomainType;
  typedef typename FunctionSpaceType::RangeType RangeType;
  typedef typename FunctionSpaceType::RangeFieldType RangeFieldType;  
  RangeFieldType factor[3];
  
public:
  LagrangeBaseFunction ( int baseNum )
    : BaseFunctionInterface< FunctionSpaceType >() 
  {
  }  

  //! evaluate the function
    virtual void evaluate ( const FieldVector<deriType, 0> &diffVariable, 
                          const DomainType & x, RangeType & phi) const 
    {
        DUNE_THROW(NotImplemented, "Second order Lagrange elements for hexahedra are not implemented yet!");
    }

  //! evaluate first derivative 
    virtual void evaluate ( const FieldVector<deriType, 1> &diffVariable, 
                          const DomainType & x, RangeType & phi) const 
    {
      DUNE_THROW(NotImplemented, "Second order Lagrange elements for hexahedra are not implemented yet!");
    }

    //! evaluate second derivative 
    virtual void evaluate ( const FieldVector<deriType, 2> &diffVariable, 
                            const DomainType & x, RangeType & phi) const 
    {
      DUNE_THROW(NotImplemented, "Second order Lagrange elements for hexahedra are not implemented yet!");
    }

};




}  // end namespace dune

#endif
