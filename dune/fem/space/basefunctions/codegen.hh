#ifndef DUNE_BASEFUNCTIONSETS_CODEGEN_HH
#define DUNE_BASEFUNCTIONSETS_CODEGEN_HH

#include <iostream>

namespace Dune
{

namespace Fem { 

  static void evaluateCodegen(std::ostream& out, 
          const int dimRange, const size_t numRows, const size_t numCols ) 
  {
    out << "template <class BaseFunctionSet>" << std::endl;
    out << "struct EvaluateRanges<BaseFunctionSet, Fem :: EmptyGeometry, " << dimRange << ", " << numRows << ", " << numCols << ">" << std::endl;
    out << "{" << std::endl;
    out << "  template< class QuadratureType,"<< std::endl;
    out << "            class RangeVectorType," << std::endl;
    out << "            class RangeFactorType," << std::endl;
    out << "            class LocalDofVectorType>" << std::endl;
    out << "  static void eval( const QuadratureType& quad," << std::endl;
    out << "                    const RangeVectorType& rangeStorage," << std::endl; 
    out << "                    const LocalDofVectorType& dofs," << std::endl;
    out << "                    RangeFactorType &rangeVector)" << std::endl;
    out << "  {" << std::endl;
    out << "    typedef typename ScalarRangeType :: field_type field_type;" << std::endl;
    out << "    typedef typename RangeVectorType :: value_type value_type;" << std::endl; 
    // make length sse conform 
    const int sseDimRange = (( dimRange % 2 ) == 0) ? dimRange : dimRange+1; 
    //out << "    typedef FieldVector< field_type, " << sseDimRange << " >  ResultType;" << std::endl;
    for( size_t row=0; row< numRows; ++row ) 
    {
      out << "    {" << std::endl;
      out << "      const value_type& rangeStorageRow = rangeStorage[ quad.cachingPoint( " << row << " ) ];" << std::endl;
      out << "      field_type result [ " << sseDimRange << " ] = { ";
      for( int r = 0; r < sseDimRange-1; ++r )  out << " 0 ,";
      out << " 0  };" << std::endl;
      for( size_t col = 0, colR = 0; col < numCols; ++col ) 
      {
        //out << "      {" << std::endl;
        out << "      const field_type phi"<< col << " = rangeStorageRow[ " << col << " ][ 0 ];" << std::endl;
        for( int r = 0; r < dimRange; ++r , ++colR ) 
        {
          out << "      result[ " << r << " ] += dofs[ " << colR << " ] * phi" << col << ";" << std::endl;
          if( sseDimRange != dimRange && r == dimRange - 1 )
          {
            //out << "      result[ " << r+1 << " ] += 0.5 * phi" << col << ";" << std::endl;
            if( (colR + 1) < numCols * dimRange )
              out << "      result[ " << r+1 << " ] += dofs[ " << colR + 1 << " ] * phi" << col << ";" << std::endl;
            else 
              out << "      result[ " << r+1 << " ] += 0.5 * phi" << col << ";" << std::endl;
          }
        }
        //out << "      }" << std::endl;
      }
      out << "      // store result in vector"  << std::endl;
      out << "      RangeType& realResult = rangeVector[ " << row << " ];"  << std::endl;
      for( int r = 0; r < dimRange; ++r) 
      {
        out << "      realResult[ " << r << " ] = result[ " << r << " ];" << std::endl;
      }
      out << "    }" << std::endl;
    }
    out << "  }" << std::endl << std::endl;
    out << "};" << std::endl;
  }

  static void axpyCodegen(std::ostream& out, 
          const int dimRange, const size_t numRows, const size_t numCols ) 
  {
    out << "template <class BaseFunctionSet>" << std::endl;
    out << "struct AxpyRanges<BaseFunctionSet, Fem :: EmptyGeometry, " << dimRange << ", " << numRows << ", " << numCols << ">" << std::endl;
    out << "{" << std::endl;
    out << "  template< class QuadratureType,"<< std::endl;
    out << "            class RangeVectorType," << std::endl;
    out << "            class RangeFactorType," << std::endl;
    out << "            class LocalDofVectorType>" << std::endl;
    out << "  static void axpy( const QuadratureType& quad," << std::endl;
    out << "                    const RangeVectorType& rangeStorage," << std::endl; 
    out << "                    const RangeFactorType &rangeFactors," << std::endl;
    out << "                    LocalDofVectorType& dofs)" << std::endl;
    out << "  {" << std::endl;
    out << "    typedef typename RangeVectorType :: value_type value_type;" << std::endl; 
    out << "    const value_type* const rangeStorageTmp[ " << numRows << " ] = {" << std::endl;
    for( size_t row = 0; row<numRows ; ++ row )
    {
      out << "       &( rangeStorage[ quad.cachingPoint( " << row << " ) ] )";
      if( row < numRows - 1 ) out << " ," << std::endl;
      else out << "  };" << std::endl;
    }
    out << std::endl;
    out << "    typedef typename ScalarRangeType :: field_type field_type;" << std::endl;
    /*
    for( int row=0; row< numRows; ++row ) 
    {
      out << "    {" << std::endl;
      out << "      const size_t baseRow = quad.cachingPoint( " << row << " );" << std::endl;
      out << "      const RangeType& factor = rangeFactors[ " << row << " ];"  << std::endl;
      out << "      const value_type& rangeStorageRow = rangeStorage[ baseRow ];" << std::endl;
      for( size_t col = 0, colR = 0; col < numCols; ++col ) 
      {
        out << "      {" << std::endl;
        out << "        const field_type phi = rangeStorageRow[ " << col << " ][ 0 ];" << std::endl;
        for( int r = 0; r < dimRange; ++r , ++colR ) 
        {
          out << "        dofs[ " << colR << " ] += phi * factor[ " << r << " ];" << std::endl;
        }
        out << "      }" << std::endl;
      }
      out << "    }" << std::endl;
    }
    */
    const int sseDimRange = (( dimRange % 2 ) == 0) ? dimRange : dimRange+1; 
    for( size_t col = 0, colR = 0; col < numCols; ++col ) 
    {
      out << "    {" << std::endl;
      out << "      field_type result [ " << sseDimRange << " ] = { ";
      for( int r = 0; r < sseDimRange-1; ++r )  out << " 0 ,";
      out << " 0  };" << std::endl << std::endl;

      out << "      const field_type phi[ " << numRows << " ] = {" << std::endl;
      for( size_t row=0; row< numRows; ++row ) 
      {
        out << "        (*(rangeStorageTmp[ " << row << " ]))[ " << col << " ][ 0 ]";
        if( row < numRows - 1) 
          out << " ," << std::endl;
        else out << "  };" << std::endl;
      }
      for( size_t row=0; row< numRows; ++row ) 
      {
        out << "      const RangeType& factor" << row  << " = rangeFactors[ " << row << " ];" << std::endl;
      }
      for( size_t row=0; row< numRows; ++row ) 
      {
        //out << "      const field_type phi" << row << " = (*(rangeStorageTmp[ " << row << " ]))[ " << col << " ][ 0 ];" << std::endl;
        for( int r = 0; r < dimRange; ++r ) 
        {
          //out << "      result[ " << r << " ]  +=  factor" << row << "[ " << r << " ] * phi" << row << ";" << std::endl; 
          //if( sseDimRange != dimRange && r == dimRange - 1 )
          //  out << "      result[ " << r+1 << " ]  +=  0.5 * phi" << row << ";" << std::endl;
          out << "      result[ " << r << " ]  +=  factor" << row << "[ " << r << " ] * phi[ " << row << " ];" << std::endl; 
          if( sseDimRange != dimRange && r == dimRange - 1 )
            out << "      result[ " << r+1 << " ]  +=  0.5 * phi[ " << row << " ];" << std::endl;
          //out << "      rangeStorage[ quad.cachingPoint( " << row << " ) ][ " << col << " ][ 0 ] * rangeFactors[ " << row << " ][ " << r << " ]"; 
          /*
          if( row < numRows - 1 ) 
            out << " + " << std::endl;
          else 
            out << " ; " << std::endl;
            */
        }
      }
      for( int r = 0; r < dimRange; ++r , ++colR ) 
        out << "      dofs[ " << colR << " ]  +=  result[ " << r << " ];" << std::endl;

      out << "    }" << std::endl;
    }
    out << "  }" << std::endl << std::endl;
    out << "};" << std::endl;
  }

  static void axpyJacobianCodegen(std::ostream& out, 
          const int dimRange, const size_t numRows, const size_t numCols ) 
  {
    out << "template <class BaseFunctionSet>" << std::endl;
    out << "struct AxpyJacobians<BaseFunctionSet, Fem :: EmptyGeometry, " << dimRange << ", " << numRows << ", " << numCols << ">" << std::endl;
    out << "{" << std::endl;
    out << "  template< class QuadratureType,"<< std::endl;
    out << "            class JacobianRangeVectorType," << std::endl;
    out << "            class JacobianRangeFactorType," << std::endl;
    out << "            class LocalDofVectorType>" << std::endl;
    out << "  static void axpy( const QuadratureType&," << std::endl;
    out << "                    const Fem :: EmptyGeometry&," << std::endl; 
    out << "                    const JacobianRangeVectorType&," << std::endl; 
    out << "                    const JacobianRangeFactorType&," << std::endl;
    out << "                    LocalDofVectorType&)" << std::endl;
    out << "  {" << std::endl;
    out << "    std::cerr << \"ERROR: wrong code generated for VectorialBaseFunctionSet::axpyJacobians\" << std::endl;" << std::endl;
    out << "    abort();" << std::endl;
    out << "  }" << std::endl;
    out << "};" << std::endl << std::endl;
    out << "template <class BaseFunctionSet, class Geometry>" << std::endl;
    out << "struct AxpyJacobians<BaseFunctionSet, Geometry, " << dimRange << ", " << numRows << ", " << numCols << ">" << std::endl;
    out << "{" << std::endl;
    out << "  template< class QuadratureType,"<< std::endl;
    out << "            class JacobianRangeVectorType," << std::endl;
    out << "            class JacobianRangeFactorType," << std::endl;
    out << "            class LocalDofVectorType>" << std::endl;
    out << "  static void axpy( const QuadratureType& quad," << std::endl;
    out << "                    const Geometry& geometry," << std::endl; 
    out << "                    const JacobianRangeVectorType& jacStorage," << std::endl; 
    out << "                    const JacobianRangeFactorType& jacFactors," << std::endl;
    out << "                    LocalDofVectorType& dofs)" << std::endl;
    out << "  {" << std::endl;
    const int dim = 3 ;
    const int sseDimRange = (( dimRange % 2 ) == 0) ? dimRange : dimRange+1; 
    out << "    typedef typename Geometry::Jacobian GeometryJacobianType;" << std::endl;
    out << "    const GeometryJacobianType gjit = geometry.jacobianInverseTransposed( quad.point( 0 ) );" << std::endl << std::endl;
    out << "    typedef typename JacobianRangeVectorType :: value_type  value_type;" << std::endl; 
    out << "    typedef typename JacobianRangeType :: field_type field_type;" << std::endl;
    const size_t dofs = sseDimRange * numCols ;
    out << "    field_type result [ " << dofs << " ] = {"; 
    for( size_t dof = 0 ; dof < dofs-1 ; ++ dof ) out << " 0,";
    out << " 0 };" << std::endl << std::endl;
    out << "    for( size_t row = 0; row < " << numRows << " ; ++ row )" << std::endl;
    out << "    {" << std::endl;
    out << "      const value_type& jacStorageRow = jacStorage[ quad.cachingPoint( row ) ];" << std::endl;
    out << "      field_type jacFactorInv[ " << dim * sseDimRange << " ];" << std::endl;
    out << "      for( int r = 0; r < " << dimRange << " ; ++r )" << std::endl;
    out << "      {"<<std::endl; 
    out << "        for( size_t i = 0; i < GeometryJacobianType :: cols; ++i )" << std::endl;
    out << "        {" << std::endl;
    out << "           const size_t ir = i * " << sseDimRange << " + r ;" << std::endl;
    out << "           jacFactorInv[ ir  ] = 0;" << std::endl;
    out << "           for( size_t j = 0; j < GeometryJacobianType :: rows; ++j )" << std::endl;
    out << "             jacFactorInv[ ir ] += gjit[ j ][ i ] * jacFactors[ row ][ r ][ j ];" << std::endl;
    out << "        }" << std::endl;
    out << "      }" << std::endl;

    for( size_t col = 0, colR = 0; col < numCols; ++col ) 
    {
      for( int d = 0, dr = 0 ; d < dim ; ++ d ) 
      {
        out << "      const field_type phi" << col << d << " = jacStorageRow[ " << col << " ][ 0 ][ " << d << " ];" << std::endl;
        for( int r = 0; r < sseDimRange; ++r, ++dr )
        {
          out << "      result[ " << colR+r << " ]  +=  phi" << col << d << " * jacFactorInv[ " << dr << " ];" << std::endl;
        }
      }
      colR += sseDimRange; 
    }
    out << "    }" << std::endl << std::endl;

    /*
    for( size_t row = 0; row<numRows ; ++ row )
    {
      out << "    {" << std::endl;
      out << "      const value_type& jacStorageRow = jacStorage[ quad.cachingPoint( " << row << " ) ];" << std::endl;
      out << "      JacobianRangeType jacFactorInv;" << std::endl;
      for( int r = 0; r < dimRange ; ++ r ) 
      {
        out << "      gjit.mtv( jacFactors[ " << row << " ][ " << r  << " ], jacFactorInv[ " << r << " ] );" << std::endl;
      }
      for( size_t col = 0, colR = 0; col < numCols; ++col ) 
      {
        for( int r = 0; r < dimRange; ++r, ++colR )
        {
          out << "      dofs[ " << colR << " ]  +=  jacStorageRow[ " << col << " ][ 0 ] * jacFactorInv[ " << r << " ];" << std::endl;
        }
      }
      out << "    }" << std::endl;
    }
    */
    for( size_t col = 0, colD = 0, colR = 0; col < numCols; ++ col ) 
    {
      for(int r = 0; r < dimRange; ++ r, ++colD, ++colR )
      {
        out << "    dofs[ " << colD << " ]  +=  result[ " << colR << " ];" << std::endl; 
      }
      if( sseDimRange != dimRange ) ++ colR ;
    }
    out << "  }" << std::endl << std::endl;
    out << "};" << std::endl;
  }

} // end namespace Fem 

} // end namespace Dune

#endif
