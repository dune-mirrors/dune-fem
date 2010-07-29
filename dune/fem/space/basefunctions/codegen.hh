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

} // end namespace Fem 

} // end namespace Dune

#endif
