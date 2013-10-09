#ifndef DUNE_FEM_BASISFUNCTIONSETS_CODEGEN_HH
#define DUNE_FEM_BASISFUNCTIONSETS_CODEGEN_HH

#include <cassert>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <set>
#include <map>

#include <dune/common/exceptions.hh>
#include <dune/fem/io/io.hh>

namespace Dune
{

  namespace Fem 
  { 

    class CodegenInfoFinished : public Dune :: Exception {};

    class CodegenInfo
    {
      std::string path_;
      std::string outFileName_;

      int nopMax_;
      int nopMin_;

      int baseMax_;
      int baseMin_;

      typedef std::set< int > DimRangeSetType;
      typedef std::set< void* > BaseSetPointerType;
      DimRangeSetType savedimRanges_; 
      mutable DimRangeSetType dimRanges_; 
      BaseSetPointerType savedBaseSets_;

      typedef void codegenfunc_t (std::ostream& out, 
                                  const int dim, 
                                  const int dimRange, 
                                  const size_t numRows, 
                                  const size_t numCols );

      typedef std::pair< std::string, int > EntryType;
      std::vector< EntryType > filenames_;

      typedef std::vector< int > KeyType;
      typedef std::set< KeyType > KeyStorageType;

      typedef std::pair< int , int > EvalPairType ;
      typedef std::set< EvalPairType > EvalSetType ;

      std::map< codegenfunc_t* , KeyStorageType > entryMap_;
      EvalSetType  evalSet_;

      CodegenInfo() 
        : path_("./autogeneratedcode"),
          outFileName_( "autogeneratedcode.hh" ),
          nopMax_(0), nopMin_(0), baseMax_(0), baseMin_(0), 
          dimRanges_(), savedBaseSets_(), filenames_(), 
          entryMap_(),
          evalSet_()
      {
      }

    public: 
      static CodegenInfo& instance() 
      {
        static CodegenInfo info;
        return info;
      }

      //! clear all registered entries 
      void clear() {
        savedBaseSets_.clear();
        filenames_.clear();
        dimRanges_.clear();
        entryMap_.clear();
        evalSet_.clear();
        nopMax_ = 0; 
        nopMin_ = 0;
        baseMax_ = 0; 
        baseMin_ = 0;
      }

      //! overwrite path 
      void setPath(const std::string& path ) { path_ = path ; }

      //! overwrite filename 
      void setFileName(const std::string& filename ) { outFileName_ = filename ; }

      template <class BaseSet>
      void addDimRange(const BaseSet* baseSet, 
                       const int dimRange) 
      {
        void* ptr = (void *) baseSet;
        if( savedBaseSets_.find( ptr ) == savedBaseSets_.end() )
        {
          savedBaseSets_.insert( ptr );
          std::cout << "Add dimRange " << dimRange << std::endl;
          dimRanges_.insert( dimRange ) ;
        }
      }

      void addEntry(const std::string& fileprefix, 
                      codegenfunc_t* codegenfunc,
                      const int dim, const int dimRange, const int quadNop, const int numBase ) 
      {
        KeyStorageType& keyMap = entryMap_[ codegenfunc ];

        typedef KeyStorageType :: iterator iterator; 
        KeyType key( 4 );
        key[ 0 ] = dim; 
        key[ 1 ] = dimRange;
        key[ 2 ] = quadNop;
        key[ 3 ] = numBase;

        // search for key, if already exists, then do nothing 
        iterator it = keyMap.find( key );
        if( it != keyMap.end() ) return ;

        // make sure dimRange was registered 
        assert( dimRanges_.find( dimRange ) != dimRanges_.end() );

        // create directory to store files 
        createDirectory( path_ );

        std::stringstream filename;
        filename << fileprefix << dimRange << "_" << quadNop << "_" << numBase << ".hh";

        // second check, if file exists, do nothing
        int pos = exists( filename.str() ); 
        if( pos >= 0 )  return  ;

        std::string filenameToWrite( path_ + "/" + filename.str() );
        std::ofstream file( filenameToWrite.c_str(), std::ios::out );

        // call code generation function 
        codegenfunc( file, dim, dimRange, quadNop, numBase );
        std::cout << "Generate code " << fileprefix << " for (" << dimRange << "," 
                  << quadNop << "," << numBase << ")" << std::endl;

        // insert evaluate pair 
        EvalPairType evalPair ( quadNop, numBase );
        evalSet_.insert( evalPair );

        if( baseMin_ == 0 ) baseMin_ = numBase;
        if( nopMin_  == 0 ) nopMin_  = quadNop;

        EntryType entry ( filename.str() , dimRange );
        filenames_.push_back( entry );
        nopMax_ = std::max( quadNop, nopMax_ );
        nopMin_ = std::min( quadNop, nopMin_ );
        baseMax_ = std::max( numBase, baseMax_ );
        baseMin_ = std::min( numBase, baseMin_ );
      }

      void finalize () const 
      {
        if( checkAbort() ) 
        {
          dumpInfo();
          std::cerr << "All automated code generated, bye, bye !! " << std::endl;
          DUNE_THROW(CodegenInfoFinished,"All automated code generated, bye, bye !!");
        }
      }


      void dumpInfo() const 
      {
        {
          std::ofstream file( outFileName_.c_str() );
          file << "#include \"" <<path_<< "/" << outFileName_ << "\"";
        }

        {
          std::string filename( path_ + "/" + outFileName_ );
          std::ofstream file( filename.c_str() );
          file << "#ifdef CODEGEN_INCLUDEMAXNUMS" << std::endl;
          file << "#ifndef CODEGEN_INCLUDEMAXNUMS_INCLUDED" << std::endl;
          file << "#define CODEGEN_INCLUDEMAXNUMS_INCLUDED" << std::endl << std::endl;
          file << "#define MAX_NUMBER_OF_QUAD_POINTS " << nopMax_ << std::endl;
          file << "#define MIN_NUMBER_OF_QUAD_POINTS " << nopMin_ << std::endl;
          file << "#define MAX_NUMBER_OF_BASE_FCT    " << baseMax_ << std::endl;
          file << "#define MIN_NUMBER_OF_BASE_FCT    " << baseMin_ << std::endl << std::endl;
          file << "/* include all headers with inner loop extern declarations */" << std::endl;
          file << "#define CODEGEN_COMPILE_INNERLOOPS 1" << std::endl;
          for( size_t i = 0; i < filenames_.size(); ++i ) 
          {
            file << "#include \""<< filenames_[ i ].first << "\"" << std::endl;
          }
          file << "#undef CODEGEN_COMPILE_INNERLOOPS" << std::endl << std::endl;
          file << "#include \"" << filename << ".c\"" <<std::endl;

          file << "#endif // CODEGEN_INCLUDEMAXNUMS_INCLUDED" << std::endl << std::endl;
          file << "#elif defined CODEGEN_INCLUDEEVALCALLERS" << std::endl;
          file << "#ifndef CODEGEN_EVALCALLERS_INCLUDED" << std::endl;
          file << "#define CODEGEN_EVALCALLERS_INCLUDED" << std::endl << std::endl;
          typedef EvalSetType :: iterator iterator ;
          const iterator endit = evalSet_.end();
          for( iterator it = evalSet_.begin(); it != endit; ++it ) 
          {
            file << "  template <class Traits>" << std::endl;
            file << "  struct EvaluateImplementation< Traits, " << it->first << " , " << it->second << " >" << std::endl;
            file << "    : public EvaluateRealImplementation< Traits, " << it->first << " , " << it->second << " >" << std::endl;
            file << "  {" << std::endl;
            file << "    typedef EvaluateRealImplementation< Traits, " << it->first << " , " << it->second << " >  BaseType;" << std::endl;
            file << "    typedef typename BaseType :: RangeVectorType  RangeVectorType;" << std::endl;
            file << "    EvaluateImplementation( const RangeVectorType& rv ) : BaseType ( rv ) {}" << std::endl;
            file << "  };"<< std::endl;
            file << std::endl;
          }
          file << "#endif // CODEGEN_EVALCALLERS_INCLUDED" << std::endl << std::endl;
          file << "#else" << std::endl << std::endl ;
          file << "#ifndef CODEGEN_INCLUDE_IMPLEMENTATION" << std::endl;
          file << "#define CODEGEN_INCLUDE_IMPLEMENTATION" << std::endl;
          file << "#undef CODEGEN_COMPILE_INNERLOOPS" << std::endl;
          for( size_t i = 0; i < filenames_.size(); ++i ) 
          {
            file << "#include \""<< filenames_[ i ].first << "\"" << std::endl;
          }
          file << "#endif  // CODEGEN_INCLUDE_IMPLEMENTATION" << std::endl << std::endl;
          file << "#endif // CODEGEN_INCLUDEMAXNUMS" << std::endl;

          // write C file with implementation of inner loop functions 
          filename += ".c";
          std::ofstream Cfile( filename.c_str() );

          Cfile << "#include <stdlib.h>" << std::endl;
          Cfile << "/* include all headers with inner loop implementation */" << std::endl;
          Cfile << "#define CODEGEN_COMPILE_INNERLOOPS 2" << std::endl;
          for( size_t i = 0; i < filenames_.size(); ++i ) 
          {
            Cfile << "#include \""<< filenames_[ i ].first << "\"" << std::endl;
          }
        }
      }
    protected:  
      int exists( const std::string& filename ) const 
      {
        for( size_t i = 0; i < filenames_.size(); ++i ) 
        {
          if( filename == filenames_[ i ].first )
          {
            return i;
          }
        }
        return -1;
      }

      bool checkAbort() const 
      {
        DimRangeSetType found ;
        bool canAbort = true ;
        for( size_t i = 0; i < filenames_.size(); ++i ) 
        {
          found.insert( filenames_[ i ].second );
          if ( filenames_[ i ].second > 0 ) 
          {
            canAbort = false ;
          }
        }
        typedef DimRangeSetType :: iterator iterator ;
        for( iterator it = found.begin(); it != found.end(); ++it ) 
        {
          dimRanges_.erase( *it );
        }

        if( canAbort && dimRanges_.size() == 0 )
          return true ; 
        else 
          return false ;
      }
    };

    /** \brief default code generator methods */
    struct DefaultCodeGenerator 
    {
      static void evaluateCodegen(std::ostream& out, 
                                  const int dim, 
                                  const int dimRange, 
                                  const size_t numRows, 
                                  const size_t numCols ) 
      {
        out << "#if ! CODEGEN_COMPILE_INNERLOOPS" << std::endl;
        out << "template <class BaseFunctionSet>" << std::endl;
        out << "struct EvaluateRanges<BaseFunctionSet, Fem :: EmptyGeometry, " << dimRange << ", " << numRows << ", " << numCols << ">" << std::endl;
        out << "{" << std::endl;
        out << "  template< class QuadratureType,"<< std::endl;
        out << "            class RangeVectorType," << std::endl;
        out << "            class RangeFactorType," << std::endl;
        out << "            class LocalDofVectorType>" << std::endl;
        out << "  static void eval( const QuadratureType& quad," << std::endl;
        out << "                    const RangeVectorType& rangeStorage," << std::endl; 
        out << "                    const LocalDofVectorType& dofStorage," << std::endl;
        out << "                    RangeFactorType &rangeVector)" << std::endl;
        out << "  {" << std::endl;
        out << "    typedef typename ScalarRangeType :: field_type field_type;" << std::endl;
        out << "    typedef typename RangeVectorType :: value_type value_type;" << std::endl; 
        // make length sse conform 
        //const int sseDimRange = dimRange + (dimRange % 2); 
        const int sseDimRange = dimRange; // + (dimRange % 2); 
        out << "    const field_type dofs[ " << numCols * dimRange << " ] = { " << std::endl;
        for( size_t col = 0, colR = 0; col < numCols; ++ col )
        {
          out << "      ";
          for( int r = 0; r < dimRange; ++ r , ++colR )
          {
            out << "dofStorage[ " << colR << " ]";    

            if( colR < (numCols * dimRange) - 1 )
            {
              out << ","; 
              if( r == dimRange - 1) out << std::endl;
            }
            else out << " };" << std::endl;
          }
        }
        out << "    for( size_t row=0; row<"<<numRows<<" ; ++row )"<<std::endl;
        out << "    {" << std::endl;
        out << "      const value_type& rangeStorageRow = rangeStorage[ quad.cachingPoint( row ) ];" << std::endl;
        out << "      field_type result [ " << sseDimRange << " ] = { ";
        for( int r = 0; r < sseDimRange-1; ++r )  out << " 0 ,";
        out << " 0  };" << std::endl;
        for( size_t col = 0, colR = 0; col < numCols; ++col ) 
        {
          out << "      const field_type phi"<< col << " = rangeStorageRow[ " << col << " ][ 0 ];" << std::endl;
          for( int r = 0; r < dimRange; ++r , ++colR ) 
          {
            out << "      result[ " << r << " ] += dofs[ " << colR << " ] * phi" << col << ";" << std::endl;
            if( sseDimRange != dimRange && r == dimRange - 1 )
            {
              if( (colR + 1) < numCols * dimRange )
                out << "      result[ " << r+1 << " ] += dofs[ " << colR + 1 << " ] * phi" << col << ";" << std::endl;
              else 
                out << "      result[ " << r+1 << " ] += 0.5 * phi" << col << ";" << std::endl;
            }
          }
        }
        out << "      // store result in vector"  << std::endl;
        out << "      RangeType& realResult = rangeVector[ row ];"  << std::endl;
        for( int r = 0; r < dimRange; ++r) 
        {
          out << "      realResult[ " << r << " ] = result[ " << r << " ];" << std::endl;
        }
        out << "    }" << std::endl;
        out << "  }" << std::endl << std::endl;
        out << "};" << std::endl;
        out << "#endif" << std::endl;
      }

      static void axpyCodegen(std::ostream& out, 
              const int dim, const int dimRange, const size_t numRows, const size_t numCols ) 
      {
        out << "#if ! CODEGEN_COMPILE_INNERLOOPS" << std::endl;
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
        for( size_t row=0; row< numRows; ++row ) 
        {
          out << "    const RangeType& factor" << row  << " = rangeFactors[ " << row << " ];" << std::endl;
        }
        out << std::endl ;

        out << "    typedef typename RangeVectorType :: value_type value_type;" << std::endl; 
        for( size_t row = 0; row<numRows ; ++ row )
        {
          out << "    const value_type& rangeStorage"<<row<<" = rangeStorage[ quad.cachingPoint( " << row << " ) ];" << std::endl;
        }
        out << std::endl;
        out << "    typedef typename ScalarRangeType :: field_type field_type;" << std::endl;
        const int sseDimRange = dimRange + (dimRange % 2);
        out << "    for( size_t col = 0, dof = 0; col <"<<numCols<< " ; ++col )" << std::endl;
        {
          out << "    {" << std::endl;
          out << "      field_type result [ " << sseDimRange << " ] = { ";
          for( int r = 0; r < sseDimRange-1; ++r )  out << " 0 ,";
          out << " 0  };" << std::endl << std::endl;

          out << "      const field_type phi[ " << numRows << " ] = {" << std::endl;
          for( size_t row=0; row< numRows; ++row ) 
          {
            out << "        rangeStorage" << row << "[ col ][ 0 ]";
            if( row < numRows - 1) 
              out << " ," << std::endl;
            else out << "  };" << std::endl;
          }
          for( size_t row=0; row< numRows; ++row ) 
          {
            for( int r = 0; r < dimRange; ++r ) 
            {
              out << "      result[ " << r << " ]  +=  factor" << row << "[ " << r << " ] * phi[ " << row << " ];" << std::endl; 
              if( sseDimRange != dimRange && r == dimRange - 1 )
                out << "      result[ " << r+1 << " ]  +=  0.5 * phi[ " << row << " ];" << std::endl;
            }
          }
          for( int r = 0; r < dimRange; ++r ) //, ++colR ) 
            out << "      dofs[ dof++ ]  +=  result[ " << r << " ];" << std::endl;

          out << "    }" << std::endl;
        }
        out << "  }" << std::endl << std::endl;
        out << "};" << std::endl;
        out << "#endif" << std::endl;
      }

      static void evaluateJacobiansCodegen(std::ostream& out, 
              const int dim, const int dimRange, const size_t numRows, const size_t numCols ) 
      {
        out << "#if ! CODEGEN_COMPILE_INNERLOOPS" << std::endl;
        out << "template <class BaseFunctionSet>" << std::endl;
        out << "struct EvaluateJacobians<BaseFunctionSet, Fem :: EmptyGeometry, " << dimRange << ", " << numRows << ", " << numCols << ">" << std::endl;
        out << "{" << std::endl;
        out << "  template< class QuadratureType,"<< std::endl;
        out << "            class JacobianRangeVectorType," << std::endl;
        out << "            class LocalDofVectorType," << std::endl;
        out << "            class JacobianRangeFactorType>" << std::endl;
        out << "  static void eval( const QuadratureType&," << std::endl;
        out << "                    const Fem :: EmptyGeometry&," << std::endl; 
        out << "                    const JacobianRangeVectorType&," << std::endl; 
        out << "                    const LocalDofVectorType&," << std::endl;
        out << "                    JacobianRangeFactorType &)" << std::endl;
        out << "  {" << std::endl;
        out << "    std::cerr << \"ERROR: wrong code generated for VectorialBaseFunctionSet::axpyJacobians\" << std::endl;" << std::endl;
        out << "    abort();" << std::endl;
        out << "  }" << std::endl;
        out << "};" << std::endl << std::endl;
        out << "template <class BaseFunctionSet, class Geometry>" << std::endl;
        out << "struct EvaluateJacobians<BaseFunctionSet, Geometry, " << dimRange << ", " << numRows << ", " << numCols << ">" << std::endl;
        out << "{" << std::endl;
        out << "  template< class QuadratureType,"<< std::endl;
        out << "            class JacobianRangeVectorType," << std::endl;
        out << "            class LocalDofVectorType," << std::endl;
        out << "            class JacobianRangeFactorType>" << std::endl;
        out << "  static void eval( const QuadratureType& quad," << std::endl;
        out << "                    const Geometry& geometry," << std::endl; 
        out << "                    const JacobianRangeVectorType& jacStorage," << std::endl; 
        out << "                    const LocalDofVectorType& dofs," << std::endl;
        out << "                    JacobianRangeFactorType& jacFactors)" << std::endl;
        out << "  {" << std::endl;
        out << "    evalJac( quad, geometry, jacStorage, dofs, jacFactors, jacFactors[ 0 ] );" << std::endl;
        out << "  }" << std::endl;
        out << "private:" << std::endl;
        out << "  template< class QuadratureType,"<< std::endl;
        out << "            class JacobianRangeVectorType," << std::endl;
        out << "            class LocalDofVectorType," << std::endl;
        out << "            class JacobianRangeFactorType," << std::endl;
        out << "            class GlobalJacobianRangeType>" << std::endl;
        out << "  static void evalJac( const QuadratureType& quad," << std::endl;
        out << "                       const Geometry& geometry," << std::endl; 
        out << "                       const JacobianRangeVectorType& jacStorage," << std::endl; 
        out << "                       const LocalDofVectorType& dofs," << std::endl;
        out << "                       JacobianRangeFactorType& jacVector," << std::endl;
        out << "                       const GlobalJacobianRangeType& )" << std::endl;
        out << "  {" << std::endl;
        out << "    typedef typename JacobianRangeVectorType :: value_type  value_type;" << std::endl; 
        out << "    typedef typename JacobianRangeType :: field_type field_type;" << std::endl;
        out << "    for( size_t row = 0; row < " << numRows << " ; ++ row )" << std::endl;
        out << "    {" << std::endl;
        out << "      const value_type& jacStorageRow = jacStorage[ quad.cachingPoint( row ) ];" << std::endl;
        out << "      typedef typename Geometry::Jacobian GeometryJacobianType;" << std::endl;
        out << "      // use reference to GeometryJacobianType to make code compile with SPGrid Geometry" << std::endl;
        out << "      const GeometryJacobianType& gjit = geometry.jacobianInverseTransposed( quad.point( row ) );" << std::endl << std::endl;
        out << "      GlobalJacobianRangeType& result = jacVector[ row ];" << std::endl;
        out << "      result = 0;" << std::endl;
        out << "      typedef typename GlobalJacobianRangeType :: row_type JacobianRangeType;" << std::endl;
        out << "      JacobianRangeType gradPhi;" << std::endl;
        for( size_t col = 0, colR = 0; col < numCols; ++col )
        { 
          out << "      gjit.mv( jacStorageRow[ "<< col << " ][ 0 ], gradPhi );" << std::endl;
          for( int r = 0; r < dimRange; ++r, ++colR )
          {
            out << "      result[ " << r << " ].axpy( dofs[ "<< colR << " ], gradPhi );" << std::endl;
          }
        }
        out << "    }" << std::endl;
        out << "  }" << std::endl;
        out << "};" << std::endl;
        out << "#endif" << std::endl;
      }

      static void axpyJacobianCodegen(std::ostream& out, 
              const int dim, const int dimRange, const size_t numRows, const size_t numCols ) 
      {
        out << "#if ! CODEGEN_COMPILE_INNERLOOPS" << std::endl;
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
        const int sseDimRange = dimRange + (dimRange % 2); 
        //const int sseDimRange = dimRange;// + (dimRange % 2); 
        out << "    typedef typename JacobianRangeVectorType :: value_type  value_type;" << std::endl; 
        out << "    typedef typename JacobianRangeType :: field_type field_type;" << std::endl;
        const size_t dofs = sseDimRange * numCols ;
        out << "    field_type result [ " << dofs << " ] = {"; 
        for( size_t dof = 0 ; dof < dofs-1 ; ++ dof ) out << " 0,";
        out << " 0 };" << std::endl << std::endl;
        out << "    for( size_t row = 0; row < " << numRows << " ; ++ row )" << std::endl;
        out << "    {" << std::endl;
        out << "      const value_type& jacStorageRow = jacStorage[ quad.cachingPoint( row ) ];" << std::endl;
        out << "      typedef typename Geometry::Jacobian GeometryJacobianType;" << std::endl;
        out << "      // use reference to GeometryJacobianType to make code compile with SPGrid Geometry" << std::endl;
        out << "      const GeometryJacobianType& gjit = geometry.jacobianInverseTransposed( quad.point( row ) );" << std::endl << std::endl;
        out << "      JacobianRangeType jacFactorTmp;" << std::endl;
        out << "      for( int r = 0; r < " << dimRange << " ; ++r )" << std::endl;
        out << "      {"<<std::endl; 
        out << "        gjit.mtv( jacFactors[ row ][ r ], jacFactorTmp[ r ] );" << std::endl;
        out << "      }" << std::endl << std::endl;
        out << "      // calculate updates" << std::endl;
        out << "      // rearrange values to have linear memory walk through" << std::endl;
        out << "      const field_type jacFactorInv[ " << dim * sseDimRange << " ] = {" << std::endl;
        std::string delimiter("        ");
        for( int i =0; i < dim; ++i )
        {
          for( int r = 0; r < dimRange; ++ r ) 
          {
            out << delimiter;
            delimiter = ", ";
            out << "jacFactorTmp[ " << r  << " ][ " << i << " ]";
          }
          if( dimRange < sseDimRange ) 
          {
            out << ", 0.5";
          }
          delimiter = ",\n        ";
        }
        out << "      };" << std::endl;
        for( size_t col = 0, colR = 0; col < numCols; ++col ) 
        {
          out << "      {" << std::endl;
          out << "        const field_type phi[ " << dim * sseDimRange << " ] = {" << std::endl;
          delimiter = std::string("        ");
          for( int d =0; d < dim; ++d )
          {
            for( int r = 0; r < dimRange; ++ r ) 
            {
              out << delimiter;
              delimiter = ", ";
              out << "    jacStorageRow[ " << col << " ][ 0 ][ " << d << " ]";
              //FactorTmp[ " << r  << " ][ " << i << " ]";
            }
            if( dimRange < sseDimRange ) 
            {
              out << ", 0.5";
            }
            delimiter = ",\n        ";
          }
          out << std::endl;
          out << "        };" << std::endl;
          for( int d = 0, dr = 0 ; d < dim ; ++ d ) 
          {
            for( int r = 0; r < sseDimRange; ++r, ++dr )
            {
              out << "        result[ " << colR+r << " ]  +=  phi[ " << dr << " ] * jacFactorInv[ " << dr << " ];" << std::endl;
            }
          }
          colR += sseDimRange; 
          out << "      }" << std::endl;
        }
        /*
        for( size_t col = 0, colR = 0; col < numCols; ++col ) 
        {
          out << "      {" << std::endl;
          for( int d = 0 ; d < dim ; ++ d ) 
          {
            out << "        const field_type& phi"<< d << " = jacStorageRow[ " << col << " ][ 0 ][ " << d << " ];" << std::endl;
          }
          for( int r = 0; r < sseDimRange; ++r )
          {
            for( int d = 0 ; d < dim ; ++ d ) 
            {
              out << "        result[ " << colR+r << " ]  +=  phi" << d << " * jacFactorTmp[ " << r << " ][ " << d << " ];" << std::endl;
            }
          }
          colR += sseDimRange; 
          out << "      }" << std::endl;
        }
        */
        out << "    }" << std::endl << std::endl;

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
        out << "#endif" << std::endl;
      }
    };

    // if this pre processor variable is defined then 
    // we assume that CODEGENERATOR_REPLACEMENT is CodeGenerator of choice 
#ifndef FEM_CODEGENERATOR_REPLACEMENT
    typedef DefaultCodeGenerator CodeGeneratorType;
#else 
    typedef FEM_CODEGENERATOR_REPLACEMENT CodeGeneratorType;  
#endif

  } // namespace Fem 

} // namespace Dune

#endif // #ifndef DUNE_FEM_BASISFUNCTIONSETS_CODEGEN_HH
