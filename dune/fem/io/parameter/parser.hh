#ifndef DUNE_FEM_IO_PARAMETER_PARSER_HH
#define DUNE_FEM_IO_PARAMETER_PARSER_HH

#include <sstream>
#include <string>
#include <type_traits>

#include <dune/common/fmatrix.hh>

namespace Dune
{

  namespace Fem
  {

    // ParameterParser
    // ---------------

    template< class T >
    struct ParameterParser
    {
      static bool parse ( const std::string &s, T &value )
      {
        std::istringstream in( s );
        in >> value;
        if( std::is_same< T, std::string >::value && s.empty() )
          return true;
        if( in.fail() )
          return false;
        char eof;
        in >> eof;
        return in.eof();
      }

      static std::string toString ( const T &value )
      {
        std::ostringstream out;
        out << value;
        return out.str();
      }
    };

    template<>
    struct ParameterParser< bool >
    {
      static bool parse ( const std::string &s, bool &value )
      {
        std::string w;
        if( ParameterParser< std::string >::parse( s, w ) )
        {
          std::transform(w.begin(), w.end(), w.begin(), ::tolower);
          if( (w == std::string( "false" )) || (w == std::string( "no" )) || (w == std::string( "0" )) )
          {
            value = false;
            return true;
          }

          if( (w == std::string( "true" )) || (w == std::string( "yes" )) || (w == std::string( "1" )) )
          {
            value = true;
            return true;
          }
        }
        return false;
      }

      static std::string toString ( const bool &value )
      {
        return std::string( value ? "true" : "false" );
      }
    };

    template< class F, int m, int n >
    struct ParameterParser< FieldMatrix< F, m, n > >
    {
      static bool parse ( const std::string &s, FieldMatrix< F, m, n > &value )
      {
        std::istringstream in( s );
        char c;
        for( int i = 0; i < m; ++i )
        {
          if( i > 0 )
          {
            in >> c;
            if( c != ',' )
              return false;
          }

          for( int j = 0; j < n; ++j )
            in >> value[ i ][ j ];
        }
        in >> c; // read eof
        return in.eof();
      }

      static std::string toString ( const FieldMatrix< F, m, n > &value )
      {
        std::ostringstream out;
        for( int i = 0; i < m; ++i )
        {
          out << (i > 0 ? "," : "");
          for( int j = 0; j< n; ++j )
            out << " " << value[ i ][ j ];
        }
        return out.str();
      }
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_IO_PARAMETER_PARSER_HH
