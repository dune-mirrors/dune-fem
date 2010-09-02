#ifndef DUNE_FEM_LATEXTABLEWRITER_HH
#define DUNE_FEM_LATEXTABLEWRITER_HH

#include <cmath>
#include <fstream>
#include <limits>
#include <iomanip>
#include <sstream>
#include <vector>

#include <dune/common/exceptions.hh>

#include <dune/fem/misc/femtuples.hh>

namespace Dune
{

  namespace Fem
  {

    // AbstractColumnWriter
    // --------------------

    template< class DataTuple >
    struct AbstractColumnWriter
    {
      enum Alignment { AlignLeft, AlignCenter, AlignRight };

      virtual ~AbstractColumnWriter () {}

      virtual Alignment alignment () const { return AlignCenter; }
      virtual std::string entry ( const DataTuple & ) const = 0;
      virtual std::string header () const = 0;
    };



    // NumberColumnWriter
    // ------------------

    template< class DataTuple, int N >
    class NumberColumnWriter
    : public AbstractColumnWriter< DataTuple >
    {
      typedef AbstractColumnWriter< DataTuple > BaseType;

    public:
      explicit NumberColumnWriter ( const std::string &header, const int decimals = 6 )
      : header_( header ),
        decimals_( decimals )
      {}

      typename BaseType::Alignment alignment () const
      {
        return BaseType::AlignRight;
      }

      std::string entry ( const DataTuple &data ) const
      {
        return toString( ElementAccess< N >::get( data ) );
      }

      std::string header () const { return header_; }

    protected:
      template< class Number >
      std::string toString ( const Number &number ) const
      {
        std::ostringstream s;
        s << std::fixed << std::setprecision( decimals_ );
        s << "$" << number << "$";
        return s.str();
      }

    private:
      std::string header_;
      int decimals_;
    };



    // EOCColumnWriter
    // ---------------

    template< class DataTuple, int Nh, int Ne >
    class EOCColumnWriter
    : public NumberColumnWriter< DataTuple, Ne >
    {
      typedef NumberColumnWriter< DataTuple, Ne > BaseType;

    public:
      explicit EOCColumnWriter ( const std::string &header, const int decimals = 2 )
      : BaseType( header, decimals ),
        hOld_( std::numeric_limits< double >::infinity() )
      {}

      std::string entry ( const DataTuple &data ) const
      {
        const double h = ElementAccess< Nh >::get( data );
        const double e = ElementAccess< Ne >::get( data );

        std::string entry = "---";
        if( hOld_ < std::numeric_limits< double >::infinity() )
        {
          const double eoc = std::log( e / eOld_ ) / std::log( h / hOld_ );
          entry = BaseType::toString( eoc );
        }
        hOld_ = h;
        eOld_ = e;
        return entry;
      }

    private:
      mutable double hOld_;
      mutable double eOld_;
    };



    // LatexTableWriter
    // ----------------

    template< class DataTuple >
    struct LatexTableWriter
    {
      typedef AbstractColumnWriter< DataTuple > ColumnWriterType;
      typedef std::vector< const ColumnWriterType * > ColumnWriterVectorType;

      LatexTableWriter ( const std::string &filename, const ColumnWriterVectorType &columnWriter );
      ~LatexTableWriter ();

      void writeRow ( const DataTuple &data );
      void writeSeparator ();

    private:
      typedef typename ColumnWriterVectorType::const_iterator ColumnWriterIteratorType;

      void cleanUp ();

      std::string preamble () const;
      std::string header () const;

      ColumnWriterVectorType columnWriters_;
      std::ofstream out_;
    };



    // Implementation of LatexTableWriter
    // ----------------------------------

    template< class DataTuple >
    inline LatexTableWriter< DataTuple >
      ::LatexTableWriter ( const std::string &filename, const ColumnWriterVectorType &columnWriters )
    : columnWriters_( columnWriters ),
      out_( filename )
    {
      if( !out_ )
      {
        cleanUp();
        DUNE_THROW( IOError, "Unable to open '" << filename << "'." );
      }
      out_ << "\\begin{tabular}{" << preamble() << "}" << std::endl;
      writeSeparator();
      out_ << header() << " \\\\"<< std::endl;
      writeSeparator();
    }


    template< class DataTuple >
    inline LatexTableWriter< DataTuple >::~LatexTableWriter ()
    {
      writeSeparator();
      out_ << "\\end{tabular}" << std::endl;
      cleanUp();
    }


    template< class DataTuple >
    inline void LatexTableWriter< DataTuple >::writeRow ( const DataTuple &data )
    {
      std::string separator = "";
      const ColumnWriterIteratorType end = columnWriters_.end();
      for( ColumnWriterIteratorType it = columnWriters_.begin(); it != end; ++it )
      {
        out_ << separator;
        separator = " & ";
        if( *it )
          out_ << (*it)->entry( data );
      }
      out_ << " \\\\" << std::endl;
    }


    template< class DataTuple >
    inline void LatexTableWriter< DataTuple >::writeSeparator ()
    {
      int b = 0, e = 0;

      const ColumnWriterIteratorType end = columnWriters_.end();
      for( ColumnWriterIteratorType it = columnWriters_.begin(); it != end; ++it, ++e )
      {
        if( !(*it) )
        {
          if( b < e )
            out_ << "\\cline{" << (b+1) << "-" << e << "}";
          b = e+1;
        }
      }
      if( b < e )
        out_ << "\\cline{" << (b+1) << "-" << e << "}";
      out_ << std::endl;
    }


    template< class DataTuple >
    inline void LatexTableWriter< DataTuple >::cleanUp ()
    {
      const ColumnWriterIteratorType end = columnWriters_.end();
      for( ColumnWriterIteratorType it = columnWriters_.begin(); it != end; ++it )
      {
        if( *it )
          delete( *it );
      }
    }


    template< class DataTuple >
    inline std::string LatexTableWriter< DataTuple >::header () const
    {
      std::string header, separator;
      const ColumnWriterIteratorType end = columnWriters_.end();
      for( ColumnWriterIteratorType it = columnWriters_.begin(); it != end; ++it )
      {
        header += separator;
        separator = " & ";
        if( *it )
          header += "\\multicolumn{1}{|c|}{" + (*it)->header() + "}";
      }
      return header;
    }


    template< class DataTuple >
    inline std::string LatexTableWriter< DataTuple >::preamble () const
    {
      const char alignment[] = { 'l', 'c', 'r' };

      std::string preamble( "|" );
      const ColumnWriterIteratorType end = columnWriters_.end();
      for( ColumnWriterIteratorType it = columnWriters_.begin(); it != end; ++it )
      {
        if( *it )
          preamble += alignment[ (*it)->alignment() ];
        else
          preamble += "@{}p{0.2em}@{}";
        preamble += "|";
      }
      return preamble;
    }

  }

}

#endif // #ifndef DUNE_FEM_LATEXTABLEWRITER_HH
