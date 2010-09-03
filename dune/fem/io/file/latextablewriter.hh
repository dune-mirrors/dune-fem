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

    /** \class AbstractColumnWriter
     *  \brief Class representing column writer in general
     *  
     *  This class represents column writer.
     *
     * \tparam DataTuple Type of the data tuple
     */
    template< class DataTuple >
    struct AbstractColumnWriter
    {
      //! The alignment for the data in this column
      enum Alignment { AlignLeft, AlignCenter, AlignRight };

      //! Destructor
      virtual ~AbstractColumnWriter () {}

      //! \return alignment
      virtual Alignment alignment () const { return AlignCenter; }
      //! \return data entry
      virtual std::string entry ( const DataTuple & ) const = 0;
      //! \return column titles as latex row
      virtual std::string header () const = 0;
    };



    // NumberColumnWriter
    // ------------------

    /** \class NumberColumnWriter
     *  \brief gets the @a N th element of a provided tuple assuming its a number
     *
     *  This class extracts the @a N th element of the tuple assuming
     *  the @a N th element is a number
     *
     *  \tparam DataTuple The type of the data tuple
     *  \tparam N Index of the entry in the data tuple tp be extracted
     */
    template< class DataTuple, int N >
    class NumberColumnWriter
    : public AbstractColumnWriter< DataTuple >
    {
      typedef AbstractColumnWriter< DataTuple > BaseType;

    public:
      /** Constructor
       *  \param[in] header Column titles in latex row format
       *  \param[in] decimals The precision of double for output to latex table
       */
      explicit NumberColumnWriter ( const std::string &header, const int decimals = 6 )
      : header_( header ),
        decimals_( decimals )
      {}

      //! set the aligment of the entries for this column in the latex table
      typename BaseType::Alignment alignment () const
      {
        return BaseType::AlignRight;
      }

      /** \brief returns @a N the element from @a data tuple
       *  \return String representing @a N th element of @a data
       *
       *  \note @a N th element of @a data is assumed to be number
       */
      std::string entry ( const DataTuple &data ) const
      {
        return toString( ElementAccess< N >::get( data ) );
      }

      //! return Column titles in latex row format
      std::string header () const { return header_; }

    protected:
      //! converts @a number to std::string
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

    /** \class EOCColumnWriter
     *  \brief calculates EOC from the @a Nh th and @a Ne th element of tuple
     *
     *  This class extracts the @a N th element of the tuple assuming
     *  the @a N th element is a number
     *
     *  \tparam DataTuple The type of the data tuple
     *  \tparam N Index of the entry in the data tuple tp be extracted
     */
    template< class DataTuple, int Nh, int Ne >
    class EOCColumnWriter
    : public NumberColumnWriter< DataTuple, Ne >
    {
      typedef NumberColumnWriter< DataTuple, Ne > BaseType;

    public:
      /** Constructor
       *  \param[in] header Column titles in latex row format
       *  \param[in] decimal Decimal placese for this column. Default 2.
       */
      explicit EOCColumnWriter ( const std::string &header, const int decimals = 2 )
      : BaseType( header, decimals ),
        hOld_( std::numeric_limits< double >::infinity() )
      {}

      /** \brief returns latex table entry 
       *
       *  \param[in] data Data tuple
       *  
       *  \return String representing EOC calculated from @a Nh th
       *          and @a Ne th elements of the @a data tuple
       *
       *  \note @a Ne th and @a Nh th elements of data are assumed to
       *        be numbers
       */
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

    /** \class LatexTableWriter
     *  \brief writes latex tables based on user-defined row structure
     *
     *  The class LatexTableWriter writes a latex table where each
     *  row corresponds to user-provided row structure of the @a DataTuple type.
     */
    template< class DataTuple >
    struct LatexTableWriter
    {
      //! Abstract column type
      typedef AbstractColumnWriter< DataTuple > ColumnWriterType;
      //! Abstract column vector type
      typedef std::vector< const ColumnWriterType * > ColumnWriterVectorType;

      /** Constructor
       *  \param[in] filename The name of the latex file
       *  \param[in] columnWriter Abstract column writer
       */
      LatexTableWriter ( const std::string &filename, const ColumnWriterVectorType &columnWriter );
      /** Destructor
       *  \brief writes "\end{tabular}" to the latex file and removes column vector
       */
      ~LatexTableWriter ();

      //! Write row to the table
      void writeRow ( const DataTuple &data );
      //! Adds extra space between two columns in the latex table
      void writeSeparator ();

    private:
      //! Type of iterator over column vector
      typedef typename ColumnWriterVectorType::const_iterator ColumnWriterIteratorType;

      //! Remove all columns in column vector
      void cleanUp ();

      //! The latex format of the table
      std::string preamble () const;
      //! The column titles
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
