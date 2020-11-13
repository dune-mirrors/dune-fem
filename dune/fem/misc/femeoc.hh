#ifndef DUNE_FEM_FEMEOC_HH
#define DUNE_FEM_FEMEOC_HH

#include <cassert>
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <tuple>

#include <dune/common/fvector.hh>

#include <dune/fem/io/io.hh>
#include <dune/fem/io/file/latextablewriter.hh>
#include <dune/fem/io/parameter.hh>
#include <dune/fem/storage/singleton.hh>

namespace Dune
{

  namespace Fem
  {

    /**
        @ingroup HelperClasses
        \brief Write a self contained tex table
        for eoc runs with timing information.

        Constructor takes base file name for the tex file and
        generates two files:
        filename_timestamp_main.tex and filename_timestamp_body.tex or
        without the time stamp (by default) in the name depending on
        the parameter fem.io.eocFileTimeStamp.
        A time stamp is added to the base file name to prevent the
        overwriting of a valuable eoc data from the previous simulation.
        The file filename_timestamp_body.tex holds the actual body
        of the eoc table which is included in filename_timestamp_main.tex
        but can also be used to combine e.g. runs with different
        parameters or for plotting using gnuplot.

        The class is singleton and thus new errors for eoc
        computations can be added in any part of the program.
        To add a new entry for eoc computations use one of the
        addEntry methods. These return a unique unsigned int
        which can be used to add error values to the table
        with the setErrors methods.
        The method write is used to write a single line
        to the eoc table.
     */
    class FemEoc
    {
      typedef std::pair< std::string, double > DoublePairType;
      typedef std::pair< std::string, int >    IntPairType;

      // level, h, size, time, counter, errors,
      // [avgTimeStep, minTimeStep, maxTimeStep],
      // [newton_iterations, ils_iterations, max_newton_iterations, max_ils_iterations]
      typedef std::tuple< int, double, double, double, int, std::vector< double >,
                          std::vector< double >, std::vector< int > >
        DataTuple;

      typedef Fem::LatexTableWriter< DataTuple > TableWriter;

      class ErrorColumnWriter;
      class EOCColumnWriter;

      TableWriter *tableWriter_;
      std::string filename_;
      int level_;
      std::vector< double > error_;
      std::vector< const EOCColumnWriter * > eocColumns_;
      std::vector< std::string > description_;
      std::vector< int > pos_;

    public:
      FemEoc ()
      : tableWriter_( 0 ),
        level_( 0 )
      {}

      ~FemEoc()
      {
        clearFile();
      }

    private:
      void clearFile()
      {
        if( tableWriter_ )
          delete tableWriter_;
        tableWriter_ = 0;
      }

      void init ( const std::string& path, const std::string &name, const std::string &descript )
      {
        if( MPIManager::rank() == 0 )
        {
          if( createDirectory( path ) )
            init( path + "/" + name, descript );
          else
            std::cerr << "Error: Unable to create path '" << path << "'" << std::endl;
        }
      }

      void init ( const std::string &filename, const std::string &descript )
      {
        if( MPIManager::rank() != 0 )
          return;

        std::string name = filename;

        // add time stamp to file name, if requested (prevents results from being overwritten)
        if( Parameter::getValue< bool >( "fem.io.eocFileTimeStamp", false ) )
        {
          time_t seconds = time(0);
          struct tm *ptm = localtime( &seconds );
          char timeString[20];
          strftime( timeString, 20, "_%d%m%Y_%H%M%S", ptm );
          name += std::string( timeString );
        }

        if( !tableWriter_ )
        {
          filename_ = name + "_body.tex";
          std::ofstream main( (name + "_main.tex").c_str() );

          if (!main) {
            std::cerr << "Could not open file : "
                      << (name+"_main.tex").c_str()
                      << " ... ABORTING" << std::endl;
            abort();
          }

          main << "\\documentclass[12pt,english]{article}\n"
               << "\\usepackage[T1]{fontenc}\n"
               << "\\usepackage[latin1]{inputenc}\n"
               << "\\usepackage{setspace}\n"
               << "\\onehalfspacing\n"
               << "\\makeatletter\n"
               << "\\providecommand{\\boldsymbol}[1]{\\mbox{\\boldmath $#1$}}\n"
               << "\\providecommand{\\tabularnewline}{\\\\}\n"
               << "\\usepackage{babel}\n"
               << "\\makeatother\n"
               << "\\begin{document}\n"
               << "\\begin{center}\\large\n"
               << "\n\\end{center}\n\n"
               << descript
               << "\\input{"
               << filename_
               << "}\n\n"
               << "\\end{document}\n" << std::endl;
          main.close();
        }
        else
        {
          std::cerr << "Could not open file : "
                    << " already opened!"
                    << " ... ABORTING" << std::endl;
          abort();
        }
      }

      template <class StrVectorType>
      size_t addentry(const StrVectorType& descript,size_t size)
      {
        if( tableWriter_ )
        {
          std::cerr << "Trying to add a new entry to FemEoc although "
                    << "entries have already been writen to disk!"
                    << " ... ABORTING" << std::endl;
          abort();
        }
        pos_.push_back(error_.size());
        for (size_t i=0;i<size;++i) {
          error_.push_back(0);
          description_.push_back(descript[i]);
        }
        return pos_.size()-1;
      }

      size_t addentry ( const std::string &descript )
      {
        if( tableWriter_ )
        {
          std::cerr << "Trying to add a new entry to FemEoc although "
                    << "entries have already been writen to disk!"
                    << " ... ABORTING" << std::endl;
          abort();
        }
        pos_.push_back(error_.size());
        error_.push_back(0);
        description_.push_back(descript);
        return pos_.size()-1;
      }

      template <class VectorType>
      void seterrors(size_t id,const VectorType& err,size_t size)
      {
        assert(id<pos_.size());
        int pos = pos_[ id ];
        assert(pos+size <= error_.size());

        for (size_t i=0; i<size; ++i)
          error_[pos+i] = err[i];
      }

      template <int SIZE>
      void seterrors(size_t id,const FieldVector<double,SIZE>& err)
      {
        seterrors(id,err,SIZE);
      }

      void seterrors(size_t id,const double& err) {
        int pos = pos_[id];
        error_[pos] = err;
      }

      void writeerr ( double h, double size, double time, int counter );
      void writeerr(double h,double size,double time,int counter,
                    const std::vector< DoublePairType>& doubleValues,
                    const std::vector< IntPairType>& intValues);

      void writeerr(double h,double size,double time,int counter,
                    double avgTimeStep,double minTimeStep,double maxTimeStep ,
                    const int newton_iterations, const int ils_iterations,
                    const int max_newton_iterations, const int max_ils_iterations);

      // do the same calculations as in write, but don't overwrite status
      void printerr(const double h,
                    const double size,
                    const double time,
                    const int counter,
                    std::ostream& out);
      void printerr(const double h,
                    const double size,
                    const double time,
                    const int counter,
                    const std::vector< DoublePairType>& doubleValues,
                    const std::vector< IntPairType>& intValues,
                    std::ostream& out);

    public:
      friend class Dune::Fem::Singleton< FemEoc >;
      static FemEoc& instance()
      {
        return Singleton< FemEoc > :: instance();
      }

      //! close file and allow FemEoc to used for a second run
      static void clear() {
        instance().clearFile();
      }
      //! open file path/name and write a description string into tex file
      static void initialize(const std::string& path, const std::string& name, const std::string& descript) {
        instance().init(path,name,descript);
      }
      //! open file name and write description string into tex file
      static void initialize(const std::string& name, const std::string& descript) {
        instance().init(name,descript);
      }
      /** \brief add a vector of new eoc values
       *
       *  \tparam  StrVectorType a vector type with operator[]
       *           returning a string (a C style array can be used)
       *           the size of the vector is given as parameter
       *  \return  a unique index used to add the error values
       */
      template <class StrVectorType>
      static size_t addEntry(const StrVectorType& descript,size_t size) {
        return instance().addentry(descript,size);
      }
      /** \brief add a vector of new eoc values
       *
       *  \tparam  StrVectorType a vector type with size() and operator[]
       *           returning a string
       *  \return  a unique index used to add the error values
       */
      template <class StrVectorType>
      static size_t addEntry(const StrVectorType& descript) {
        return instance().addentry(descript,descript.size());
      }
      /** \brief add a single new eoc output
       *
       *  \return  a unique index used to add the error values
       */
      static size_t addEntry(const std::string& descript) {
        return instance().addentry(descript);
      }
      /** \brief add a single new eoc output
       *
       *  \return  a unique index used to add the error values
       */
      static size_t addEntry(const char* descript) {
        return addEntry(std::string(descript));
      }
      /** \brief add a vector of error values for the given id (returned by
       *         addEntry)
       *  \tparam  VectorType a vector type with an operator[]
       *           returning a double (C style array can be used)
       */
      template <class VectorType>
      static void setErrors(size_t id,const VectorType& err,int size)
      {
        instance().seterrors(id,err,size);
      }
      /** \brief add a vector of error values for the given id (returned by
       *         addEntry)
       *  \tparam  VectorType a vector type with a size() and an operator[]
       *           returning a double
       */
      template <class VectorType>
      static void setErrors(size_t id,const VectorType& err) {
        instance().seterrors(id,err,err.size());
      }
      /** \brief add a vector in a FieldVector of error values for the given id (returned by
       *         addEntry)
       */
      template <int SIZE>
      static void setErrors(size_t id,const FieldVector<double,SIZE>& err) {
        instance().seterrors(id,err);
      }
      /** \brief add a single error value for the given id (returned by
       *         addEntry)
       */
      static void setErrors(size_t id,const double& err) {
        instance().seterrors(id,err);
      }
      /** \brief commit a line to the eoc file
       *
       *  \param h grid width (e.g. given by GridWith utitlity class)
       *  \param size number of elements in the grid or number of dofs...
       *  \param time computational time
       *  \param counter number of timesteps or iterations for a solver...
       */
      static void write(double h,double size,double time,int counter)
      {
        instance().writeerr(h,size,time,counter);
      }

      /** \brief commit a line to the eoc file
       *
       *  \param h grid width (e.g. given by GridWith utitlity class)
       *  \param size number of elements in the grid or number of dofs...
       *  \param time computational time
       *  \param counter number of timesteps or iterations for a solver...
       *  \param out std::ostream to print data to (e.g. std::cout)
       */
      static void write(const double h,
                        const double size,
                        const double time,
                        const int counter,
                        std::ostream& out)
      {
        // print last line to out
        instance().printerr( h, size, time, counter, out );

        // now write to file
        instance().writeerr(h,size,time,counter);
      }

      /** \brief commit a line to the eoc file
       *
       *  \param  h                      grid width (e.g. given by GridWith utility class)
       *  \param  size                   number of grid elements
       *  \param  time                   computational time
       *  \param  counter                number of time steps
       *  \param  avgTimeStep            average time step for the ODE solver
       *  \param  minTimeStep            minimal time step for the ODE solver
       *  \param  maxTimeStep            maximal time step for the ODE solver
       *  \param  newton_iterations      number of newton iterations
       *  \param  ils_iterations         number of iteration of the iterative linear solver
       *  \param  max_newton_iterations  maximal number of newton iterations
       *  \param  max_ils_iterations     maximal number of iteration of the iterative linear solver
       */
      static void write(const double h,
                        const double size,
                        const double time,
                        const int counter,
                        const double avgTimeStep,
                        const double minTimeStep,
                        const double maxTimeStep,
                        const int newton_iterations,
                        const int ils_iterations,
                        const int max_newton_iterations,
                        const int max_ils_iterations)
      {
        std::vector< DoublePairType > doubleValues;
        doubleValues.push_back( DoublePairType( "avg dt", avgTimeStep ) );
        doubleValues.push_back( DoublePairType( "min dt", minTimeStep ) );
        doubleValues.push_back( DoublePairType( "max dt", maxTimeStep ) );

        std::vector< IntPairType > intValues;
        intValues.push_back( IntPairType( "Newton", newton_iterations ) );
        intValues.push_back( IntPairType( "ILS", ils_iterations ) );
        intValues.push_back( IntPairType( "max{Newton/linS}", max_newton_iterations ) );
        intValues.push_back( IntPairType( "max{ILS/linS}", max_ils_iterations ) );

        // now write to file
        instance().writeerr(h,size,time,counter, doubleValues, intValues );
      }

      /** \brief commit a line to the eoc file
       *
       *  \param  h                      grid width (e.g. given by GridWith utility class)
       *  \param  size                   number of grid elements
       *  \param  time                   computational time
       *  \param  counter                number of time steps
       *  \param  doubleValues           list of extra double values to be written
       *  \param  intValues              list of extra int values to be written
       */
      static void write(const double h,
                        const double size,
                        const double time,
                        const int counter,
                        const std::vector< DoublePairType >& doubleValues,
                        const std::vector< IntPairType >& intValues )
      {
        // now write to file
        instance().writeerr(h,size,time,counter, doubleValues, intValues );
      }

      /** \brief commit a line to the eoc file
       *
       *  \param  h                      grid width (e.g. given by GridWith utility class)
       *  \param  size                   number of grid elements
       *  \param  time                   computational time
       *  \param  counter                number of time steps
       *  \param  avgTimeStep            average time step for the ODE solver
       *  \param  minTimeStep            minimal time step for the ODE solver
       *  \param  maxTimeStep            maximal time step for the ODE solver
       *  \param  newton_iterations      number of newton iterations
       *  \param  ils_iterations         number of iteration of the iterative linear solver
       *  \param  max_newton_iterations  maximal number of newton iterations
       *  \param  max_ils_iterations     maximal number of iteration of the iterative linear solver
       */
      static void write(const double h,
                        const double size,
                        const double time,
                        const int counter,
                        const double avgTimeStep,
                        const double minTimeStep,
                        const double maxTimeStep,
                        const int newton_iterations,
                        const int ils_iterations,
                        const int max_newton_iterations,
                        const int max_ils_iterations,
                        std::ostream& out)
      {
        std::vector< DoublePairType > doubleValues;
        doubleValues.push_back( DoublePairType( "avg dt", avgTimeStep ) );
        doubleValues.push_back( DoublePairType( "min dt", minTimeStep ) );
        doubleValues.push_back( DoublePairType( "max dt", maxTimeStep ) );

        std::vector< IntPairType > intValues;
        intValues.push_back( IntPairType( "Newton", newton_iterations ) );
        intValues.push_back( IntPairType( "ILS", ils_iterations ) );
        intValues.push_back( IntPairType( "max{Newton/linS}", max_newton_iterations ) );
        intValues.push_back( IntPairType( "max{ILS/linS}", max_ils_iterations ) );

        // print last line to out
        instance().printerr( h, size, time, counter, doubleValues, intValues, out );

        // now write to file
        instance().writeerr(h,size,time,counter, doubleValues, intValues );
      }

      /** \brief commit a line to the eoc file
       *
       *  \param  h                      grid width (e.g. given by GridWith utility class)
       *  \param  size                   number of grid elements
       *  \param  time                   computational time
       *  \param  counter                number of time steps
       *  \param  avgTimeStep            average time step for the ODE solver
       *  \param  minTimeStep            minimal time step for the ODE solver
       *  \param  maxTimeStep            maximal time step for the ODE solver
       *  \param  newton_iterations      number of newton iterations
       *  \param  ils_iterations         number of iteration of the iterative linear solver
       *  \param  max_newton_iterations  maximal number of newton iterations
       *  \param  max_ils_iterations     maximal number of iteration of the iterative linear solver
       */
      static void write(const double h,
                        const double size,
                        const double time,
                        const int counter,
                        const std::vector< DoublePairType >& doubleValues,
                        const std::vector< IntPairType >& intValues,
                        std::ostream& out)
      {
        // print last line to out
        instance().printerr( h, size, time, counter, doubleValues, intValues, out );

        // now write to file
        instance().writeerr(h,size,time,counter, doubleValues, intValues );
      }

    }; // end class FemEoc


    class FemEoc::ErrorColumnWriter
    : public Fem::AbstractColumnWriter< FemEoc::DataTuple >
    {
      typedef Fem::AbstractColumnWriter< FemEoc::DataTuple > BaseType;

    public:
      ErrorColumnWriter ( const std::string &header, const int index )
      : header_( header ),
        index_( index )
      {}

      std::string entry ( const FemEoc::DataTuple &data ) const
      {
        return toString( error( data ) );
      }

      std::string header () const { return header_; }

    protected:
      double error ( const FemEoc::DataTuple &data ) const
      {
        return std::get< 5 >( data )[ index_ ];
      }

      std::string toString ( const double &error ) const
      {
        std::ostringstream s;
        s << "$" << error << "$";
        // s << " " << error << " ";
        return s.str();
      }

    private:
      std::string header_;
      int index_;
    };


    class FemEoc::EOCColumnWriter
    : public FemEoc::ErrorColumnWriter
    {
      typedef FemEoc::ErrorColumnWriter BaseType;

    public:
      explicit EOCColumnWriter ( const int index )
      : BaseType( "EOC", index ),
        hOld_( std::numeric_limits< double >::infinity() )
      {}

      std::string entry ( const DataTuple &data ) const
      {
        const double h = std::get< 1 >( data );
        const double e = BaseType::error( data );

        std::string entry = "---";
        if( hOld_ < std::numeric_limits< double >::infinity() )
          entry = BaseType::toString( eoc( h, e ) );
        hOld_ = h;
        eOld_ = e;
        return entry;
      }

      double eoc ( const double h, const double e ) const
      {
        return std::log( e / eOld_ ) / std::log( h / hOld_ );
      }

    private:
      mutable double hOld_;
      mutable double eOld_;
    };



    inline void FemEoc
      ::writeerr ( double h, double size, double time, int counter )
    {
      std::vector< DoublePairType > doubleValues;
      std::vector< IntPairType > intValues;
      writeerr( h, size, time, counter, doubleValues, intValues);
    }


    inline void FemEoc
      ::writeerr(double h,double size,double time,int counter,
                 const std::vector< DoublePairType >& doubleValues,
                 const std::vector< IntPairType >& intValues )
    {
      if( MPIManager::rank() != 0 )
        return;

      if( !tableWriter_ )
      {
        TableWriter::ColumnWriterVectorType columns;
        columns.push_back( new Fem::NumberColumnWriter< DataTuple, Fem::TupleDataSource< 0 > >( "level" ) );
        columns.push_back( new Fem::NumberColumnWriter< DataTuple, Fem::TupleDataSource< 1 > >( "h" ) );
        columns.push_back( new Fem::NumberColumnWriter< DataTuple, Fem::TupleDataSource< 2 > >( "size" ) );
        columns.push_back( new Fem::NumberColumnWriter< DataTuple, Fem::TupleDataSource< 3 > >( "CPU-time" ) );
        columns.push_back( new Fem::NumberColumnWriter< DataTuple, Fem::TupleDataSource< 4 > >( "counter" ) );
        columns.push_back( (const TableWriter::ColumnWriterType *)0 );

        typedef Fem::ArrayDataSource< Fem::TupleDataSource< 6 > > DoubleValueSource;
        for( unsigned int i = 0; i < doubleValues.size(); ++i )
        {
          columns.push_back( new Fem::NumberColumnWriter< DataTuple, DoubleValueSource >( doubleValues[ i ].first, DoubleValueSource( i ) ) );
        }

        typedef Fem::ArrayDataSource< Fem::TupleDataSource< 7 > > IntValueSource;
        for( unsigned int i = 0; i < intValues.size(); ++i )
        {
          columns.push_back( new Fem::NumberColumnWriter< DataTuple, IntValueSource >( intValues[ i ].first, IntValueSource( i ) ) );
        }

        eocColumns_.resize( error_.size(), (const EOCColumnWriter *)0 );
        for( unsigned int i = 0; i < error_.size(); ++i )
        {
          columns.push_back( (const TableWriter::ColumnWriterType *)0 );
          columns.push_back( new ErrorColumnWriter( description_[ i ], i ) );
          eocColumns_[ i ] = new EOCColumnWriter( i );
          columns.push_back( eocColumns_[ i ] );
        }

        tableWriter_ = new TableWriter( filename_, columns );
      }

      std::vector< double > doubleVals( doubleValues.size() );
      for( unsigned int i=0; i<doubleValues.size(); ++i )
        doubleVals[ i ] =  doubleValues[ i ].second;

      std::vector< int > intVals( intValues.size() );
      for( unsigned int i=0; i<intValues.size(); ++i )
        intVals[ i ] =  intValues[ i ].second;

      DataTuple data( level_, h, size, time, counter, error_, doubleVals, intVals );
      tableWriter_->writeRow( data );
      ++level_;
    }


    inline void FemEoc
      ::printerr(const double h,
                  const double size,
                  const double time,
                  const int counter,
                  std::ostream& out)
    {
      std::vector< DoublePairType > doubleValues;
      std::vector< IntPairType >    intValues;
      printerr( h, size, time, counter, doubleValues, intValues, out );
    }

    inline void FemEoc
      ::printerr(const double h,
                  const double size,
                  const double time,
                  const int counter,
                  const std::vector< DoublePairType >& doubleValues,
                  const std::vector< IntPairType >& intValues,
                  std::ostream& out)
    {
      if (!Parameter::verbose()) return;

      out << "level:   " << level_  << std::endl;
      out << "h        " << h << std::endl;
      out << "size:    " << size << std::endl;
      out << "time:    " << time << " sec. " << std::endl;
      out << "counter: " << counter << std::endl;
      for( unsigned int i=0; i<doubleValues.size(); ++i )
      {
        out << doubleValues[ i ].first << ": " << doubleValues[ i ].second << std::endl;
      }
      for( unsigned int i=0; i<intValues.size(); ++i )
      {
        out << intValues[ i ].first << ": " << intValues[ i ].second << std::endl;
      }

      for (unsigned int i=0;i<error_.size();++i)
      {
        out << description_[i] << ":       " << error_[i] << std::endl;
        if( tableWriter_ )
        {
          const double eoc = eocColumns_[ i ]->eoc( h, error_[ i ] );
          out << "EOC (" <<description_[i] << "): " << eoc << std::endl;
        }
        out << std::endl;
      }
    }

  } // end namespace Fem

} // end namespace Dune

#endif // #ifndef DUNE_FEM_FEMEOC_HH
