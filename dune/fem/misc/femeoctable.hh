#ifndef DUNE_FEM_FEMEOCTABLE_HH
#define DUNE_FEM_FEMEOCTABLE_HH

#include <cassert>
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>

#include <dune/fem/io/io.hh>
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
        generate the files:
        filename.dat

        The class is singleton but severeal numbers of tables can be written with different
        errors. Each time the class is initialized a new table will be written.
        Since the class is singelton new errors for eoc
        computations can be added in any part of the program for each table, specified by
        the table Id returned by the initializer.
        To add a new entry for eoc computations use one of the
        addEntry methods. These return a unique unsigned int
        which can be used to add error values to the table
        with the setErrors methods.
        The method write is used to write a single line
        to the eoc table.
     */

    /*
    //! Base Class of all eoc calculator types
    class BaseEocCalculator
    {
      public:
      static double calculate (double,double,double,double)=0;
    }
    */

    //! Default eoc calculator, it can be re-implemented if another Eoc calulator is needed
    class DefaultEocCalculator
    {
      public:
      static double calculate (double &eold, double &enew, double &hold, double &hnew )
      {
        double ret=0;
        ret = log(eold/enew)/log(hold/hnew);
        return ret;
      }
    };



    //! The Fem Eoc Table writer
    class FemEocTable
    {
      int nrOfTabs_;
      std::vector< std::stringstream* > outputFile_;
      std::vector< std::string > fileNames_;
      std::vector< int > level_;
      std::vector< std::vector< double > > prevError_;
      std::vector< std::vector< double > > error_;
      std::vector< std::vector< std::string > > description_;
      std::vector< double > prevh_;
      std::vector< bool > initial_;
      std::vector< std::vector< int > > pos_;

    public:
      FemEocTable() :
        nrOfTabs_(0),
        outputFile_(),
        fileNames_(),
        level_(),
        prevError_(),
        error_(),
        description_(),
        prevh_(),
        initial_(),
        pos_()
      {
      }

      ~FemEocTable() {
        std::ofstream filewriter;

        for(size_t k=0;k<outputFile_.size();++k)
        {
          std::stringstream filename;
          filename << fileNames_[k] <<".dat";
          filewriter.open( filename.str().c_str() );
          filewriter <<  outputFile_[k]->str();
          filewriter.close();
          delete outputFile_[k];
          outputFile_[k] = 0;
        }
      }

    private:
      int init(const std::string& path,
               const std::string& name, const std::string& descript)
      {
        if (MPIManager::rank() != 0) return -1;
        if( !createDirectory( path ) )
          DUNE_THROW( IOError, "Failed to create path `" << path << "'." );
        return init(path+"/"+name,descript);
      }

      int init(const std::string& filename, const std::string& descript)
      {
        if (MPIManager::rank() != 0) return -1;

        //! create new vector components
        outputFile_.push_back((std::stringstream *) 0);

        fileNames_.push_back("");
        level_.push_back(0);
        prevError_.push_back( std::vector< double >() );
        error_.push_back( std::vector< double >() );
        description_.push_back( std::vector< std::string >() );
        prevh_.push_back(0);
        initial_.push_back(1);
        pos_.push_back( std::vector< int >() );

        const int tabId = outputFile_.size() -1;

        fileNames_[tabId] = filename;
        outputFile_[tabId] = new std::stringstream;

        // write together with table the description
        // of the scheme and problem in the simulation
        *outputFile_[tabId] << descript << "\n";
        nrOfTabs_ ++;

        return tabId;
      }

      void checkTabId (const int tabId)
      {
        if (tabId > nrOfTabs_ || tabId <0)
        {
          std::cout<<"No table with id:"<<tabId<<" existing!"<<std::endl;
            abort();
        }
      }

      template <class StrVectorType>
      size_t addentry(const int tabId, const StrVectorType& descript,size_t size)
      {
        checkTabId(tabId);

        if (!initial_[tabId])
          abort();
        pos_[tabId].push_back(error_[tabId].size());
        for (size_t i=0;i<size;++i) {
          error_[tabId].push_back(0);
          prevError_[tabId].push_back(0);
          description_[tabId].push_back(descript[i]);
        }
        return pos_[tabId].size()-1;
      }

      size_t addentry(const int tabId, const std::string& descript) {

        checkTabId(tabId);

        if (!initial_[tabId])
          abort();
        pos_[tabId].push_back(error_[tabId].size());
        error_[tabId].push_back(0);
        prevError_[tabId].push_back(0);
        description_[tabId].push_back(descript);
        return pos_[tabId].size()-1;
      }

      template <class VectorType>
      void seterrors(const int tabId, size_t id,const VectorType& err,size_t size)
      {
        checkTabId(tabId);

        assert(id<pos_[tabId].size());
        int pos = pos_[tabId][ id ];
        assert(pos+size <= error_[tabId].size());

        for (size_t i=0; i<size; ++i)
          error_[tabId][pos+i] = err[i];
      }

      template <int SIZE>
      void seterrors(const int tabId, size_t id,const FieldVector<double,SIZE>& err)
      {
        seterrors(tabId,id,err,SIZE);
      }

      void seterrors(const int tabId, size_t id,const double& err) {
        checkTabId(tabId);

        int pos = pos_[tabId][id];
        error_[tabId][pos] = err;
      }

      //! write data to file
      template<class EocCalculator>
      void writeerr( const int tabId,
                     std::vector<double> &vals,
                     std::vector<std::string> &descriptions,
                     std::string &delimiter,
                     std::string &terminatingChar,
                     std::string &header,
                     std::string &tableSpacer,
                     std::string &footer)

      {

        typedef EocCalculator EocCalculatorType;

        checkTabId(tabId);

        assert(vals.size() == descriptions.size());
        if (MPIManager::rank() != 0) return;

        if (initial_[tabId]) {
          *outputFile_[tabId] << header;
          for(unsigned int k=0;k<descriptions.size();++k)
          {
           *outputFile_[tabId] <<  descriptions[k] << delimiter;
          }
          for (unsigned int i=0;i<error_[tabId].size();i++)
          {
            *outputFile_[tabId] <<  description_[tabId][i] <<  delimiter << "EOC" << delimiter;
          }
          *outputFile_[tabId] << terminatingChar << "\n" << tableSpacer <<"\n";
        }

        *outputFile_[tabId] << level_[tabId] << delimiter;
        for(unsigned int k =0; k<vals.size(); ++k)
          *outputFile_[tabId] << vals[k] << delimiter;

        for (unsigned int i=0;i<error_[tabId].size();++i) {
          *outputFile_[tabId] <<delimiter << error_[tabId][i] << delimiter;
          if (initial_[tabId]) {
            *outputFile_[tabId] << " --- ";
          }
          else {
           *outputFile_[tabId] <<  EocCalculatorType :: calculate(prevError_[tabId][i], error_[tabId][i], prevh_[tabId], vals[0] );
          }
          prevError_[tabId][i]=error_[tabId][i];
          error_[tabId][i] = -1;  // uninitialized
        }
        *outputFile_[tabId] << terminatingChar<<"\n" << footer;

        //! set back the put point in the stream to prohibit writing the footer each time in
        // the file
        outputFile_[tabId]->seekp(0,std::ios::end);
        int length = outputFile_[tabId]->tellp();
        length -= footer.length();
        outputFile_[tabId]->seekp(length, std::ios::beg);


        prevh_[tabId] = vals[0];
        level_[tabId] ++;
        initial_[tabId] = false;
      }

      template<class EocCalculator>
      void printerr(const int tabId,
                    std::vector<double> vals,
                    std::vector<std::string> descriptions,
                    std::ostream& out)
      {
        typedef EocCalculator EocCalculatorType;

        checkTabId(tabId);

        assert(descriptions.size() == vals.size());

        if (!Parameter::verbose()) return;

        out << "level:   " << level_[tabId]  << std::endl;
        for(unsigned int k =0 ;k< vals.size();++k)
          out << descriptions[k]<<":       " << vals[k] << std::endl;

        for (unsigned int i=0;i<error_[tabId].size();++i)
        {

          out << description_[tabId][i] << ":       " << error_[tabId][i] << std::endl;
          if (! initial_[tabId])
          {
            const double eoc = EocCalculatorType :: calculate( prevError_[tabId][i], error_[tabId][i],  prevh_[tabId], vals[0]);

            out << "EOC (" <<description_[tabId][i] << "): " << eoc << std::endl;
          }
          out << std::endl;
        }
      }
     public:
      friend class Dune::Fem::Singleton< FemEocTable >;

      static FemEocTable& instance()
      {
        return Singleton< FemEocTable >::instance();
      }

      //! creates a new table and opens the corresponding file path/name
      //! returns the Id of the created table.
      static int initialize(const std::string& path, const std::string& name, const std::string& descript) {
        return instance().init(path,name,descript);
      }
      //! creates a new table and opens file name
      //! as above returns the Id of the table opened.
      static int initialize(const std::string& name, const std::string& descript) {
        return instance().init(name,descript);
      }

      /** \brief add a vector of new eoc values
       *
       *  \param tabId Id of the table we inserte a value
       *  \tparam  StrVectorType a vector type with operator[]
       *           returning a string (a C style array can be used)
       *           the size of the vector is given as parameter
       *  \param descript vector with entry description
       *  \param size length of description
       *  \return  a unique index used to add the error values
       */
      template <class StrVectorType>
      static size_t addEntry(const int tabId, const StrVectorType& descript,size_t size) {
        return instance().addentry(tabId,descript,size);
      }

      template <class StrVectorType>
      static size_t addEntry(const StrVectorType& descript,size_t size) {
        return instance().addentry(0,descript,size);
      }

      /** \brief add a vector of new eoc values
       *
       *  \param tabId Id of the table we inserte a value
       *  \tparam  StrVectorType a vector type with size() and operator[]
       *           returning a string
       *  \param descript vector with entry description
       *  \return  a unique index used to add the error values
       */
      template <class StrVectorType>
      static size_t addEntry(const int tabId,const StrVectorType& descript) {
        return instance().addentry(tabId,descript,descript.size());
      }

      template <class StrVectorType>
      static size_t addEntry(const StrVectorType& descript) {
        return instance().addentry(0,descript,descript.size());
      }
      /** \brief add a single new eoc output
       *
       *  \param tabId Id of the table we want to add an entry
       *  \param descript vector with entry description
       *  \return  a unique index used to add the error values
       */
      static size_t addEntry(const int tabId, const std::string& descript) {
        return instance().addentry(tabId, descript);
      }

      static size_t addEntry(const std::string& descript) {
        return instance().addentry(0, descript);
      }

      /** \brief add a single new eoc output
       *
       *  \param tabId Id of the table we want to add an entry
       *  \param descript vector with entry description
       *  \return  a unique index used to add the error values
       */
      static size_t addEntry(const int tabId, const char* descript) {
        return addEntry(tabId,std::string(descript));
      }

      static size_t addEntry(const char* descript) {
        return addEntry(0,std::string(descript));
      }


      /** \brief add a vector of error values for the given id (returned by
       *         addEntry)
       *  \param tabId Id of the table we want to set errors
       *  \param id    Id of the error
       *  \param err   Vector containing the error
       *  \param size  Size of error Vector
       *  \tparam  VectorType a vector type with an operator[]
       *           returning a double (C style array can be used)
       */
      template <class VectorType>
      static void setErrors(const int tabId, size_t id,const VectorType& err,int size)
      {
        instance().seterrors(tabId,id,err,size);
      }

      template <class VectorType>
      static void setErrors(size_t id,const VectorType& err,int size)
      {
        instance().seterrors(0,id,err,size);
      }


      /** \brief add a vector of error values for the given id (returned by
       *         addEntry)
       *  \param tabId Id of the table we want to set errors
       *  \param id    Id of the error
       *  \param err   Vector containing the error
       *  \tparam  VectorType a vector type with a size() and an operator[]
       *           returning a double
       */
      template <class VectorType>
      static void setErrors(const int tabId, size_t id,const VectorType& err) {
        instance().seterrors(tabId,id,err,err.size());
      }

      template <class VectorType>
      static void setErrors(size_t id,const VectorType& err) {
        instance().seterrors(0,id,err,err.size());
      }


      /** \brief add a vector in a FieldVector of error values for the given id (returned by
       *         addEntry)
       *  \param tabId Id of the table we want to set errors
       *  \param id    Id of the error
       *  \param err   Vector containing the error
       */
      template <int SIZE>
      static void setErrors(const int tabId, size_t id,const FieldVector<double,SIZE>& err) {
        instance().seterrors(tabId,id,err);
      }

      template <int SIZE>
      static void setErrors(size_t id,const FieldVector<double,SIZE>& err) {
        instance().seterrors(0,id,err);
      }
      /** \brief add a single error value for the given id (returned by
       *         addEntry)
       *  \param tabId Id of the table we want to set errors
       *  \param id    Id of the error
       *  \param err   Vector containing the error
       */
      static void setErrors(const int tabId,size_t id,const double& err) {
        instance().seterrors(tabId,id,err);
      }

      static void setErrors(size_t id,const double& err) {
        instance().seterrors(0,id,err);
      }

      /** \brief commit a line to the eoc file
       *
       *
       *  \param tabId table Id returned by the initial function
       *  \param vals  std::vector of vals that should appear in the EOC table, vals[0]  is
       *         expected to be a charateristical value.
       *  \param descriptions std::vector with descriptions of the values that should appear
       *  \param terminatingChar  char which ends an entry, default =  " "
       *  \param header header string for Latex output, default = ""
       *  \param tableSpacer spacer for empty columns in the table, default = ""
       *  \param footer footer string for Latex output, default = ""
       *  \param delimiter spacer between the entries, default =" "
       *
       *
       */
      static void write(const int tabId,
                        std::vector<double> &vals,
                        std::vector<std::string> &descriptions,
                        std::string delimiter = " ",
                        std::string terminatingChar = "",
                        std::string header ="",
                        std::string tableSpacer ="",
                        std::string footer ="" )
      {
        instance().writeerr<DefaultEocCalculator> (tabId, vals, descriptions, delimiter, terminatingChar, header, tableSpacer, footer);
      }

      static void write(std::vector<double> &vals,
                        std::vector<std::string> &descriptions,
                        std::string delimiter = " ",
                        std::string terminatingChar ="",
                        std::string header ="",
                        std::string tableSpacer ="",
                        std::string footer =""
                        )
      {
        instance().writeerr<DefaultEocCalculator> (0, vals, descriptions, delimiter, terminatingChar, header, tableSpacer, footer);
      }


     /** \brief commit a line to the eoc file, using EocCalculatorType to calculate the eoc.
      *
      *
      *  \param tabId table Id returned by the initial function
      *  \param vals  std::vector of vals that should appear in the EOC table, vals[0]  is
      *         expected to be a charateristical value.
      *  \param descriptions std::vector with descriptions of the values that should appear
      *  \param terminatingChar  char which ends an entry, default =  " "
      *  \param header header string for Latex output, default = ""
      *  \param tableSpacer spacer for empty columns in the table, default = ""
      *  \param footer footer string for Latex output, default = ""
      *  \param delimiter spacer between the entries, default =" "
      *
      */
      template<class EocCalculatorType>
      static void write(const int tabId,
                        std::vector<double> &vals,
                        std::vector<std::string> &descriptions,
                        std::string delimiter = " ",
                        std::string terminatingChar ="",
                        std::string header ="",
                        std::string tableSpacer ="",
                        std::string footer =""
                        )
      {
        instance().template writeerr<EocCalculatorType>(tabId, vals, descriptions, delimiter, terminatingChar, header, tableSpacer, footer);
      }

      template<class EocCalculatorType>
      static void write(std::vector<double> &vals,
                        std::vector<std::string> &descriptions,
                        std::string delimiter = " ",
                        std::string terminatingChar ="",
                        std::string header ="",
                        std::string tableSpacer ="",
                        std::string footer =""
                        )
      {
        instance().template writeerr<EocCalculatorType>(0, vals, descriptions, delimiter, terminatingChar, header, tableSpacer, footer );
      }

      /** \brief commit a line to the eoc file
       *
       *  \param tabId table Id returned by the initial function
       *  \param vals  std::vector of vals that should appear in the EOC table
       *  \param descriptions std::vector with descriptions of the values that should appear
       *  \param terminatingChar  char which ends an entry, default =  " "
       *  \param out std::ostream to print data to (e.g. std::cout)
       *  \param header header string for Latex output, default = " "
       *  \param tableSpacer spacer for empty columns in the table, default = " "
       *  \param footer footer string for Latex output, default = " "
       *  \param delimiter spacer between the entries, default " "
       *
       */
      static void write(const int tabId,
                        std::vector<double> &vals,
                        std::vector<std::string> &descriptions,
                        std::ostream& out,
                        std::string delimiter = " ",
                        std::string terminatingChar ="",
                        std::string header = "",
                        std::string tableSpacer = "",
                        std::string footer = "")
      {
        // print last line to out
        instance().printerr<DefaultEocCalculator>( tabId, vals, descriptions, out );

        // now write to file
        instance().writeerr<DefaultEocCalculator>(tabId, vals, descriptions, delimiter, terminatingChar, header, tableSpacer, footer);
      }

      static void write(std::vector<double> &vals,
                        std::vector<std::string> &descriptions,
                        std::ostream& out,
                        std::string delimiter = " ",
                        std::string terminatingChar ="",
                        std::string header = "",
                        std::string tableSpacer = "",
                        std::string footer = "")
      {
        // print last line to out
        instance().printerr<DefaultEocCalculator>( 0, vals, descriptions, out );

        // now write to file
        instance().writeerr<DefaultEocCalculator>(0, vals, descriptions, delimiter, terminatingChar, header, tableSpacer, footer);
      }


      /** \brief commit a line to the eoc file, using EocCalculatorType for non standart Eoc calculations.
       *
       *  \param tabId table Id returned by the initial function
       *  \param vals  std::vector of vals that should appear in the EOC table
       *  \param descriptions std::vector with descriptions of the values that should appear
       *  \param terminatingChar  char which ends an entry, default =  " "
       *  \param out std::ostream to print data to (e.g. std::cout)
       *  \param header header string for Latex output, default = " "
       *  \param tableSpacer spacer for empty columns in the table, default = " "
       *  \param footer footer string for Latex output, default = " "
       *  \param delimiter spacer between the entries, default " "
       *
       */
      template <class EocCalculatorType>
      static void write(const int tabId,
                        std::vector<double> &vals,
                        std::vector<std::string> &descriptions,
                        std::ostream& out,
                        std::string delimiter = " ",
                        std::string terminatingChar ="",
                        std::string header = "",
                        std::string tableSpacer = "",
                        std::string footer = "")
      {
        // print last line to out
        instance().template printerr<EocCalculatorType>( tabId, vals, descriptions, out );

        // now write to file
        instance().template writeerr<EocCalculatorType>(tabId, vals, descriptions, delimiter, terminatingChar, header, tableSpacer, footer);
      }

      template <class EocCalculatorType>
      static void write(std::vector<double> &vals,
                        std::vector<std::string> &descriptions,
                        std::ostream& out,
                        std::string delimiter = " ",
                        std::string terminatingChar ="",
                        std::string header = "",
                        std::string tableSpacer = "",
                        std::string footer = "")
      {
        // print last line to out
        instance().template printerr<EocCalculatorType>(0, vals, descriptions, out );

        // now write to file
        instance().template writeerr<EocCalculatorType>(0, vals, descriptions,delimiter, terminatingChar, header, tableSpacer, footer);
      }

    }; // class FemEocTable

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_FEMEOCTABLE_HH
