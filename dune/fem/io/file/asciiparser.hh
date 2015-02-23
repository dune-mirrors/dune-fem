#ifndef DUNE_FEM_ASCIIPARSER_HH
#define DUNE_FEM_ASCIIPARSER_HH

#include <iostream>
#include <fstream>
#include <string>

namespace Dune
{

  namespace Fem
  {

    //! reads data folowing the given keyword
    //! if verbose is true then an output of what was read is given
    //! the token '%' or '#' stands for comment
    template< class T >
    static inline bool readParameter( std::istream& file,
                        const std::string keyword,
                        T & data,
                        bool verbose = true,
                        bool warn = true )
    {
      bool readData = false;
      while (! file.eof() )
      {
        std::string  keyHelp;
        file >> keyHelp;

        // % or # means comment
        if((keyHelp[0] == '%') || (keyHelp[0] == '#'))
        {
          std::string tmp;
          std::getline(file,tmp);
        }
        else
        {
          // copy only keyword size and compare
          int pos = keyHelp.find_first_of(':');
          int pos1 = keyHelp.find_first_of(' ');

          if (pos > 0)
          {
            if(pos1 > 0)
              pos -=  pos1;
            else
              pos = 1;
          }
          else
            pos = 0;

          std::string key (keyHelp,0,keyHelp.size() - pos);
          if(key == keyword)
          {
            file >> data;
            readData = true;
            break;
          }
        }
      }

      // if data was read sucessfully
      if(readData)
      {
        if(verbose)
        {
          static const int MAXTAB = 30;
          int length = MAXTAB - keyword.size();
          std::cout << "Reading " << keyword;
          for(int i=0; i<length; i++) std::cout << ".";
          std::cout << " " << data << std::endl;
        }
      }
      else
      {
        if( warn )
          std::cerr << "WARNING: couldn't read " << keyword << std::endl;
      }

      return readData;
    }

    //! reads data folowing the given keyword
    //! if verbose is true then an output of what was read is given
    //! the token '%' or '#' stands for comment
    template< class T >
    static inline bool readParameter( const std::string filename,
                        const std::string keyword,
                        T & data,
                        bool verbose = true,
                        bool warn = true )
    {
      std::ifstream file (filename.c_str());
      if( !file.is_open() )
      {
        if( warn )
          std::cerr << "WARNING: couldn't open file '" << filename << "' in " <<  __FILE__<< " line: " << __LINE__ << std::endl;
        return false;
      }
      const bool result = readParameter( file, keyword, data, verbose, warn );

      // close file
      file.close();

      return result;
    }

  } // namespace Fem

} // namespace  Dune
#endif // #ifndef DUNE_FEM_ASCIIPARSER_HH
