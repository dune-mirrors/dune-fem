#ifndef DUNE_FEM_IO_HH
#define DUNE_FEM_IO_HH

#include <iostream>

namespace Dune
{

  /** \brief create a directory
   *
   *  \param[in]  name  name of the directory to create
   *
   *  \returns whether the directory has been successfully created
   */
  bool createDirectory ( const std::string &name );

  /** \brief check whether a directory exists 
   *
   *  \param[in]  name  name of the directory to create
   *
   *  \returns true if directory exists, false otherwise 
   */
  bool directoryExists ( const std::string &name );


  /** \brief executes a command and return the output  
   *
   *  \param[in]  command   command to execute
   *
   *  \returns the output of the executed command
   *
   *  \note  This command throws an exception if the command 
   *         returns unsuccsessfully.
   */
  std::string executeCommand ( const std::string &command );
}

#endif // #ifndef DUNE_FEM_IO_HH
