#ifndef DUNE_EOCOUTPUT_HH
#define DUNE_EOCOUTPUT_HH

namespace Dune{

/*======================================================================*/
/*!  @ingroup HelperClasses
 *   Write a self contained tex table and a gnuplot
 *   file for eoc runs with timing information.
 *   Constructor takes base name (filename) of file and
 *   generates two files:
 *   filename.tex and filename.gnu
 */
/*======================================================================*/
class EocOutput {
  std::string outputFile;
  int level;
  std::vector<double> prevError;
  bool initial;
  bool finalized;
 public:
  EocOutput(std::string name) :
    outputFile(name)
  , level(0)
  , prevError(0)
  , initial(true)
  , finalized(false) {
    {
      std::ostringstream filestream;
      filestream << outputFile << ".tex";
      std::ofstream ofs(filestream.str().c_str(), std::ios::out);

      ofs << "\\documentclass[12pt,english]{article}\n"
         	<< "\\usepackage[T1]{fontenc}\n"
	        << "\\usepackage[latin1]{inputenc}\n"
	        << "\\usepackage{setspace}\n"
	        << "\\onehalfspacing\n"
	        << "\\makeatletter\n"
	        << "\\providecommand{\\boldsymbol}[1]{\\mbox{\\boldmath $#1$}}\n"
	        << "\\providecommand{\\tabularnewline}{\\\\}\n"
	        << "\\usepackage{babel}\n"
	        << "\\makeatother\n"
	        << "\\begin{document}\n";
      ofs.close();	
    }
    {
      std::ostringstream filestream;
      filestream << outputFile << ".gnu";
      std::ofstream ofs(filestream.str().c_str(), std::ios::out);
      ofs.close();
    }
  }
  ~EocOutput() {
    if (!finalized) {
      std::ostringstream filestream;
      filestream << outputFile << ".tex";
      std::ofstream ofs(filestream.str().c_str(), std::ios::app);
      ofs  << "\\end{tabular}\\\\\n\n"
	         << "\\end{document}\n" << std::endl;
      ofs.close();
    }
  }
  void printTexEnd(double totaltime) {
    std::ostringstream filestream;
    filestream << outputFile << ".tex";
    std::ofstream ofs(filestream.str().c_str(), std::ios::app);
    ofs  << "\\end{tabular}\\\\\n\n"
	       << "Total time: " << totaltime << "\n"
	       << "\\end{document}\n" << std::endl;
    ofs.close();
    finalized = true;
  }
	
  template <class GridType,class ProblemType,
            class SolverType>
  void printInput(const ProblemType& problem, 
                  const GridType& grid,
                  const SolverType& solver) {
    std::ostringstream filestream;
    filestream << outputFile << ".tex";
    std::ofstream ofs(filestream.str().c_str(), std::ios::app);
    filestream << outputFile << ".gnu";
    std::ofstream ofsGnu(filestream.str().c_str(), std::ios::app);

    ofs  << "Grid: "      << grid.name() << "\\\\"
	       << "\n";
	  solver.printTexInfo(ofs);
    problem.printTexInfo(ofs);
    ofs.close();
  }
  void printTexAddError(std::vector<double> error, 
                        double time, 
                        int counter) {
    std::ostringstream filestream;
    filestream << outputFile << ".tex";
    std::ofstream ofs(filestream.str().c_str(), std::ios::app);
    std::ostringstream filestreamGnu;
    filestreamGnu << outputFile << ".gnu";
    std::ofstream ofsGnu(filestreamGnu.str().c_str(), std::ios::app);
		
    if (initial) {
	    ofs << "\\begin{tabular}{|c|c|c|";
      for (int i=0;i<error.size();i++) {
        ofs << "|cc|";
      }
      ofs << "}\n"
	        << "\\hline \n"
          << "level & CPU-time & counter";
      for (int i=0;i<error.size();i++) {
        ofs << " & error & eoc ";
      }
	    ofs << "\n \\tabularnewline\n"
	        << "\\hline\n"
	        << "\\hline\n";
    } else {
      assert(error.size()==prevError.size());
    }
	  ofs <<  "\\hline \n"
	      << level << " & "
        << time << " & " 
        << counter;
	  ofsGnu << level << "  "
           << time << "  " 
           << counter;
    for (int i=0;i<error.size();++i) {
      ofs << " & " << error[i] << " & ";
      if (initial)
        ofs << " --- ";
      else 
        ofs << log(prevError[i]/error[i])/M_LN2;
      ofsGnu << "    " << error[i] << " ";
      if (initial)
        ofsGnu << " 1/0 ";
      else 
        ofsGnu << log(prevError[i]/error[i])/M_LN2;
    }
	  ofs << "\n"
        << "\\tabularnewline\n"
	      << "\\hline \n";
    ofsGnu << "\n";
    prevError = error;
    level++;
    initial = false;
    ofs.close();
    ofsGnu.close();
  }

};

}
#endif
