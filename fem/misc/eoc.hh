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
    outputFile = name;
    std::ostringstream filestream;
    filestream << outputFile;
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
  ~EocOutput() {
    if (!finalized) {
      std::ostringstream filestream;
      filestream << outputFile;
      std::ofstream ofs(filestream.str().c_str(), std::ios::app);
		
      ofs  << "\\end{tabular}\\\\\n\n"
	         << "\\end{document}\n" << std::endl;
		
      ofs.close();
    }
  }
  void printTexEnd(double totaltime) {
    std::ostringstream filestream;
    filestream << outputFile;

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
    filestream << outputFile;

    std::ofstream ofs(filestream.str().c_str(), std::ios::app);

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
    filestream << outputFile;
    std::ofstream ofs(filestream.str().c_str(), std::ios::app);
		
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
    for (int i=0;i<error.size();++i) {
      ofs << " & " << error[i] << " & ";
      if (initial)
        ofs << " --- ";
      else 
        ofs << log(prevError[i]/error[i])/M_LN2;
    }
	  ofs << "\n"
        << "\\tabularnewline\n"
	      << "\\hline \n";
    prevError = error;
    level++;
    initial = false;
    ofs.close();
  }

};

