#ifndef DUNE_FEMEOC_HH
#define DUNE_FEMEOC_HH

#include <cassert>
#include <sstream>
#include <fstream>
#include <vector>

#include <dune/common/fvector.hh>
#include <dune/fem/io/file/iointerface.hh>

namespace Dune
{

/**  
    @ingroup HelperClasses
    \brief Write a self contained tex table and a gnuplot
    file for eoc runs with timing information.
    
    Constructor takes base name (filename) of file and
    generates two files:
    filename.tex, filename_body.tex and filename.gnu
    The file filename_body.tex hold the actual body
    of the tex file which is included in filename.tex
    but can also be used to combine runs with different
    parameters.
 */
class FemEoc
{
  std::ofstream outputFile_;
  int level_;
  std::vector<double> prevError_;
  std::vector<double> error_;
  std::vector<std::string> description_;
  double prevh_;
  bool initial_;
  std::vector<int> pos_;
  FemEoc() :
    outputFile_()
  , level_(0)
  , prevError_(0)
  , error_(0)
  , description_(0)
  , prevh_(0)
  , initial_(true)
  , pos_(0)
  {
  }
  ~FemEoc() {
    outputFile_.close();
  }
  void init(const std::string& path,
            const std::string& name, const std::string& descript) {
    IOInterface::createPath(path);
    init(path+"/"+name,descript);
  }
  void init(const std::string& name, const std::string& descript) {
    if (!outputFile_.is_open()) {
      std::ofstream main((name+"_main.tex").c_str());
      if (!main) {
        std::cerr << "Could not open file : "
                  << (name+"_main.tex").c_str() 
                  << " ... ABORTING" << std::endl;
        abort();
      }
      std::ostringstream filestreamBody;
      filestreamBody << name << "_body.tex";
      outputFile_.open(filestreamBody.str().c_str(), std::ios::out);
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
           << descript 
           << "\n\\end{center}\n\n"
           << "\\input{"
           << filestreamBody.str() 
           << "}\n";
      main << "\\end{tabular}\\\\\n\n"
	         << "\\end{document}\n" << std::endl;
      main.close();	
    } else {
      abort();
    }
  }
  template <class StrVectorType>
  size_t addentry(const StrVectorType& descript,size_t size) {
    if (!initial_) 
      abort();
    pos_.push_back(error_.size());
    for (size_t i=0;i<size;++i) {
      error_.push_back(0);
      prevError_.push_back(0);
      description_.push_back(descript[i]);  
    }
    return pos_.size()-1;
  }
  size_t addentry(const std::string& descript) {
    if (!initial_) 
      abort();
    pos_.push_back(error_.size());
    error_.push_back(0);
    prevError_.push_back(0);
    description_.push_back(descript);  
    return pos_.size()-1;
  }
  template <class VectorType>
  void seterrors(size_t id,const VectorType& err,size_t size) {
    int pos = pos_[id];
    for (size_t i=0;i<size;++i)
      error_[pos+i] = err[i];
  }
  template <int SIZE>
  void seterrors(size_t id,const FieldVector<double,SIZE>& err) {
    int pos = pos_[id];
    for (int i=0;i<SIZE;++i)
      error_[pos+i] = err[i];
  }
  void seterrors(size_t id,const double& err) {
    int pos = pos_[id];
    error_[pos] = err;
  }
  void writeerr(double h,double size,double time,int counter) {
    if (initial_) {
	    outputFile_ << "\\begin{tabular}{|c|c|c|c|c|";
      for (unsigned int i=0;i<error_.size();i++) {
        outputFile_ << "|cc|";
      }
      outputFile_ << "}\n"
	        << "\\hline \n"
          << "level & h & size & CPU-time & counter";
      for (unsigned int i=0;i<error_.size();i++) {
        outputFile_ << " & " << description_[i]
                    << " & EOC ";
      }
	    outputFile_ << "\n \\tabularnewline\n"
	                << "\\hline\n"
	                << "\\hline\n";
    }
	  outputFile_ <<  "\\hline \n"
	              << level_ << " & "
                << h      << " & "
                << size   << " & "
                << time   << " & " 
                << counter;
    for (unsigned int i=0;i<error_.size();++i) {
      outputFile_ << " & " << error_[i] << " & ";
      if (initial_) {
        outputFile_ << " --- ";
      }
      else {
        double factor = prevh_/h;
        outputFile_ << log(prevError_[i]/error_[i])/log(factor);
      }
      prevError_[i]=error_[i];
    }
	  outputFile_ << "\n"
                << "\\tabularnewline\n"
	              << "\\hline \n";
    outputFile_.flush();
    prevh_ = h;
    level_++;
    initial_ = false;
  }
 public:
  static FemEoc& instance() {
    static FemEoc instance_;
    return instance_;
  }
  static void initialize(const std::string& path, const std::string& name, const std::string& descript) {
    instance().init(path,name,descript);
  }
  static void initialize(const std::string& name, const std::string& descript) {
    instance().init(name,descript);
  }
  template <class StrVectorType>
  static size_t addEntry(const StrVectorType& descript,size_t size) {
    return instance().addentry(descript,size);
  }
  template <class StrVectorType>
  static size_t addEntry(const StrVectorType& descript) {
    return instance().addentry(descript,descript.size());
  }
  static size_t addEntry(const std::string& descript) {
    return instance().addentry(descript);
  }
  template <class VectorType>
  static void setErrors(size_t id,const VectorType& err,int size) {
    instance().seterrors(id,err,size);
  }
  template <class VectorType>
  static void setErrors(size_t id,const VectorType& err) {
    instance().seterrors(id,err,err.size());
  }
  template <int SIZE>
  static void setErrors(size_t id,const FieldVector<double,SIZE>& err) {
    instance().seterrors(id,err);
  }
  static void setErrors(size_t id,const double& err) {
    instance().seterrors(id,err);
  }
  static void write(double h,double size,double time,int counter) {
    instance().writeerr(h,size,time,counter);
  }
};

}
#endif
