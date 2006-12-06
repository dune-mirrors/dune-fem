template <class DiscreteFunctionType, class FunctionType, int polOrd>
class L2Projection
{
  typedef typename DiscreteFunctionType::FunctionSpaceType FunctionSpaceType;
  
 public:
  static void project (const FunctionType &f, DiscreteFunctionType &discFunc) {
    typedef typename DiscreteFunctionType::Traits::DiscreteFunctionSpaceType FunctionSpaceType;
    typedef typename FunctionSpaceType::Traits::GridPartType GridPartType;
    typedef typename FunctionSpaceType::Traits::GridType GridType;
    typedef typename FunctionSpaceType::Traits::IteratorType Iterator;
    
    const FunctionSpaceType& space =  discFunc.space();
    
    discFunc.clear();
    
    typedef typename DiscreteFunctionType::LocalFunctionType LocalFuncType;
    
    typename FunctionSpaceType::RangeType ret (0.0);
    typename FunctionSpaceType::RangeType phi (0.0);
    
    Iterator it = space.begin();
    Iterator endit = space.end();
    
    // Get quadrature rule
    CachingQuadrature<GridPartType,0> quad(*it, 2*polOrd+1);
    
    for( ; it != endit ; ++it) {
      LocalFuncType lf = discFunc.localFunction(*it);
      const typename FunctionSpaceType::BaseFunctionSetType & set =
	lf.baseFunctionSet();
      for(int i=0; i<lf.numDofs(); i++) {
        for(int qP = 0; qP < quad.nop(); qP++) {
	  f.evaluate((*it).geometry().global(quad.point(qP)), ret);
	  set.eval(i,quad,qP,phi);
	  lf[i] += quad.weight(qP) * (ret * phi) ;
        }
      }
    }
  }
};
struct SStruct {
  SStruct(int n1, int n2, double lx, double ly, double hx, double hy)
  {
    n_[0] = n1;
    n_[1] = n2;
    n_[2] = n2;
    l_[0] = lx;
    l_[1] = ly;
    l_[2] = ly;
    h_[0] = hx;
    h_[1] = hy;
    h_[2] = hy;
  }

  SStruct(int n, double h) {
    n_[0] = n;
    n_[1] = n;
    n_[2] = n;
    l_[0] = .0;
    l_[1] = .0; // h/2.0;
    l_[2] = .0; // h/2.0;
    h_[0] = 2.0;
    h_[1] = 2.0; // h/2.0;
    h_[2] = 2.0; // h/2.0;
  }
  SStruct(int n) {
    n_[0] = n;
    n_[1] = n;
    n_[2] = n;
    l_[0] = .0;
    l_[1] = .0; // 0.5/double(n);
    l_[2] = .0; // 0.5/double(n);
    h_[0] = 2.0;
    h_[1] = 2.0; // 0.5/double(n);
    h_[2] = 2.0; // 0.5/double(n);
  }

  int n_[3];
  double l_[3];
  double h_[3];
};
/*
template <class Geometry>
void midPoint(const Geometry& geo, FieldVector<double, 3>& result)
{
  result *= 0.0;
  for (int i = 0; i < geo.corners(); ++i) {
    result += geo[i];
  }

  result /= static_cast<double>(geo.corners());
}
template <class Geometry>
void midPoint(const Geometry& geo, FieldVector<double, 2>& result)
{
  result *= 0.0;
  for (int i = 0; i < geo.corners(); ++i) {
    result += geo[i];
  }

  result /= static_cast<double>(geo.corners());
}
template <class Geometry>
void midPoint(const Geometry& geo, FieldVector<double, 1>& result)
{
  result *= 0.0;
  for (int i = 0; i < geo.corners(); ++i) {
    result += geo[i];
  }

  result /= static_cast<double>(geo.corners());
}
*/
template <class Geometry,int n>
void midPoint(const Geometry& geo, FieldVector<double, n>& result)
{
  result *= 0.0;
  for (int i = 0; i < geo.corners(); ++i) {
    result += geo[i];
  }

  result /= static_cast<double>(geo.corners());
}
template <class Sol, class SpaceType>
void printSGrid(double time, int timestep, const SpaceType& space, const Sol& sol)
{
  typedef typename SpaceType::IteratorType Iterator;
  typedef typename Sol::DiscreteFunctionType DiscreteFunctionType;
  typedef typename DiscreteFunctionType::LocalFunctionType LocalFunctionType;

  std::ostringstream filestream;
  filestream << "sgrid" << timestep;

  std::ofstream ofs(filestream.str().c_str(), std::ios::out);

  ofs << "# " << time << std::endl;

  typename SpaceType::DomainType mid (0.0);
  typename SpaceType::DomainType localMid (0.5);
  typename SpaceType::RangeType result (0.0);

  Iterator endit = space.end();
  for (Iterator it = space.begin(); it != endit; ++it) {
    midPoint(it->geometry(), mid);
    {
      LocalFunctionType lf = sol.localFunction(*it);
      lf.evaluate(localMid, result);

      ofs << mid << " " << result << "\n";
    }
  }
  ofs << std::endl;
}


template <class StupidFunction,class DFType>
void initialize(const StupidFunction& f,DFType& df)
{
  //- Typedefs and enums
  typedef typename DFType::Traits::DiscreteFunctionSpaceType SpaceType;
  typedef typename SpaceType::Traits::IteratorType Iterator;
  typedef typename DFType::Traits::LocalFunctionType LocalFunction;
  typedef typename SpaceType::Traits::GridType Grid;
  typedef typename Grid::template Codim<0>::Entity::Geometry Geometry;

  typedef typename SpaceType::Traits::RangeType RangeType;
  typedef typename SpaceType::Traits::DomainType DomainType;

  enum { dim = Grid::dimension };

  typedef FieldVector<double, dim> Coordinate;

  //- Actual method
  L2Projection<DFType, StupidFunction, 2>::project(f, df);

  typedef typename DFType::DofIteratorType DofIterator;
  /*for (DofIterator it = df.dbegin(); it != df.dend(); ++it) {
    std::cout << *it << std::endl;
    }*/
}

template <class DFType>
void printIt(DFType& df)
{
  typedef typename DFType::DofIteratorType DofIt;

  std::cout << "print it\n";
  for (DofIt it = df.dbegin(); it != df.dend(); ++it) {
    std::cout << *it << std::endl;
  }

}


class EocOutput {

  string outputFile;
	
 public:
		
  EocOutput(string name)
  {
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
	
  void printTexEnd(double totaltime)
  {
    std::ostringstream filestream;
    filestream << outputFile;

    std::ofstream ofs(filestream.str().c_str(), std::ios::app);
		
    ofs  << "\\end{tabular}\\\\\n\n"
	 << "Total time: " << totaltime << "\n"
	 << "\\end{document}\n" << std::endl;
		
    ofs.close();
  }
	
  void printTexAddError(double error, double prevError, double time, int level, int counter,double averagedt)
  {
    std::ostringstream filestream;
    filestream << outputFile;

    std::ofstream ofs(filestream.str().c_str(), std::ios::app);
		
    if(prevError > 0.0)
      {	       
	ofs <<  "\\hline \n"
	    << level << " & " << error << " & " << log(prevError/error)/M_LN2 << " & " << time << " & " << counter << " & " << averagedt << "\n"
	    << "\\tabularnewline\n"
	    << "\\hline \n";
      }
    else
      {	       
	ofs << "\\begin{tabular}{|c|c|c|c|c|c|}\n"
	    << "\\hline \n"
	    << "Size & $\\left\\Vert u-u_{h}\\right\\Vert _{L_{2}}$ & EOC & CPU & \\#Iterations & a-dt\n"
	    << "\\tabularnewline\n"
	    << "\\hline\n"
	    << "\\hline\n"
	    << level << " & " << error << " & " << "---" << " & " << time << " & " << counter << " & " << averagedt << "\n"
	    << "\\tabularnewline\n"
	    << "\\hline \n";
      }
		
    ofs.close();
  }

		
  void printInput(InitialDataType& u0, GridType& grid,ODEType& ode, 
		  char *arg)
  {
    std::ostringstream filestream;
    filestream << outputFile;

    std::ofstream ofs(filestream.str().c_str(), std::ios::app);

    ofs  << "Grid: " << grid.name() << "\n\n"
	 << "Macrogrid: " << arg << "\\\\\n\n";
		
    ofs.close();
    ode.printmyInfo(outputFile);
    u0.printmyInfo(outputFile);		
  }
};


//Old version of latex output
/*
  void printError(double totaltime, double time[10], double error[10], \
  int maxit, InitialDataType u0, GridType *grid)
  {
  std::ostringstream filestream;
  filestream << "eoc.tex";

  std::ofstream ofs(filestream.str().c_str(), std::ios::out);
	
  ofs << "\\documentclass[12pt,english]{article}\n"
  <<	"\\usepackage[T1]{fontenc}\n"
  <<	"\\usepackage[latin1]{inputenc}\n"
  <<	"\\usepackage{setspace}\n"
  <<	"\\onehalfspacing\n"
  <<	"\\makeatletter\n"
  <<	"\\providecommand{\\boldsymbol}[1]{\\mbox{\\boldmath $#1$}}\n"
  <<	"\\providecommand{\\tabularnewline}{\\\\}\n"
  <<	"\\usepackage{babel}\n"
  <<	"\\makeatother\n"
  <<	"\\begin{document}\n"
  <<	"\\input{eoc" << u0.myName << ".tex}\n\n"
  <<	"Grid: " << transformToGridName(grid->type()) << "\n\n"
  << "Total Time: " << totaltime << "\n\n"
  << "Runge-Kutta Steps: " << rksteps << "\n\n"
  <<	"Polynomial Order: " << order << "\\\\\n\n"
  <<	"\\begin{tabular}{|c|c|c|c|}\n"
  <<	"\\hline \n"
  <<	"h & $\\left\\Vert u-u_{h}\\right\\Vert _{L_{2}}$ & EOC & CPU\n"
  <<	"\\tabularnewline\n"
  <<	"\\hline\n"
  <<	"\\hline \n"
  <<	" & " << error[0] << " & --- & " << time[0] << "\n"
  <<	"\\tabularnewline\n"
  <<	"\\hline \n";
			 
  for(int i=1;i<maxit;++i)
  { 
  ofs <<	"\\hline \n"
  <<	" & " << error[i] << " & " << log( error[i-1]/error[i])/M_LN2 << " & " << time[i] <<"\n"
  <<	"\\tabularnewline\n"
  <<	"\\hline \n";	
  }
				 
  ofs	 <<	"\\end{tabular}\n"
  <<	"\\end{document}\n" << std::endl;
  }
*/
