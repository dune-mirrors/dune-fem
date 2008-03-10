template <class DiscreteFunctionType, class FunctionType, int polOrd>
class L2ProjectionLocal
{
  typedef typename DiscreteFunctionType::FunctionSpaceType FunctionSpaceType;
  
 public:
  static void project (const FunctionType &f, DiscreteFunctionType &discFunc) 
  {
    typedef typename DiscreteFunctionType::Traits::DiscreteFunctionSpaceType FunctionSpaceType;
    typedef typename FunctionSpaceType::Traits::GridPartType GridPartType;
    typedef typename FunctionSpaceType::Traits::GridType GridType;
    typedef typename GridType :: template Codim<0> :: Geometry Geometry;
    typedef typename FunctionSpaceType::Traits::IteratorType Iterator;
    
    const FunctionSpaceType& space =  discFunc.space();
    
    discFunc.clear();
    
    typedef typename DiscreteFunctionType::LocalFunctionType LocalFuncType;
    
    typename FunctionSpaceType::RangeType ret (0.0);
    typename FunctionSpaceType::RangeType phi (0.0);
    
    Iterator it = space.begin();
    Iterator endit = space.end();

    // if empty grid, do nothing 
    if( it == endit ) return ;
    
    // Get quadrature rule
    CachingQuadrature<GridPartType,0> quad(*it, 2*polOrd+1);
    
    const int quadNop = quad.nop();
    for( ; it != endit ; ++it) 
    {
      LocalFuncType lf = discFunc.localFunction(*it);

      const Geometry& geo = (*it).geometry();
      for(int qP = 0; qP < quadNop; ++qP) 
      {
	      f.evaluate(geo.global(quad.point(qP)), ret);
        ret *=  quad.weight(qP);

        lf.axpy( quad[ qP ], ret );
      }
    }
  }
};
template <class Sol, class SpaceType>
void printSGrid(double time, int timestep, const SpaceType& space, Sol& sol)
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
  L2ProjectionLocal<DFType, StupidFunction, 2>::project(f, df);

  typedef typename DFType::DofIteratorType DofIterator;
  /*for (DofIterator it = df.dbegin(); it != df.dend(); ++it) {
    std::cout << *it << std::endl;
    }*/
}

class EocOutput {

	std::string outputFile;
	
 public:
		
  EocOutput(std::string name)
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

