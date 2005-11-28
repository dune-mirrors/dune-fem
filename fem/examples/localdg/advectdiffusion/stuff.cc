template <class DiscreteFunctionType, class FunctionType, int polOrd>
class L2Projection
{
  typedef typename DiscreteFunctionType::FunctionSpaceType FunctionSpaceType;
  
 public:
  static void project (FunctionType &f, DiscreteFunctionType &discFunc) {
    typedef typename DiscreteFunctionType::Traits::DiscreteFunctionSpaceType FunctionSpaceType;
    typedef typename FunctionSpaceType::Traits::GridType GridType;
    typedef typename FunctionSpaceType::Traits::IteratorType Iterator;
    
    const FunctionSpaceType& space =  discFunc.getFunctionSpace();
    
    discFunc.clear();
    
    typedef typename DiscreteFunctionType::LocalFunctionType LocalFuncType;
    
    const int dim = GridType::dimension;
    
    typename FunctionSpaceType::RangeType ret (0.0);
    typename FunctionSpaceType::RangeType phi (0.0);
    
    Iterator it = space.begin();
    Iterator endit = space.end();
    
    // Get quadrature rule
    const QuadratureRule<double, dim>& quad =
      QuadratureRules<double, dim>::rule(it->geometry().type(), polOrd);
    
    for( ; it != endit ; ++it) {
      LocalFuncType lf = discFunc.localFunction(*it);
      
      const typename FunctionSpaceType::BaseFunctionSetType & set =
	space.getBaseFunctionSet(*it);
      
      for(int i=0; i<lf.numDofs(); i++) {
        for(int qP = 0; qP < quad.size(); qP++) {
	  double det =
	    (*it).geometry().integrationElement(quad[qP].position());
	  f.evaluate((*it).geometry().global( quad[qP].position() ), ret);
	  set.eval(i,quad[qP].position(),phi);
          lf[i] += quad[qP].weight() * (ret * phi) ;
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
    l_[0] = lx;
    l_[1] = ly;
    h_[0] = hx;
    h_[1] = hy;
  }

  SStruct(int n, double h) {
    n_[0] = n;
    n_[1] = 1;
    l_[0] = -1.0;
    l_[1] = -h/2.0;
    h_[0] = 1.5;
    h_[1] = h/2.0;
  }
  SStruct(int n) {
    n_[0] = n;
    n_[1] = 1;
    l_[0] = -1.0;
    l_[1] = -0.5/double(n);
    h_[0] = 1.5;
    h_[1] = 0.5/double(n);
  }

  int n_[2];
  double l_[2];
  double h_[2];
};
template <class Geometry>
void midPoint(const Geometry& geo, FieldVector<double, 2>& result)
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
    if (mid[1] < 0.1 && mid[1] > -0.1) {
      LocalFunctionType lf = sol.localFunction(*it);
      lf.evaluateLocal(*it, localMid, result);
      ofs << mid[0] << " " << result << "\n";
    }
  }
  ofs << std::endl;
}


template <class StupidFunction,class DFType>
void initialize(DFType& df)
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
  StupidFunction f;
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
