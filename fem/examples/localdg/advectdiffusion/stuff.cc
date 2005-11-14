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
template <class Loop, class SpaceType>
void printSGrid(double time, int timestep, SpaceType& space, Loop& loop)
{
  typedef typename SpaceType::IteratorType Iterator;
  typedef typename Loop::DiscreteFunctionType DiscreteFunctionType;
  typedef typename DiscreteFunctionType::LocalFunctionType LocalFunctionType;

  std::ostringstream filestream;
  filestream << "sgrid" << timestep;

  std::ofstream ofs(filestream.str().c_str(), std::ios::out);

  FieldVector<double, 2> mid(0.0);
  FieldVector<double, 2> localMid(0.5);

  FieldVector<double, 1> result;

  DiscreteFunctionType& sol = loop;

  Iterator endit = space.end();
  for (Iterator it = space.begin(); it != endit; ++it) {
    midPoint(it->geometry(), mid);
    if (mid[1] < 0.1 && mid[1] > -0.1) {
      LocalFunctionType lf = sol.localFunction(*it);
      lf.evaluateLocal(*it, localMid, result);
      ofs << mid[0] << " " << result[0] << "\n";
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
