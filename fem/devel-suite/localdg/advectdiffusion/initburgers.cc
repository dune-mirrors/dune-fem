
template <class GridType>
class U0 {
  double startTime_;
public:
  enum { ConstantVelocity = true };
  enum { dimDomain = GridType::dimensionworld };  
  typedef FieldVector<double,dimDomain> DomainType;
  typedef FieldVector<double,1> RangeType;
  U0() :
    startTime_(Parameter::getValue<double>("fem.localdg.starttime",0.0)) { 
    velocity(0) {  
      Parameter::get("fem.localdg.epsilon", epsilon); 
      myName = "Burgers-Diffusion";
    }
  double endtime() {
    return 0.4;
  }
  void evaluate(const DomainType& arg, RangeType& res) const {
    evaluate(startTimne_,arg,res);
  }

  void evaluate(const DomainType& arg, double t, RangeType& res) const
  {
    evaluate(t, arg, res);
    return;
  }

  void evaluate(double t,const DomainType& arg, RangeType& res) const {
    if (std::abs(arg[0])<0.5) 
      res = 0.5*(cos(2.*arg[0]*M_PI)+1.);
    else
      res = 0.;
    return;
    if(epsilon < 1e-9)
      res = ((2.0*arg[0]-1.0)<0.0)? 1.0:-1.0;
    else
      res = -tanh((2.0*arg[0]-1.0)/(4.*epsilon));
  }
  
  void printmyInfo(string filename)
  {
    std::ostringstream filestream;
    filestream << filename;
    
    std::ofstream ofs(filestream.str().c_str(), std::ios::app);
    
    ofs << "Problem: " << myName << "\n\n"
	<< "Epsilon = " << epsilon << "\n\n"
	<< "Exact solution: $u(x,y,z,t):=\\displaystyle{-\\tanh\\left( \\frac{2x-1}{4\\varepsilon} \\right) }$\\\\\n\n";
    
    ofs.close();
  }
  
  
  DomainType velocity;
  double epsilon;
  string myName;
};
