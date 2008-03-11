template <class GridType>
class U0 {
  double startTime_;
 public:
  enum { ConstantVelocity = true };
  enum { dimDomain = GridType::dimensionworld };  
  typedef FieldVector<double,dimDomain> DomainType;
  typedef FieldVector<double,1> RangeType;
  U0() : startTime_(Parameter::getValue<double>("fem.localdg.starttime",0.0)),
	 velocity_(0) { 
    Parameter::get("fem.localdg.epsilon",epsilon);
    velocity_[0]=0.8;
    velocity_[1]=0.8;
    
    max_n_of_coefs = 2;
    
    //x coordinate
    common_coef_x[0] = 2.0;
    sin_coef_x[0] = 0.8;
    cos_coef_x[0] = 0.6;
    
    //x coordinate
    common_coef_x[1] = 0.7;
    sin_coef_x[1] = 0.2;
    cos_coef_x[1] = 0.9;
    
    //y coordinate    
    common_coef_y[0] = 1.0;
    sin_coef_y[0] = 0.4;
    cos_coef_y[0] = 1.2;
    
    
    //y coordinate    
    common_coef_y[1] = 0.5;
    sin_coef_y[1] = 0.1;
    cos_coef_y[1] = 0.3;
    
    
    //z coordinate
    common_coef_z[0] = 1.3;
    sin_coef_z[0] = -0.4;
    cos_coef_z[0] = 0.1;
    
    //z coordinate    
    common_coef_z[1] = 0.1;
    sin_coef_z[1] = 0.2;
    cos_coef_z[1] = -0.3;
    
    myName = "AdvectDiff";
  }
  
  void velocity(const DomainType& x, DomainType& v) const
  {
    v = velocity_;
  }
  
  double endtime() {
    return 0.1;
  }
    
  void evaluate(const DomainType& arg, RangeType& res) const {
    evaluate(startTime_,arg,res);
  }
  
  
  void evaluate(const DomainType& arg, double t, RangeType& res) const 
	{
	  evaluate(t, arg, res);
		return;
	}
	
  void evaluate(double t,const DomainType& arg, RangeType& res) const {
    
    res = 0.;

    for(int i=0;i<max_n_of_coefs;++i)
      {
  if(dimDomain == 1)
    res += exp(-epsilon*t*(SQR(common_coef_x[i]*M_PI)))\
      *((cos_coef_x[i]*cos(common_coef_x[i]*M_PI*(arg[0]-velocity_[0]*t))\
         +  sin_coef_x[i]*sin(common_coef_x[i]*M_PI*(arg[0]-velocity_[0]*t))));
  else if(dimDomain == 2)
    res += exp(-epsilon*t*(SQR(common_coef_x[i]*M_PI)+SQR(common_coef_y[i]*M_PI)))\
      *((cos_coef_x[i]*cos(common_coef_x[i]*M_PI*(arg[0]-velocity_[0]*t))\
         +  sin_coef_x[i]*sin(common_coef_x[i]*M_PI*(arg[0]-velocity_[0]*t)))\
        \
        *(cos_coef_y[i]*cos(common_coef_y[i]*M_PI*(arg[1]-velocity_[1]*t))\
    + sin_coef_y[i]*sin(common_coef_y[i]*M_PI*(arg[1]-velocity_[1]*t))));
  else if(dimDomain == 3)
    res += exp(-epsilon*t*(SQR(common_coef_x[i]*M_PI)+SQR(common_coef_y[i]*M_PI)+SQR(common_coef_z[i]*M_PI)))\
      *((cos_coef_x[i]*cos(common_coef_x[i]*M_PI*(arg[0]-velocity_[0]*t))\
         +  sin_coef_x[i]*sin(common_coef_x[i]*M_PI*(arg[0]-velocity_[0]*t)))\
        \
        *(cos_coef_y[i]*cos(common_coef_y[i]*M_PI*(arg[1]-velocity_[1]*t))\
    + sin_coef_y[i]*sin(common_coef_y[i]*M_PI*(arg[1]-velocity_[1]*t)))\
        \
        *(cos_coef_z[i]*cos(common_coef_z[i]*M_PI*(arg[2]-velocity_[2]*t))\
    + sin_coef_z[i]*sin(common_coef_z[i]*M_PI*(arg[2]-velocity_[2]*t))));
      }
  }
  
  void printmyInfo(std::string filename)
  {
    std::ostringstream filestream;
    filestream << filename;

    std::ofstream ofs(filestream.str().c_str(), std::ios::app);
  
    ofs << "Problem: " << myName << "\n\n"
  << "Epsilon = " << epsilon << "\n\n"
  << "Exact solution: $u(x,y,z,t):=\\displaystyle{\\sum_{i=0}^{" << max_n_of_coefs-1 << "}} v_i(t) \\cdot \\mu_i(x) \\cdot \\nu_i(y) \\cdot \\omega_i(z)$\n\n";
  
    for(int i=0;i<max_n_of_coefs;++i)
      { 
  ofs << "$v_" << i << "(t):=e^{-\\varepsilon t \\pi^2 (" << common_coef_x[i] << "^2 + " << common_coef_y[i] << "^2 )} $\n\n"
    << "$\\mu_" << i << "(x):=" << cos_coef_x[i] << "\\cdot \\cos(" << common_coef_x[i] << "\\pi (x-at)) + " << sin_coef_x[i] 
    << "\\cdot \\sin(" << common_coef_x[i] << "\\pi (x-at)) $";
  if(dimDomain > 1)
    {
      ofs << "\n\n"
    << "$\\nu_" << i << "(y):=" << cos_coef_y[i] << "\\cdot \\cos(" << common_coef_y[i] << "\\pi (y-at)) + " << sin_coef_y[i] 
          << "\\cdot \\sin(" << common_coef_y[i] << "\\pi (y-at)) $";
    }
  if(dimDomain >2)
    {
      ofs << "\n\n"
    << "$\\omega_" << i << "(z):=" << cos_coef_z[i] << "\\cdot \\cos(" << common_coef_z[i] << "\\pi (z-at)) + " << sin_coef_z[i] 
          << "\\cdot \\sin(" << common_coef_z[i] << "\\pi (z-at)) $";
    }
  ofs << "\n\n";
      }
    ofs << "\n\n";
  
    ofs.close();
      
  }
    
  DomainType velocity_;
  double epsilon;
	std::string myName;
  int max_n_of_coefs;
  double common_coef_x[2];
  double sin_coef_x[2];
  double cos_coef_x[2];
  double common_coef_y[2];
  double sin_coef_y[2];
  double cos_coef_y[2];
  double common_coef_z[2];
  double sin_coef_z[2];
  double cos_coef_z[2];
};
template <class GridType>
class U0Disc : public U0<GridType> {
  double startTime_;
 public:
  enum { ConstantVelocity = true };
  typedef U0<GridType> BaseType;
  enum { dimDomain = GridType::dimensionworld };  
  typedef FieldVector<double,dimDomain> DomainType;
  typedef FieldVector<double,1> RangeType;
  U0Disc() : BaseType(),
	     startTime_(Parameter::getValue<double>("fem.localdg.starttime",0.0)) { 
    this->myName = "Discontinuous AdvectDiff";
    this->velocity_[0] = 1.0;
    this->velocity_[1] = 0.1;
  }

  double endtime() {
    return 0.5;
  }
    
  void evaluate(const DomainType& arg, RangeType& res) const {
    evaluate(startTime_,arg,res);
  }
	       
  void evaluate(const DomainType& arg, double t, RangeType& res) const {
    evaluate(t, arg, res);
    return;
  }

  void evaluate(double t,const DomainType& arg, RangeType& res) const 
  {
    //BaseType::evaluate(t,arg,res); 
    res = 1.0;
    DomainType x(this->velocity_);
    x *= -t;
    x += arg-DomainType(0.2);
    if (x[0]*x[1]>0.)
      res = 2.0;
  }
  
  void printmyInfo(std::string filename)
  {
    std::ostringstream filestream;
    filestream << filename;

    std::ofstream ofs(filestream.str().c_str(), std::ios::app);
  
    ofs << "Problem: " << this->myName << "\n\n"
  << "Exact solution: $u(x,y,z,t):=\\displaystyle{{\\rm sign}(x-\\frac{1}{2}-" << this->velocity_[0] << "t)\\sum_{i=0}^{" << this->max_n_of_coefs-1 << "}} \\mu_i(x) \\cdot \\nu_i(y) \\cdot \\omega_i(z)$\n\n";
  
    for(int i=0;i<this->max_n_of_coefs;++i)
      { 
  ofs << "$\\mu_" << i << "(x):=" << this->cos_coef_x[i] << "\\cdot \\cos(" << this->common_coef_x[i] << "\\pi (x-at)) + " << this->sin_coef_x[i] 
    << "\\cdot \\sin(" << this->common_coef_x[i] << "\\pi (x-at)) $";
  if(dimDomain > 1)
    {
      ofs << "\n\n"
    << "$\\nu_" << i << "(y):=" << this->cos_coef_y[i] << "\\cdot \\cos(" << this->common_coef_y[i] << "\\pi (y-at)) + " << this->sin_coef_y[i] 
          << "\\cdot \\sin(" << this->common_coef_y[i] << "\\pi (y-at)) $";
    }
  if(dimDomain >2)
    {
      ofs << "\n\n"
    << "$\\omega_" << i << "(z):=" << this->cos_coef_z[i] << "\\cdot \\cos(" << this->common_coef_z[i] << "\\pi (z-at)) + " << this->sin_coef_z[i] 
          << "\\cdot \\sin(" << this->common_coef_z[i] << "\\pi (z-at)) $";
    }
  ofs << "\n\n";
      }
    ofs << "\n\n";
  
    ofs.close();
      
  }
};
