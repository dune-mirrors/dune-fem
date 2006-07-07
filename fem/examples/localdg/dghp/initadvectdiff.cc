template <class GridType>
class U0 {
 public:
  enum { dimDomain = GridType::dimensionworld };  
  typedef FieldVector<double,dimDomain> DomainType;
  typedef FieldVector<double,1> RangeType;
  U0(double eps,bool diff_timestep=true) :
    velocity_(0.1), epsilon(eps), diff_tstep(diff_timestep) {
      velocity_[0]=-1.5;
 						
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

  double endtime() {
    return 0.1;
  }
	 
  void velocity(double t, const DomainType &x, DomainType &res) const{
    res = velocity_;
  }
	
  void evaluate(const DomainType& arg, RangeType& res) const {
    evaluate(0,arg,res);
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
	
  void printmyInfo(string filename)
  {
    std::ostringstream filestream;
    filestream << filename;

    std::ofstream ofs(filestream.str().c_str(), std::ios::app);
	
    ofs << "Problem: " << myName << "\n\n"
	<< "Epsilon = " << epsilon << "\n\n"
	<< "Exact solution: $u(x,y,z,t):=\\displaystyle{\\sum_{i=0}^{" << max_n_of_coefs-1 << "}} v_i(t) \\cdot \\mu_i(x) \\cdot \\nu_i(y) \\cdot \\omega_i(z)$\n\n";
	
    for(int i=0;i<max_n_of_coefs;++i)
      { 
	ofs	<< "$v_" << i << "(t):=e^{-\\varepsilon t \\pi^2 (" << common_coef_x[i] << "^2 + " << common_coef_y[i] << "^2 )} $\n\n"
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
	ofs	<< "\n\n";
      }
    ofs	<< "\n\n";
	
    ofs.close();
			
  }
		
  DomainType velocity_;
  double epsilon;
  bool diff_tstep;
  string myName;
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
 public:
  typedef U0<GridType> BaseType;
  enum { dimDomain = GridType::dimensionworld };  
  typedef FieldVector<double,dimDomain> DomainType;
  typedef FieldVector<double,1> RangeType;
  U0Disc(double eps,bool diff_timestep=true) : BaseType(0.0,diff_timestep) {
    this->myName = "Discontinuous AdvectDiff";
  }

  double endtime() {
    return 0.1;
  }
      
	
  void evaluate(const DomainType& arg, RangeType& res) const {
    evaluate(0,arg,res);
  }
	
	
  void evaluate(double t,const DomainType& arg, RangeType& res) const {
    const double x0=0.5;
    const double x1=1.3177653851118265;
    BaseType::evaluate(t,arg,res); 
    if (arg[0]-this->velocity_[0]*t<x0 || arg[0]-this->velocity_[0]*t>x1)
      res *= -1.;
  }
	
  void printmyInfo(string filename)
  {
    std::ostringstream filestream;
    filestream << filename;

    std::ofstream ofs(filestream.str().c_str(), std::ios::app);
	
    ofs << "Problem: " << this->myName << "\n\n"
	<< "Exact solution: $u(x,y,z,t):=\\displaystyle{{\\rm sign}(x-\\frac{1}{2}-" << this->velocity_[0] << "t)\\sum_{i=0}^{" << this->max_n_of_coefs-1 << "}} \\mu_i(x) \\cdot \\nu_i(y) \\cdot \\omega_i(z)$\n\n";
	
    for(int i=0;i<this->max_n_of_coefs;++i)
      { 
	ofs	<< "$\\mu_" << i << "(x):=" << this->cos_coef_x[i] << "\\cdot \\cos(" << this->common_coef_x[i] << "\\pi (x-at)) + " << this->sin_coef_x[i] 
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
	ofs	<< "\n\n";
      }
    ofs	<< "\n\n";
	
    ofs.close();
			
  }
};



template <class GridType>
class U0RotCone {
 public:
  enum { dimDomain = GridType::dimensionworld };  
  typedef FieldVector<double,dimDomain> DomainType;
  typedef FieldVector<double,1> RangeType;
  U0RotCone(double eps,bool diff_timestep=true) :
    velocity_(1.0), epsilon(eps), diff_tstep(diff_timestep)  {
    center_ = 1.0;
    center_[0] -= 0.5;
    radius_ = 0.2;

    myName = "Rotating Cone";
  }

  double endtime() {
    return 3.14159;
  }
		
  double dist (const DomainType& x, const DomainType& y) const{
    double ret = 0.0;
    for(int i =0 ; i<dimDomain; ++i)
      ret += (x[i]-y[i])*(x[i]-y[i]);

    return(sqrt(ret));
  }

  void velocity(double t, const DomainType &x, DomainType &res) const{
    
    res = 0.0;

    if(dimDomain == 1)
      res[0] = sin(x[0]);
    else{
      res[0] = - (x[1]-1.0);
      res[1] =   (x[0]-1.0);
    }

  }

  void evaluate(const DomainType& arg, RangeType& res) const {
    
    double r = dist(arg, center_);

    res = 0.0;
    if ( r < radius_ ){
      //res = 1.0 - r/radius_;
	res = 1.0;
      //res = 1.0 - (r*r/(radius_*radius_));
    }
  }
	
	
  void evaluate(double t,const DomainType& arg, RangeType& res) const {
    DomainType x;
    DomainType vel;

    velocity(t,arg,vel);

    x = vel;
    x *= -1.0;
    x += arg;

    evaluate(x,res);
  }
	
  void printmyInfo(string filename)
  {
    std::ostringstream filestream;
    filestream << filename;

    std::ofstream ofs(filestream.str().c_str(), std::ios::app);
	
    ofs << "Problem: " << myName << "\n\n";
    ofs	<< "\n\n";
	
    ofs.close();
			
  }

  DomainType velocity_;
  double epsilon;
  bool diff_tstep;	
  string myName;
  DomainType center_;
  double radius_;
};


template <class GridType>
class U0BuckLev {
 public:
  enum { dimDomain = GridType::dimensionworld };  
  typedef FieldVector<double,dimDomain> DomainType;
  typedef FieldVector<double,1> RangeType;
  U0BuckLev(double eps,bool diff_timestep=true) :
    velocity_(0.0), epsilon(eps), diff_tstep(diff_timestep)  {
    velocity_[0] = 1.0;

    myName = "Buckley Leverett";
  }

  double endtime() {
    return 0.4;
  }

  double f(double u) const {
    if (u<0.0)
      return 0.;
    if (u>1.0)
      return 1.; 
    return u*u/(u*u+0.5*(1.-u)*(1.-u));
  }

  double f1(double u) const {
    if (u<0.0)
      return 0.;
    if (u>1.0)
      return 0.; 
    double d=3.*u*u+1.-2.*u;
    return -4.*u*(-1.+u)/d/d;
  }

  double f2(double u) const {
    double d=3.*u*u+1.-2.*u;
    return 4.*(-9.*u*u+6.*u*u*u+1.)/d/d/d;
  }

  double Newton(double u0,double xi) const {
    double u=u0;
    while (fabs(f1(u)-xi)>1e-8) {
      double u1=u;
      double df2=f2(u1);
      assert(fabs(df2)>1e-10);
      u=u1-(f1(u1)-xi)/df2;
      // cerr << u1 << " " << u << " " << f1(u)-xi << endl;
    }
    return u;
  }

  double rp_sol(double ul,double ur,double t,double x) const {
    double u1=(3.-sqrt(6.))/3.;
    double u2=sqrt(3.)/3.;
    double xi=x/t;
    if (ul<ur) {
      if (ul>u1) {
	double fl=f(ul);
	double fr=f(ur);
	double s=(fl-fr)/(ul-ur);
	if (xi<s) return ul;
	else return ur;
      } else if (ur<u1) {
	double dfl=f1(ul);
	double dfr=f1(ur);
	if (xi<dfl)
	  return ul;
	else if (xi>dfr) 
	  return ur;
	else 
	  return Newton(0.5*(ul+ur),xi);
      } else {
	double dfl=f1(ul);
	double df1=f1(u1);
	assert(dfl<df1);
	if (xi<dfl)
	  return ul;
	else if (xi>df1) {
	  double f1=f(u1);
	  double fr=f(ur);
	  double s=(f1-fr)/(u1-ur);
	  if (xi<s) return u1;
	  else return ur;
	}
	else 
	  return Newton(0.5*(ul+u1),xi);      
      }
    }
    else if (ul>ur) {
      if (ul<u2) {
	double fl=f(ul);
	double fr=f(ur);
	double s=(fl-fr)/(ul-ur);
	if (xi<s) return ul;
	else return ur;
      } else if (ur>u2) {
	double dfl=f1(ul);
	double dfr=f1(ur);
	assert(dfl<dfr);
	if (xi<dfl)
	  return ul;
	else if (xi>dfr) 
	  return ur;
	else
	  return Newton(0.5*(ul+ur),xi);
      } else {
	double dfl=f1(ul);
	double df2=f1(u2);
	assert(dfl<df2);
	if (xi<dfl)
	  return ul;
	else if (xi>df2) {
	  double f2=f(u2);
	  double fr=f(ur);
	  double s=(f2-fr)/(u2-ur);
	  if (xi<s) return u2;
	  else return ur;
	}
	else 
	  return Newton(0.5*(ul+u2),xi);      
      }
    } 
    else return ul;
  }


  void velocity(double t, const DomainType &x, DomainType &res) const{
    res = 0.0;
    res[0] = 1.0;
  }

  void evaluate(const DomainType& arg, RangeType& res) const {
    evaluate(0.0,arg,res); 
  }
	
	
  void evaluate(double t,const DomainType& arg, RangeType& res) const {
    double x = arg[0];
    if (x>1.1)
      res[0]=rp_sol(0.,1.,t,x-1.2);
    else
      res[0]=rp_sol(1.,0.,t,x-0.4);
  }
	
  void printmyInfo(string filename)
  {
    std::ostringstream filestream;
    filestream << filename;

    std::ofstream ofs(filestream.str().c_str(), std::ios::app);
	
    ofs << "Problem: " << myName << "\n\n";
    ofs	<< "\n\n";
	
    ofs.close();
			
  }

  DomainType velocity_;
  double epsilon;
  bool diff_tstep;	
  string myName;
};

