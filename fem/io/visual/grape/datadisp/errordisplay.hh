template <class ConsType,class InitialDataType>
class DisplayErrorFunction {
  struct Error {
    typedef typename ConsType::DiscreteFunctionSpaceType::GridPartType 
      GridPartType;
    typedef typename GridPartType :: GridType :: 
      template Codim<0> :: Entity EntityType;
    typedef typename EntityType :: Geometry GeometryImp;
    typedef typename ConsType :: DiscreteFunctionSpaceType :: FunctionSpaceType
                     ConsFunctionSpaceType;
    typedef typename ConsFunctionSpaceType :: DomainType ConsDomainType;
    typedef typename ConsFunctionSpaceType :: RangeType ConsRangeType;
    typedef FunctionSpace<double,double,
			  ConsDomainType::dimension,ConsRangeType::dimension>
			  FunctionSpaceType;
    typedef typename FunctionSpaceType :: DomainType DomainType;
    typedef typename FunctionSpaceType :: RangeType RangeType;
    typedef typename FunctionSpaceType :: JacobianRangeType JacobianRangeType;
    Error(const ConsType& Uh,double time) : 
      lUh_(Uh), 
      initU0_(),
      time_(time),
      geometry_(0),
      initialized_(false) {
    }
    ~Error() {
    }
    template< class PointType >
    void evaluate ( const PointType &x, RangeType& ret) const {
      assert(initialized_);
      ConsRangeType u;
      lUh_.evaluate(x,ret);
      ConsRangeType u0;
      DomainType global = geometry_->global( coordinate( x ) );
      initU0_.evaluate(time_,global,u0);
      ret -= u0;
    }
    template< class PointType >
    void jacobian ( const PointType &x, JacobianRangeType& ret) const {
      abort();
    }
    //! init local function
    void init(const EntityType& en)
    {
      lUh_.init(en);
      geometry_ = &(en.geometry());
      initialized_ = true;
    }
  private:
    typedef typename ConsType::LocalFunctionType ConsLocalFunctionType;
    ConsLocalFunctionType lUh_;
    InitialDataType initU0_;
    double time_;
    const GeometryImp* geometry_;
    bool initialized_;
  };
public:
  // Add the error function to the display 
  template <class GrapeDispType, 
          class GR_GridType,
          class DestinationType>
  static void apply(GrapeDispType& disp,
		    const GR_GridType& grid,
		    const double time,
		    const DestinationType& Uh)
  {
    typedef LocalFunctionAdapter<Error> ErrorFunction;
    Error* evalError = new Error(Uh,time);
    ErrorFunction* error = new 
      ErrorFunction("error",*evalError,Uh.space().gridPart()); 
    disp.addData(*error,"error",time);  
  }
};
