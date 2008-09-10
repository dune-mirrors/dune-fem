# don't follow the full GNU-standard; we definitely need automake >= 1.5
AUTOMAKE_OPTIONS = foreign 1.5

SUBDIRS = common 0order 1order 2order matrix model projection 

operatordir =  $(includedir)/dune/fem/operator
operator_HEADERS = artificialdiffusion.hh \
		   cachedlinearoperator.hh \
		   discreteoperatorimp.hh \
		   diffusionoperator.hh \
                   dirichletboundaryoperator.hh \
		   elementintegrators.hh \
		   feop.hh \
		   identityoperator.hh \
		   integrationoperator.hh \
		   inverseoperators.hh \
		   lagrangeinterpolation.hh \
		   lagrangeinterpolation_inline.hh \
		   matrixadapter.hh \
		   minimizationoperator.hh \
		   sourceprojection.hh

include $(top_srcdir)/am/global-rules

