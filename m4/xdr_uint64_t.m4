# searches for xdr_uint64_t or xdr_u_int64_t 

AC_DEFUN([DUNE_PATH_XDR_UINT64_T],[
  AC_REQUIRE([DUNE_PATH_XDR])
  AC_CHECK_FUNC(xdr_uint64_t,
                [AC_DEFINE(XDR_UINT64_FUNC, [xdr_uint64_t], [xdr_uint64_t routine])],
                [AC_CHECK_FUNC(xdr_u_int64_t,
                               [AC_DEFINE(XDR_UINT64_FUNC, [xdr_u_int64_t],[xdr_uint64_t routine])])
                ]
               )
])
