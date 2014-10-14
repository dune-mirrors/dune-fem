
#make libdunefem known locally
set(LOCAL_LIBS "${PROJECT_BINARY_DIR}/lib/libdunefem.a"
		CACHE STRING "path to local libs in dune-fem" )

#find endian headers
find_package(Endian)
add_definitions("-DSYSTEM_ENDIAN_HEADER=${SYSTEM_ENDIAN_HEADER}")

message(AUTHOR_WARNING "TODO. Improve module test.")

