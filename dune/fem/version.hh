#include <dune/common/version.hh>

// include DUNE_DEPRECATED header 
#include <dune/common/deprecated.hh>

#define DUNE_VERSION_DEPRECATED_1_4 DUNE_DEPRECATED 
#define DUNE_VERSION_DEPRECATED_1_5 DUNE_DEPRECATED
#define DUNE_VERSION_DEPRECATED_1_6 DUNE_DEPRECATED
#define DUNE_VERSION_DEPRECATED_1_7 DUNE_DEPRECATED
#define DUNE_VERSION_DEPRECATED_1_8 DUNE_DEPRECATED
#define DUNE_VERSION_DEPRECATED_1_9 DUNE_DEPRECATED

#define DUNE_REMOVE_DEPREACTED dune_remove_deprecated

// current version should fail 
#define DUNE_VERSION_DEPRECATED_1_3 DUNE_REMOVE_DEPREACTED

#define DUNE_VERSION_DEPRECATED(major,minor) \
 DUNE_VERSION_DEPRECATED_##major##_##minor
