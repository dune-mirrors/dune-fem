#include <dune/common/version.hh>

// include DUNE_DEPRECATED header 
#include <dune/common/deprecated.hh>

// for all future versions use the usual DUNE_DEPRECATED
#define DUNE_VERSION_DEPRECATED_1_3(m) DUNE_DEPRECATED
#define DUNE_VERSION_DEPRECATED_1_4(m) DUNE_DEPRECATED
#define DUNE_VERSION_DEPRECATED_1_5(m) DUNE_DEPRECATED
#define DUNE_VERSION_DEPRECATED_1_6(m) DUNE_DEPRECATED
#define DUNE_VERSION_DEPRECATED_1_7(m) DUNE_DEPRECATED
#define DUNE_VERSION_DEPRECATED_1_8(m) DUNE_DEPRECATED
#define DUNE_VERSION_DEPRECATED_1_9(m) DUNE_DEPRECATED

// current version should fail 
#define DUNE_VERSION_DEPRECATED_1_2(newmethod) \
  dune_remove_deprecated_method__use_##newmethod##_instead

#define DUNE_VERSION_DEPRECATED(major,minor,newmethod) \
 DUNE_VERSION_DEPRECATED_##major##_##minor(newmethod)
