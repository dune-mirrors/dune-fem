#include <dune/common/version.hh>
#include <dune/common/deprecated.hh>

// for all future versions use the usual DUNE_DEPRECATED
#define DUNE_VERSION_DEPRECATED_3_0(m) DUNE_DEPRECATED

// for all future versions use the usual DUNE_DEPRECATED
#define DUNE_VERSION_DEPRECATED_2_7(m) DUNE_DEPRECATED

// for all future versions use the usual DUNE_DEPRECATED
#define DUNE_VERSION_DEPRECATED_2_6(m) DUNE_DEPRECATED

// current version should fail
#define DUNE_VERSION_DEPRECATED_2_5(newmethod) \
  dune_remove_deprecated_method__use_##newmethod##_instead

// current version should fail
#define DUNE_VERSION_DEPRECATED_2_4(newmethod) \
  dune_remove_deprecated_method__use_##newmethod##_instead

#define DUNE_VERSION_DEPRECATED(major,minor,newmethod) \
 DUNE_VERSION_DEPRECATED_##major##_##minor(newmethod)
