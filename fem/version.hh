#ifndef DUNE_VERSION_HH
#define DUNE_VERSION_HH

#ifndef DUNE_COMMON_VERSION_MAJOR
#include "version.inc"
#endif

#define DUNE_VERSION_JOIN(module,type) module##_VERSION_##type

#define DUNE_VERSION_EQUAL(module,major,minor) \
    ((DUNE_VERSION_JOIN(module,MAJOR) == major) && \
     (DUNE_VERSION_JOIN(module,MINOR) == minor))

#define DUNE_VERSION_EQUAL_REV(module,major,minor,revision) \
    ( DUNE_VERSION_EQUAL(module,major,minor) && \
     (DUNE_VERSION_JOIN(module,REVISION) == revision))

#define DUNE_VERSION_NEWER(module,major,minor,revision) \
  ((DUNE_VERSION_JOIN(module,MAJOR) > major) \
   || ((DUNE_VERSION_JOIN(module,MAJOR) == major) && (DUNE_VERSION_JOIN(module,MINOR) >= minor)))

#define DUNE_VERSION_NEWER_REV(module,major,minor,revision) \
  ((DUNE_VERSION_JOIN(module,MAJOR) > major) \
   || ((DUNE_VERSION_JOIN(module,MAJOR) == major) && (DUNE_VERSION_JOIN(module,MINOR) > minor)) \
   || ((DUNE_VERSION_JOIN(module,MAJOR) == major) && (DUNE_VERSION_JOIN(module,MINOR) == minor) \
       && (DUNE_VERSION_JOIN(module,REVISION) >= revision)))

#define DUNE_VERSION_ID(major,minor,revision) \
  ((unsigned int)((major << 24) + (minor << 16) + revision ))

#define DUNE_MODULE_VERSION_ID(module) \
  DUNE_VERSION_ID( DUNE_VERSION_JOIN(module,MAJOR), DUNE_VERSION_JOIN(module,MINOR), DUNE_VERSION_JOIN(module,REVISION) )

#endif
