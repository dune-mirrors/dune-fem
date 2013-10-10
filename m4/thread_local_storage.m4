AC_DEFUN([THREAD_LOCAL_STORAGE],[
  AC_REQUIRE([ACX_PTHREAD])

  AC_MSG_CHECKING(for working __thread)
  # store old values
  ac_save_LDFLAGS="$LDFLAGS"
  ac_save_CPPFLAGS="$CPPFLAGS"
  ac_save_LIBS="$LIBS"

  CPPFLAGS="$CPPFLAGS $PTHREAD_CFLAGS"
  LDFLAGS="$LDFLAGS $PTHREAD_LDFLAGS"
  LIBS="$LIBS $PTHREAD_LIBS"

  AC_LANG_PUSH([C++])
  AC_TRY_RUN([#include <pthread.h>
              static const int maxThreads = 2 ;
              static int& threadNumber () {
                  static __thread int tn;
                  return tn;
              }
              static int threadnumber[ maxThreads ];
              pthread_barrier_t barrier;
              static void* setThreadNumber(void *number) {
                const int num = *((int*) number);
                threadNumber () = num ;
                //sleep (1);
                threadnumber[ num ] = threadNumber ();
                pthread_barrier_wait( &barrier );
                return 0;
              }
              int main () {
                pthread_t thread[ maxThreads ];
                int num[ maxThreads ]; 
                for (int i=0; i<maxThreads; ++i ) 
                {
                  num[ i ] = i;
                  threadnumber[ i ] = -1;  
                }
                pthread_barrier_init( &barrier, 0, maxThreads );
                for( int i=0; i<maxThreads; ++ i) 
                  pthread_create (&thread[ i ], 0, &setThreadNumber, (void *) &num[ i ]);

                for( int i=0; i<maxThreads; ++ i) 
                  pthread_join (thread[ i ], 0);
                int result = 0;
                for( int i=0; i<maxThreads; ++i )
                  if( threadnumber[ i ] != i ) result = 1;
                return result;
              } ],
              [AC_MSG_RESULT(yes)
               AC_DEFINE(HAVE_PTHREAD_TLS, 1,[This is true if the keyword __thread can be used])
               ],
              [AC_MSG_RESULT(no)])

  AC_LANG_POP

  # reset old values
  LIBS="$ac_save_LIBS"
  CPPFLAGS="$ac_save_CPPFLAGS"
  LDFLAGS="$ac_save_LDFLAGS"
])
