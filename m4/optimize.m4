AC_DEFUN([AC_C_EXTERN_INLINE],
[
AC_CACHE_CHECK([for extern inline],
[ac_cv_c_extern_inline],
[ac_cv_c_extern_inline=no
for ac_kw in "extern inline" "extern __inline__" "extern __inline"; do
  AC_COMPILE_IFELSE([AC_LANG_SOURCE([
  #define inline inline
  #define __inline __inline
  #define __inline__ __inline__
  #ifndef __cplusplus
  typedef int foo_t;
  $ac_kw foo_t foo () { return 0; }
  #endif
  ])],
    [ac_cv_c_extern_inline="$ac_kw"; break])
done
]
AH_VERBATIM([__EXTERN_INLINE__],
[/* Define to `__inline__' or `__inline' if that's what the C compiler
 calls it, or to nothing if 'inline' is not supported under any name.  */
#ifndef __cplusplus
#undef __EXTERN_INLINE__
#endif])
dnl block 2
    [
    case "$ac_cv_c_extern_inline" in
      "extern inline" | yes) ac_val="$ac_cv_c_extern_inline";;
      *)
      case $ac_cv_c_extern_inline in
           no) ac_val=;;
            *) ac_val="$ac_cv_c_extern_inline";;
      esac
    esac
cat >>confdefs.h <<_ACEOF
#ifndef __cplusplus
#define __EXTERN_INLINE__ $ac_val
#endif
_ACEOF
]
)])

AC_DEFUN([AC_PROG_CC_BIGOPTS],
[ac_test_CFLAGS=${CFLAGS+set}
ac_save_CFLAGS=$CFLAGS
for ac_opt in O3 xO5; do
  CFLAGS="-$ac_opt"
  AC_CACHE_CHECK([whether $CC accepts -$ac_opt], ac_cv_prog_cc_bigopts,
  [AC_COMPILE_IFELSE([AC_LANG_PROGRAM()], [
   AC_LINK_IFELSE([AC_LANG_PROGRAM()],
    [ac_cv_prog_cc_bigopts=$ac_opt],
    [])
   ],
    [])])
  if test "x$ac_cv_prog_cc_bigopts" = "x$ac_opt"; then
    CFLAGS="`echo $ac_save_CFLAGS | sed -e 's/O[0-9]*/$ac_opt/'` -$ac_opt"
    break
  fi
done
])
