# AC_REQUIRE_LIB(LIBRARY, FUNCTIONS)
# ------------------------------------------------------
#
# This demands that the library be present with the given functions compilable.
# A friendly error message is displayed and configure halts if these conditions
# are not met.  Note that unlike AC_CHECK_LIB, you can supply a list of 
# functions to the second argument of this macro, not just one.
AC_DEFUN([AC_REQUIRE_LIB], [
  AC_CHECK_LIB([$1],[main],[ 
    AH_TEMPLATE(AS_TR_CPP([HAVE_LIB$1]),
	     [Define to 1 if you have the `$1' library (-l$1).])
    AC_DEFINE_UNQUOTED(AS_TR_CPP(HAVE_LIB$1))
    LIBS="-l$1 $LIBS"

    for ac_func in $2
    do
    AC_CHECK_FUNC($ac_func,
            [AC_DEFINE_UNQUOTED([AS_TR_CPP([HAVE_$ac_func])])],
            [AC_MSG_FAILURE([
You appear to be missing the function $ac_func expected to be in library 
$1.  Perhaps you have an outdated version of this library, or the library 
is not properly compiled.])
    ])dnl
    done
   ],[
    AC_MSG_FAILURE([
You appear to be missing the library $1.  It may not be in your
normal library path.  You can fix this by
  a. Running configure with CFLAGS='-L/path/to/lib' ./configure [options]
  b. Setting LD_LIBRARY_PATH or DYLD_LIBRARY_PATH to include /path/to/lib])
  ])
]) # AC_REQUIRE_LIB
