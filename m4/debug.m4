AC_DEFUN([AC_CHECK_DEBUG],
[
DEBUG_CFLAGS="-O0 -g -pg --pedantic" #-Wall 
DEBUG_LDFLAGS="-g -pg"
CFLAGS="${CFLAGS=}"
AC_MSG_CHECKING([if we're enabling debugging])
AC_ARG_ENABLE(debug, , [
if test "x$enable_debug" != "xno"; then
  CFLAGS="$CFLAGS $DEBUG_CFLAGS -DDEBUGGING=$enable_debug"
  LDFLAGS="$LDFLAGS $DEBUG_LDFLAGS"
  AC_MSG_RESULT([yes])
else
  dnl cc on solaris doesn't let you compile -g with -xO*; so 
  dnl disable this when not debugging
  case $host_os in
    *solaris*) CFLAGS=`echo $CFLAGS | sed -e 's/-g//g'`;;
    *) ;;
  esac
fi],[
  AC_MSG_RESULT([no])
])
])
