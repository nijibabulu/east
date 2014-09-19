# AC_CHECK_DOXYGEN([PATH_TO_LEAVE_DOCS])
# 
# Check for doxygen, and if so, substitute the variable WITH_APIDOCS_TARGET
# to the passed parameter, or $(top_builddir)/docs by default.
#
AC_DEFUN([AC_CHECK_DOXYGEN],
[
AC_PATH_PROG(__DOXYGEN, doxygen, no, $PATH)
withval=auto
AC_ARG_WITH(apidocs, [  --with-apidocs          build API docs ])

if test $withval = auto -a $__DOXYGEN != no ; then
  withval=yes
  WITH_APIDOCS_TARGET='$(top_builddir)/docs'
  if test  "x$1" != "x"; then
    WITH_APIDOCS_TARGET=$1
  fi
  AC_SUBST(WITH_APIDOCS_TARGET)
elif test $withval = yes -a $__DOXYGEN = no ; then
  AC_MSG_ERROR(--> docs needs doxygen in PATH)
fi
])
