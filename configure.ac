#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.64])
AC_INIT([east], [0.1], [robert.zimmermann@univie.ac.at])
AC_CONFIG_SRCDIR([src/matrix.h])
AM_INIT_AUTOMAKE([foreign])
AC_CONFIG_HEADERS([config.h])

# Checks for programs.
AC_PROG_CC_C99

# Checks for libraries.
AC_CHECK_LIB(m,sqrt,[],[AC_MSG_ERROR("math lib required")])
AC_CHECK_LIB(m,exp,[],[AC_MSG_ERROR("exp() function not found in the math library")])
AC_CHECK_LIB(m,log,[],[AC_MSG_ERROR("log() function not found in the math library")])


# Checks for header files.
AC_CHECK_HEADERS([stdlib.h string.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_INLINE
AC_TYPE_SIZE_T

# Checks for user options
AC_CHECK_DEBUG

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_CHECK_FUNCS([floor strtol])

AC_CONFIG_FILES([ \
                 Makefile
                 src/Makefile
                 ])
AC_OUTPUT
