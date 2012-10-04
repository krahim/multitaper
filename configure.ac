#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.61)
AC_INIT(FULL-PACKAGE-NAME, VERSION, BUG-REPORT-ADDRESS)
AC_CONFIG_SRCDIR([src/dpss1.c])
#AC_CONFIG_HEADER([config.h])

# Checks for programs.


: ${R_HOME=`R RHOME`}
if test -z "${R_HOME}"; then
   echo "could not determine R_HOME"
        exit 1
fi

CC=`"${R_HOME}/bin/R" CMD config CC`
AC_PROG_CC

# Checks for libraries.

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([stdlib.h math.h complex.h])

# Checks for typedefs, structures, and compiler characteristics.

# Checks for library functions.
AC_FUNC_MALLOC
##AC_CHECK_FUNCS([pow])
AC_CONFIG_FILES([src/Makevars])
AC_OUTPUT
