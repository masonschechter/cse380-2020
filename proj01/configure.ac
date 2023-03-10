#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([Project_1], [1.0], [masonschechter@utexas.edu])
AC_CONFIG_SRCDIR([src/])
AC_CONFIG_AUX_DIR([build-aux])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_MACRO_DIR([m4])
AM_INIT_AUTOMAKE


# Checks for programs.
AC_PROG_CXX

# Checks for libraries.
#AX_PATH_MASA(0.30, [],AC_MSG_ERROR([Unable to find a suitable version of MASA]))
#AX_PATH_GRVY(0.32, [],AC_MSG_ERROR([Unable to find a suitable version of GRVY]))
AX_PATH_MASA([0.30], [yes])
AX_PATH_GRVY([0.32], [yes])

# Checks for header files.

# Checks for typedefs, structures, and compiler characteristics.

# Checks for library functions.

#AC_CONFIG_FILES([Makefile 
#		src/Makefile])

AC_OUTPUT(Makefile src/Makefile test/Makefile)
