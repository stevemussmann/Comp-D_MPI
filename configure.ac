# Initialize AC: ACINIT( package_name, version, bug_report_address)
AC_INIT([compDmpi], [1.0], [mussmann@email.uark.edu])
AM_INIT_AUTOMAKE
AC_PROG_RANLIB
AC_LANG(C++)
AC_PROG_CXX([mpic++ mpicxx],,[AC_MSG_ERROR([mpi is required but was not found on your system])])
AC_CONFIG_FILES([Makefile src/Makefile])
AC_OUTPUT
AC_CONFIG_MACRO_DIR([m4])
AX_BOOST_BASE([1.53],,[AC_MSG_ERROR([Boost >= 1.53 required, but it was not found on your system])])
#AX_BOOST_SYSTEM
AX_BOOST_PROGRAM_OPTIONS
