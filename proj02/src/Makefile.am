bin_PROGRAMS = heat_equation_solver

AM_CXXFLAGS = $(MASA_CXXFLAGS) $(GRVY_CFLAGS) $(HDF5_CFLAGS)
LIBS = $(MASA_LIBS) $(GRVY_LIBS) $(HDF5_LIBS)

#if PETSC_ENABLED
#    include $(PETSC_DIR)/lib/petsc/conf/variables
#    AM_CXXFLAGS += -DINCLUDE_PETSC $(PETSC_CC_INCLUDES)
#endif
#
#if PETSC_ENABLED
#    heat_equation_solver_LDADD = ${PETSC_LIB}
#endif

heat_equation_solver_SOURCES = main.cpp
