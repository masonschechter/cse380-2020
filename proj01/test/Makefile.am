AM_CPPFLAGS = $(MASA_CXXFLAGS) $(GRVY_CFLAGS)
LIBS = $(MASA_LIBS) $(GRVY_LIBS)

TESTS = 1D_linear_system_check
check_PROGRAMS = 1D_linear_system_check
1D_linear_system_check_SOURCES = 1D_check.cpp

TESTS += 2D_linear_system_check
check_PROGRAMS += 2D_linear_system_check
2D_linear_system_check_SOURCES = 2D_check.cpp

TESTS += jacobi_test_1D
check_PROGRAMS += jacobi_test_1D
jacobi_test_1D_SOURCES = jacobi_test_1D.cpp

TESTS += jacobi_test_2D
check_PROGRAMS += jacobi_test_2D
jacobi_test_2D_SOURCES = jacobi_test_2D.cpp

TESTS += gauss_seidel_test_1D
check_PROGRAMS += gauss_seidel_test_1D
gauss_seidel_test_1D_SOURCES = gauss_seidel_test_1D.cpp

TESTS += gauss_seidel_test_2D
check_PROGRAMS += gauss_seidel_test_2D
gauss_seidel_test_2D_SOURCES = gauss_seidel_test_2D.cpp

TESTS += tests.sh
