# Files
EXEC := integrate
SRC := $(wildcard *.cpp)
OBJ := $(patsubst %.cpp,%.o,$(SRC))

# Options
CC := g++
CFLAGS := -I$(TACC_GSL_INC)
LDLIBS := -lgsl -lgslcblas
LDFLAGS := -L$(TACC_GSL_LIB)

# Rules
$(EXEC): $(OBJ)
	$(CC) $(LDFLAGS) $(LDLIBS) -o $@ $^
%.o: %.cpp
	$(CC) $(CFLAGS) -c $<

# Phony targets
.PHONY: clobber clean

clobber: clean
	$(RM) $(EXEC)

clean:
	$(RM) $(OBJ)
