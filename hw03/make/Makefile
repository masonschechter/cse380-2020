# Files
EXEC := tardis
SRC := $(wildcard *.c)
OBJ := $(patsubst %.c,%.o,$(SRC))

# Options
CC := icc
CFLAGS := -I$(TACC_GSL_INC)
LDLIBS := -lgsl -lgslcblas
LDFLAGS := -L$(TACC_GSL_LIB)

# Rules
$(EXEC): $(OBJ)
	$(CC) $(LDFLAGS) $(LDLIBS) -o $@ $^
%.o: %.c
	$(CC) $(CFLAGS) -c $<

# Phony targets
.PHONY: clobber clean

clobber: clean
	$(RM) $(EXEC)

clean:
	$(RM) $(OBJ)
