FC = gfortran-mp-8
EXTRA_COMPILE_FLAGS =-ffree-line-length-none -O0 -g
EXTRA_LINK_FLAGS =  
LIBSTELL_DIR=/Users/mattland/stellopt_github/develop/STELLOPT/LIBSTELL/Release

# End of system-dependent variable assignments

TARGET = scale_vmec

.PHONY: all clean

all: $(TARGET)

%.o: %.f90
	$(FC) $(EXTRA_COMPILE_FLAGS) -I $(LIBSTELL_DIR) -c $<

$(TARGET): $(TARGET).o
	$(FC) -o $(TARGET) $(TARGET).o $(LIBSTELL_DIR)/libstell.a $(EXTRA_LINK_FLAGS)

clean:
	rm -f *.o *.mod *.MOD *~ $(TARGET) *.a

test_make:
	@echo HOSTNAME is $(HOSTNAME)
	@echo FC is $(FC)
	@echo LIBSTELL_DIR is $(LIBSTELL_DIR)
	@echo EXTRA_COMPILE_FLAGS is $(EXTRA_COMPILE_FLAGS)
	@echo EXTRA_LINK_FLAGS is $(EXTRA_LINK_FLAGS)
