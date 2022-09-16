# For compiling setting, Makefile will point to the uncommented file below:
include ../src/Makefile_include_local.txt
#include Makefile_include_mac.txt
#include Makefile_include_snellius.txt

OBJ_PATH = test_obj
MOD_PATH = mod
UFEMISM_OBJ_PATH = obj

UNIT_TESTS=$(wildcard test_*.f90)
UNIT_TESTS_OBJ_1=$(UNIT_TESTS:%=$(OBJ_PATH)/%)
UNIT_TESTS_OBJ=$(UNIT_TESTS_OBJ_1:.f90=.o)

$(OBJ_PATH)/%.o : %.f90 Makefile UFEMISM.a
	$(F90) $(F90FLAGS) -J${MOD_PATH} -fbackslash -c -o $@ $<


UFEMISM_unit_tests: ${OBJ_PATH}/main.o $(UNIT_TESTS_OBJ) UFEMISM.a
	$(F90) $(F90FLAGS) -fprofile-arcs -ftest-coverage -o $@ $^ $(LDFLAGS)

main.f90: create_test_main.sh ${UNIT_TESTS} Makefile
	./create_test_main.sh

${OBJ_PATH}/main.o: $(UNIT_TESTS_OBJ)

# Requested, but not used, just to signal completion of compilation
# Always done, to make it up to date
.PHONY: UFEMISM.a
UFEMISM.a:
	cd ../src && $(MAKE) unit_tests_lib

clean:
	rm -fv */*.o */*.gcda */*.gcno */*.mod UFEMISM.a main.f90 UFEMISM_unit_tests coverage.xml html/* coverage.csv

.PHONY: run_tests
run: UFEMISM_unit_tests
	./UFEMISM_unit_tests

coverage_html:
	gcovr -r ../src --object-directory obj --html-details -o html/index.html

coverage:
	gcovr -r ../src --object-directory obj --csv=coverage.csv -s
