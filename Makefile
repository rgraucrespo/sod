## This makefile must be executed with gmake (gnu make).

f90comp = gfortran
PREFIX   = /usr/local
UNAME := $(shell uname)
ifeq ($(UNAME), Darwin)
  LDFLAGS = -isysroot $(shell xcrun --show-sdk-path)
else
  LDFLAGS =
endif
FFLAGS = -Jbin -ffree-line-length-none

.PHONY: all clean test install docs

all:
	$(f90comp) $(FFLAGS) $(LDFLAGS) -o bin/combsod src/insod_reader.f90 src/ensemble_io.f90 src/factorials.f90 src/bubble.f90  src/ksubset.f90  src/member.f90  src/cell.f90 src/ccf.f90 src/combsod.f90
	$(f90comp) $(FFLAGS) $(LDFLAGS) -o bin/genersod src/insod_reader.f90 src/ensemble_io.f90 src/member.f90 src/cell.f90 src/genersod.f90
	$(f90comp) $(FFLAGS) $(LDFLAGS) -o bin/pmesod src/insod_reader.f90 src/ensemble_io.f90 src/structwriters.f90 src/pmemod.f90 src/pmesod.f90
	$(f90comp) $(FFLAGS) $(LDFLAGS) -o bin/mcsod src/insod_reader.f90 src/ensemble_io.f90 src/structwriters.f90 src/pmemod.f90 src/mcsod.f90
	$(f90comp) $(FFLAGS) $(LDFLAGS) -o bin/mcstatsod src/insod_reader.f90 src/ensemble_io.f90 src/structwriters.f90 src/pmemod.f90 src/mcstatsod.f90
	$(f90comp) $(FFLAGS) $(LDFLAGS) -o bin/invertENSEMBLE src/invertENSEMBLE.f90
	$(f90comp) $(FFLAGS) $(LDFLAGS) -o bin/statsod  src/ensemble_io.f90 src/statsod.f90
	$(f90comp) $(FFLAGS) $(LDFLAGS) -o bin/gcstatsod  src/ensemble_io.f90 src/factorials.f90 src/momenta.f90 src/gcstatsod.f90
	$(f90comp) $(FFLAGS) $(LDFLAGS) -o bin/peaks2spec  src/peaks2spec.f90
	$(f90comp) $(FFLAGS) $(LDFLAGS) -o bin/sqssod  src/ksubset.f90 src/cell.f90 src/sqssod.f90
	$(f90comp) $(FFLAGS) $(LDFLAGS) -o bin/gqssod src/ensemble_io.f90 src/ksubset.f90 src/cell.f90 src/gqssod.f90
	rm -f *.o

clean:
	rm -f bin/combsod bin/invertENSEMBLE bin/statsod bin/gcstatsod bin/genersod bin/pmesod bin/mcsod bin/mcstatsod bin/peaks2spec bin/sqssod bin/gqssod
	rm -f *.o bin/*.mod

test:
	bin/sod_run_tests.sh

install:
	install -d $(PREFIX)/bin
	install -m 755 bin/combsod bin/genersod bin/statsod bin/gcstatsod \
	    bin/invertENSEMBLE bin/pmesod bin/mcsod bin/mcstatsod bin/sqssod bin/gqssod \
	    bin/peaks2spec $(PREFIX)/bin
	install -m 755 bin/*.sh $(PREFIX)/bin

docs:
	$(MAKE) -C docs html
