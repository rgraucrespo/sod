## This makefile must be executed with gmake (gnu make).

f90comp = gfortran
UNAME := $(shell uname)
ifeq ($(UNAME), Darwin)
  LDFLAGS = -isysroot $(shell xcrun --show-sdk-path)
else
  LDFLAGS =
endif

all:
	$(f90comp) $(LDFLAGS) -o bin/combsod src/factorials.f90 src/bubble.f90  src/ksubset.f90  src/member.f90  src/cell.f90 src/ccf.f90 src/combsod.f90
	$(f90comp) $(LDFLAGS) -o bin/genersod src/member.f90 src/cell.f90 src/genersod.f90
	$(f90comp) $(LDFLAGS) -o bin/spbesod src/spbesod.f90
	$(f90comp) $(LDFLAGS) -o bin/invertOUTSOD src/invertOUTSOD.f90
	$(f90comp) $(LDFLAGS) -o bin/statsod  src/statsod.f90
	$(f90comp) $(LDFLAGS) -o bin/gcstatsod  src/factorials.f90 src/momenta.f90 src/gcstatsod.f90
	$(f90comp) $(LDFLAGS) -o bin/peaks2spec  src/peaks2spec.f90

clean:
	rm -f bin/combsod bin/invertOUTSOD bin/statsod bin/gcstatsod bin/genersod bin/spbesod bin/peaks2spec
	rm -f src/*.o src/*.mod


