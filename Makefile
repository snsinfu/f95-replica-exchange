MPIFC = mpifort
MPIFFLAGS =

OBJECTS = main.o rex.o
ARTIFACTS = main rex.mod $(OBJECTS)

.PHONY: run clean
.SUFFIXES: .f95 .mod

run: main
	mpirun -np 3 ./main

clean:
	rm -f $(ARTIFACTS)

main: $(OBJECTS)
	$(MPIFC) $(MPIFFLAGS) -o $@ $(OBJECTS) $(LDFLAGS)

.f95.o:
	$(MPIFC) $(MPIFFLAGS) -c -o $@ $<

main.o: rex.o
