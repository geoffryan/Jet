MAKEFILE_MACH = $(PWD)/Makefile.mach
MAKEFILE_IN = $(PWD)/Makefile.in

include $(MAKEFILE_MACH)
include $(MAKEFILE_IN)

FLAGS = -O3 -Wall -g

INC = -I$(H55)/include
LIB = -L$(H55)/lib -lm -lhdf5

OBJ = main.o readpar.o onestep.o exchange.o plm.o domain.o faces.o cooling.o nozzle.o gravity.o misc.o mpisetup.o gridsetup.o $(RIEMANN).o $(TIMESTEP).o $(INITIAL).o $(HYDRO).o $(GEOMETRY).o $(BOUNDARY).o $(OUTPUT).o snapshot.o report.o $(RESTART).o

default: jet

%.o: %.c paul.h
	$(CC) $(FLAGS) $(LOCAL_CFLAGS) $(INC) -c $<

$(RIEMANN).o: Riemann/$(RIEMANN).c paul.h
	$(CC) $(FLAGS) $(LOCAL_CFLAGS) $(INC) -c Riemann/$(RIEMANN).c

$(TIMESTEP).o: Timestep/$(TIMESTEP).c paul.h
	$(CC) $(FLAGS) $(LOCAL_CFLAGS) $(INC) -c Timestep/$(TIMESTEP).c

$(INITIAL).o : Initial/$(INITIAL).c paul.h
	$(CC) $(FLAGS) $(LOCAL_CFLAGS) $(INC) -c Initial/$(INITIAL).c

$(HYDRO).o : Hydro/$(HYDRO).c paul.h
	$(CC) $(FLAGS) $(LOCAL_CFLAGS) $(INC) -c Hydro/$(HYDRO).c

$(GEOMETRY).o : Geometry/$(GEOMETRY).c paul.h
	$(CC) $(FLAGS) $(LOCAL_CFLAGS) $(INC) -c Geometry/$(GEOMETRY).c

$(BOUNDARY).o : Boundary/$(BOUNDARY).c paul.h
	$(CC) $(FLAGS) $(LOCAL_CFLAGS) $(INC) -c Boundary/$(BOUNDARY).c

$(OUTPUT).o : Output/$(OUTPUT).c paul.h
	$(CC) $(FLAGS) $(LOCAL_CFLAGS) $(INC) -c Output/$(OUTPUT).c

$(RESTART).o : Restart/$(RESTART).c paul.h
	$(CC) $(FLAGS) $(LOCAL_CFLAGS) $(INC) -c Restart/$(RESTART).c

jet: $(OBJ) paul.h
	$(CC) $(FLAGS) $(LOCAL_LDFLAGS) $(LIB) -o jet $(OBJ)

clean:
	rm -f *.o jet
