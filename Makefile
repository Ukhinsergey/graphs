CC=g++
CXX=g++
MPICXX=mpicxx
CXXFLAGS=-O3 -Wall -std=c++0x
LDFLAGS=-lm -lrt 

all:  compare gen_RMAT graphs_reference  gen_random  sssp_shmem

gen_RMAT_serial: gen_RMAT.o
	$(CXX) $^ -o $@ $(LDFLAGS)

gen_random_serial: gen_random.o
	$(CXX) $^ -o $@ $(LDFLAGS)

gen_RMAT_mpi: gen_RMAT_mpi.cpp
	$(MPICXX) $^ -o $@ -lrt


graphs_reference: main.o sssp_serial_reference.o
	$(CXX) $^ -o $@ $(LDFLAGS)

sssp_shmem: sssp_shmem.cpp sssp_serial_reference.cpp
	acxx $(CXXFLAGS) -o sssp_shmem sssp_shmem.cpp sssp_serial_reference.cpp
# write your own file
#sssp: main_mpi.o <file>.o
#	$(CXX) $^ -o $@ $(LDFLAGS)

compare: compare.o
	$(CXX) $^ -o $@ 

.cpp.o:
	$(CXX) $(CXXFLAGS) -o $@ -c $<

clean:
	rm -rf *.o gen_RMAT gen_RMAT_mpi gen_random gen_random_mpi graphs_reference compare

