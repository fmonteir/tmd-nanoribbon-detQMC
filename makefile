# tmd-nanoribbon-detQMC
#
# Determinant Quantum Monte Carlo for a minimal multi-orbital Hubbard model
# of a TMD nanoribbon
#
# Created by Francisco Brito (6 May 2019)
#
#	The program is based on
#
# 	"Numerical Methods for
# 	Quantum Monte Carlo Simulations of the Hubbard Model" by Zhaojun Bai,
# 	Wenbin Chen, Richard Scalettar, and Ichitaro Yamazaki (2009)
#
#		"Stable solutions of linear systems involving
#   long chain of matrix multiplications" by Zhaojun Bai, Che-Rung Lee,
#		Ren-Cang Li and Shufang Xu (2010)
#
#   Liu et al., Phys Rev B 88, 085433, 2013
#
# It uses the first principles GCA parameters of the last paper in the
#	multi-orbital tight-binding part of the TMD nanoribbon Hamiltonian.

#	Input the number corresponding to the desired TMD
#	(1) MoS2
#	(2) WS2
#	(3) MoSe2
#	(4) WSe2
#	(5) MoTe2
#	(6) WTe2

#	DEFAULT PARAMETERS

# Number of threads
nthreads=4
# Length of the ribbon
nx=10
# Width of the ribbon
ny=5
# Number of orbitals in the model
norb=3
# Trotter error
dt_inv=8
# Inverse temperature
beta=8
# Frequency of recomputing G
green_afresh_freq=4
# Toggle prints
verbose=1
# Choose source file
source=simulation_eq_time

# Set parameters of the simulation here.
CXX = g++-8 -DNX=$(nx) -DNY=$(ny) -DNORB=$(norb) -DDT_INV=$(dt_inv)\
 -DBETA=$(beta) -DGREEN_AFRESH_FREQ=$(green_afresh_freq) -DVERBOSE=$(verbose)\
 -DEIGEN_DONT_PARALLELIZE -DNTHREADS=$(nthreads) -fopenmp

include_dir=./includes

CXXFLAGS = -Wall -g -O3 -std=c++11 -I$(include_dir)

simulation: src/$(source).o
ifeq ($(verbose),1)
	@echo ""
	@echo "		tmd-nanoribbon-detQMC - Determinant Quantum Monte Carlo for a \
	minimal multi-orbital Hubbard model of a TMD nanoribbon"
	@echo ""
	@echo "			Created by Francisco Brito (6 May 2019)"
	@echo ""
	@echo "The code has compiled successfully. To change the number of threads, the system's dimensions,\
	inverse Trotter error, inverse temperature, the frequency of recomputing \
	the Green's functions, whether or not information about the run is output,\
	or the compiled source file type:"
	@echo ""
	@echo "make clean"
	@echo ""
	@echo "make nx=<Number of sites along x> ny=<Number of sites along y> \
	dt_inv=<Inverse Trotter Error> beta=<Inverse Temperature> \
	green_afresh_freq=<Frequency of Recomputing G>\
	verbose=<0:silent, 1: outputs information about the run>\
	source=<Source File> nthreads=<Number of threads>"
	@echo ""
	@echo "To run a simulation, type ./simulation followed by its arguments:"
	@echo ""
	@echo "./simulation <tmd> <U> <mu> \
	 <Total Number of Sweeps (Space-Time)>\
	 <Number of Warm-up Sweeps (Space-Time)> \
	 <Number of Auto-correlation Sweeps (Space-Time)>"
	@echo ""
	@echo "where U is the on-site interaction (in units of t0),\
	mu is the chemical potential (in particle-hole symmetric form, i.e.\
	\mu -> mu + U /2 )"
	@echo "Input the number corresponding to the desired TMD."
	@echo ""
	@echo "(1) MoS2"
	@echo ""
	@echo "(2) WS2"
	@echo ""
	@echo "(3) MoSe2"
	@echo ""
	@echo "(4) WSe2"
	@echo ""
	@echo "(5) MoTe2"
	@echo ""
	@echo "(6) WTe2"
	@echo ""
endif

	$(CXX) $(CXXFLAGS) -o simulation src/$(source).o


src/$(source).o: src/$(source).cpp $(include_dir)/model.hpp $(include_dir)/aux.hpp\
	 $(include_dir)/green.hpp $(include_dir)/QDT.hpp $(include_dir)/TDQ.hpp

clean:
	rm -f simulation src/*.o
	rm -f ./tmp/*
