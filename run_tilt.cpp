#include <iostream> 
#include <fstream>
#include <random>
#include <algorithm>
#include <math.h>
#include <ctime>
#include <string>
#include <vector>
#include <mpi.h>
#include <chrono>
#include "saruprng.hpp"
#include "dimers_tilt.hpp"
using namespace std;

int main(int argc, char* argv[]) {
    // Initialize MPI
    MPI_Init(&argc, &argv);
    DimerTilt system;
    system.GetParams("param", 0);
    system.Equilibriate(system.cycles_equil);
    ofstream myfile_equil;
    myfile_equil.precision(10);
    myfile_equil.open("config_equil.xyz", std::ios_base::app);
    system.DumpXYZ(myfile_equil);
    system.Simulate(system.cycles);
    system.DumpPhi();
    system.DumpBond();
    system.RDFAnalyzer();
    system.DumpVoronoi();
    // Finalize the MPI environment
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
}
