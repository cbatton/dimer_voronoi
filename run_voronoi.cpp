#include <iostream> 
#include <fstream>
#include <random>
#include <algorithm>
#include <math.h>
#include <ctime>
#include <string>
#include <vector>
#include <chrono>
#include "saruprng.hpp"
#include "dimers_voronoi.hpp"
using namespace std;

int main(int argc, char* argv[]) {
    Dimer system;
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
    //system.DumpStates();
}
