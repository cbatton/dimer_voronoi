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
#ifndef DIMER_T_
#define DIMER_T_

class DimerTilt : public Dimer {
    public:
        void Simulate(int);
        void FreeEnergies();
        void ReactionRate();
        vector<float> free_energies;
};

#endif
