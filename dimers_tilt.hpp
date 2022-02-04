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
        override void GetParams(string, int, char*);
        override void Simulate(int);
        void FreeEnergies();
        void ReactionRate();
        vector<float> free_energies;
        string free_energies_txt;
        vector<float> free_energies_ref;

        // MPI variables
        int world_size=1;
        int world_rank=0;
        string local_path;
        vector<string> rank_paths;
        string output;
        string output_path;
        ofstream my_cout;

        // Have to get all stored config files from the other runs
        vector<vector<vector<vector<float>>>> state_database;

};

#endif
