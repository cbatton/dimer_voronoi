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
        void GetParams(string, int&, char**);
        void ReentryProbs(vector<float>&, vector<vector<int>>&, vector<vector<float>>&, vector<vector<float>>&);
        void ReentryProbs(vector<float>&, vector<vector<float>>&, vector<vector<float>>&, vector<vector<float>>&);
        void InitializeState();
        void Simulate(int);
        void FreeEnergies(vector<vector<float>>&);
        void ReactionRate();
        double time_counter; // total simulation time
        vector<float> free_energies; // more aptly the probabilities of states
        string free_energies_txt; // file to read free energies from
        double k_hits_ref_time; // time for reference fluxes
        string k_hits_txt; // file to read references fluxes from
        vector<float> free_energies_ref; // reference free energies
        vector<vector<float>> k_hits_ref; // reference fluxes
        vector<vector<float>> reentry_probs; 
        vector<vector<float>> reentry_cdf;

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
