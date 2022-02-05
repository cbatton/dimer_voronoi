#include <iostream> 
#include <fstream>
#include <stdlib.h>
#include <random>
#include <algorithm>
#include <math.h>
#include <ctime>
#include <string>
#include <vector>
#include <mpi.h>
#include <chrono>
#include "saruprng.hpp"
#include "dimers_voronoi.hpp"
#include "dimers_tilt.hpp"
using namespace std;

// LAPACK
extern "C" {
  void dgesv_(int *n, int *nrhs, double *a, int *lda, int *ipiv, double *b, int *ldb, int *info);
}

void DimerTilt::GetParams(string name, int& argc, char* argv[]) {
    // MPI housekeeping
    // Initialize MPI
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    // Get local path
    string path_file(argv[1]);
    ifstream input_0;
    input_0.open(path_file);
    input_0 >> local_path;
    input_0.close();
    // Get path for each processor
    string path_files(argv[2]);
    rank_paths.resize(world_size);
    ifstream input_1;
    string line;
    input_1.open(path_files);
    for(int i=0; i<world_size; i++) {
        input_1 >> rank_paths[i];
        getline(input_1,line);
    }
    input_1.close();
    // Print off a hello world message from all processors
    string output = local_path+rank_paths[world_rank]+"/out";
    output_path = local_path+rank_paths[world_rank]+"/";
    my_cout.open(output, std::ios_base::app);

    // Now get parameters
    ifstream input;
    input.open(output_path+name);
    if(input.fail()) {
        my_cout << "No input file" << endl;
    }
    else {
        string line;
        //my_cout << "Param file detected. Changing values." << endl;
        input >> line >> temp;
        //my_cout << "temp is now " << temp << endl;
        getline(input, line);
        input >> line >> mass;
        //my_cout << "mass is now " << mass << endl;
        getline(input, line);
        input >> line >> gamma;
        //my_cout << "gamma is now " << gamma << endl;
        getline(input, line);
        input >> line >> dt;
        //my_cout << "dt is now " << dt << endl;
        getline(input, line);
        input >> line >> box[0] >> box[1] >> box[2];
        //my_cout << "box is now " << box[0] << " " << box[1] << " " << box[2] << endl;
        getline(input, line);
        input >> line >> height;
        //my_cout << "height is now " << height << endl;
        getline(input, line);
        input >> line >> r_0;
        //my_cout << "r_0 is now " << r_0 << endl;
        getline(input, line);
        input >> line >> width;
        //my_cout << "width is now " << width << endl;
        getline(input, line);
        input >> line >> dist_init;
        //my_cout << "dist_init is now " << dist_init << endl;
        getline(input, line);
        input >> line >> cycles >> storage_time;
        //my_cout << "Cycles " << cycles << " storage_time " << storage_time << endl;
        getline(input, line);
        input >> line >> cycles_equil >> max_step;
        //my_cout << "cycles_equil " << cycles_equil << " max_step " << max_step << endl;
        getline(input, line);
        input >> line >> seed_base >> count_step >> frame_time >> check_time;
        //my_cout << "seed_base " << seed_base << " count_step " << count_step << " frame_time " << frame_time << " check_time " << check_time << endl;
        getline(input, line);
        input >> line >> num_solv;
        //my_cout << "num_solv is now " << num_solv << endl;
        getline(input, line);
        input >> line >> epsilon;
        //my_cout << "epsilon is now " << epsilon << endl;
        getline(input, line);
        input >> line >> dr >> gr_time;
        //my_cout << "dr is now " << dr << " gr_time " << gr_time << endl;
        getline(input, line);
        input >> line >> voronoi_num;
        //my_cout << "voronoi_num is now " << voronoi_num << endl;
        getline(input, line);
        input >> line >> voronoi_txt;
        //my_cout << "voronoi_txt is now " << voronoi_txt << endl;
        getline(input, line);
        input >> line >> cell_tar;
        //my_cout << "cell_tar is now " << cell_tar << endl;
        getline(input, line);
        input >> line >> k_umb >> bond_umb;
        //my_cout << "k_umb " << k_umb << " bond_umb " << bond_umb << endl;
        getline(input, line);
        input >> line >> free_energies_txt;
        //my_cout << "free_energies_txt is now " << free_energies_txt << endl;
        getline(input, line);
        input >> line >> k_hits_ref_time;
        //my_cout << "k_hits_ref_time is now " << k_hits_ref_time << endl;
        getline(input, line);
        input >> line >> k_hits_txt;
        //my_cout << "k_hits_txt is now " << k_hits_txt << endl;
        getline(input, line);

        // Initialize system
        // Initialize particles such that they have distance of dist_init
        // Please don't make dist_init greater than box[2]
        num_particles = num_solv+2;
        state.resize(num_particles, vector<float>(3,0));
        // Put particles on an incomplete cubic lattice
        int num_spacing = ceil(pow(num_particles,1.0/3.0));
        float spacing_x = box[0]/num_spacing;
        float spacing_y = box[1]/num_spacing;
        float spacing_z = box[2]/num_spacing;
        int count = 0;
        int id_x = 0;
        int id_y = 0;
        int id_z = 0;
        while((num_particles)>count) {
            state[id_z+id_y*num_spacing+id_x*num_spacing*num_spacing][0] = spacing_x*id_x-0.5*box[0];
            state[id_z+id_y*num_spacing+id_x*num_spacing*num_spacing][1] = spacing_y*id_y-0.5*box[1];
            state[id_z+id_y*num_spacing+id_x*num_spacing*num_spacing][2] = spacing_z*id_z-0.5*box[2];
            count++;
            id_z++;
            if(id_z==num_spacing) {
                id_z = 0;
                id_y++;
            }
            if(id_y==num_spacing) {
                id_y = 0;
                id_x++;
            }
        }
        // Read voronoi list first
        voronoi.resize(voronoi_num, 0);
        voronoi_boundaries = vector<vector<float>>(voronoi_num,vector<float>(2,0));
        ifstream input_voronoi;
        input_voronoi.open(voronoi_txt);
        for(int i=0; i<voronoi_num; i++) {
            input_voronoi >> voronoi[i];
            //my_cout << i << " " << voronoi[i] << endl;
            getline(input_voronoi, line);
        }
        voronoi_boundaries[0][0] = 0;
        voronoi_boundaries[0][1] = 0.5*(voronoi[0]+voronoi[1]);
        voronoi_boundaries[voronoi_num-1][0] = 0.5*(voronoi[voronoi_num-2]+voronoi[voronoi_num-1]);
        voronoi_boundaries[voronoi_num-1][1] = 0.5*sqrt(3)*box[0];
        for(int i=1; i<(voronoi_num-1); i++) {
            voronoi_boundaries[i][0] = 0.5*(voronoi[i-1]+voronoi[i]);
            voronoi_boundaries[i][1] = 0.5*(voronoi[i]+voronoi[i+1]);
        }
        bond_umb = voronoi[cell_tar];
        // Now read in reference free energy
        free_energies_ref.resize(voronoi_num, 0);
        ifstream input_free_energies;
        input_free_energies.open(free_energies_txt);
        for(int i=0; i<voronoi_num; i++) {
            input_free_energies >> free_energies_ref[i];
            //my_cout << i << " " << free_energies[i] << endl;
            getline(input_free_energies, line);
        }
        // Set free_energies equal to reference to start
        free_energies_ref = free_energies;
        //By convention, first two particles are the dimer
        float phi_bond = 0;
        float phi_wca = 0;
        Energy(phi_bond,phi_wca);
        phi = phi_bond+phi_wca;
        phi_storage = vector<vector<float>>(cycles/storage_time,vector<float>(2,0));
        bond_storage = vector<float>(cycles/storage_time,0.0);
        state_storage = vector<vector<vector<float>>>(cycles/storage_time, vector<vector<float>>(num_particles, vector<float>(3,0)));
        // Hash seed_base
        seed_base = seed_base*0x12345677 + 0x12345;
        seed_base = seed_base^(seed_base>>16);
        seed_base = seed_base*0x45679;
        generator = Saru(seed_base, count_step);
        // Prepare g_r
        num_bins_gr = int(box[0]*0.5*sqrt(3)/dr); 
        g_r_storage = vector<vector<float>>(4,vector<float>(num_bins_gr,0));
        // Prepare k_hits
        k_hits.resize(voronoi_num, vector<int>(voronoi_num,0));
        // Now read in reference k_hits
        k_hits_ref.resize(voronoi_num, vector<float>(voronoi_num,0));
        ifstream input_k_hits;
        input_k_hits.open(k_hits_txt);
        for(int i=0; i<voronoi_num; i++) {
            for(int j=0; j<voronoi_num; j++) {
                input_k_hits >> k_hits_ref[i][j];
                k_hits_ref[i][j] /= k_hits_ref_time;
            }
            getline(input_k_hits, line);
        }
        // Calculate initial reentry probabilities
        reentry_probs.resize(voronoi_num, vector<float>(voronoi_num,0));
        ReentryProbs(free_energies_ref, k_hits_ref, reentry_probs, reentry_cdf);

        // Prepare state database
        vector<vector<vector<float>>> state_database_i;
        vector<vector<float>> states;
        vector<float> pos(3,0);
        state_database.resize(voronoi_num,state_database_i);
        vector<int> config_file_sizes(voronoi_num,0);
        for(int i=0; i<voronoi_num; i++) {
            ifstream input_config;
            // Count number of lines
            input_config.open(output_path+"config_"+to_string(i)+".xyz");
            while(getline(input_config, line)) {
                config_file_sizes[i]++; 
            }
            input_config.close();
            // Now actually read them if it has lines
            if(config_file_sizes[i] > 2) {
                input_config.open(output_path+"config_"+to_string(i)+".xyz");
                int num_configs = config_file_sizes[i]/(num_particles+2);
                state_database[i].resize(num_configs,states);
                for(int num=0; num<num_configs; num++) {
                    state_database[i][num].resize(32, pos);
                    // Skip first two lines, then read in everything
                    getline(input_config, line);
                    getline(input_config, line);
                    for(int part=0; part<num_particles; part++) {
                        input_config >> line >> state_database[i][num][part][0] >> state_database[i][num][part][1] >> state_database[i][num][part][2];
                        getline(input_config, line);
                    }
                }
            }
        }
        // Initialize configuration from state database
        InitializeState();
    }
    // also modify config path
    config_file.open("string_"+to_string(world_rank)+"_config.xyz", std::ios_base::app);

}

void DimerTilt::ReentryProbs(vector<float>& fe_, vector<vector<int>>& k_, vector<vector<float>>& prob_, vector<vector<float>>& cdf_) {
    // Calculate reentry probabilities given free energies fe_ and hitting estimates k_
    // First evaluate normalizing factor
    vector<float> norm(voronoi_num, 0);
    for(int i=0; i<voronoi_num; i++) {
        for(int j=0; j<voronoi_num; j++) {
            norm[i] = fe_[j]*k_[j][i];
        }
    }
    // Now evaluate probabilities
    for(int i=0; i<voronoi_num; i++) {
        for(int j=0; j<voronoi_num; j++) {
            prob_[i][j] = fe_[j]*k_[j][i]/norm[i];
            cdf_[i][j] = 0;
        }
    }
    for(int i=0; i<voronoi_num; i++) {
        for(int j=0; j<voronoi_num; j++) {
            cdf_[j][i] += cdf_[j][i]+prob_[j][i];
        }
    }
}

void DimerTilt::ReentryProbs(vector<float>& fe_, vector<vector<float>>& k_, vector<vector<float>>& prob_, vector<vector<float>>& cdf_) {
    // Calculate reentry probabilities given free energies fe_ and hitting estimates k_
    // First evaluate normalizing factor
    vector<float> norm(voronoi_num, 0);
    for(int i=0; i<voronoi_num; i++) {
        for(int j=0; j<voronoi_num; j++) {
            norm[i] = fe_[j]*k_[j][i];
        }
    }
    // Now evaluate probabilities
    for(int i=0; i<voronoi_num; i++) {
        for(int j=0; j<voronoi_num; j++) {
            prob_[i][j] = fe_[j]*k_[j][i]/norm[i];
            cdf_[i][j] = 0;
        }
    }
    for(int i=0; i<voronoi_num; i++) {
        for(int j=0; j<voronoi_num; j++) {
            cdf_[j][i] += cdf_[j][i]+prob_[j][i];
        }
    }
}

void DimerTilt::InitializeState() {
    // Reinitialize state from database
    float chance_state = generator.f();
    // Figure out what what state to choose
    int alpha_init = -1;
    for(int i=0; i<voronoi_num; i++) {
        if(chance_state <= reentry_cdf[i][cell_tar]) {
            alpha_init = i; 
        }
    }
    // Now have state, determine what configuration to select
    int database_size = state_database[alpha_init].size();
    if(database_size == 0) {
        // Look to left and right to see if there is another option
        if((alpha_init > 0) && ( alpha_init < (voronoi_num-1))) {
            int database_size_left = state_database[alpha_init-1].size(); 
            int database_size_right = state_database[alpha_init+1].size(); 
            if(database_size_left > database_size_right) {
                alpha_init = alpha_init-1;
            }
            else {
                alpha_init = alpha_init+1;
            }
        }
    }
    int database_config = generator.rand_select(database_size-1);
    state = state_database[alpha_init][database_config];
}

void DimerTilt::Simulate(int steps) {
    // Run simulation
    ofstream config_file_2;
    config_file_2.precision(10);
    config_file_2.open(config, std::ios_base::app);
    // Dump configurations that cross Voronoi cells
    vector<ofstream> config_files;
    config_files.resize(voronoi_num);
    for(int i=0; i<voronoi_num; i++) {
        config_files[i].precision(10);
        config_files[i].open("config_"+to_string(i)+".xyz", std::ios_base::app);
    }
    double time_counter = 0;
    // Run things in stages where we update the free energies
    for(int stage=0; stage<steps/1000; stage++) {
        if(stage > 0) {
            // Re-evaluate free energies
            // Use MPI to gather all k_hits values using reduce
            MPI_Barrier(MPI_COMM_WORLD);
            vector<vector<int>> k_hits_local;
            k_hits_local.resize(voronoi_num, vector<int>(voronoi_num,0));
            for(int i=0; i<voronoi_num; i++) {
                MPI_Allreduce(k_hits[i].data(), k_hits_local.data(), voronoi_num, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
            }
            // Modify k_hits[i][voronoi-1] to 0
            for(int i=0; i<voronoi_num; i++) {
                k_hits_local[i][voronoi_num-1] = 0;
            }
            // Now process as needed for free energies
            vector<vector<float>> k_hits_local_2;
            k_hits_local_2.resize(voronoi_num, vector<float>(voronoi_num,0));
            for(int i=0; i<voronoi_num; i++) {
                for(int j=0; j<voronoi_num; j++) {
                    k_hits_local_2[i][j] = k_hits_local[i][j]/time_counter;
                }
            }
            // Modify k_hits[0] to match reference
            for(int i=0; i<voronoi_num; i++) {
                k_hits_local_2[0][i] = k_hits_ref[0][i];
            }
            // Now calculate free energies, reentry probabilities
            FreeEnergies(k_hits_local_2);
            ReentryProbs(free_energies, k_hits_local_2, reentry_probs, reentry_cdf);
        }
        for(int i=0; i<steps; i++) {
            generator = Saru(seed_base, count_step++);
            BDStep();
            // Check to see if in Voronoi cell
            float bond_len = BondLength();
            time_counter += dt;
            // Check to see if we crossed a voronoi
            int voronoi_check = VoronoiIndex(bond_len);
            if(voronoi_check != cell_tar) {
                // Reset
                k_hits[cell_tar][voronoi_check] += 1;
                InitializeState();
            }
            if(count_step%check_time==0) {
                float phi_bond = 0;
                float phi_wca = 0;
                Energy(phi_bond,phi_wca);
                phi = phi_bond+phi_wca;
                my_cout << "Cycle " << i << " phi_bond " << phi_bond << " phi_wca " << phi_wca << endl;
            }
            if(count_step%storage_time==0) {
                float phi_bond = 0;
                float phi_wca = 0;
                Energy(phi_bond,phi_wca);
                float bond_len = BondLength();
                phi_storage[i/storage_time][0] = phi_bond;
                phi_storage[i/storage_time][1] = phi_wca;
                bond_storage[i/storage_time] = bond_len;
                state_storage[i/storage_time]= state;
            }
            if(count_step%frame_time==0) {
                DumpXYZ(config_file_2);
            }
        }
    }
}

void DimerTilt::FreeEnergies(vector<vector<float>>& k_) {
    // Evaluate free energies
    // Have to solve a matrix equation to do so
    // Construct those terms 
    // \sum_b fe_b k_b,a = \sum_b fe_a k_a,b
    // with fe_0 = fe_{0,ref}, fe_{voronoi_num-1} = 0
    // Solve for the rest of them
    vector<double> k_matrix((voronoi_num-2)*(voronoi_num-2));
    vector<double> b(voronoi_num-2);
    // Diagonal terms first
    for(int i=1; i<voronoi_num-1; i++) {
        for(int j=1; j<voronoi_num-1; j++) {
            k_matrix[(i-1)*voronoi_num] -= k_[i][j];
        }
    }
    // Now off-diagonal
    for(int i=1; i<voronoi_num-1; i++) {
        for(int j=1; j<voronoi_num-1; j++) {
            if(i != j) {
                k_matrix[(i-1)+(j-1)*voronoi_num] += k_[j][i];
            }
        }
    }
    // Now terms associated with a = 0
    for(int i=1; i<voronoi_num-1; i++) {
        b[i-1] = -k_hits[0][i]*free_energies_ref[0];
        k_matrix[(i-1)*voronoi_num] -= k_[i][0];
    }
    // a = voronoi_num-1 has zero free energy, so those are handled automatically
    // solving time
    int dim = voronoi_num-2;
    int info;
    int one = 1;
    char N = 'N';
    vector<int> ipiv(voronoi_num-2);
    dgesv_(&dim, &one, k_matrix.data(), &dim, ipiv.data(), b.data(), &dim, &info);
    // Free energies are stored in b, so extract them
    for(int i=1; i<voronoi_num-1; i++) {
        free_energies[i] = b[i-1];
    }
}
