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
#include "dimers_voronoi.hpp"
#include "dimers_tilt.hpp"
using namespace std;

void DimerTilt::GetParams(string name, int argc, char* argv[]) {
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
        mycout << "No input file" << endl;
    }
    else {
        string line;
        //mycout << "Param file detected. Changing values." << endl;
        input >> line >> temp;
        //mycout << "temp is now " << temp << endl;
        getline(input, line);
        input >> line >> mass;
        //mycout << "mass is now " << mass << endl;
        getline(input, line);
        input >> line >> gamma;
        //mycout << "gamma is now " << gamma << endl;
        getline(input, line);
        input >> line >> dt;
        //mycout << "dt is now " << dt << endl;
        getline(input, line);
        input >> line >> box[0] >> box[1] >> box[2];
        //mycout << "box is now " << box[0] << " " << box[1] << " " << box[2] << endl;
        getline(input, line);
        input >> line >> height;
        //mycout << "height is now " << height << endl;
        getline(input, line);
        input >> line >> r_0;
        //mycout << "r_0 is now " << r_0 << endl;
        getline(input, line);
        input >> line >> width;
        //mycout << "width is now " << width << endl;
        getline(input, line);
        input >> line >> dist_init;
        //mycout << "dist_init is now " << dist_init << endl;
        getline(input, line);
        input >> line >> cycles >> storage_time;
        //mycout << "Cycles " << cycles << " storage_time " << storage_time << endl;
        getline(input, line);
        input >> line >> cycles_equil >> max_step;
        //mycout << "cycles_equil " << cycles_equil << " max_step " << max_step << endl;
        getline(input, line);
        input >> line >> seed_base >> count_step >> frame_time >> check_time;
        //mycout << "seed_base " << seed_base << " count_step " << count_step << " frame_time " << frame_time << " check_time " << check_time << endl;
        getline(input, line);
        input >> line >> num_solv;
        //mycout << "num_solv is now " << num_solv << endl;
        getline(input, line);
        input >> line >> epsilon;
        //mycout << "epsilon is now " << epsilon << endl;
        getline(input, line);
        input >> line >> dr >> gr_time;
        //mycout << "dr is now " << dr << " gr_time " << gr_time << endl;
        getline(input, line);
        input >> line >> voronoi_num;
        //mycout << "voronoi_num is now " << voronoi_num << endl;
        getline(input, line);
        input >> line >> voronoi_txt;
        //mycout << "voronoi_txt is now " << voronoi_txt << endl;
        getline(input, line);
        input >> line >> cell_tar;
        //mycout << "cell_tar is now " << cell_tar << endl;
        getline(input, line);
        input >> line >> k_umb >> bond_umb;
        //mycout << "k_umb " << k_umb << " bond_umb " << bond_umb << endl;
        getline(input, line);
        input >> line >> free_energies_txt;
        //mycout << "free_energies_txt is now " << free_energies_txt << endl;
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
        // Initialize particles such that they are in middle of cell_tar
        // Read voronoi list first
        voronoi.resize(voronoi_num, 0);
        voronoi_boundaries = vector<vector<float>>(voronoi_num,vector<float>(2,0));
        ifstream input_voronoi;
        input_voronoi.open(voronoi_txt);
        for(int i=0; i<voronoi_num; i++) {
            input_voronoi >> voronoi[i];
            //mycout << i << " " << voronoi[i] << endl;
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
        for(int i=0; i<free_energies_num; i++) {
            input_free_energies >> free_energies_ref[i];
            //mycout << i << " " << free_energies[i] << endl;
            getline(input_free_energies, line);
        }
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
    }
    // also modify config path
    config_file.open("string_"+to_string(rank_in)+"_config.xyz", std::ios_base::app);

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
    for(int i=0; i<steps; i++) {
        generator = Saru(seed_base, count_step++);
        vector<vector<float>> state_old(state);
        BDStep();
        // Check to see if in Voronoi cell
        float bond_len = BondLength();
        time_counter += dt;
        // Check to see if we crossed a voronoi
        int voronoi_check = VoronoiIndex(bond_len);
        if(voronoi_check != cell_tar) {
            // Reset
            k_hits[cell_tar][voronoi_check] += 1;
            if(k_hits[cell_tar][voronoi_check]%1000==0) { 
                DumpXYZ(config_files[voronoi_check]);
            }
            state = state_old;
        }
        if(i%check_time==0) {
            float phi_bond = 0;
            float phi_wca = 0;
            Energy(phi_bond,phi_wca);
            phi = phi_bond+phi_wca;
            mycout << "Cycle " << i << " phi_bond " << phi_bond << " phi_wca " << phi_wca << endl;
        }
        if(i%storage_time==0) {
            float phi_bond = 0;
            float phi_wca = 0;
            Energy(phi_bond,phi_wca);
            float bond_len = BondLength();
            phi_storage[i/storage_time][0] = phi_bond;
            phi_storage[i/storage_time][1] = phi_wca;
            bond_storage[i/storage_time] = bond_len;
            state_storage[i/storage_time]= state;
        }
        if(i%frame_time==0) {
            DumpXYZ(config_file_2);
        }
        if(i%gr_time==0) {
            RDFSample();
        }
    }
}
