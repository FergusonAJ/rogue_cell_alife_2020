// Standard 
#include <iostream>
#include <algorithm>
#include <random>
#include <iterator>
#include <fstream>
// Empirical
#include "config/ArgManager.h"
#include "config/command_line.h"
#include "config/command_line.h"
#include "tools/Random.h"
#include "tools/random_utils.h"
#include "base/vector.h"
// Local
#include "./rogue_cell_mc_config.h"

// Main goal: See just how far from the threshold the germs will evolve
// Approach: Abstract out the innards of the multicell, using replication times gathered
    // from empiricial data

struct Multicell{
    size_t genome;
    size_t update_birth;
    size_t update_repro;
};

class Experiment{
private:
    // Config vars
    MulticellEvolutionConfig config;
    std::string config_filename = "experiment_config.cfg";
    int random_seed;
    size_t num_updates;
    size_t num_multicells;
    bool print_every_update;
    size_t genome_length;
    size_t restraint_style;  
    std::string profiling_data_dir;
    bool write_output_file;
    std::string output_filename;
    size_t output_interval;
    size_t min_num_ones;
    size_t max_num_ones;
    // Population
    emp::vector<Multicell> multicell_vec;
    // Other evolution vars
    size_t update_cur;
    emp::vector<emp::vector<size_t>> profiling_data_vec;
    size_t multicell_side_length;
    double mutation_rate;    
    emp::vector<size_t> tmp_genome;
    //Bookkeeping
    size_t num_attempted_mutations = 0;
    size_t num_actual_mutations = 0;
    std::ofstream fp_out;
    size_t next_output_update;
    // Misc.
    emp::Random random;
         
    void LoadConfig(int argc, char* argv[]){
        // Use command line args from Empirical
        auto args = emp::cl::ArgManager(argc, argv);
        config.Read(config_filename);
        if (args.ProcessConfigOptions(config, std::cout, config_filename, "mc_evo_config-macros.h")
                == false)
            exit(0);
        if (args.TestUnknown() == false)
            exit(0); // If there are leftover args, throw an error. 
        // Slurp in variables from config file
        random_seed =           (int)           config.RANDOM_SEED();
        num_updates =           (size_t)        config.NUM_UPDATES();
        num_multicells =        (size_t)        config.NUM_MULTICELLS();
        print_every_update =    (bool)          config.PRINT_EVERY_UPDATE();
        genome_length =         (size_t)        config.GENOME_LENGTH();
        profiling_data_dir =    (std::string)   config.PROFILING_DATA_DIR();
        multicell_side_length = (size_t)        config.MULTICELL_SIDE_LENGTH();
        mutation_rate =         (double)        config.MUTATION_RATE();
        restraint_style =       (size_t)        config.RESTRAINT_STYLE();
        write_output_file =     (bool)          config.WRITE_OUTPUT_FILE();
        output_filename =       (std::string)   config.OUTPUT_FILENAME();
        output_interval =       (size_t)        config.OUTPUT_INTERVAL();
        min_num_ones =          (size_t)        config.MIN_NUM_ONES();
        max_num_ones =          (size_t)        config.MAX_NUM_ONES();
    }
    void RandomizeGenome(Multicell& multicell){
        multicell.genome = genome_length / 2;//random.GetUInt(genome_length + 1);
    }
    void SetMulticellTimes(Multicell& multicell){
        size_t alignment_proxy = multicell.genome;
        if(multicell.genome < min_num_ones)
            alignment_proxy = min_num_ones;
        else if(multicell.genome > max_num_ones)
            alignment_proxy = max_num_ones;
        multicell.update_birth = update_cur;
        multicell.update_repro = update_cur 
            + profiling_data_vec[alignment_proxy][
                random.GetUInt(profiling_data_vec[alignment_proxy].size())
            ]; 
    }
    void InitializePopulation(){
        multicell_vec.resize(num_multicells);
        for(size_t mc_idx = 0; mc_idx < num_multicells; ++mc_idx){
            RandomizeGenome(multicell_vec[mc_idx]);
            SetMulticellTimes(multicell_vec[mc_idx]);
        }
        //std::sort(multicell_vec.begin(), multicell_vec.end(), 
        //        [](const Multicell& mc_a, const Multicell& mc_b){
        //            return mc_a.update_repro < mc_b.update_repro;
        //        });
        //for(size_t mc_idx = 0; mc_idx < num_multicells; ++mc_idx){
        //    std::cout 
        //        << multicell_vec[mc_idx].genome << " "
        //        << multicell_vec[mc_idx].update_birth << " "
        //        << multicell_vec[mc_idx].update_repro << std::endl;
        //}
    }
    void LoadProfilingData(){
        profiling_data_vec.resize(genome_length + 1); // +1 to account for zero ones
        std::stringstream ss; 
        std::ifstream fp;
        std::string line;
        size_t num_lines = 0;
        for(size_t i = 0; i <= genome_length; i++){
            emp::vector<size_t>& vec = profiling_data_vec[i];
            if(i < min_num_ones || i > max_num_ones){
                std::cout << i << " outside range!" << std::endl;
                continue;
            }
            double evil_pct = ((double)i) / genome_length;
            std::cout << "Loading data for " << evil_pct << std::endl;
            ss.str("");
            ss << profiling_data_dir << "profiling_data" 
                << "__MCX_" << multicell_side_length
                << "__MUT_" << mutation_rate 
                << "__RES_" << restraint_style 
                << "__WEP_" << evil_pct 
                << ".csv";
            fp.open(ss.str(), std::ios::in);
            // If we failed to open file, error out!
            if(!fp.is_open()){
                std::cerr << "Error! Unable to load profiling data file: " 
                    << ss.str() << std::endl;
                exit(-1);
            }
            // Count the number of lines
            num_lines = 0;
            while(std::getline(fp, line))
                ++num_lines;
            std::cout << "\tNumber of lines: " << num_lines << std::endl;
            vec.resize(num_lines);
            fp.clear();
            fp.seekg(0, fp.beg);
            for(size_t val_idx = 0; val_idx < num_lines; ++val_idx){
                fp >> vec[val_idx];
            }
            fp.close();
        }
        std::cout << "Done loading profiling data!" << std::endl;
    }
    void MutateGenome(Multicell& child, Multicell& parent){
        ++num_attempted_mutations;
        size_t num_ones = 0;
        for(size_t gene_idx = 0; gene_idx < genome_length; ++gene_idx){
            if(gene_idx < parent.genome)
                tmp_genome[gene_idx] = 1;
            else
                tmp_genome[gene_idx] = 0;
            if(random.GetDouble(1.0) < mutation_rate)
                tmp_genome[gene_idx] = (tmp_genome[gene_idx] + 1) % 2;
            if(tmp_genome[gene_idx] == 1)
                ++num_ones;
        }
        child.genome = num_ones;
        if(child.genome != parent.genome)
            ++num_actual_mutations;
    }
    void ReplicateMulticell(size_t mc_idx){
        size_t child_idx = random.GetUInt(num_multicells);
        while(child_idx == mc_idx)
            child_idx = random.GetUInt(num_multicells);
        Multicell& mc_child = multicell_vec[child_idx];
        Multicell& mc_parent = multicell_vec[mc_idx];
        // Add offspring
        MutateGenome(mc_child, mc_parent);
        SetMulticellTimes(mc_child);
        // Reset self
        MutateGenome(mc_parent, mc_parent);
        SetMulticellTimes(mc_parent);
    }
    double GetAlignmentMean(){
        double running_sum = 0;
        for(size_t mc_idx = 0; mc_idx < num_multicells; ++mc_idx){
            running_sum += multicell_vec[mc_idx].genome;
        }
        return running_sum / num_multicells;
    }
    void OpenOutputFile(){
        if(fp_out.is_open()){
            std::cerr << "Error! Attempted to open output file twice!" << std::endl;
            exit(-1);
        }
        fp_out.open(output_filename, std::ios::out);
        if(!fp_out.is_open()){
            std::cerr << "Error! Could not open output file: " << output_filename << std::endl;
            exit(-1);
        }
        fp_out << "update,mean_germ_alignment\n";
    }
    void AddOutputLine(){
        if(update_cur >= next_output_update){
            fp_out << update_cur << "," << GetAlignmentMean() << "\n"; 
            next_output_update = ((update_cur / output_interval) + 1) * output_interval;
        }
    }
public:
    Experiment(int argc, char* argv[]){
        LoadConfig(argc, argv);
        WriteConfig(std::cout);
        random.ResetSeed(random_seed);
        tmp_genome.resize(genome_length, 0);
    }
    void WriteConfig(std::ostream& oss){
        // Write to screen how the experiment is configured
        oss << "==============================" << std::endl;
        oss << "|    Current configuration   |" << std::endl;
        oss << "==============================" << std::endl;
        config.Write(oss);
        oss << "==============================\n" << std::endl;
    }
    void Run(){
        update_cur = 0;
        next_output_update = 0;
        LoadProfilingData();
        InitializePopulation();
        if(write_output_file){
            OpenOutputFile();
        }
        size_t next_repro_update;
        while(update_cur <= num_updates){
            // We can fast forward to the next replication event, we just need to find it!
            next_repro_update = -1; // MAX
            for(size_t mc_idx = 0; mc_idx < num_multicells; ++mc_idx){
                if(multicell_vec[mc_idx].update_repro < next_repro_update)
                    next_repro_update = multicell_vec[mc_idx].update_repro;
            }
            update_cur = next_repro_update;
            // Shuffle the vector to make sure we treat multicells fairly
            emp::Shuffle(random, multicell_vec); 
            for(size_t mc_idx = 0; mc_idx < num_multicells; ++mc_idx){
                if(multicell_vec[mc_idx].update_repro == update_cur){
                    ReplicateMulticell(mc_idx);
                }
            }
            if(print_every_update){
                std::cout << "Finshed update " << update_cur 
                    << " mean alignment: " << GetAlignmentMean() << std::endl;
            }
            AddOutputLine();
        }
        std::cout << "Final update: " << update_cur 
            << " mean alignment: " << GetAlignmentMean() << std::endl;
        if(write_output_file){
            std::cout << "Ouput file saved to: " << output_filename << std::endl;
            fp_out.close();
        }
        //std::cout << "Realized *multicell* mutation rate: " 
        //    << ((double)num_actual_mutations) / num_attempted_mutations
        //    << std::endl;
    }
};

int main(int argc, char* argv[])
{
    Experiment experiment(argc, argv);
    experiment.Run();
    return 0;
}
