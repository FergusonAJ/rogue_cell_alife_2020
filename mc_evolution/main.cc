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
    size_t generation;
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
    int max_generations;
    // Population
    emp::vector<Multicell> multicell_vec;
    emp::vector<size_t> replicating_idx_vec;
    size_t num_replicating;
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
    size_t last_output_update;
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
        max_generations =       (int)           config.MAX_GENERATIONS();
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
            multicell_vec[mc_idx].generation = 0;
        }
    }
    void LoadProfilingData(){
        profiling_data_vec.resize(genome_length + 1); // +1 to account for zero ones
        std::stringstream ss; 
        std::ifstream fp;
        std::string line;
        size_t num_lines = 0;
        for(size_t i = 0; i <= genome_length; i++){
            emp::vector<size_t>& vec = profiling_data_vec[i];
            // Only look for alignments that are between the specified parameters
            if(i < min_num_ones || i > max_num_ones){
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
        ++mc_child.generation;
        // Reset self
        MutateGenome(mc_parent, mc_parent);
        SetMulticellTimes(mc_parent);
        ++mc_parent.generation;
    }
    double GetAlignmentMean(){
        double running_sum = 0;
        for(size_t mc_idx = 0; mc_idx < num_multicells; ++mc_idx){
            running_sum += multicell_vec[mc_idx].genome;
        }
        return running_sum / num_multicells;
    }
    double GetGenerationMean(){
        double running_sum = 0;
        for(size_t mc_idx = 0; mc_idx < num_multicells; ++mc_idx){
            running_sum += multicell_vec[mc_idx].generation;
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
        fp_out << "update,mean_germ_alignment,mean_generation\n";
    }
    void AddOutputLine(){
        fp_out << update_cur << "," << GetAlignmentMean() << "," << GetGenerationMean() << "\n"; 
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
        last_output_update = 0;
        LoadProfilingData();
        InitializePopulation();
        if(write_output_file){ // Create output file and save off update 0
            OpenOutputFile();
            AddOutputLine();
        }
        replicating_idx_vec.resize(num_multicells, 0);
        num_replicating = 0;
        size_t next_repro_update;
        Multicell& cur_mc = multicell_vec[0];
        while(1){
            // We can fast forward to the next replication event, we just need to find it!
            next_repro_update = -1; // MAX
            for(size_t mc_idx = 0; mc_idx < num_multicells; ++mc_idx){
                cur_mc = multicell_vec[mc_idx];
                if(cur_mc.update_repro > next_repro_update)
                    continue;
                if(cur_mc.update_repro == next_repro_update)
                    replicating_idx_vec[num_replicating++] = mc_idx;
                else if(cur_mc.update_repro < next_repro_update){
                    next_repro_update = cur_mc.update_repro;
                    num_replicating = 0;
                    replicating_idx_vec[num_replicating++] = mc_idx;
                }
            }
            if(write_output_file){
                if(next_repro_update - last_output_update > output_interval){
                    // Jump forward as many intervals as we can (leverage integer division)
                    update_cur = last_output_update + output_interval * 
                        ((next_repro_update - last_output_update) / output_interval);
                    AddOutputLine();
                    last_output_update = update_cur;
                }
            }
            update_cur = next_repro_update;
            if(update_cur > num_updates)
                break;
            if(max_generations != -1){
                if(GetGenerationMean() > max_generations)
                    break;
            }
            // Shuffle the vector to make sure we treat multicells fairly
            //Modified emp::shuffle to only shuffle first n elements
            if(num_replicating > 1){
                for (size_t i = 0; i < num_replicating; i++) {
                  const size_t pos = random.GetUInt(i, num_replicating);
                  if (pos == i) continue;
                  std::swap(replicating_idx_vec[i], replicating_idx_vec[pos]);
                }
            }
            for(size_t idx = 0; idx < num_replicating; ++idx){
                if(multicell_vec[replicating_idx_vec[idx]].update_repro == update_cur){
                    ReplicateMulticell(replicating_idx_vec[idx]);
                }
            }
            if(print_every_update){
                std::cout << "Finshed update " << update_cur 
                    << " mean alignment: " << GetAlignmentMean() << std::endl; 
            }
        }
        std::cout << "Final update: " << update_cur 
            << " mean alignment: " << GetAlignmentMean()
            << " mean generations: " << GetGenerationMean() << std::endl;
        if(write_output_file){
            if(update_cur != last_output_update)
                AddOutputLine();
            std::cout << "Ouput file saved to: " << output_filename << std::endl;
            fp_out.close();
        }
    }
};

int main(int argc, char* argv[])
{
    Experiment experiment(argc, argv);
    experiment.Run();
    return 0;
}
