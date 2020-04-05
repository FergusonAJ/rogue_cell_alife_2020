#ifndef ROGUE_MC_CONFIG_H
#define ROGUE_MC_CONFIG_H

#include "config/config.h"
#include <string>

EMP_BUILD_CONFIG(MulticellEvolutionConfig,
    // General Group 
    GROUP(GENERAL, "General settings"),
    VALUE(RANDOM_SEED,              int,    -1,     "Random number seed (-1 to use current time)."),
    VALUE(NUM_UPDATES,              size_t, 1000,   "Number of updates in one expeirment."),
    VALUE(NUM_MULTICELLS,           size_t, 1000,   "How many mutlicells to simulate?"),
    VALUE(PRINT_EVERY_UPDATE,       bool,   false,  "Print the result of every update?"),
    VALUE(GENOME_LENGTH,            size_t, 100,    "Number of bits in the genome."),
    VALUE(PROFILING_DATA_DIR,       std::string, "","Path to the directory containing all "
                                                    "pre-generated replication time data."),
    VALUE(MULTICELL_SIDE_LENGTH,    size_t, 16,     "Multicells are square with this value squared "
                                                    "cells"),
    VALUE(MUTATION_RATE,            double, 0.01,   "Per-bit rate at which genomes are mutated"),
    VALUE(RESTRAINT_STYLE,          size_t, 1,      "0 for basic restraint, 1 for enhanced"),
    VALUE(WRITE_OUTPUT_FILE,        bool,   true,   "If true, write to a summary file every update"),
    VALUE(OUTPUT_FILENAME,          std::string, "output.csv", "Path to write summary data to if "
                                                                " WRITE_OUTPUT_FILE is 1"),
    VALUE(OUTPUT_INTERVAL,          size_t, 1000,  "Only write output roughly this often (updates)"),
    VALUE(MIN_NUM_ONES,             size_t, 40,     "If germ has fewer than this many ones, treat "
                                                        "it as if it has this many"),
    VALUE(MAX_NUM_ONES,             size_t, 60,     "If germ has more than this many ones, treat "
                                                        "it as if it has this many")
)
#endif
