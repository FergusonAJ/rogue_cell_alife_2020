num_updates = 4000000
num_multicells = 1000000
seed_offset = 12000
output_interval = 10000
replicates_per_condition = 100
genome_length = 100
min_num_ones = 40
max_num_ones = 60
profiling_data_dir = '/mnt/home/fergu358/research/rogue_cell/rogue_cell_profiling/data/100_bits/'
executable_path = '/mnt/home/fergu358/research/rogue_cell/mc_evolution/rogue_cell_alife_2020/mc_evolution/multicell_evolution'
output_dir = '/mnt/gs18/scratch/users/fergu358/rogue_cell/mc_evolution/100_bits_large_pop/'

#mut_rate_L = ['0.05', '0.01', '0.002'] 
mut_rate_L = ['0.002'] 
restraint_type_L = ['0', '1'] 
#side_length_L = ['16', '32', '64']
side_length_L = ['16', '32', '64', '128']
with open('./jobs/run_all_jobs.sh', 'w') as fp_run_all:
    fp_run_all.write('#!/bin/bash\n\n')
    for mut_rate in mut_rate_L:
        for restraint_style in restraint_type_L:
            for side_length in side_length_L: 
                filename_job = \
                        'MUT_' + str(mut_rate) + \
                        '__MCX_' + str(side_length) + \
                        '__RES_' + str(restraint_style) + \
                        '.sb'
                fp_run_all.write('sbatch ' + filename_job + '\n')
                with open('./jobs/' + filename_job, 'w') as fp_job:
                    fp_job.write('#!/bin/bash --login\n')
                    fp_job.write('\n\n')
                    fp_job.write('#SBATCH --time=003:00:00\n')
                    fp_job.write('#SBATCH --nodes=1\n')
                    fp_job.write('#SBATCH --ntasks=1\n')
                    fp_job.write('#SBATCH --cpus-per-task=1\n')
                    fp_job.write('#SBATCH --mem-per-cpu=100m\n')
                    fp_job.write('#SBATCH --job-name rc_mc_evo\n')
                    fp_job.write('#SBATCH --array=1-' + str(replicates_per_condition) + '\n')
                    fp_job.write('\n\n')
                    fp_job.write('module purge\n')
                    fp_job.write('module load GCC/9.1.0-2.32\n')
                    
                    fp_job.write('SEED_OFFSET=' + str(seed_offset - 1) + '\n')
                    fp_job.write('SEED=$((SLURM_ARRAY_TASK_ID + SEED_OFFSET))\n')       
             
                    fp_job.write('\n\n')
                    fp_job.write(executable_path + ' ' + \
                        '-RANDOM_SEED ${SEED} ' + \
                        '-NUM_UPDATES ' + str(num_updates) + ' ' + \
                        '-NUM_MULTICELLS ' + str(num_multicells) + ' ' +\
                        '-PRINT_EVERY_UPDATE 0 ' + \
                        '-GENOME_LENGTH ' + str(genome_length) + ' ' + \
                        '-PROFILING_DATA_DIR ' + str(profiling_data_dir) + ' ' + \
                        '-MULTICELL_SIDE_LENGTH ' + side_length + ' ' + \
                        '-MUTATION_RATE ' + mut_rate + ' ' + \
                        '-RESTRAINT_STYLE ' + restraint_style + ' ' + \
                        '-WRITE_OUTPUT_FILE 1 ' + \
                        '-OUTPUT_FILENAME ' + output_dir + \
                            'germ_alignment__' + \
                            'MUT_' + str(mut_rate) + \
                            '__MCX_' + str(side_length) + \
                            '__RES_' + str(restraint_style) + \
                            '__SEED_${SEED}' + \
                            '.csv' + ' ' + \
                        '-OUTPUT_INTERVAL ' + str(output_interval) + ' ' +\
                        '-MIN_NUM_ONES ' + str(min_num_ones) + ' ' +\
                        '-MAX_NUM_ONES ' + str(max_num_ones) + ' ' +\
                        '\n')
