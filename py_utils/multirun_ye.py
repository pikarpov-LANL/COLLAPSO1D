# The script prepares and launches multiple COLLPASO1D runs,
# each running independently on the cores provided. 
# Run initialization is performed in serial, but compilation
# and execution is spread between all available cores via MPI.

# -pikarpov

import numpy as np
import time
from mpi4py import MPI
from setup_run_mpi import multirun

def main():
    
    comm = MPI.COMM_WORLD
    size = comm.Get_size()
    rank = comm.Get_rank()    
    
    suffixs              = ['_g9k_c8.4k_p0.3k_ye','_g9k_c8.4k_p0.3k_ye']#,'_g10k_c9.4k_p0.3k_1.3P']
    masses               = [13.0,15.0]
    enclosed_mass_cutoff = [1.61,1.52]
    
    pns_cutoff      = [1.25 for i in masses]    
    pns_grid_goal   = [300 for i in masses]
    conv_grid_goals = [8400,8400,9000]
    grid_goals      = [9000,9000,10000]
    maxrads         = [1.5e9,1.5e9,2e9] # 1e9 for 9.0 and 10.0 #2e9 for 19.0
    dataset         = 'sukhbold2016'
    base_path       = '/home/pkarpov/runs/ye'
    template_path   = f'{base_path}/template_Ye'
    output_path     = '/home/pkarpov/scratch/1dccsn/sfho_s/encm_tuned/ye'
    eos_table_path  = '/home/pkarpov/COLLAPSO1D/project/1dmlmix/Hempel_SFHoEOS_rho222_temp180_ye60_version_1.3_20190605.h5'
    mlmodel         = 'None'
    read_dump       = 0
    dump_interval   = 5e-4
    restart         = True
    
    mr = multirun(suffixs,masses,enclosed_mass_cutoff,pns_cutoff,
                  dataset,base_path,template_path,output_path,eos_table_path,
                  pns_grid_goal, conv_grid_goals, grid_goals,maxrads, 
                  mlmodel, read_dump, dump_interval, restart)
    
    if rank == 0:        
        if len(mr.masses) != size:
            print(f'Ranks are not distributed well!\nRank size {size} for {len(mr.masses)} datasets')
            comm.Abort()
        if not restart: mr.initialize(rank) # will initialize in serial
            
    if restart: mr.initialize(rank) # processes existing data in parallel
    
    comm.Barrier()
    time.sleep(0.1)
    
    mr.run_name         = f's{mr.masses[rank]}{mr.suffixs[rank]}'
    mr.run_path         = f'{mr.base_path}/{mr.run_name}'
    mr.sim_path         = f'{mr.run_path}/project/1dmlmix'
    mr.full_output_path = f'{mr.output_path}/{mr.run_name}'

    mr.run(rank)              
        
if __name__ == '__main__':
    main()