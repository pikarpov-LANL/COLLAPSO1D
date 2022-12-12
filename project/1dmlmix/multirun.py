# The script prepares and launches multiple COLLPASO1D runs,
# each running independently on the cores provided. 
# Run initialization is performed in serial, but compilation
# and execution is spread between all available cores via MPI.

# -pikarpov

import numpy as np
import os
import sys
import shutil
import time
from subprocess import Popen, PIPE
from mpi4py import MPI

def main():
    
    comm = MPI.COMM_WORLD
    size = comm.Get_size()
    rank = comm.Get_rank()    
    
    suffix = '_g9k_c8.4k_p_0.3k'
    masses = [#11.0,]
              12.0,13.0,14.0,15.0,
              16.0,17.0,18.0]#,19.0,20.0]
    #enclosed_mass_cutoff = [1.5 for i in range(len(masses))]
    #enclosed_mass_cutoff = [1.31,1.31,1.33,1.35,1.35,1.37,1.37,1.35,1.37,1.39]
    #enclosed_mass_cutoff = [1.4,1.42,1.46,1.46,1.46,1.46,1.47,1.48,1.49,1.5]
    enclosed_mass_cutoff = [#1.48,
                            1.49,1.61,1.61,1.52,
                            1.55,1.57,1.55]#,1.63,1.8]
    # failed: 15,16,17,18,19,20
    #enclosed_mass_cutoff = [1.3,1.3] # for 9.0 and 10.0
    #pns_cutoff     = [i-0.15 for i in enclosed_mass_cutoff]     
    #pns_cutoff[-1]-= 0.10
    pns_cutoff     = [1.25 for i in enclosed_mass_cutoff]
    pns_grid_goal  = 300
    conv_grid_goal = 8400
    grid_goal      = 9000
    maxrad         = 1.5e9 # 1e9 for 9.0 and 10.0
    dataset        = 'sukhbold2016'
    base_path      = '/home/pkarpov/runs'
    output_path    = '/home/pkarpov/scratch/1dccsn/sfho_s/encm_tuned'
    eos_table_path = '/home/pkarpov/COLLAPSO1D/project/1dmlmix/Hempel_SFHoEOS_rho222_temp180_ye60_version_1.3_20190605.h5'
    mlmodel        = 'None'
    read_dump      = 0
    dump_interval  = 5e-4
    restart        = True
    
    mr = multirun(suffix,masses,enclosed_mass_cutoff,pns_cutoff,
                  dataset,base_path,output_path,eos_table_path,
                  pns_grid_goal, conv_grid_goal, grid_goal,maxrad, 
                  mlmodel, read_dump, dump_interval, restart)
    
    if rank == 0:        
        if len(mr.masses) != size:
            print(f'Ranks are not distributed well!\nRank size {size} for {len(mr.masses)} datasets')
            comm.Abort()
        else: mr.initialize()
    
    comm.Barrier()
    
    mr.run_name = f's{mr.masses[rank]}{mr.suffix}'
    mr.run_path = f'{mr.base_path}/{mr.run_name}'
    mr.sim_path = f'{mr.run_path}/project/1dmlmix'
    mr.full_output_path = f'{mr.output_path}/{mr.run_name}'

    mr.run(rank)
    
class multirun:
    def __init__(self, suffix, masses, enclmass_conv_cutoff,pns_cutoff,
                 dataset,base_path, output_path, eos_table_path, 
                 pns_grid_goal, conv_grid_goal, grid_goal,
                 maxrad=5e9, mlmodel='None',
                 read_dump=0, dump_interval=1e-3, restart=False,
                 data_names=None, maxtime=0.5):
        
        self.suffix         = suffix
        self.masses         = masses
        self.enclmass_conv_cutoff = enclmass_conv_cutoff
        self.pns_cutoff     = pns_cutoff
        self.dataset        = dataset
        self.base_path      = base_path
        self.output_path    = output_path
        self.data_names     = data_names
        self.template_path  = f'{self.base_path}/template'
        self.data_path      = self.template_path#f'{self.base_path}/template'#produce_data'
        self.eos_table_path = eos_table_path
        self.maxrad         = maxrad
        self.mlmodel        = mlmodel
        self.pns_grid_goal  = pns_grid_goal
        self.conv_grid_goal = conv_grid_goal
        self.grid_goal      = grid_goal
        self.read_dump      = read_dump
        self.dump_interval  = dump_interval
        self.maxtime        = maxtime
        self.restart        = restart
        
    def run(self, rank):
        
        if rank == 0: print('Compiling...')        
                    
        # run 'make project'
        os.chdir(f'{self.run_path}')

        if not self.restart:
            p = Popen('make eos', shell=True, stdin=PIPE, stdout=PIPE, stderr=PIPE)
            output = p.stdout.read()
            p.stdout.close()
            
            p = Popen('make project', shell=True, stdin=PIPE, stdout=PIPE, stderr=PIPE)
            output = p.stdout.read()
            p.stdout.close()        
                
        # run the simulation
        os.chdir(self.sim_path)
        
        #stdout = f'{self.full_output_path}/stdout'
        #if os.path.exists(stdout): os.remove(stdout) 
        stdout = '/dev/null'
        stderr = f'{self.full_output_path}/stderr'
        if os.path.exists(stderr): os.remove(stderr)
                
        counter = 0
        while not os.path.isfile('1dmlmix'):
            time.sleep(1)
            counter += 1
            if counter == 120: sys.exit("ERROR: executable '1dmlmix' not found (waited 2 mins).")
        print(f'rank {rank} compiled project {self.run_name}; running...')
        
        p = Popen(f'time ./1dmlmix > {stdout} 2> {stderr}', shell=True, stdin=PIPE, stdout=PIPE, stderr=PIPE)
        output = p.stdout.read()
        p.stdout.close()
        print('--------------------------------------------')
        print(f'rank {rank} finished running {self.run_name}')
        print('--------------------------------------------')
        
            
    def initialize(self):
        for i in range(len(self.masses)):
            self.mass          = self.masses[i]
            self.enclmass_conv = self.enclmass_conv_cutoff[i]
            self.enclmass_pns  = self.pns_cutoff[i]
            
            if self.data_names != None: self.data_in = self.data_names[i]
            
            self.run_name         = f's{self.mass}{self.suffix}'
            self.run_path         = f'{self.base_path}/{self.run_name}'
            self.full_output_path = f'{self.output_path}/{self.run_name}'
                                
            print(f'---------------- {self.run_name} ---------------')
            
            print(f'Mass:                   {self.mass}')
            print(f'Enclosed PNS  Cutoff:   {self.enclmass_pns}')
            print(f'Enclosed Conv Cutoff:   {self.enclmass_conv}')            
                                
            # check if run_folder exists; create and copy template if not
            self.edit_copy_template()
                        
            if self.restart: 
                self.find_last_dump()
                print(f'Restart from dump:      {self.read_dump}')
            else:
                # prep data and copy it to the run folder
                if self.data_names==None: self.data_in = 'Data'
                self.prep_data()                
                self.data_out = f'{self.full_output_path}/DataOut'
            
            # edit 'setup' to include unique output path
            self.setup()
            
            # edit 'setup_readout' to include unique output path
            self.setup_readout()     
            
        print('--------- Initialization Completed ---------')                                       


    def edit_copy_template(self):     
        if not os.path.exists(self.run_path): 
            # Edit Makefile
            run_project_dir = f'PROJECT_DIR:={self.run_path}/project\n'
            filepath = f'{self.data_path}/Makefile' 
            with open(filepath, 'r') as file:    
                data = file.readlines()            
                for i, line in enumerate(data):                
                    if 'PROJECT_DIR:=' in line:
                        data[i] = run_project_dir
                    
            self.write_data(filepath, data)
            
            shutil.copytree(self.template_path, self.run_path, ignore = shutil.ignore_patterns("*presn*", '.git', 'docs', 
                                                                                               'mkdocs', 'examples', 'legacy',
                                                                                               'papers'))
            
    def check_path(self, path):
        if not os.path.exists(path):
            os.makedirs(path)
            
    def write_data(self, filepath, data):
        with open(filepath, 'w') as file:   
            file.writelines(data)     
            
    def find_last_dump(self):
        outfiles = [filename for filename in os.listdir(f'{self.full_output_path}') if "restart" in filename]
        if any("restart" in file for file in outfiles):             
            last_num   = max([int(filename.split('_')[-1]) for filename in outfiles])            
            input_name = f"DataOut_restart_{last_num}" 
        else:
            last_num   = 0
            input_name = "DataOut"             
            
        next_num  = last_num + 1
                
        read      = Readout(self.full_output_path, base_file = 'DataOut_read', outfile=input_name)
        last_dump = read.run_readable()
        
        self.data_in   = f'{self.full_output_path}/{input_name}'
        self.data_out  = f'{self.full_output_path}/DataOut_restart_{next_num}'
        self.read_dump = last_dump
                
        return 
            
    def prep_data(self):
        
        # Edit setup_prep
        filepath = f'{self.data_path}/prep_data/setup_prep'
        with open(filepath, 'r') as file:    
            data = file.readlines()            
            for i, line in enumerate(data):                
                if 'Input File' in line:
                    data[i+1] = f'{self.dataset}/s{self.mass}_presn\n'
                elif 'Output File' in line: 
                    data[i+1] = f'{self.data_in}\n'  
                elif 'Goal Size of the PNS' in line: 
                    data[i+1] = f'{self.pns_grid_goal}\n'
                elif 'Goal Size of Convective Grid' in line: 
                    data[i+1] = f'{self.conv_grid_goal}\n'
                elif 'Goal Total Resolution' in line: 
                    data[i+1] = f'{self.grid_goal}\n'                    
                elif 'Enclosed Mass Cutoff for the PNS' in line: 
                    data[i+1] = f'{self.enclmass_pns}\n'                                              
                elif 'Enclosed Mass Cutoff for Convective Region' in line: 
                    data[i+1] = f'{self.enclmass_conv}\n'    
                elif 'Maximum Radius of the Grid' in line: 
                    data[i+1] = f'{self.maxrad}\n'                                      
                elif 'EOS Table Path' in line: 
                    data[i+1] = f'{self.eos_table_path}\n'                                                
                                                           
        self.write_data(filepath, data)       
        
        # run 'make data'
        os.chdir(self.data_path)
        
        if not os.path.exists(f'{self.data_path}/prep_data/eosmodule.mod'):
            p = Popen('make eos', shell=True, stdin=PIPE, stdout=PIPE, stderr=PIPE)
            output = p.stdout.read()
            p.stdout.close()            
        
        p = Popen('make data', shell=True, stdin=PIPE, stdout=PIPE, stderr=PIPE)
        output = p.stdout.read()
        p.stdout.close()
        print('Data prepared')
        
    def setup(self):
                
        # Edit setup
        filepath = f'{self.run_path}/project/1dmlmix/setup' 
        with open(filepath, 'r') as file:    
            data = file.readlines()            
            for i, line in enumerate(data):    
                if 'Input File' in line: 
                    data[i+1] = f'{self.data_in}\n'                            
                elif 'Output File' in line: 
                    data[i+1] = f'{self.data_out}\n'  
                elif 'Input PyTorch Model' in line: 
                    data[i+1] = f'{self.mlmodel}\n'    
                elif 'Dump # to read' in line: 
                    data[i+1] = f'{self.read_dump}\n' 
                elif 'Dump time intervals' in line: 
                    data[i+1] = f'{self.dump_interval}\n'  
                elif 'Max time' in line: 
                    data[i+1] = f'{self.maxtime}\n'                                                           
                elif 'EOS Table Path' in line: 
                    data[i+1] = f'{self.eos_table_path}\n'                                                       
                    
        self.write_data(filepath, data)            
            
        # check output if output folder exists
        self.check_path(self.full_output_path)     
        
    def setup_readout(self):            
        
        # Edit setup_readout
        filepath = f'{self.run_path}/project/1dmlmix/setup_readout' 
        with open(filepath, 'r') as file:    
            data = file.readlines()            
            for i, line in enumerate(data):                
                if 'Data File Name' in line: 
                    data[i+1] = f'{self.full_output_path}/DataOut\n' 
                if 'Output Basename' in line:                     
                    data[i+1] = f'{self.full_output_path}/DataOut_read\n'                    
                    
        self.write_data(filepath, data)            
            
        # check output if output folder exists
        self.check_path(self.full_output_path)     
        
        
class Readout:
    def __init__(self, full_output_path, base_file, outfile):
        self.base_file        = base_file
        self.outfile          = outfile
        self.full_output_path = full_output_path

    def run_readable(self):
        
        self.setup_readout()
        
        p = Popen('./readout', shell=True, stdin=PIPE, stdout=PIPE, stderr=PIPE)
        output = p.stdout.read()
        p.stdout.close()
         
        return self.get_lastdump()
    
    def get_lastdump(self):
        alldumps = [filename for filename in os.listdir(f'{self.full_output_path}') if self.base_file in filename]
        alldumps = [int(filename.split('.')[-1]) for filename in alldumps]

        return max(alldumps)
        
    def setup_readout(self):                    
        # Edit setup
        filepath = f'setup_readout' 
        with open(filepath, 'r') as file:    
            data = file.readlines()            
            for i, line in enumerate(data):                
                if 'Data File Name' in line: 
                    data[i+1] = f'{self.full_output_path}/{self.outfile}\n' 
                if 'Output Basename' in line:                     
                    data[i+1] = f'{self.full_output_path}/{self.base_file}\n'     
                if 'Number of dumps' in line:                     
                    data[i+1] = f'10000\n'                                      
                    
        self.write_data(filepath, data)               
    
    def write_data(self, filepath, data):
        with open(filepath, 'w') as file:   
            file.writelines(data)                  
        
if __name__ == '__main__':
    main()