# The script prepares and launches multiple COLLPASO1D runs,
# each running independently on the cores provided. 
# Run initialization is performed in serial, but compilation
# and execution is spread between all available cores via MPI.

# -pikarpov

import numpy as np
import os
import shutil
from subprocess import Popen, PIPE
from mpi4py import MPI

def main():
    
    comm = MPI.COMM_WORLD
    size = comm.Get_size()
    rank = comm.Get_rank()    
    
    mr = multirun(suffix = '_4k',
                  masses = [19.0,20.0],
                  deltam = [9.2e-5,1.05e-4],
                  dataset = 'sukhbold2016',
                  base_path = '/home/pkarpov/runs',
                  output_path = '/home/pkarpov/scratch/1dccsn',
                  )
    
    if rank == 0:
        mr.initialize()
        if len(mr.masses) != size:
            print(f'Ranks are not distributed well!\nRank size {size} for {len(mr.masses)} datasets')
            comm.Abort()
    
    comm.Barrier()
    
    mr.run_name = f's{mr.masses[rank]}{mr.suffix}'
    mr.run_path = f'{mr.base_path}/{mr.run_name}'
    mr.sim_path = f'{mr.run_path}/project/1dmlmix'
    mr.full_output_path = f'{mr.output_path}/{mr.run_name}'
    
    mr.run(rank)
    
class multirun:
    def __init__(self, suffix, masses, deltam, dataset, 
                 base_path, output_path, data_names=None):
        self.suffix = suffix
        self.masses = masses
        self.deltam = deltam
        self.dataset = dataset
        self.base_path = base_path
        self.output_path = output_path
        self.data_names = data_names
        self.template_path = f'{self.base_path}/template'
        self.data_path = f'{self.base_path}/produce_data'
        
    def run(self, rank):
        
        # run 'make project'
        os.chdir(f'{self.run_path}')
        
        p = Popen('make project', shell=True, stdin=PIPE, stdout=PIPE, stderr=PIPE)
        output = p.stdout.read()
        p.stdout.close()
        print(f'rank {rank} compiled project {self.run_name}; running...')
        
        
        # run the simulation
        os.chdir(self.sim_path)
        
        stdout = f'{self.full_output_path}/stdout'
        if os.path.exists(stdout): os.remove(stdout) 
        
        p = Popen(f'time ./1dmlmix >> {stdout}', shell=True, stdin=PIPE, stdout=PIPE, stderr=PIPE)
        output = p.stdout.read()
        p.stdout.close()
        print('--------------------------------------------')
        print(f'rank {rank} finished running {self.run_name}')
        print('--------------------------------------------')
        
            
    def initialize(self):
        for i in range(len(self.masses)):
            self.mass = self.masses[i]
            self.dm = self.deltam[i]
            if self.data_names != None: self.data_name = self.data_names[i]
            
            self.run_name = f's{self.mass}{self.suffix}'
            self.run_path = f'{self.base_path}/{self.run_name}'
            self.full_output_path = f'{self.output_path}/{self.run_name}'
                                
            print(f'---------------- {self.run_name} ---------------')
            
            print(f'Mass:   {self.mass}')
            print(f'Deltam: {self.dm:e}')
                                
            # check if run_folder exists; create and copy template if not
            self.copy_template()
            
            # prep data and copy it to the run folder
            self.prep_data()
            
            # edit 'setup' to include unique output path
            self.setup()
            
            # edit 'setup_readout' to include unique output path
            self.setup_readout()     
            
        print('--------- Initialization Completed ---------')                                       


    def copy_template(self):    
        if not os.path.exists(self.run_path): 
            shutil.copytree(self.template_path, self.run_path)
            
    def check_path(self, path):
        if not os.path.exists(path):
            os.makedirs(path)
            
    def write_data(self, filepath, data):
        with open(filepath, 'w') as file:   
            file.writelines(data)         
            
    def prep_data(self):
        
        # Edit setup_prep
        filepath = f'{self.data_path}/prep_data/setup_prep'
        with open(filepath, 'r') as file:    
            data = file.readlines()            
            for i, line in enumerate(data):                
                if 'Input File' in line:
                    data[i+1] = f'{self.dataset}/s{self.mass}_presn\n'
                if 'Output File' in line and self.data_names!=None: 
                    data[i+1] = f'{self.data_name}\n'  
                if 'Initial Cell Mass' in line: 
                    data[i+1] = f'{self.dm:e}\n' 
                                                           
        self.write_data(filepath, data)
                    
        # Edit Makefile
        run_project_dir = f'PROJECT_DIR:={self.run_path}/project\n'
        filepath = f'{self.data_path}/Makefile' 
        with open(filepath, 'r') as file:    
            data = file.readlines()            
            for i, line in enumerate(data):                
                if 'PROJECT_DIR:=' in line:
                    data[i] = run_project_dir
                
        self.write_data(filepath, data)        
        
        # run 'make data'
        os.chdir(self.data_path)
        
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
                if 'Output File' in line: 
                    data[i+1] = f'{self.full_output_path}/DataOut\n'                    
                    
        self.write_data(filepath, data)            
            
        # check output if output folder exists
        self.check_path(self.full_output_path)     
        
    def setup_readout(self):            
        
        # Edit setup
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
        
if __name__ == '__main__':
    main()