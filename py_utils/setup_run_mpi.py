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

class multirun:
    def __init__(self, suffixs, masses, enclmass_conv_cutoff,pns_cutoff,
                 dataset,base_path,template_path,output_path, eos_table_path, 
                 pns_grid_goals, conv_grid_goals, grid_goals,
                 maxrads, mlmodel='None',
                 read_dump=0, dump_interval=1e-3, restart=False,
                 data_names=None, maxtime=0.5):
        
        self.suffixs         = suffixs
        self.masses          = masses
        self.enclmass_conv_cutoff = enclmass_conv_cutoff
        self.pns_cutoff      = pns_cutoff
        self.dataset         = dataset
        self.base_path       = base_path
        self.output_path     = output_path
        self.data_names      = data_names
        self.template_path   = template_path
        self.data_path       = self.template_path#f'{self.base_path}/template'#produce_data'
        self.eos_table_path  = eos_table_path
        self.maxrads         = maxrads
        self.mlmodel         = mlmodel
        self.pns_grid_goals  = pns_grid_goals
        self.conv_grid_goals = conv_grid_goals
        self.grid_goals      = grid_goals
        self.read_dump       = read_dump
        self.dump_interval   = dump_interval
        self.maxtime         = maxtime
        self.restart         = restart
        
    def setup_pars(self, i):
        self.mass          = self.masses[i]
        self.enclmass_conv = self.enclmass_conv_cutoff[i]
        self.enclmass_pns  = self.pns_cutoff[i]

        self.pns_grid_goal  = self.pns_grid_goals[i]                                        
        self.conv_grid_goal = self.conv_grid_goals[i]
        self.grid_goal      = self.grid_goals[i]
        self.maxrad         = self.maxrads[i]
        self.suffix         = self.suffixs[i]
        
        if self.data_names != None: self.data_in = self.data_names[i]
        
        self.run_name         = f's{self.mass}{self.suffix}'
        self.run_path         = f'{self.base_path}/{self.run_name}'
        self.full_output_path = f'{self.output_path}/{self.run_name}'   
             
    def run(self, rank):                     
        
        if rank == 0: colored.head('\n<<<<<<<< Running Simulations >>>>>>>>')   
        # run 'make project'
        os.chdir(f'{self.run_path}')

        if not self.restart:
            
            print(f'Rank {rank}: Compiling...')   
            
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
            if counter == 120: colored.error("executable '1dmlmix' not found (waited 2 mins)")
        print(f'rank {rank} prepared {self.run_name}; running...')
        
        p = Popen(f'time ./1dmlmix > {stdout} 2> {stderr}', shell=True, stdin=PIPE, stdout=PIPE, stderr=PIPE)
        output = p.stdout.read()
        p.stdout.close()
        
        colored.subhead('--------------------------------------------')
        colored.subhead(f'rank {rank} finished running {self.run_name}')
        colored.subhead('--------------------------------------------')
        
            
    def initialize(self, rank):
        
        if self.restart:
            # Restarts in parallel utilizing all given ranks
            if rank == 0: colored.head('<<< Converting Binary to Readable >>>') 
            self.setup_pars(rank)            
            self.find_last_dump()
            self.setup()
            self.setup_readout()
            print(f'Rank',f'{rank}'.ljust(2, ' '),
                  f'{self.run_name} restarts from dump: {self.read_dump}')
        else:
            # Initializes all fresh runs in serial
            for i in range(len(self.masses)):
                self.setup_pars(i)
                                    
                colored.subhead(f'--- {self.run_name} ---')
                
                if os.path.exists(self.full_output_path): 
                    valid_input = False
                    yes, no     = ['y', 'yes'], ['n', 'no']
                    while not valid_input:
                        overwrite = input("Output exists. Overwrite? [Y/N]")
                        if overwrite.lower() in yes+no: valid_input = True
                        else: print(f"Invalid input: '{overwrite}'")
                        
                    if overwrite.lower() in no: 
                        colored.warn(f'Aborted: {self.run_name}')
                        continue
                
                print(f'Mass:                   {self.mass}')
                print(f'Enclosed PNS  Cutoff:   {self.enclmass_pns}')
                print(f'Enclosed Conv Cutoff:   {self.enclmass_conv}')                                                    
                            
                # check if run_folder exists; create and copy template if not            
                # also prep data and copy it to the run folder
                if self.data_names==None: self.data_in = 'Data'
                self.prep_data()                
                self.data_out = f'{self.full_output_path}/DataOut'
                
                # edit 'setup' to include unique output path
                self.setup()
                
                # edit 'setup_readout' to include unique output path
                self.setup_readout()     
                                                
            colored.head('<<<< Initialization Completed >>>>')                           

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
        
        # check if run_folder exists; create and copy template if not
        self.edit_copy_template()      
        
        # run 'make data'
        os.chdir(self.data_path)
        
        if not os.path.exists(f'{self.data_path}/prep_data/eosmodule.mod'):
            p = Popen('make eos', shell=True, stdin=PIPE, stdout=PIPE, stderr=PIPE)
            output = p.stdout.read()
            p.stdout.close()            
        
        p = Popen('make data', shell=True, stdin=PIPE, stdout=PIPE, stderr=PIPE)
        output = p.stdout.read()
        p.stdout.close()
        
        p = Popen('make clean', shell=True, stdin=PIPE, stdout=PIPE, stderr=PIPE)
        output = p.stdout.read()
        p.stdout.close()
        
        print('Data prepared')

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
                        
        read      = Readout(self.run_path, self.full_output_path,
                            base_file = 'DataOut_read', outfile=input_name)
        last_dump = read.run_readable()
        
        self.data_in   = f'{self.full_output_path}/{input_name}'
        self.data_out  = f'{self.full_output_path}/DataOut_restart_{next_num}'
        self.read_dump = last_dump
                
        return             
        
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
    def __init__(self, run_path, full_output_path, base_file, outfile):
        self.base_file        = base_file
        self.outfile          = outfile
        self.run_path         = run_path
        self.sim_path         = f'{self.run_path}/project/1dmlmix'
        self.full_output_path = full_output_path

    def run_readable(self):
        
        self.setup_readout()
        
        os.chdir(f'{self.sim_path}')
        
        if not os.path.isfile('readout'): colored.error("readout executable doesn't exist")
        
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
        filepath = f'{self.run_path}/project/1dmlmix/setup_readout' 
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
            
class colored:
    RED    = '\033[31m'
    GREEN  = '\033[32m'
    YELLOW = '\033[33m'
    ORANGE = '\033[34m'
    PURPLE = '\033[35m' 
    CYAN   = "\033[36m"
    RESET  = "\033[0m"
    
    @classmethod
    def head(cls, message): print(cls.CYAN+f"{message}"+cls.RESET) 
    
    @classmethod
    def subhead(cls, message): print(cls.PURPLE+f"{message}"+cls.RESET)   
    
    @classmethod
    def warn(cls, message): print(cls.YELLOW+f"WARNING: {message}"+cls.RESET)

    @classmethod        
    def error(cls, message): sys.exit(cls.RED+f"ERROR: {message}"+cls.RESET)            