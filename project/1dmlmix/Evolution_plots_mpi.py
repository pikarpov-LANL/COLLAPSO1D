# This is an analysis script that can:
#     1. convert binary output to readable [.txt]
#     2. produce time evolution plots for selected variables
#     3. produce summary plots of 
#       - the convection region grid size post-bounce
#       - time evolution of electron neutrinos flux
#     4. combines all of the plots into movies with ffmpeg
#
# and all of this in parllel with MPI. For examples, to run on 4 cores:
#
# mpirun -n 4 python Evolution_plots_mpi.py
#    
# -pikarpov

import os
import sys
from subprocess import Popen, PIPE
from mpi4py import MPI

sys.path.append("/home/pkarpov/Sapsan")

import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
from sapsan.utils import line_plot, plot_params

def main():
    
    # --- Datasets and values to plot ---        
    vals             = [
                        'rho', 
                        # 'v',
                        # 'P',
                        # 'T',
                        # 'enclmass',
                        # 'vsound',
                       ]        
    versus           = 'r'   # options are either 'r' or 'encm' for enclosed mass
    masses           = [11.0,12.0,13.0,14.0,15.0,
                        16.0,17.0,18.0,19.0,20.0]        
    #masses           = [12.0,16.0,19.0,20.0]

    # --- Paths & Names ---
    datasets         = [f's{m}_g2k_c1k_p0.4k' for m in masses]
    base_file        = f'DataOut_read'
    base_path        = '/home/pkarpov/scratch/1dccsn/sfho_s/'
    #base_path        = '/home/pkarpov/scratch/1dccsn/sleos/funiek/'
    #base_path        = '/home/pkarpov/COLLAPSO1D/project/1dmlmix/output/'
    save_name_amend  = ''      # add a custom index to the saved plot names
    
    # --- Extra ---
    convert2read     = False#True    # convert binary to readable (really only needed to be done once) 
    only_post_bounce = True   # only produce plots after the bounce    
    
    # --- Compute Bounce Time, PNS & Shock Positions ---
    compute          = False
    rho_threshold    = 1e13    # for the PNS radius - above density is considered a part of the Proto-Neutron Star

    # --- Plots & Movie Parameters ---
    dpi              = 60      # increase for production plots
    make_movies      = False#True 
    fps              = 10   
    save_plot        = True    # shouldn't be touched for the mpi routine
    
    # === No need to go beyond this point ===========================

    # --- MPI setup ---
    comm = MPI.COMM_WORLD
    size = comm.Get_size()
    rank = comm.Get_rank()
    
    for dataset in datasets:
        
        # --- Create initial assignment ---
        if rank == 0:                              
            print(f'<<<<<<<<< {dataset} >>>>>>>>>')
            
            numfiles = get_numfiles(base_path, dataset, base_file)
            
            if convert2read: numfiles = Readout(base_path, dataset, base_file).run_readable()

            print( '\n--------- Summary ---------')
            print(f'Total number of files:  {numfiles}') 
            
            last_file = numfiles 
                    
            pf = Profiles(rank = rank, numfiles = numfiles, 
                          base_path = base_path, base_file = base_file, dataset = dataset,
                          save_name_amend=save_name_amend, only_post_bounce = only_post_bounce)
            
            bounce_files = pf.check_bounce(compute=compute)
            bounce_shift = pf.bounce_ind
            if only_post_bounce: numfiles = bounce_files
            
            print(f'Bounce at file:         {bounce_shift+1}')
            print(f'Post bounce files:      {bounce_files}')
            
            interval_size = int(numfiles/size)
            leftover = numfiles%size
            
            print(f'Average interval size:  {interval_size}')
            print(f'Leftover to distribute: {leftover}')
            print( '\n-------- Intervals --------')

            interval = np.array([[i*interval_size,i*interval_size+interval_size] for i in range(size)])
            
            if leftover != 0: interval = spread_leftovers(interval, leftover, size)
            
            if only_post_bounce: interval+=bounce_shift
            
            # creates (if needed) directories to store all plots
            if save_plot: [pf.set_paths(val, versus, check_path=True) for val in vals]
            
        else:
            last_file    = 0
            interval     = 0   
            bounce_shift = 0   

        numfiles = comm.bcast(last_file, root=0)
        interval = comm.scatter(interval, root=0)                

        pf = Profiles(rank = rank, numfiles = numfiles, 
                      base_path = base_path, base_file = base_file, dataset = dataset,
                      save_name_amend=save_name_amend, only_post_bounce = only_post_bounce,                    
                      interval = interval, dpi = dpi)
        
        pf.bounce_ind = comm.bcast(bounce_shift, root=0)

        print(f'rank {rank}, interval {interval}')

        comm.Barrier()
        if rank == 0: print( '\n-------- Progress ---------', flush=True)

        # --- Main Parallel Loop ---
        for i in range(interval[0], interval[1]):  
            pf.plot_profile(i             = i, 
                            vals          = vals, 
                            versus        = versus, 
                            show_plot     = False, 
                            save_plot     = save_plot,
                            compute       = compute,
                            rho_threshold = rho_threshold
                           )     

        pf.progress_bar(i+1, 'Done!', done = True)   

        gather_pns_ind    = comm.gather(pf.pns_ind_ar,    root=0)
        gather_pns_x      = comm.gather(pf.pns_x_ar,      root=0)
        gather_pns_encm   = comm.gather(pf.pns_encm_ar,   root=0)
        gather_shock_ind  = comm.gather(pf.shock_ind_ar,  root=0)
        gather_shock_x    = comm.gather(pf.shock_x_ar,    root=0)
        gather_shock_encm = comm.gather(pf.shock_encm_ar, root=0)
        gather_lumnue     = comm.gather(pf.lumnue,        root=0)
        
        # --- Back to Rank 0 to produce Summary Plots & Movies ---
        if rank == 0: 
            if save_plot:
                print( '\n--------- Plot Path --------', flush=True)
                print(pf.base_save_path) 
                    
            pf.pns_ind_ar    = sum(gather_pns_ind)
            pf.pns_x_ar      = sum(gather_pns_x)
            pf.pns_encm_ar   = sum(gather_pns_encm)
            pf.shock_ind_ar  = sum(gather_shock_ind)
            pf.shock_x_ar    = sum(gather_shock_x)
            pf.shock_encm_ar = sum(gather_shock_encm)
            pf.lumnue        = sum(gather_lumnue)            
            pf.bounce_ind    = bounce_shift
            
            ax = pf.plot_convection()
            ax = pf.plot_lumnue()
            ax = pf.plot_pns_shock()
                    
            if make_movies:
                print( '\n---------- Movies ----------')
                for val in vals: pf.movie(val,fps=fps)
                 
            print(f'<<<<<<<<<<< Done >>>>>>>>>>>\n')                

# === Backend ===============================================
                        
def spread_leftovers(interval, leftover, size):    
    shift = 0
    for i in range(size):        
        interval[i,0] += shift
        if leftover != 0:
            shift    += 1
            leftover -= 1             
        interval[i,1] += shift
    return interval

def get_numfiles(base_path, dataset, basefile):
    numfiles = len([filename for filename in os.listdir(f'{base_path}{dataset}') if basefile in filename])
    return numfiles

class Readout:
    def __init__(self, base_path, dataset, base_file):
        self.base_path        = base_path
        self.dataset          = dataset
        self.base_file        = base_file
        self.full_output_path = f'{self.base_path}{self.dataset}'

    def run_readable(self):
        print( '---- Converting Binary ----')           
        
        self.setup_readout()
        
        p = Popen('./readout', shell=True, stdin=PIPE, stdout=PIPE, stderr=PIPE)
        output = p.stdout.read()
        p.stdout.close()
        print(f'Files are now readable for {self.dataset}!')
                
        return get_numfiles(self.base_path, self.dataset, self.base_file)
        
    def setup_readout(self):            
        
        # Edit setup
        filepath = f'setup_readout' 
        with open(filepath, 'r') as file:    
            data = file.readlines()            
            for i, line in enumerate(data):                
                if 'Data File Name' in line: 
                    data[i+1] = f'{self.full_output_path}/DataOut\n' 
                if 'Output Basename' in line:                     
                    data[i+1] = f'{self.full_output_path}/{self.base_file}\n'     
                if 'Number of dumps' in line:                     
                    data[i+1] = f'10000\n'                                      
                    
        self.write_data(filepath, data)               
    
    def write_data(self, filepath, data):
        with open(filepath, 'w') as file:   
            file.writelines(data)                  
        
   
class ComputeRoutines:
    #
    # Routines to calculate PNS radius and shock position
    #
    def __init__(self, x, rho, v, vsound=None):
        self.x          = x
        self.rho        = rho
        self.v          = v  
        self.vsound     = vsound      
        self.shock_ind  = 0
        self.shock_x    = 0
        self.pns_ind    = 0
        self.pns_x      = 0        
        
    def shock_radius(self):        
                        
        mach = abs(self.v/self.vsound)
        
        for i in range(np.argmin(self.v),-1,-1):            
            if mach[i] < 1:
                self.shock_ind = i
                self.shock_x   = self.x[self.shock_ind]
                break
    
        #print('shock position: %.2e'%self.shock_x, self.shock_ind)
        return self.shock_ind, self.shock_x
        
    def pns_radius(self, rho_threshold = 2e11):      
          
        for i in range(len(self.rho)):            
            if self.rho[i] > rho_threshold:
                self.pns_ind = i
                self.pns_x   = self.x[i]
                
        return self.pns_ind, self.pns_x        

class Profiles:
    #
    # All things plotting related (+ bounce check)
    #
    def __init__(self, rank, numfiles, base_path, base_file, dataset, 
                 save_name_amend='', only_post_bounce = False, interval=[0,0], dpi=60):
        self.numfiles         = numfiles
        self.lumnue           = np.zeros((self.numfiles))
        self.times            = np.zeros((self.numfiles))
        self.base_path        = base_path
        self.dataset          = dataset
        self.base_save_path   = f'{self.base_path}{self.dataset}/plots/'
        self.movie_save_path  = f'{self.base_path}{self.dataset}/movies/'
        self.save_name_amend  = save_name_amend
        self.only_post_bounce = only_post_bounce
        self.interval         = interval
        self.dpi              = dpi
        self.cm2km            = 1e-5
        #self.progress_bar(0)
        
        self.pns_ind_ar       = np.zeros(self.numfiles)
        self.pns_x_ar         = np.zeros(self.numfiles)
        self.pns_encm_ar      = np.zeros(self.numfiles)        
        self.shock_ind_ar     = np.zeros(self.numfiles)
        self.shock_x_ar       = np.zeros(self.numfiles)
        self.shock_encm_ar    = np.zeros(self.numfiles)
        self.time_ar          = np.zeros(self.numfiles)
        self.ind_ar           = np.arange(1, self.numfiles+1)
        self.bounce_ind       = 0
        self.rank             = rank
        self.base_file        = base_file 
        
    def progress_bar(self, current, val='', done = False, bar_length=20):
        current -= self.interval[0]
        lastfile = self.interval[1]-self.interval[0]
        fraction = current / lastfile

        arrow       = int(fraction * bar_length - 1) * '-' + '>'
        padding     = int(bar_length - len(arrow)) * ' '
        padding_val = int(8-len(val)) * ' '

        ending = '\n' if done == True else '\r'

        print(f'rank {self.rank}: [{arrow}{padding}] {current}/{lastfile}, val: {val}{padding_val}', end=ending) 

    def set_paths(self, val, versus, check_path=False):
        self.plot_path = f'{self.base_save_path}{val}'
        
        if check_path: 
            if not os.path.exists(self.plot_path): os.makedirs(self.plot_path)
            
        if   versus == 'r'   : self.versus_name = '_r'
        elif versus == 'encm': self.versus_name = '_encm'    
        
        self.plot_file = f'{self.plot_path}/{val}{self.versus_name}'

    def check_bounce(self, compute=False):     
        self.bounce_ind = -1
        bounced = 0           
        for i in range(self.numfiles):
            file   = f'{self.base_file}.{i+1}'
            file1d = f'{self.base_path}{self.dataset}/{file}' 
                      
            with open(file1d, "r") as file:
                line        = file.readline()        
                header_vals = file.readline()
                vals_strip  = header_vals[:-1].split(' ')        
                time1d, bounce_time, pns_ind, pns_x, shock_ind, shock_x, rlumnue = [float(x) for x in vals_strip if x!='']        
                
            pns_ind = int(pns_ind)-1              
            
            if compute: 

                ps = np.genfromtxt(file1d, skip_header=3)                
                ps = np.moveaxis(ps,0,1) 
                rho = ps[3]
                if np.amax(rho) > 2e14:
                    if bounced==0: bounced = time1d 
                    if (time1d-bounced) > 2e-3:
                        self.bounce_ind = i
                        return self.numfiles - self.bounce_ind 
                     
            elif shock_ind > 0: 
                self.bounce_ind = i
                return self.numfiles - self.bounce_ind
                                                                                
        print("WARNING: Bounce has not been found :(")
        return -1        
        
    def plot_format(self, series, xlabel, ylabel, title, save_path, 
                          show_plot, plot_style = 'plot',
                          bounce_lim=False, label=None):                       
        
        style = 'tableau-colorblind10'
        mpl.style.use(style)
        mpl.rcParams.update(plot_params())  
        
        if not label: label = ["Data {}".format(i) for i in range(len(series))]
                    
        fig = plt.figure(figsize=(10,6), dpi=self.dpi)
        ax  = fig.add_subplot(111)        
        for idx, data in enumerate(series): 
            if   plot_style == 'plot'    : plot_func = ax.plot
            elif plot_style == 'semilogx': plot_func = ax.semilogx
            elif plot_style == 'semilogy': plot_func = ax.semilogy
            elif plot_style == 'loglog'  : plot_func = ax.loglog  
                      
            plot_func(data[0], data[1], linewidth=1.5, marker='.', label=label[idx])
        
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        ax.xaxis.set_major_locator(mpl.ticker.MaxNLocator(integer=True))
        
        if bounce_lim:
            if self.only_post_bounce: ax.set_xlim(self.bounce_ind,self.numfiles)
            
        plt.tight_layout()
        
        plt.savefig(save_path)        
        if not show_plot: plt.close()        
        print(save_path, flush=True)
         
        return ax            
        
    def plot_convection(self, show_plot=False):        
        print( '\n-------- Convection --------')  
                 
        save_path = f'{self.base_save_path}{self.save_name_amend}convgrid.png'
        
        ax = self.plot_format(series     = [[self.ind_ar, self.shock_ind_ar-self.pns_ind_ar]],
                              xlabel     = 'index', 
                              ylabel     = 'Convection Grid Size',
                              title      = f'Bounce index = {self.bounce_ind+1}', 
                              save_path  = save_path, 
                              show_plot  = show_plot,
                              bounce_lim = True)                        
        return ax
    
    def plot_lumnue(self, show_plot=False):                  
        print( '\n---------- lumnue ----------')                   
        
        if self.only_post_bounce: name_amend = '_bounce'            
        else: name_amend = ''
        
        save_path = f'{self.base_save_path}{self.save_name_amend}lumnue{name_amend}.png'
        
        ax = self.plot_format(series     = [[self.ind_ar, self.lumnue]],
                              xlabel     = 'index', 
                              ylabel     = r'$F_{\nu_{e}} \; [foe/s]$', 
                              title      = f'Bounce index = {self.bounce_ind+1}',
                              save_path  = save_path, 
                              show_plot  = show_plot,
                              plot_style = 'semilogy',
                              bounce_lim = True)          
        return ax   
    
    def plot_pns_shock(self, show_plot=False):                  
        print( '\n--------- pns_shock --------')  
                
        plot_data = [[self.pns_x_ar[self.bounce_ind:]*self.cm2km,
                      self.pns_encm_ar[self.bounce_ind:]         ],
                     [self.shock_x_ar[self.bounce_ind:]*self.cm2km, 
                      self.shock_encm_ar[self.bounce_ind:]       ]]
        
        labels    = ['pns', 'shock']   
        
        save_path = f'{self.base_save_path}{self.save_name_amend}pns_shock.png'
        
        ax = self.plot_format(series    = plot_data,
                              xlabel    = r'$Radius \; [km]$', 
                              ylabel    = r'$M_{enc} \; [M_{sol}]$', 
                              title     = f'Bounce index = {self.bounce_ind+1}',
                              save_path = save_path, 
                              show_plot = show_plot,
                              label     = labels) 
        ax.legend(loc=0) 
                   
        return ax        
                    
    def plot_profile(self, i, vals, versus,
                     show_plot=False, save_plot=False, 
                     compute=False, rho_threshold = 2e11):
        file = f'{self.base_file}.{i+1}'
        file1d = f'{self.base_path}{self.dataset}/{file}'                        

        with open(file1d, "r") as file:
            line = file.readline()        
            header_vals = file.readline()
            vals_strip  = header_vals[:-1].split(' ')        
            time1d, bounce_time, pns_ind, pns_x, shock_ind, shock_x, rlumnue = [float(x) for x in vals_strip if x!='']        
            self.lumnue[i] = rlumnue
            self.times[i]  = time1d    
            #print(file.readline())
            
        pns_ind   = int(pns_ind)-1
        shock_ind = int(shock_ind)-1
                        
        #print('Time %.2f ms'%(float(time1d)*1e3))

        ps = np.genfromtxt(file1d, skip_header=3)
        ps = np.moveaxis(ps,0,1)    
        
        # Columns in the readable files:
        # Cell  M_enclosed [M_sol]  Radius [cm]  Rho [g/cm^3]  Velocity [cm/s]  Ye  Pressure [g/cm/s^2]  Temperature [K]  Sound [cm/s] 
        ncell   = ps[0]
        encm    = ps[1]
        r       = ps[2]
        rho     = ps[3]
        v       = ps[4]
        ye      = ps[5]
        P       = ps[6]
        T       = ps[7]
        vsound  = ps[8]
        
        # Print out enclosed mass at 200km without plotting anything
        # if i == self.numfiles-1:         
        #     for j in range(len(cell)):
        #         if r[j] >= 2e7:
        #             print(f'at {r[j]:.2e}, mass is {encm[j]:.3e}')
        #             return
        # else: return
                    
        for val in vals:
            if versus=='encm':
                x         = encm
                xlabel    = r'$M_{enc} \; [M_{sol}]$'
                xlim      = (0,3) #or (0,4)
                plot_type = 'semilogy'
                unit      = 1        
            elif versus=='r':
                x         = r
                xlabel    = r'$Radius \; [km]$'
                xlim      = (1e0,1e5)
                plot_type = 'loglog'
                unit      = self.cm2km
            elif versus=='Placeholder':
                x         = r
                xlabel    = r'$(V/V_{sound})^2$'
                xlim      = None
                plot_type = 'plot'
                unit      = 1
            else: print("ERROR: unknown 'versus' {versus}, trying to exit"); sys.exit()
            
            if val == 'rho': 
                y      = rho
                ylabel = r'Density $[g/cm^3]$'
                ylim   = (1e4,1e15)
                loc    = 1
            elif val == 'v':
                y      = v
                ylabel = r'Velocity $[cm/s]$'
                ylim   = (-1e10, 1e9)
                loc    = 4
                if versus == 'r': plot_type = 'semilogx'
                elif versus == 'encm': plot_type = 'plot'                
            elif val == 'vsound':
                y      = vsound
                ylabel = r'$V_{sound}$ $[cm/s]$'
                ylim   = (0, 1.4e10)
                loc    = 1
                if versus == 'r': plot_type = 'semilogx'
                elif versus == 'encm': plot_type = 'plot'                                
            elif val == 'P':
                y      = P
                ylabel = r'$P_{gas} \; [\frac{g}{cm\;s^2}]$'
                ylim   = (1e20,1e36)   
                loc    = 1
            elif val == 'T':
                y      = T
                ylabel = r'Temperature $[K]$'
                ylim   = (1e7,4e11)
                loc    = 1
            elif val == 'enclmass':
                y      = encm
                ylabel = r'Enclosed Mass $[M_{\odot}]$'
                ylim   = None
                loc    = 4
            else: sys.exit(f"ERROR: unknown val {val}, trying to exit")            

            #print('diff shock', shock_ind, np.argmin(ps[4]), shock_ind-np.argmin(ps[4]))
            #shock_ind = np.argmin(ps[4])

            ax = line_plot([[x*unit, y],],
                            #[ps[2,:504]*1e0, ps[ps_ind,:504]],
                           plot_type = plot_type,
                           label     = [f'ind     {i+1}'],
                           linestyle = ['-','--','-','--'],               
                           figsize   = (10,6), 
                           dpi=self.dpi
                           )    
            
            # check if after bounce                
            if i >= self.bounce_ind:                
                
                if compute and vals.index(val)==0:
                    rt = ComputeRoutines(x, rho=rho, v=v, vsound=vsound)
                    shock_ind, shock_x = rt.shock_radius()
                    pns_ind,   pns_x   = rt.pns_radius(rho_threshold = rho_threshold)   
                    
                #print(f'shock {ps[2,int(shock_ind)]:.3e}, {shock_x:.3e}')
                
                pns_edge    = pns_x*unit
                shock_front = x[int(shock_ind)]*unit
                if versus == 'r'   : line_label = '%.2e km'          
                if versus == 'encm': line_label = '%.3f $M_{sol}$'      
                         
                ax.axvline(x=pns_edge,linestyle='-',color='r',linewidth=1,
                           label=f'PNS    {line_label%pns_edge}')
                ax.axvline(x=shock_front,linestyle='--',color='r',linewidth=1,
                           label=f'shock {line_label%shock_front}')
                
            ax.set_xlabel(xlabel)
            ax.set_ylabel(ylabel)
            ax.set_xlim(xlim)
            ax.set_ylim(ylim)
            ax.set_title('$t_{1d}=$%.2f ms'%(float(time1d)*1e3))

            plt.legend(loc=loc)
            plt.tight_layout()
            if save_plot:        
                self.set_paths(val, versus)                
                plt.savefig(f'{self.plot_file}{self.save_name_amend}_{i+1}.png')
                
            if not show_plot: plt.close()
            
            if vals.index(val)==0:
                self.pns_ind_ar[i]    = pns_ind
                self.pns_x_ar[i]      = pns_x
                self.pns_encm_ar[i]   = encm[pns_ind]                
                self.shock_ind_ar[i]  = shock_ind
                self.shock_x_ar[i]    = shock_x
                self.shock_encm_ar[i] = encm[shock_ind]
                self.time_ar[i]       = time1d 
            
            done=True if (i==self.numfiles and vals.index(val)==(len(vals)-1)) else False
            if self.rank == 0: self.progress_bar(i+1, val, done = done)          
        return
    
    def movie(self, val, versus='r', fps=15, start=0, printout=False):      
        self.set_paths(val,versus)
        padding_val = int(8-len(val)) * ' '
        
        if self.only_post_bounce:
            name_amend = '_bounce'
            start      = self.bounce_ind
        else: name_amend = ''
            
        name = f'{self.plot_file}{self.save_name_amend}'
        
        movie_name = f'{self.movie_save_path}{val}{self.versus_name}{self.save_name_amend}'
        if not os.path.exists(self.movie_save_path): os.makedirs(self.movie_save_path)
            
        print(f'{val}{padding_val}: {movie_name}{name_amend}.mp4')
        
        result = Popen(['ffmpeg', '-r', f'{fps}', '-start_number', f'{start}',
                        '-i', f'{name}_%d.png', 
                        '-vcodec', 'libx264', f'{movie_name}{name_amend}.mp4', '-y'],
                        stdin=PIPE, stdout=PIPE, stderr=PIPE)   
        
        output, error = result.communicate()
        if printout: print(output, error)  
        return    
    
if __name__=='__main__':
    main()
