# This is an analysis script that can:
#     1. convert binary output to readable .txt
#     2. produce time evolution plots for selected variables
#     3. produce a plot of the convection region grid size post-bounce
#     4. combines all of the plots into movies with ffmpeg
    
# -pikarpov

import os
import sys
from subprocess import Popen, PIPE
from mpi4py import MPI

sys.path.append("/home/pkarpov/Sapsan")

import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
from sapsan.utils import line_plot, slice_plot, pdf_plot, cdf_plot

def main():
    vals = [
            'rho', 
            'v',
            'P',
            'T',
            #'mach'
           ]
    versus = 'r'
    #versus = '(v/vsound)^2'
    base_path = '/home/pkarpov/scratch/1dccsn/'
    #base_path = '/home/pkarpov/COLLAPSO1D/project/1dmlmix/output/'
    dataset = 's19.0_4k'
    base_file = f'DataOut_read'
    save_name_amend = ''
    convert2read = False
    only_post_bounce = True
    save_plot = True,
    compute = False,
    rho_threshold = 1e13
    make_movies = True    
        
    # -----------------------------------------------------------------------------------

    comm = MPI.COMM_WORLD
    size = comm.Get_size()
    rank = comm.Get_rank()
        
    if rank == 0:                              
        
        numfiles = get_numfiles(base_path, dataset, base_file)
        
        if convert2read: numfiles = readout(base_path, dataset, base_file).run_readable()

        print( '---------- Files ----------')             
        print(f'Total number of files:  {numfiles}') 
        
        last_file = numfiles 
                   
        pf = Profiles(rank = rank, numfiles = numfiles, 
                        base_path = base_path, base_file = base_file, dataset = dataset,
                        save_name_amend=save_name_amend, only_post_bounce = only_post_bounce)
        
        bounce_files = pf.check_bounce()
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
        
    else:
        last_file = 0
        interval = 0        
          
    numfiles = comm.bcast(last_file, root=0)
    interval = comm.scatter(interval, root=0)        

    pf = Profiles(rank = rank, numfiles = numfiles, 
                  base_path = base_path, base_file = base_file, dataset = dataset,
                  save_name_amend=save_name_amend, only_post_bounce = only_post_bounce,                    
                  interval = interval)

    print(f'rank {rank}, interval {interval}')

    comm.Barrier()
    if rank == 0: print( '\n-------- Progress ---------', flush=True)

    for i in range(interval[0], interval[1]):  
        #i = 247 #hrg looks good! from 50 cells to 300+ while res diffs from ~1800 to ~3600
        #versus can be either 'r' or 'encm'
        pf.plot_profile(i = i, vals = vals, 
                        versus = versus, 
                        show_plot=False, 
                        save_plot=save_plot,
                        compute = compute,
                        rho_threshold = rho_threshold
                       )     

    pf.progress_bar(i+1, 'Done!', done = True)   
        
    gather_pns = comm.gather(pf.pns_ind_ar, root=0)
    gather_shock = comm.gather(pf.shock_ind_ar, root=0)
    
    if rank == 0: 
        if save_plot:
            print( '\n--------- Plot Path --------', flush=True)
            print(pf.base_save_path) 
                  
        pf.pns_ind_ar = sum(gather_pns)
        pf.shock_ind_ar = sum(gather_shock)
        
        pf.bounce_ind = bounce_shift
        ax = pf.plot_convection()
                
        if make_movies:
            print( '\n---------- Movies ----------')
            for val in vals: pf.movie(val,fps=2) 
                    
def spread_leftovers(interval, leftover, size):    
    shift = 0
    for i in range(size):        
        interval[i,0] += shift
        if leftover != 0:
            shift+=1
            leftover -=1             
        interval[i,1] += shift
    return interval

def get_numfiles(base_path, dataset, basefile):
    numfiles = len([filename for filename in os.listdir(f'{base_path}{dataset}') if basefile in filename])
    return numfiles

class readout:
    def __init__(self, base_path, dataset, base_file):
        self.base_path = base_path
        self.dataset = dataset
        self.base_file = base_file
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
                    
        self.write_data(filepath, data)               
    
    def write_data(self, filepath, data):
        with open(filepath, 'w') as file:   
            file.writelines(data)                  
        
   
class routines:
    def __init__(self, x, rho, v):
        self.x = x
        self.rho = rho
        self.v = v        
        self.shock_ind = 0
        self.shock_x = 0
        self.pns_ind = 0
        self.pns_x = 0        
        
    def shock_radius(self):        
        self.shock_ind = np.argmin(self.v)
        self.shock_x = self.x[self.shock_ind]
    
        #print('shock position: %.2e'%self.shock_x, self.shock_ind)
        return self.shock_ind, self.shock_x
        
    def pns_radius(self, rho_threshold = 2e11):        
        for i in range(len(self.rho)):
            if self.rho[i] > rho_threshold:
                self.pns_ind = i
                self.pns_x = self.x[i]
        return self.pns_ind, self.pns_x        

class Profiles:
    def __init__(self, rank, numfiles, base_path, base_file, dataset, 
                 save_name_amend='', only_post_bounce = False, interval=[0,0]):
        self.numfiles = numfiles
        self.lumnue = np.zeros((self.numfiles))
        self.times = np.zeros((self.numfiles))
        self.base_path = base_path
        self.dataset = dataset
        self.base_save_path = f'{self.base_path}{self.dataset}/plots/'
        self.save_name_amend = save_name_amend
        self.only_post_bounce = only_post_bounce
        self.interval = interval
        #self.progress_bar(0)
        
        self.shock_ind_ar = np.zeros(self.numfiles)
        self.shock_x_ar = np.zeros(self.numfiles)
        self.pns_ind_ar = np.zeros(self.numfiles)
        self.pns_x_ar = np.zeros(self.numfiles)
        self.time_ar = np.zeros(self.numfiles)
        self.ind_ar = np.arange(1, self.numfiles+1)
        self.bounce_ind = 0
        self.rank = rank
        self.base_file = base_file 
        
    def progress_bar(self, current, val='', done = False, bar_length=20):
        current -= self.interval[0]
        lastfile = self.interval[1]-self.interval[0]
        fraction = current / lastfile

        arrow = int(fraction * bar_length - 1) * '-' + '>'
        padding = int(bar_length - len(arrow)) * ' '
        padding_val = int(7-len(val)) * ' '

        ending = '\n' if done == True else '\r'

        print(f'rank {self.rank}: [{arrow}{padding}] {current}/{lastfile}, val: {val}{padding_val}', end=ending) 

    def set_paths(self, val, versus):
        self.plot_path = f'{self.base_save_path}{val}'
        
        if not os.path.exists(self.plot_path): os.makedirs(self.plot_path)
        if versus=='r': self.plot_file = f'{self.plot_path}/{val}_r'
        elif versus=='encm': self.plot_file = f'{self.plot_path}/{val}_encm'    

    def movie(self, val, versus='r', fps=15, start=0, printout=False):      
        self.set_paths(val,versus)
        padding_val = int(7-len(val)) * ' '
        
        if self.only_post_bounce:
            name_amend = '_bounce'
            start = self.bounce_ind
        else: name_amend = ''
            
        name = f'{self.plot_file}{self.save_name_amend}'
            
        print(f'{val}{padding_val}: {name}{name_amend}.mp4')
        result = Popen(['ffmpeg', '-r', f'{fps}', '-start_number', f'{start}',
                        '-i', f'{name}_%d.png', 
                        '-vcodec', 'libx264', f'{name}{name_amend}.mp4', '-y'],
                        stdin=PIPE, stdout=PIPE, stderr=PIPE)   
        output, error = result.communicate()
        if printout: print(output, error)       
        
    def plot_convection(self):
        ax = line_plot([[self.ind_ar, self.shock_ind_ar-self.pns_ind_ar]])
        ax.get_legend().remove()        
        ax.set_xlim(self.bounce_ind+1,self.numfiles)
        ax.set_xlabel('index')
        ax.set_ylabel('Convection Grid Size')
        ax.set_title(f'Bounce index = {self.bounce_ind+1}')
        ax.xaxis.set_major_locator(mpl.ticker.MaxNLocator(integer=True))
        plt.tight_layout()
        save_convection_path = f'{self.base_save_path}{self.save_name_amend}convgrid.png'
        plt.savefig(save_convection_path)
        
        print( '\n-------- Convection --------')
        print(save_convection_path, flush=True) 
        
        return ax

    def check_bounce(self):                
        for i in range(self.numfiles):
            file = f'{self.base_file}.{i+1}'
            file1d = f'{self.base_path}{self.dataset}/{file}' 

            with open(file1d, "r") as file:
                line = file.readline()        
                header_vals = file.readline()
                vals_strip = header_vals[:-1].split(' ')        
                time1d, bounce_time, pns_ind, pns_x, shock_ind, shock_x, rlumnue = [float(x) for x in vals_strip if x!='']        
                
            pns_ind = int(pns_ind)-1            
            
            if shock_ind > 0: 
                self.bounce_ind = i
                return self.numfiles - self.bounce_ind
        
        
        sys.exit("ERROR: Bounce has not been found :(")
                
    def plot_profile(self, i, vals, versus,
                     show_plot=True, save_plot=False, 
                     compute=False, rho_threshold = 2e11):
        file = f'{self.base_file}.{i+1}'
        file1d = f'{self.base_path}{self.dataset}/{file}'                        

        with open(file1d, "r") as file:
            line = file.readline()        
            header_vals = file.readline()
            vals_strip = header_vals[:-1].split(' ')        
            time1d, bounce_time, pns_ind, pns_x, shock_ind, shock_x, rlumnue = [float(x) for x in vals_strip if x!='']        
            self.lumnue[i] = rlumnue
            self.times[i] = time1d    
            #print(file.readline())
            
        pns_ind = int(pns_ind)-1
        shock_ind = int(shock_ind)-1
                        
        #print('Time %.2f ms'%(float(time1d)*1e3))

        ps = np.loadtxt(file1d, skiprows=3)
        ps = np.moveaxis(ps,0,1)        
        
        for val in vals:
            if versus=='encm':
                x = ps[1]
                xlabel = r'$M_{enc} \; [M_{sol}]$'
                xlim = (0,3) #or (0,4)
                plot_type = 'semilogy'
                unit = 1        
            elif versus=='r':
                x = ps[2]
                xlabel = r'$Radius \; [km]$'
                xlim = (1e0,1e5)
                plot_type = 'loglog'
                unit = 1e-5
            elif versus=='Placeholder':
                x = ps[2]
                xlabel = r'$(V/V_{sound})^2$'
                xlim = None
                plot_type = 'plot'
                unit = 1
            else: print("ERROR: unknown 'versus' {versus}, trying to exit"); sys.exit()

            #Cell M_enclosed Position Rho V Ye Pressure Temperature
            if val == 'rho': 
                y = ps[3]
                ylabel = r'Density $[g/cm^3]$'
                ylim = (1e4,1e15)
            elif val == 'v':
                y = ps[4]
                ylabel = r'Velocity $[cm/s]$'
                ylim = (-6e9, 2.5e9)
                if versus == 'r': plot_type = 'semilogx'
                elif versus == 'encm': plot_type = 'plot'
            elif val == 'P':
                y = ps[6]
                ylabel = r'$P_{gas} \; [\frac{g}{cm\;s^2}]$'
                ylim = (1e20,1e36)   
            elif val == 'T':
                y = ps[7]
                ylabel = r'Temperature $[K]$'
                ylim = (1e5,4e11)
            elif val == 'Placeholder':
                y = ps[7]
                ylabel = r'$P_{turb}/P_{gas}$'
                ylim = None
                sys.exit('Nope')
            else: sys.exit(f"ERROR: unknown val {val}, trying to exit")            

            #print('diff shock', shock_ind, np.argmin(ps[4]), shock_ind-np.argmin(ps[4]))
            #shock_ind = np.argmin(ps[4])

            ax = line_plot([[x*unit, y],],
                            #[ps[2,:504]*1e0, ps[ps_ind,:504]],
                           plot_type = plot_type,
                           label = [f'ind     {i+1}'],
                           linestyle=['-','--','-','--'],               
                           figsize=(10,6))    
            
            # check if after bounce                
            if shock_ind > 0:
                if self.bounce_ind == 0: self.bounce_ind = i                
                
                if compute and vals.index(val)==0:
                    rt = routines(x, ps[3], ps[4])
                    shock_ind, shock_x = rt.shock_radius()
                    pns_ind, pns_x = rt.pns_radius(rho_threshold = rho_threshold)   
                    
                #print(f'shock {ps[2,int(shock_ind)]:.3e}, {shock_x:.3e}')
                
                pns_edge = pns_x*unit
                shock_front = x[int(shock_ind)]*unit
                if versus == 'r': line_label = '%.2e km'          
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

            plt.legend()
            plt.tight_layout()
            if save_plot:
                self.set_paths(val, versus)                
                plt.savefig(f'{self.plot_file}{self.save_name_amend}_{i+1}.png')
                
            if not show_plot: plt.close()
            
            if vals.index(val)==0:
                self.shock_ind_ar[i] = shock_ind
                self.shock_x_ar[i] = shock_x
                self.pns_ind_ar[i] = pns_ind
                self.pns_x_ar[i] = pns_x
                self.time_ar[i] = time1d 
            
            done=True if (i==self.numfiles and vals.index(val)==(len(vals)-1)) else False
            if self.rank == 0: self.progress_bar(i+1, val, done = done)          
        return
    
if __name__=='__main__':
    main()