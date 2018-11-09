import PAA_LISA
  
input_param = {
        'plot_on':True, #If plots will be made
        'dir_savefig': False, # The directory where the figures will be saved. If False, it will be in the current working directory
        'noise_check':False,
        'home':'/home/ester/git/synthlisa/', # Home directory
        'directory_imp': False,
        'num_back': 0,
        'dir_orbits': '/home/ester/git/synthlisa/orbits/', # Folder with orbit files
        'length_calc':20, # Length of number of imported datapoints of orbit files. 'all' is also possible
        'dir_extr': 'new_1_synthlisa_armcalc', # This will be added to the folder name of the figures
        'timeunit':'Default', # The timeunit of the plots (['minutes'],['days']['years'])
        'LISA_opt':True, # If a LISA object from syntheticLISA will be used for further calculations (not sure if it works properly if this False)
        'arm_influence': True, # Set True to consider the travel time of the photons when calculating the nominal armlengths
        'tstep':False,
        'delay':True, #'Not ahead' or False
        'method':'fsolve', # Method used to solve the equation for the photon traveling time
        'valorfunc':'Function', #
        'select':'Hallion' # Select which orbit files will be imported ('all' is all)
        }

data_all,PAA_res = PAA_LISA.runfile.do_run(input_param)
