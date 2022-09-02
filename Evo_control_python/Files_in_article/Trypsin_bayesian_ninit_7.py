#We import the necessary packages

import csv
import time
import random
import ProcessOptimizer as po # used to suggest next experiment
from ProcessOptimizer.utils import cook_estimator #Used to customize the algorithm
from ProcessOptimizer.learning.gaussian_process.kernels import Matern #A specific "brain" in the algorithm
from ProcessOptimizer import expected_minimum
import glob
import os
import pandas as pd


#Here we specify the path for files and parameters for the bayesian optimization algorithm
Tecan_running_file = 'C:\\Users\\Tecan\\Desktop\\Variable.txt' #file acting as communication between evoware and this script
recipe_folder = 'C:\\Users\\Public\\Documents\\Tecan\\' #folder where csv recipies for evoware are put
log_file = 'C:\\temp\\Trypsin_10buffers_n20_w_control_log_of_experimets.csv' #log file specifying performed experiments and results
log_file_predictions='C:\\temp\\Trypsin_10buffers_n20_w_control_log_of_p½redictions.csv' #log file specifying the predicted best recipe at each step
log_file_predictions_results='C:\\temp\\Test_log_of_prediction_results.csv' #log file specifying the predicted best result and the obtained result testing the acording recipe
batch_size = 1 #no of experiments per iteration
comparator_delta = 48 #the offset to the reference wells
space = [(0, 10)]*10 #space for input parameters
n_objectives = 1 
n_initial=7 #number of initial experiments performed by random or latin hyper cube before gausian process is used


def get_latest_file(path, fileextension):
    '''
    A function to open a path and find the newest file in the folder.
    If fileextension (e.g. 'csv') is given, then only specified files are
    searched
    Return: the absolute path to the newest file
    '''
    list_of_files = glob.glob(path + '**\\' + '*.' + fileextension, recursive=True)
    latest_file = max(list_of_files, key=os.path.getmtime)
    return latest_file

def humanize(well):
    """
    Inspired from: https://github.com/CyrusK/96-Well-Plate.py
    Given a number, return its human form
    """
    ROW = 8
    COL = 12

    # offset = (well-1)//(ROW)*i
    if well % ROW != 0:
        offset = well % ROW - 1
        colIndex = well // ROW + 1
    else:
        offset = ROW - 1
        colIndex = well // ROW

    rowIndex = chr(65 + offset)

    return rowIndex, str(colIndex)

def extract_SPARK_reading(experiment_number, latest_file, comparator_delta=None):
    '''
    A function to find, open and extract an absorbance reading from the latest
    experiment

    Parameters
    ----------
    experiment_number : Int
        Helps pinpoint position on plate to read number.
    comparator_delta : Int
        If given, we assume that a blank was run in a well. This well is assumed
        positioned as the experiment number + this given delta.
    latest_file : Str
        Absolute path to the result file.

    Returns
    -------
    Absorbance read in given well. If comparator is given, then blank reading
    is subtracted from reading

    '''

    df1 = pd.read_excel(latest_file)
    rowstoskip = df1[df1.iloc[:, 0] == "<>"].index.to_list()[0]
    footertoskip=int(df1.shape[0])-df1[df1.iloc[:,0] == "H"].index.to_list()[0]
    df = pd.read_excel(latest_file,
                       skiprows=int(rowstoskip)+1,
                       skipfooter=int(footertoskip)-1,
                       index_col='<>')
    sample_row, sample_column = humanize(experiment_number)
    absorbance_sample = df.loc[sample_row, str(sample_column)]
    if comparator_delta:
        blank_row, blank_column = humanize(experiment_number + comparator_delta)
        absorbance_blank = df.loc[blank_row, str(blank_column)]
        return absorbance_sample - absorbance_blank
    else:
        return absorbance_sample

def get_experiment_numbers(iteration, batch_size):
    '''
    Extracts a list of experiments

    Parameters
    ----------
    iteration : Int
        DESCRIPTION.
    batch_size : Ing
        DESCRIPTION.

    Returns
    -------
    A list of experiment numbers.

    '''
    start = (iteration - 1) * batch_size + 1
    stop = iteration * batch_size + 1
    return list(range(start, stop))

def read_tecan_status(link_to_control_file): # this function reads the control file
    with open(link_to_control_file) as file:
        line = file.readlines()[0]
    return line


def make_recipe_file(link_to_recipe_folder, iteration, recipe):
    
    '''
    Based on suggestion of next experiment 
    given by bayesian algorithm, 
    this function writes a CSV for evoware to use.

    Parameters
    ----------
    link_to_recipe_folder : str
        path to folder where recipe is read by evoware
    iteration : Int
        Iteration number used to specify where the reagents are pipetted to
    recipe : list
        list of sugested levels for each parameter in the input space

    Returns
    -------
    Writes a CSV that can be used by evoware to the given path

    '''    
    
    next_x = recipe

    if batch_size == 1:
        recipe_file_first_add_falcon = os.path.join(link_to_recipe_folder, 'TRYPSIN_first_adds.csv')
        with open(recipe_file_first_add_falcon, mode='w', newline='') as file:
            writer = csv.writer(file, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
            for i in range(10):
                if next_x[i] != 0:
                    writer.writerow(['SOURCE4',str(i+1),'DEST1', str(iteration), str(next_x[i])]) #Buffers 1-10 in the "real" experiment
            for i in range(10):
                if next_x[i] != 0:
                    writer.writerow(['SOURCE4',str(i+1),'DEST1', str(iteration+comparator_delta), str(next_x[i])]) #Buffers 1-10 in the reference experiment
            writer.writerow(['SOURCE4','21','DEST1',str(iteration),str(100-(sum(next_x)))]) #Vol(water until 100 uL)
            writer.writerow(['SOURCE4','23','DEST1',str(iteration),str(10)]) #Add Trypsin solution until 38 µg/mL
            writer.writerow(['SOURCE4','24','DEST1',str(iteration),str(20)]) #Add substrate until 0.5 mM
            #reference samples
            writer.writerow(['SOURCE4','22','DEST1',str(iteration+comparator_delta),str(100-(sum(next_x))+10)]) #Vol(water until 110uL because trypsin is not added to reference)
            writer.writerow(['SOURCE4','24','DEST1',str(iteration+comparator_delta),str(20)]) #Add substrate until 0.5 mM



def make_log(link_to_log_file, iteration, latest_exp, list_of_transformed_meassurements):
    ''' 
    This function takes the latest recipe and the acording result and appends this to a log-file

    Parameters
    ----------
    link_to_log_file : str
        path to log file
    iteration : int
        iteration number in experiment
    latest_exp : list 
        list of latest performed recipes with entries for each input parameter
    list_of_transformed_meassurements : list
        list of transformed meassurements from the latest experiment
    ''' 
    with open(link_to_log_file, mode='a', newline='') as file:
        writer = csv.writer(file, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)

        if iteration == 1:
            writer.writerow([str(i) for i in list(range(len(space) + n_objectives))])
        elif iteration > 1 and batch_size == 1:
            writer.writerow(latest_exp+list_of_transformed_meassurements)
        elif iteration > 1 and batch_size > 1:
            for i in range(batch_size):
                temp_list = latest_exp[i][:]
                temp_list.append(list_of_transformed_meassurements[i])
                writer.writerow(temp_list)
        else:
            print('Error # 42')


def update_tecan_status(link_to_control_file, status): #
    '''
     this function is used to update the control file that is used for communication between evoware and script

    Parameters 
    ----------
    link_to_control_file : str
        path to the control file
    status: str
        The new string to put in the control file
    '''
    with open(link_to_control_file, 'w') as file:
        file.write(str(status))

def make_optimizer():

        
    '''
    This function builds an instance of a po.Optimizer with a bayesian optimization gaussian process using the matern kernel
    The hyperparameters of the kernel and optimizer can be tuned from here

    Returns
    -------
    An instance of the optimizer
    '''    
    space = space
    lhs = True
    n_objectives = n_objectives
    kernel_intern = 1 ** 2 * Matern(length_scale=[1]*10, length_scale_bounds=[(0.05, 10)], nu=2.5)
    base_est = cook_estimator(base_estimator='GP', space=space, normalize_y=True, noise='gaussian',
                              kernel=kernel_intern)
    # Instantiate the optimizer, we are building:
    optimizer = po.Optimizer(dimensions=space,
                                  acq_func='EI',
                                  n_initial_points=n_initial,
                                  lhs=lhs,
                                  n_objectives=n_objectives,
                                  base_estimator=base_est,
                                  acq_func_kwargs={'xi': 0.1})
    return optimizer


#Now that the functions are defined, the actual optimization loop can begin.
#The loop runs for 100000 iterations (forever) but the actual number of runs is governed by the control file


update_tecan_status(Tecan_running_file, '"standby"') #the value of the control file is set to standby as we want to make the recipe for the inital experiment

iteration = 0 # the iteration variable counts the iteration number, which is used to specify where to pipet to

for i in range(100000):

    tecan_running = read_tecan_status(Tecan_running_file) # the value of the control file is read

    print(f'Read Status of tecan_running-file. it is {tecan_running} and of type {type(tecan_running)}') # debugging line
    if tecan_running == '"running"':
        time.sleep(60)
        print('tecan must still be running. I\'ll keep waiting') # if the value is "running" the tecan is running and the script will wait with performing the next predictions
    elif tecan_running[:9] == '"standby"': #if the value is "standby" the script will perform a new prediction
        print(f'tecan is done running, now it\'s my turn') #vocal script :-)
        if iteration == 0: #in the first iteration no data needs to be extracted
            iteration = 1
            print(f'Running on iteration {iteration}') #vocal script :-)
            well_numbers = get_experiment_numbers(iteration, batch_size) #current well-numbers are estimated
            print(f'Just found wellnumbers to be {well_numbers}') #vocal script :-)

            opt = make_optimizer() #optimizer is made
            next_x = opt.ask(batch_size) #this line asks the optimizer for prediction of the next result to perform

            make_recipe_file(recipe_folder, iteration, next_x) #this function takes the next_x and writes it to a recipe CSV that the evoware can read

            update_tecan_status(Tecan_running_file, '"running"')  #The script is now done and the control file is changed to "running"
            time.sleep(5)

        else:
            old_well_numbers = get_experiment_numbers(iteration, batch_size) #this gets the wellnumbers form the latest iteration
            latest_file = get_latest_file('C:\\Users\\Public\\Documents\\Tecan\\SparkControl\\Workspaces\\',
                                               'xlsx') # this finds the xlsx data file outputtet from the plate reader
            list_of_meassurements = [
                extract_SPARK_reading(well_number, latest_file, comparator_delta=comparator_delta) for
                well_number in old_well_numbers] # this function extracts the relevant data from the data file (since all wells are read each time)
            list_of_transformed_meassurements = [1 - float(meassurement) for meassurement in list_of_meassurements] #this line transforms the meassurements to 1-meassurement since the algorithm is preset for minimization
            temporary_model=opt.tell([latest_exp], list_of_transformed_meassurements) #the model is fed the experiment(s) tried (latest_exp) and the acording result (list_of_transformed_meassurements)
            
            next_x = opt.ask(batch_size) #the optimiser predicts the next best experiments to perform
        
            iteration += 1  # the iteration is incremented by one since we are now ready to perform the next experiment
            print(f'Running on iteration {iteration}, because we have {len(opt.Xi)} datas') # debugging
            well_numbers = get_experiment_numbers(iteration, batch_size) #the well numbers for the experiment to be done are made
            print(f'Just found wellnumbers to be {well_numbers}') #debugging
            make_log(log_file, iteration, latest_exp, list_of_transformed_meassurements) #the latest experiment and the acordin results are written to the log file

                
            make_recipe_file(recipe_folder, iteration, next_x) #the experiments suggested by the optimiser are written to a CSV form for evoware to read

            update_tecan_status(Tecan_running_file, '"running"') #the control file status is set to running
            time.sleep(5)

        latest_exp = next_x[:] #latest experiment is updated to be the latest experiment outside the if - else statement

    elif tecan_running[:7] == '"final"': #when the evoware script finished the n number of iteartions it changes the control file to "final" amd the script stops!
        print('Finished')
        time.sleep(5)
        break
    else:
        print(f'tecan_running was {tecan_running}, I had expected it to be in [0,1]')

