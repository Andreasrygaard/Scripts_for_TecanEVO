import csv
import time
import random
import ProcessOptimizer as po # used to suggest next experiment
from ProcessOptimizer.utils import cook_estimator #Used to customize the algorithm
from ProcessOptimizer.learning.gaussian_process.kernels import Matern #A specific "brain" in the algorithm
import glob
import os
import pandas as pd


Tecan_running_file = 'C:\\Users\\Tecan\\Desktop\\Variable.txt'
recipe_folder = 'C:\\Users\\Public\\Documents\\Tecan\\'
log_file = 'C:\\temp\\compressedLog_20220318.csv'
batch_size = 1
comparator_delta = None
space = [(10, 36), (10, 36)]
n_objectives = 1


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

def read_tecan_status(link_to_control_file):
    with open(link_to_control_file) as file:
        line = file.readlines()[0]
    return line


def make_recipe_file(link_to_recipe_folder, iteration, recipe):
    next_x = recipe

    if batch_size == 1:
        recipe_file_first_add_falcon = os.path.join(link_to_recipe_folder, 'TRYPSIN_first_adds.csv')
        with open(recipe_file_first_add_falcon, mode='w', newline='') as file:
            writer = csv.writer(file, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
            if next_x[0] != 0:
                writer.writerow(['SOURCE4','1','DEST1', str(iteration), str(next_x[0])]) #Acidic buffer
            if next_x[1] != 0:
                writer.writerow(['SOURCE4','2','DEST1', str(iteration), str(next_x[1])]) #Basic buffer
            writer.writerow(['SOURCE4','3','DEST1',str(iteration),str(90-(next_x[0]+next_x[1]))]) #Vol(water until 225uL)
            writer.writerow(['SOURCE4','4','DEST1',str(iteration),str(8)]) ##Add Trypsin solution until 10 ug/mL
            writer.writerow(['SOURCE4','5','DEST1',str(iteration),str(16)]) #Add substrate until 250uM
            #antal rækker tilsvarende antal komponenter i buffer
    else:
        recipe_file_first_add_falcon = os.path.join(link_to_recipe_folder, 'LOCI_first_adds_FALCON1.csv')
        recipe_file_first_add_Amber = os.path.join(link_to_recipe_folder, 'LOCI_first_adds_AMBER1.csv')
        recipe_file_second_add_falcon = os.path.join(link_to_recipe_folder, 'LOCI_second_adds_FALCON1.csv')
        with open(recipe_file_first_add_falcon, mode='w', newline='') as file:
            writer = csv.writer(file, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
        with open(recipe_file_first_add_Amber, mode='w', newline='') as file:
            writer = csv.writer(file, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
        with open(recipe_file_second_add_falcon, mode='w', newline='') as file:
            writer = csv.writer(file, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
        print(f'Ordering experiments for wellnumbers {well_numbers}')

        for i in range(batch_size):
            with open(recipe_file_first_add_falcon, mode='a', newline='') as file:
                writer = csv.writer(file, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
                writer.writerow(['FALCON1', '1', '384WellMP', str(well_numbers[i]), str(next_x[i][0])])  # ACC
                writer.writerow(['FALCON1', '2', '384WellMP', str(well_numbers[i]),
                                 str((15 - next_x[i][0]) * next_x[i][1])])  # Bio
                writer.writerow(['FALCON1', '3', '384WellMP', str(well_numbers[i]), str(20 - (
                            next_x[i][0] + ((15 - next_x[i][0]) * next_x[i][1]) + next_x[i][2]))])  # Vol(buffer), as a remainder of the total volume to 20 uL
                writer.writerow(['FALCON1', '1', '384WellMP', str(well_numbers[i] + comparator_delta),
                                 str(next_x[i][0])])  # ACC
                writer.writerow(['FALCON1', '2', '384WellMP', str(well_numbers[i] + comparator_delta),
                                 str((15 - next_x[i][0]) * next_x[i][1])])  # Bio
                writer.writerow(['FALCON1', '3', '384WellMP', str(well_numbers[i] + comparator_delta), str(20 - (
                            next_x[i][0] + ((15 - next_x[i][0]) * next_x[i][1])))])  # Vol(buffer), as a remainder of the total volume to 20 uL

            with open(recipe_file_first_add_Amber, mode='a', newline='') as file:
                writer = csv.writer(file, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
                writer.writerow(['AMBER1', '1', '384WellMP', str(well_numbers[i]), str(next_x[i][2])])  # Sample

            with open(recipe_file_second_add_falcon, mode='a', newline='') as file:
                writer = csv.writer(file, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
                writer.writerow(['FALCON1', '4', '384WellMP', str(well_numbers[i]), str(next_x[i][3])])  # Donor
                writer.writerow(['FALCON1', '5', '384WellMP', str(well_numbers[i]), str(30 - next_x[i][3])])
                writer.writerow(['FALCON1', '4', '384WellMP', str(well_numbers[i] + comparator_delta),
                                 str(next_x[i][3])])  # Donor
                writer.writerow(['FALCON1', '5', '384WellMP', str(well_numbers[i] + comparator_delta),
                                 str(30 - next_x[i][3])])




def make_log(link_to_log_file, iteration, latest_exp, list_of_transformed_meassurements):
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

def update_tecan_status(link_to_control_file, status):
    with open(link_to_control_file, 'w') as file:
        file.write(str(status))

def make_optimizer():
    space = [(10,36),(10,36)] #(0,10) x antal solvent
    lhs = True
    n_objectives = 1
    kernel_intern = 1 ** 2 * Matern(length_scale=[1,1], length_scale_bounds=[(0.05, 10)], nu=2.5)
    base_est = cook_estimator(base_estimator='GP', space=space, normalize_y=True, noise='gaussian',
                              kernel=kernel_intern)
    # Instantiate the optimizer, we are building:
    optimizer = po.Optimizer(dimensions=space,
                                  acq_func='EI',
                                  n_initial_points=3,
                                  lhs=lhs,
                                  n_objectives=n_objectives,
                                  base_estimator=base_est,
                                  acq_func_kwargs={'xi': 0.1})
    return optimizer



update_tecan_status(Tecan_running_file, '"standby"')


iteration = 0



for i in range(100000):

    tecan_running = read_tecan_status(Tecan_running_file)

    print(f'Read Status of tecan_running-file. it is {tecan_running} and of type {type(tecan_running)}')
    if tecan_running == '"running"':
        time.sleep(60)
        print('tecan must still be running. I\'ll keep waiting')
    elif tecan_running[:9] == '"standby"':
        print(f'tecan is done running, now it\'s my turn')
        if iteration == 0:
            iteration = 1
            print(f'Running on iteration {iteration}')
            well_numbers = get_experiment_numbers(iteration, batch_size)
            print(f'Just found wellnumbers to be {well_numbers}')

            opt = make_optimizer()
            next_x = opt.ask(batch_size)

            make_recipe_file(recipe_folder, iteration, next_x)

            update_tecan_status(Tecan_running_file, '"running"')
            time.sleep(5)

        else:
            old_well_numbers = get_experiment_numbers(iteration, batch_size)
            latest_file = get_latest_file('C:\\Users\\Public\\Documents\\Tecan\\SparkControl\\Workspaces\\',
                                               'xlsx')
            list_of_meassurements = [
                extract_SPARK_reading(well_number, latest_file, comparator_delta=comparator_delta) for
                well_number in old_well_numbers]
            list_of_transformed_meassurements = [1 - float(meassurement) for meassurement in list_of_meassurements]
            opt.tell([latest_exp], list_of_transformed_meassurements)
            next_x = opt.ask(batch_size)
            iteration += 1
            print(f'Running on iteration {iteration}, because we have {len(opt.Xi)} datas')
            well_numbers = get_experiment_numbers(iteration, batch_size)
            print(f'Just found wellnumbers to be {well_numbers}')
            make_log(log_file, iteration, latest_exp, list_of_transformed_meassurements)

            make_recipe_file(recipe_folder, iteration, next_x)

            update_tecan_status(Tecan_running_file, '"running"')
            time.sleep(5)

        latest_exp = next_x[:]

    elif tecan_running[:7] == '"final"':
        print('Finished')
        break
    else:
        print(f'tecan_running was {tecan_running}, I had expected it to be in [0,1]')
