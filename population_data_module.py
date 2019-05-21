from classes_module import Target, PopulationData
import numpy as np
import os
def load_population_data_source(name, file_path, info_file_path):
    time_multiplier, density_multiplier,bin_size, max_population, max_for_uninhabited, is_active, ascending_time,likelihood_parameters = read_population_data_info(info_file_path)

    filename=file_path
    data=np.load(filename)
    lats_txt=data['lats_txt']
    lons_txt=data['lons_txt']
    ts_txt=data['ts_txt']
    dens_txt=data['dens_txt']
    new_population_data = PopulationData(name, is_active, lats_txt, lons_txt, ts_txt, dens_txt,  time_multiplier, density_multiplier,bin_size, max_population, max_for_uninhabited, ascending_time,likelihood_parameters)
    return new_population_data
    

def load_population_data(base_path, population_data_sources):
    
    pop_data_path = os.path.join(base_path, "population_data")

    # Synthetic
    synthetic_binary_path = os.path.join(pop_data_path, "synth2.npz")
    synthetic_info_path = os.path.join(pop_data_path,'synth2_info.txt')
    synthetic = load_population_data_source("Synthetic", synthetic_binary_path, synthetic_info_path)
    population_data_sources.append(synthetic)

    # Eriksson
    eriksson_binary_path = os.path.join(pop_data_path, 'eriksson.npz')
    eriksson_info_path = os.path.join(pop_data_path, 'eriksson_info.txt')
    eriksson = load_population_data_source("Eriksson", eriksson_binary_path, eriksson_info_path)
    population_data_sources.append(eriksson)

    # Timmermann
    tim_fn=os.path.join(pop_data_path, 'timmermann.npz')
    data=np.load(tim_fn)
    lats_ncf=data['lats_ncf']
    lons_ncf=data['lons_ncf']
    ts_ncf=data['ts_ncf']
    dens_ncf=data['dens_ncf']

    tim_info_path = os.path.join(pop_data_path, 'timmermann_info.txt')
    time_multiplier, density_multiplier,bin_size, max_population, max_for_uninhabited, is_active, ascending_time, likelihood_parameters = read_population_data_info(tim_info_path)
     # why do we use a different method to load timmermann data?
    timmermann = PopulationData("Timmermann", is_active, lats_ncf, lons_ncf, ts_ncf, dens_ncf, time_multiplier, density_multiplier,bin_size, max_population, max_for_uninhabited, ascending_time, likelihood_parameters)
    population_data_sources.append(timmermann)

    return population_data_sources

def read_population_data_info(file_path):
    info_file = open(file_path)

    # time multiplier
    line = info_file.readline()
    data = line.split("\t")
    time_multiplier = float(data[1])
    
    # density multiplier
    
    line = info_file.readline()
    data = line.split("\t")
    density_multiplier = float(data[1])
    

    # bin size
    line = info_file.readline()
    data = line.split("\t")
    bin_size = int(data[1])

    # max population
    line = info_file.readline()
    data = line.split("\t")
    max_population = float(data[1])

    # max for uninhabited
    line = info_file.readline()
    data = line.split("\t")
    max_for_uninhabited = float(data[1])

    # is active
    line = info_file.readline()
    data = line.strip().split("\t")
    is_active = data[1] == 'True'

    # ascending time
    line = info_file.readline()
    data = line.strip().split("\t")
    ascending_time = data[1] == 'True'
    
    # parameters for likelihood estimate - this is very ugly
    line = info_file.readline()
    data = line.strip().split("\t")
    lambda_start=float(data[1])
    line = info_file.readline()
    data = line.strip().split("\t")
    lambda_end=float(data[1])
    line = info_file.readline()
    data = line.strip().split("\t")
    zetta_start=float(data[1])
    line = info_file.readline()
    data = line.strip().split("\t")
    zetta_end=float(data[1])
    line = info_file.readline()
    data = line.strip().split("\t")
    eps_start=float(data[1])
    line = info_file.readline()
    data = line.strip().split("\t")
    eps_end=float(data[1])
    line = info_file.readline()
    data = line.strip().split("\t")
    y_acc_start=float(data[1])
    line = info_file.readline()
    data = line.strip().split("\t")
    y_acc_end=float(data[1])
    likelihood_parameters=np.array((lambda_start,lambda_end, zetta_start,zetta_end, eps_start,eps_end,y_acc_start, y_acc_end))
    

    return time_multiplier, density_multiplier, bin_size, max_population, max_for_uninhabited, is_active, ascending_time, likelihood_parameters