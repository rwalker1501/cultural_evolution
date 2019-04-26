from classes_module import Target, PopulationData
import numpy as np
import os
import json

def load_population_data_source(name, file_path, info_file_path, label='txt'):
    

    filename=file_path
    data=np.load(filename)
    lats_txt=data['lats_' + label]
    lons_txt=data['lons_' + label]
    ts_txt=data['ts_' + label]
    dens_txt=data['dens_' + label]
    
    print(info_file_path)
    population_data_info = json.load(open(info_file_path));
    new_population_data = PopulationData(name, lats_txt, lons_txt, ts_txt, dens_txt, population_data_info)
    return new_population_data
    

def load_population_data(base_path, population_data_sources):
    
    pop_data_path = os.path.join(base_path, "population_data")

    # Eriksson
    eriksson_binary_path = os.path.join(pop_data_path, 'eriksson.npz')
    eriksson_info_path = os.path.join(pop_data_path, 'eriksson_info.txt')
    eriksson = load_population_data_source("Eriksson", eriksson_binary_path, eriksson_info_path)
    population_data_sources['Eriksson'] = eriksson;

    # Timmermann
    timmermann_binary_path = os.path.join(pop_data_path, 'timmermann.npz')
    timmermann_info_path = os.path.join(pop_data_path, 'timmermann_info.txt')
    timmermann = load_population_data_source("Timmermann", timmermann_binary_path, timmermann_info_path, label='ncf')
    population_data_sources['Timmermann'] = timmermann;
    
    return population_data_sources
