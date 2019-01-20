"""Convert specific txt files to netcdf format"""

# pylint: disable=R0914
# THIS SEEMS LIKE UP TO DATE VERSION
# from __future__ import division
import numpy as np
import os
import sys
import random
import target_module as tam
import population_data_module as pdm
import stats_module as stm
import plot_module as plm
import write_module as wrm
import pandas as pd
import math
from os import listdir
from os.path import isfile, join
from copy import deepcopy
from scipy.stats import binom_test,wilcoxon,linregress,ks_2samp
from scipy.optimize import curve_fit
from sklearn.metrics import r2_score
from sklearn.model_selection import RepeatedKFold
from datetime import date
from classes_module import Target, PopulationData
from datetime import datetime
from mpmath import mpf
pd.options.mode.chained_assignment = None 



#==============================================================================
# 'plt.style.use('ggplot')
#==============================================================================

class MainProgram:


    def __init__(self):
        self.base_path = os.getcwd();
        self.population_data_sources = []
        self.target_list=[]
        self.dataframe = pd.DataFrame()
        self.globals_dataframe = pd.DataFrame()

        parameters = {};
        parameters['filters_applied'] = "";
        parameters['dataframe_loaded'] = False;
        parameters['globals_type'] = "All Except Equatorials";
        parameters['date_window'] = 24
        parameters['date_lag'] = 0;
        parameters['user_max_for_uninhabited'] = 1;
        parameters['default_max_for_uninhabited'] = False;
        parameters['max_lat'] = 60;
        parameters['min_lat'] = -40;
        parameters['max_date'] = 50000;
        parameters['min_date'] = 0;
        parameters['min_globals'] = 1;
        parameters['min_p'] = 0.01;

        parameters['clustering_on'] = False;
        parameters['critical_distance'] = 0;
        parameters['critical_time'] = 10000;

        self.parameters = parameters;

    ##################
    # Setter Methods #
    ##################


    def set_population_data_active(self, is_active, index):
        self.population_data_sources[index].is_active = is_active;
    

    def set_target_list(self, some_target_list):
        self.target_list = some_target_list;

    def set_dataframe(self, dataframe, globals_dataframe):
        self.dataframe = dataframe;
        self.globals_dataframe = globals_dataframe;
        self.dataframe_loaded = True;

    def set_user_max_for_uninhabited(self, mfi):
        self.parameters['user_max_for_uninhabited'] = mfi;
        self.parameters['default_max_for_uninhabited'] = False;

    def set_parameter(self, parameter_name, value):
        self.parameters[parameter_name] = value;

    ##################
    # Getter Methods #
    ##################

    def get_base_path(self):
        return self.base_path;

    def get_population_data(self):
        return self.population_data_sources;

    def get_current_target_list(self):
        return self.target_list;

    def get_parameters(self):
        return self.parameters;


    #############################
    # Population Data Functions #
    #############################

    def load_population_data(self):
        self.population_data_sources = pdm.load_population_data(self.base_path, self.population_data_sources)

    def add_population_data(self, name, binary_path, info_path):
        new_population_data = pdm.load_population_data_source(name, binary_path, info_path)
        self.population_data_sources.append(new_population_data)

    #########################
    # Target List Functions #
    #########################

    def read_target_list(self, filename):
        new_list=tam.read_target_list_from_csv(filename)
        self.filters_applied=""
        self.dataframe = pd.DataFrame()
        self.dataframe_loaded = False;
        return new_list

    def save_target_list(self, filename, some_target_list):
        tests_path=os.path.join(self.base_path,"targets")
        tam.save_target_list_to_csv(some_target_list, tests_path, filename)
        self.dataframe = pd.DataFrame()
        self.dataframe_loaded = False;

    ##################
    # Plot Functions #
    ##################

    def adjust_time_value(self, time, time_multiplier):
        if  time % time_multiplier != 0:
            time = time + (time_multiplier - time % time_multiplier)
        return time

    def plot_population_by_time(self, population_data, time):
        name = population_data.name
        time = self.adjust_time_value(time, population_data.time_multiplier)


        print("Plotting map...")
        print("Time: " + str(time))
        print("Source: " + name)
        plm.plot_densities_on_map_by_time_point(population_data, time)

    def plot_population_by_time_range(self, population_data, start_time, end_time):
        name = population_data.name
        time_multiplier = population_data.time_multiplier
        start_time = self.adjust_time_value(start_time, population_data.time_multiplier)
        end_time = self.adjust_time_value(end_time, population_data.time_multiplier)

        print("Plotting map...")
        print("Start Time: " + str(start_time))
        print("End Time: " + str(end_time))
        print("Source: " + name)
        plm.plot_densities_on_map_by_time_range(population_data, start_time, end_time)

    def plot_population_by_range(self, population_data, min_density, max_density, start_time, end_time):
        name = population_data.name
        time_multiplier = population_data.time_multiplier
        start_time = self.adjust_time_value(start_time, population_data.time_multiplier)
        end_time = self.adjust_time_value(end_time, population_data.time_multiplier)

        print("Plotting map...")
        print("Min Density: " + str(min_density))
        print("Max Density: " + str(max_density))
        print("Start Time: " + str(start_time))
        print("End Time: " + str(end_time))
        print("Source: " + name)
        plm.plot_densities_on_map_by_range(population_data, min_density, max_density, start_time, end_time);

     

    #############################
    # Generate Results Function #
    #############################

    def generate_results(self, population_data, original_target_list, base_path, directory):


        # Check if a target list has been loaded, otherwise abort
        if len(original_target_list) == 0:
            return "Load target list before generating results"

        # Use pre-set max_for_uninhabited per population data by default
        # If user set a max_for_uninhabited, use that value
        if self.parameters['default_max_for_uninhabited']:
            max_for_uninhabited = population_data.max_for_uninhabited
        else:
            max_for_uninhabited = self.parameters['user_max_for_uninhabited'];

        #####################################
        # Create directory and results file #
        #####################################
        results_path = os.path.join(base_path, "results")
        if not os.path.exists(results_path):
            os.makedirs(results_path)

        new_path = os.path.join(results_path, directory)
        if not os.path.exists(new_path):
            os.makedirs(new_path)
        results_filename= os.path.join(new_path, directory + "_results.csv") 

        print("Results: " + results_filename)
        f2= open(results_filename, 'w')


        ############################
        # Write header information #
        ############################
        dateTime=str(datetime.now())
        f2.write('Date: '+dateTime)
        f2.write('\n')

        wrm.write_parameters(f2, self.parameters);

        #######################
        # Process target list #
        #######################
        #   - clusters targets and returns 2D array of targets grouped in clusters
        #   - If dataframe has not been loaded (through load processed targets), extracts dataframe and saves it.
        #   - dataframe: contains all locations and population densities in the population data that is relevant to the target list
        clustered_target_list, self.dataframe, self.globals_dataframe = tam.process_targets(self.base_path, population_data, original_target_list, self.dataframe, self.globals_dataframe, self.parameters, max_for_uninhabited, directory)
        if self.dataframe.empty or len(clustered_target_list)<10:
            f2.write("Not enough sites in Target Areas")
            f2.close()
            return "Not enough sites in target area"

        print("Processing sites and controls dataframe...")
        self.dataframe = stm.process_dataframe(self.dataframe)

        print("Saving merged sites and globals dataframes...")
        tam.generate_merged_dataframe(base_path, directory, self.dataframe, self.globals_dataframe);
       
        ########################################
        # Write filtered clustered target list #
        ########################################
        print("Writing target list...")
        wrm.write_label(f2, "Target list after clustering")
        wrm.write_cluster_table(f2, self.dataframe, self.parameters)
        
        #################
        # Generate bins #
        #################
        # - extracts bin values
        # - write bin values to file
        print("Writing bins...")
        bin_values_df = stm.generate_bin_values_dataframe(self.dataframe, self.globals_dataframe, population_data, self.parameters['min_globals'])
        wrm.write_bin_table(f2, bin_values_df, self.parameters['min_globals'])


        #################
        # Stat Analysis #
        #################
        # - n, median, mean, std of samples and globals
        # - K-S2 test for bin distribution of samples (p_samples) vs. bin distribution of globals (p_globals)
        print("Calculating statistics...")
        stat_dictionary, trimmed_bin_values_df = stm.generate_statistics(self.dataframe, self.globals_dataframe, bin_values_df, self.parameters['min_globals'])
        
        if stat_dictionary is None:
            f2.write('insufficient non-zero bins for analysis');
            return 'insufficient non-zero bins for analysis';
        
        print("Writing analysis...")
        wrm.write_analysis(f2, stat_dictionary, self.parameters['min_p']);

        ###############
        # Plot Graphs #
        ###############

        print("Generating graphs...")
        plm.plot_stat_graphs(stat_dictionary, trimmed_bin_values_df, population_data, directory, new_path);

        # plots targets and globals on a map
        plm.plot_targets_on_map(self.dataframe, self.globals_dataframe, new_path, directory)

        f2.close();

        return "Generated results"

def run_experiment(results_path, target_list_file, output_directory, population_data_name="Eriksson", globals_type="All Except Equatorials", date_window=24, date_lag = 0, user_max_for_uninhabited=-1, clustering_on = False, critical_distance=0, critical_time=10000, filter_date_before=-1, filter_not_direct=False, filter_not_exact=False, filter_not_figurative=False, filter_not_controversial = False,  filter_min_date=-1, filter_max_date=-1, filter_min_lat=-1, filter_max_lat=-1, processed_targets=False):
  # Note: current setting of minimum_globals is overwritten in stats_modulegenera
    mp = MainProgram()
    base_path = mp.get_base_path()
    pop_data_path = os.path.join(base_path, "population_data")
    if population_data_name == "Eriksson":
        eriksson_binary_path = os.path.join(pop_data_path, "eriksson.npz") #base_path+'/population_data/eriksson.npz'
        eriksson_info_path = os.path.join(pop_data_path, "eriksson_info.txt") #base_path+'/population_data/eriksson_info.txt'
        population_data = pdm.load_population_data_source("Eriksson", eriksson_binary_path, eriksson_info_path)
    else:
        tim_fn= os.path.join(pop_data_path, "timmermann.npz") #base_path+'population_data/timmermann.npz'
        data=np.load(tim_fn)
        lats_ncf=data['lats_ncf']
        lons_ncf=data['lons_ncf']
        ts_ncf=data['ts_ncf']
        dens_ncf=data['dens_ncf']

        tim_info_path = os.path.join(pop_data_path, "timmermann_info.txt") #base_path+'/population_data/timmermann_info.txt'
        # time_likelihood_ratio should be time_multiplier
        time_likelihood_ratio, density_multiplier,bin_size, max_population, max_for_uninhabited, is_active, ascending_time = pdm.read_population_data_info(tim_info_path)
        population_data = PopulationData("Timmermann", is_active, lats_ncf, lons_ncf, ts_ncf, dens_ncf, time_likelihood_ratio, density_multiplier,bin_size, max_population, max_for_uninhabited, ascending_time)
            

    mp.set_parameter('date_window', date_window)
    mp.set_parameter('date_lag', date_lag)
    if user_max_for_uninhabited != -1:
        mp.set_parameter('user_max_for_uninhabited', user_max_for_uninhabited)
    mp.set_parameter('clustering_on', clustering_on)
    mp.set_parameter('critical_distance', critical_distance)
    mp.set_parameter('critical_time', critical_time)
    mp.set_parameter('globals_type', globals_type)
    if processed_targets:
        target_list, dataframe, globals_dataframe = tam.load_processed_targets(results_path, output_directory)
        mp.set_target_list(target_list)
        mp.set_dataframe(dataframe, globals_dataframe)
    else:
        filters_applied = ""
        target_list = tam.read_target_list_from_csv(base_path+"/targets/" + target_list_file)


        if filter_date_before != -1:
            target_list, filters_applied = tam.filter_targets_for_date_before(target_list, filter_date_before, filters_applied)
        if filter_not_direct:
            target_list, filters_applied = tam.filter_targets_for_not_direct(target_list, filters_applied)
        if filter_not_exact:
            target_list, filters_applied = tam.filter_targets_for_not_exact_age(target_list, filters_applied)
        if filter_not_figurative:
            target_list, filters_applied = tam.filter_targets_for_not_figurative(target_list, filters_applied)
        if filter_not_controversial:
            target_list, filters_applied = tam.filter_targets_for_not_controversial(target_list, filters_applied)
        if filter_min_date != -1:
            target_list, filters_applied = tam.filter_targets_for_date(target_list, filter_min_date, filter_max_date, filters_applied)
        if filter_min_lat != -1:
            target_list, filters_applied = tam.filter_targets_for_latitude(target_list, filter_min_lat, filter_max_lat, filters_applied)
        mp.set_parameter('filters_applied', filters_applied)

    # mp.set_perform_cross_validation(perform_cross_validation)
    # if perform_cross_validation:
    #     mp.set_number_of_kfolds(number_of_kfolds)
    #     mp.set_minimum_likelihood_ratio(minimum_likelihood_ratio)
    mp.generate_results(population_data, target_list, results_path, output_directory)

