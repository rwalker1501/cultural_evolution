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
from os import listdir
from os.path import isfile, join
from copy import deepcopy
from scipy.stats import binom_test,wilcoxon,linregress
from scipy.optimize import curve_fit
from sklearn.metrics import r2_score
from sklearn.model_selection import RepeatedKFold
from datetime import date
from classes_module import Target, PopulationData
from datetime import datetime

pd.options.mode.chained_assignment = None 



#==============================================================================
# 'plt.style.use('ggplot')
#==============================================================================

class MainProgram:


    def __init__(self):
        self.base_path = os.getcwd()
        self.filters_applied = ""
        self.population_data_sources = []
        self.target_list=[]
        self.dataframe_loaded = False
        self.dataframe = pd.DataFrame()
        self.controls = "No Empty Lats"
        self.controls_dataframe = pd.DataFrame()
        self.clustering_on=False
        self.critical_distance=250
        self.date_window=1500
        self.number_of_kfolds = 10
        self.minimum_likelihood_ratio = 0
        self.perform_cross_validation = False
        self.user_max_for_uninhabited = 1000
        self.default_mfu = True
        self.min_date_window = 0
        self.critical_time = 10000
        self.minimum_controls = 385 # check statistics power if correct

    def set_critical_time(self, critical_time):
        self.critical_time = critical_time

    def set_clustering(self, clustering_on):
        self.clustering_on = clustering_on

    def set_date_window(self, date_window):
        self.date_window = date_window
        
    def set_critical_distance(self, critical_distance):
        self.critical_distance = critical_distance

    def set_perform_cross_validation(self, perform_cv):
        self.perform_cross_validation = perform_cv

    def set_number_of_kfolds(self, kfolds):
        self.number_of_kfolds = kfolds

    def set_minimum_likelihood_ratio(self, min_mult):
        self.set_minimum_likelihood_ratio = min_mult

    def set_user_max_for_uninhabited(self, mfi):
        self.user_max_for_uninhabited = mfi
        self.default_mfu = False

    def set_default_mfu(self, d_mfi):
        self.default_mfu = d_mfi

    def set_filters_applied(self, filters_applied):
        self.filters_applied = filters_applied

    def set_target_list(self, some_target_list):
        self.target_list = some_target_list

    def set_dataframe(self, dataframe, controls_dataframe):
        self.dataframe = dataframe
        self.controls_dataframe = controls_dataframe
        self.dataframe_loaded = True

    def set_controls(self, controls):
        self.controls = controls

    def set_perform_cross_validation(self, perform_cv):
        self.perform_cross_validation = perform_cv

    def set_number_of_kfolds(self, num_kfolds):
        self.number_of_kfolds = num_kfolds

    def set_population_data_active(self, is_active, index):
        self.population_data_sources[index].is_active = is_active

    def set_min_date_window(self, min_date_window):
        self.min_date_window = min_date_window

    def set_minimum_controls(self, minimum_controls):
        self.minimum_controls = minimum_controls

    def get_clustering(self):
        return self.clustering_on

    def get_critical_distance(self):
        return self.critical_distance

    def get_critical_time(self):
        return self.critical_time

    def get_current_target_list(self):
        return self.target_list

    def get_population_data(self):
        return self.population_data_sources

    def get_date_window(self):
        return self.date_window

    def get_dataframe_loaded(self):
        return self.dataframe_loaded

    def get_filters_applied(self):
        return self.filters_applied

    def get_base_path(self):
        return self.base_path

    def get_user_max_for_uninhabited(self):
        return self.user_max_for_uninhabited

    def get_default_mfu(self):
        return self.default_mfu

    def get_minimum_likelihood_ratio(self):
        return self.minimum_likelihood_ratio

    def get_perform_cross_validation(self):
        return self.perform_cross_validation

    def get_number_of_kfolds(self):
        return self.number_of_kfolds

    def get_controls(self):
        return self.controls

    def get_minimum_controls(self):
        return self.minimum_controls

    def load_population_data(self):
        self.population_data_sources = pdm.load_population_data(self.base_path, self.population_data_sources)

    def plot_population(self, time):
        plm.plot_densities_on_map_by_time(self.population_data_sources[1], time)
                
    def read_target_list(self, filename):
        new_list=tam.read_target_list_from_csv(filename)
        filters_applied=""
        dataframe = pd.DataFrame()
        dataframe_loaded = False;
        return new_list

    def save_target_list(self, filename, some_target_list):
        tests_path=os.path.join(self.base_path,"targets")
        tam.save_target_list_to_csv(some_target_list, tests_path, filename)
        dataframe = pd.DataFrame()
        dataframe_loaded = False;

    def add_population_data(self, name, binary_path, info_path):

        new_population_data = pdm.load_population_data_source(name, binary_path, info_path)
        self.population_data_sources.append(new_population_data)

    def generate_results(self, population_data, original_target_list, base_path, directory):


        # Check if a target list has been loaded, otherwise abort
        if len(original_target_list) == 0:
            return "Load target list before generating results"

        # Use pre-set max_for_uninhabited per population data by default
        # If user set a max_for_uninhabited, use that value
        if self.default_mfu:
            max_for_uninhabited = population_data.max_for_uninhabited
        else:
            max_for_uninhabited = self.user_max_for_uninhabited

        print("Date window: " + str(self.date_window))
        print("Max for uninhabited: " + str(max_for_uninhabited))
        print("Clustering: " + str(self.clustering_on))

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

        header_labels = ['population_data', 'max_for_uninhabited', 'date_window', 'critical_distance', "filters_applied", "minimum_controls"]
        header_values = [population_data.name, max_for_uninhabited, self.date_window, self.critical_distance, self.filters_applied, self.minimum_controls]
        wrm.write_information(f2, header_labels, header_values, ", ")

        ##############################
        # Write original target list #
        ##############################
        wrm.write_label(f2, "Target list")
        f2.write('\n')
        wrm.write_target_table(f2, original_target_list, self.date_window)
        f2.write('\n\n')

        #######################
        # Process target list #
        #######################
        #   - clusters targets and returns 2D array of targets grouped in clusters
        #   - If dataframe has not been loaded (through load processed targets), extracts dataframe and saves it.
        #   - dataframe: contains all locations and population densities in the population data that is relevant to the target list
        clustered_target_list, self.dataframe, self.controls_dataframe = tam.process_targets(self.base_path, population_data, original_target_list, self.dataframe, self.controls_dataframe, self.controls, self.dataframe_loaded, self.clustering_on, self.date_window, self.critical_distance, self.critical_time, directory, self.min_date_window)

        if self.dataframe.empty:
            f2.write("No Geographic Points Fall in Target Areas")
            f2.close()
            return "No Geographic Points Fall in Target Areas"

        #####################
        # Process dataframe #
        #####################
        # - gets statistics (means, growth coefficients) of dataframe
        # - filters target list/dataframe by removing clusters with 0 sample means 
        all_sample_means, all_control_means, growth_coefficients, samples_gt_controls, n_targets_gt_0, self.dataframe = stm.process_dataframe(self.dataframe, max_for_uninhabited)

        ########################################
        # Write filtered clustered target list #
        ########################################
        wrm.write_label(f2, "Target list after clustering")
        wrm.write_cluster_table(f2, self.dataframe, growth_coefficients)

        # plm.plot_confounder_likelihood_ratio(self.dataframe, original_target_list)
        ################################
        # Compute and write statistics #
        ################################
        # - binomial test
        # - wilcoxon
        wrm.write_label(f2, "Statistics")
        p_binomial=binom_test(samples_gt_controls,n_targets_gt_0,0.5)
        stats_header_labels = ["Number of successes", "Number of targets", "pBinomial"]
        stats_header_values = [samples_gt_controls, n_targets_gt_0, p_binomial]
        wrm.write_information(f2, stats_header_labels, stats_header_values, "   ")

        t_wilcox,p_wilcox=wilcoxon(all_sample_means,all_control_means)
        f2.write( 'Wilcoxon stat for sample vs controls means whole period:'+str(float(t_wilcox))+ '   p='+str(float(p_wilcox))+'\n')



        #################
        # Generate bins #
        #################
        # - extracts bin values
        # - write bin values to file
        bin_array, sample_counts, control_counts, likelihood_ratios, p_samples, p_controls, p_likelihood_ratios = stm.generate_bin_values(self.dataframe, self.controls_dataframe, population_data, max_for_uninhabited, self.minimum_controls)
        wrm.write_bin_table(f2, bin_array, sample_counts, control_counts, likelihood_ratios, p_samples, p_controls, p_likelihood_ratios)
        

        ##################################
        # Generate graphs and statistics #
        ##################################
        # - plots likelihood ratios and sites graphs
        stm.generate_stats_for_likelihood_ratio(bin_array, likelihood_ratios, sample_counts, control_counts, p_likelihood_ratios, p_samples, p_controls, f2, new_path, directory)

        # - plots p_graphs and write statistics (binomial and wilcoxon)
        threshold_binomial, threshold, threshold_success_count, threshold_trial_count, threshold_samples, threshold_controls = stm.generate_p_threshold_and_binomial(p_samples, p_controls, bin_array)
        plm.plot_p_graphs(bin_array, p_samples, p_controls, threshold, directory, new_path)

        t_threshold_wilcoxon, p_threshold_wilcoxon = wilcoxon(threshold_controls, threshold_samples)

        wrm.write_label(f2, "Statistics for threshold bins")
        f2.write("Threshold: " + str(threshold))
        stats_header_labels = ["Number of successes", "Number of targets", "pBinomial"]
        stats_header_values = [threshold_success_count, threshold_trial_count, threshold_binomial]
        wrm.write_information(f2, stats_header_labels, stats_header_values, "   ")
        f2.write( 'Wilcoxon stat for pControls vs pSamples:'+str(float(t_threshold_wilcoxon))+ '   p='+str(float(p_threshold_wilcoxon))+'\n')


        # - plots growth coefficients
        plm.plot_growth_coefficient_boxplot(growth_coefficients, new_path)

        # - plots targets and controls on a map
        plm.plot_targets_on_map(self.dataframe, self.controls_dataframe, new_path, directory)


        plm.plot_densities_on_map_by_range(population_data, 5700, 5800, 75, 44000)
        f2.close()

        return "Generated results."


def run_experiment(results_path, target_list_file, output_directory, population_data_name="Eriksson", controls="All", date_window=1500, user_max_for_uninhabited=-1, clustering_on = False, critical_distance=1, filter_date_before=-1, filter_not_direct=False, filter_not_figurative=False, filter_not_controversial = False, perform_cross_validation=False, number_of_kfolds = 100, minimum_controls=385, min_date_window=0, critical_time=10000, filter_min_date=-1, filter_max_date=-1, filter_min_lat=-1, filter_max_lat=-1, processed_targets=False):
    
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
        time_likelihood_ratio, bin_size, max_population, max_for_uninhabited, is_active, ascending_time = pdm.read_population_data_info(tim_info_path)
        population_data = PopulationData("Timmermann", is_active, lats_ncf, lons_ncf, ts_ncf, dens_ncf, time_likelihood_ratio, bin_size, max_population, max_for_uninhabited, ascending_time)
            

    mp.set_date_window(date_window)
    mp.set_min_date_window(min_date_window)
    if user_max_for_uninhabited != -1:
        mp.set_user_max_for_uninhabited(user_max_for_uninhabited)
    mp.set_clustering(clustering_on)
    mp.set_critical_distance(critical_distance)
    mp.set_critical_time(critical_time)
    mp.set_controls(controls)
    mp.set_minimum_controls(minimum_controls)

    if processed_targets:
        target_list, dataframe, controls_dataframe = tam.load_processed_targets(results_path, output_directory)
        mp.set_target_list(target_list)
        mp.set_dataframe(dataframe, controls_dataframe)
    else:
        filters_applied = ""
        target_list = tam.read_target_list_from_csv(base_path+"/targets/" + target_list_file)


        if filter_date_before != -1:
            target_list, filters_applied = tam.filter_targets_for_date_before(target_list, filter_date_before, filters_applied)
        if filter_not_direct:
            target_list, filters_applied = tam.filter_targets_for_not_direct(target_list, filters_applied)
        if filter_not_figurative:
            target_list, filters_applied = tam.filter_targets_for_not_figurative(target_list, filters_applied)
        if filter_not_controversial:
            target_list, filters_applied = tam.filter_targets_for_not_controversial(target_list, filters_applied)
        if filter_min_date != -1:
            target_list, filters_applied = tam.filter_targets_for_date(target_list, filter_min_date, filter_max_date, filters_applied)
        if filter_min_lat != -1:
            target_list, filters_applied = tam.filter_targets_for_latitude(target_list, filter_min_lat, filter_max_lat, filters_applied)
        mp.set_filters_applied(filters_applied)

    # mp.set_perform_cross_validation(perform_cross_validation)
    # if perform_cross_validation:
    #     mp.set_number_of_kfolds(number_of_kfolds)
    #     mp.set_minimum_likelihood_ratio(minimum_likelihood_ratio)

    print(len(target_list))

    return mp.generate_results(population_data, target_list, results_path, output_directory)
