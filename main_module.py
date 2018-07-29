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
from scipy.stats import binom_test,wilcoxon,linregress,ks_2samp
from scipy.optimize import curve_fit
from sklearn.metrics import r2_score
from sklearn.model_selection import RepeatedKFold
from datetime import date
from classes_module import Target, PopulationData
from datetime import datetime
import math

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
        self.globals = "All"
        # self.globals = "No Empty Lats"
        self.globals_dataframe = pd.DataFrame()
        self.clustering_on=False
        self.critical_distance=0
        self.date_window=50
        self.number_of_kfolds = 10
        self.minimum_likelihood_ratio = 0
        self.perform_cross_validation = False
        self.user_max_for_uninhabited = 1
        self.default_mfu = True
        self.min_date_window = 0
        self.critical_time = 10000
        self.max_lat = 60;
        self.min_lat = -40;
        self.max_date = 50000;
        self.min_date = 0;
        self.min_globals=100 #this is (currently abitrary) value for min number of globals in a valid bin)
        self.min_p=0.01 #this is (currently abitrary) min p value for qualitative tests of match to model)
        self.reweighting=False
        

    def set_max_lat(self, value):
        self.max_lat = value;
    def set_min_lat(self, value):
        self.min_lat = value;
    def set_max_date(self, value):
        self.max_date = value;
    def set_min_date(self, value):
        self.min_date = value;

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

    def set_dataframe(self, dataframe, globals_dataframe):
        self.dataframe = dataframe
        self.globals_dataframe = globals_dataframe
        self.dataframe_loaded = True

    def set_globals(self, globals):
        self.globals = globals

    def set_perform_cross_validation(self, perform_cv):
        self.perform_cross_validation = perform_cv

    def set_number_of_kfolds(self, num_kfolds):
        self.number_of_kfolds = num_kfolds

    def set_population_data_active(self, is_active, index):
        self.population_data_sources[index].is_active = is_active

    def set_min_date_window(self, min_date_window):
        self.min_date_window = min_date_window

    def set_minimum_globals(self, minimum_globals):
        self.minimum_globals = minimum_globals

    def get_max_lat(self):
        return self.max_lat;
    def get_min_lat(self):
        return self.min_lat;
    def get_max_date(self):
        return self.max_date;
    def get_min_date(self):
        return self.min_date;

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

    def get_globals(self):
        return self.globals

    def load_population_data(self):
        self.population_data_sources = pdm.load_population_data(self.base_path, self.population_data_sources)

    def plot_population(self, time):
        plm.plot_densities_on_map_by_time(self.population_data_sources[1], time)

    def plot_population_by_time(self, population_data, time):
        name = population_data.name
        time_multiplier = population_data.time_multiplier
        print(time)
        if  time % time_multiplier != 0:
            # get next value divisible by time_multiplier
            # For example: 
            #   time_multiplier = 25
            #   date_to = 20
            #   time = 20 + (25 - 20) = 25
            print(time_multiplier)
            time = time + (time_multiplier - time % time_multiplier)
            print(time)

        print("Getting map...")
        print("Time: " + str(time))
        print("Source: " + name)
        plm.plot_densities_on_map_by_time(population_data, time)

    def plot_min_densities_in_time_range(self, population_data, time_from, time_to, min_density):
        plm.plot_min_densities_in_time_range(population_data, time_from, time_to, min_density)
                
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

    def generate_reweighted_data(self, population_data, original_target_list, base_path, directory):
        test_labels=[]
        test_results=[]
        tests_passed_array=[]
        results_path = os.path.join(base_path, "results")
        if not os.path.exists(results_path):
            os.makedirs(results_path)

        new_path = os.path.join(results_path, directory)
        if not os.path.exists(new_path):
            os.makedirs(new_path)

        clustered_target_list, self.dataframe, self.globals_dataframe = tam.process_targets(self.base_path, population_data, original_target_list, self.dataframe, self.globals_dataframe, self.globals, self.dataframe_loaded, self.clustering_on, self.date_window, self.critical_distance, self.critical_time, directory, self.min_date_window, self.min_lat, self.max_lat, self.min_date, self.max_date)

        columns = ["period", "latitude"];
        

        for column in columns:
            results_directory = directory + "_" + column + "_results.csv"; #This is not a directory - needs to be cleaned up
            test_results_filename = os.path.join(new_path, results_directory);
            test_results_file = open(test_results_filename, 'w');

            unique_keys = []
            min_val = 0
            max_val = 15000
            interval = 10
            if column == "period":
                unique_keys = self.dataframe.period.unique();
                min_val = self.min_date
                max_val = self.max_date
                interval = 10000
            elif column == "latitude":
                unique_keys = self.dataframe.latitude.unique();
                min_val = self.min_lat
                max_val = self.max_lat
                interval = 10
            else:
                return ("Invalid column: " + column);

            for x in range(0, 2): 

                length = int(math.ceil((max_val - min_val)/float(interval)))


                ranges = [max_val-interval*i for i in range(0, length)]
                ranges.append(min_val);
                ranges.reverse()

                # Generate random values
                base_random_values = np.random.uniform(1,5,[1,length])[0];

                # Normalize random values
                random_values_sum = sum(base_random_values)
                random_values = base_random_values/random_values_sum

                # Reweighting per key
                dataframe = self.dataframe[self.dataframe.type=='s'];
                globals_dataframe = self.globals_dataframe

                dataframe['multiplier'] = pd.cut(dataframe[column], bins=ranges, labels=random_values, include_lowest=True).astype(float);
                globals_dataframe['multiplier'] = pd.cut(dataframe[column], bins=ranges, labels=random_values, include_lowest=True).astype(float);

                print(dataframe[dataframe.multiplier.isnull()])


                bin_size = population_data.bin_size
                max_population = population_data.max_population
                # minimum_bin=max_for_uninhabited
                minimum_bin = 0
                bins_to_omit=int(minimum_bin/bin_size)

                ######################
                # Create bin columns #
                ######################
                # creating bins according to density of the row

                # main dataframe
                dataframe['bin_index'] = (dataframe.density/bin_size)-bins_to_omit
                dataframe['bin_index'] = dataframe.bin_index.astype(int)
                dataframe = dataframe[dataframe.bin_index >= 0]
                dataframe['bin'] = dataframe.bin_index*bin_size+minimum_bin

                # globals dataframe
                globals_dataframe['bin_index'] = (globals_dataframe.density/bin_size)-bins_to_omit
                globals_dataframe['bin_index'] = globals_dataframe.bin_index.astype(int)
                globals_dataframe['bin'] = globals_dataframe.bin_index*bin_size+minimum_bin

                bin_array = [minimum_bin+bin_size*i for i in range(0, int(max_population/bin_size) + 1)]
                if bin_array[-1] != max_population:
                    bin_array.append(max_population);

                dataframe['m_sample'] = dataframe['contribution'].groupby(dataframe['bin']).transform('sum')*dataframe['multiplier']
                dataframe['sample_count'] = dataframe['m_sample'].groupby(dataframe['bin']).transform('sum')

                globals_dataframe['m_global'] = globals_dataframe['density'].groupby(globals_dataframe['bin']).transform('count')*globals_dataframe['multiplier']
                globals_dataframe['global_count'] = globals_dataframe['m_global'].groupby(globals_dataframe['bin']).transform('sum')
                label="Reweighting"
                identifier=directory + "_" + column + "_" + str(x); #x is the number of the attempt
                results_filename=os.path.join(new_path,str(identifier) + "_" + label +".csv")   
                print "results_filename=",results_filename
                print "path=",new_path
                results_file=open(results_filename, 'w')
                tests_passed=stm.write_results(results_file,identifier,new_path,self.dataframe, self.globals_dataframe, population_data,self.min_globals, self.min_p)
                test_results_file.write(str(x)+";"+str(tests_passed) +"\n")

# =============================================================================
#                 Richard cut this but it probably is necessary - if so the previous write results instruction comes too early                
# 
#                 temp_df = globals_dataframe.groupby('bin').first().reset_index()
#                 temp_global_counts = temp_df['global_count'].values
#                 g_bin_array = temp_df['bin'].values;
#                 extra_bins = len(bin_array) - len(g_bin_array)
#                 if extra_bins > 0:
#                     g_bin_array = np.concatenate((g_bin_array, [0 for i in range(0, extra_bins)]))
# 
#                 temp_df = dataframe.groupby('bin').first().reset_index()
#                 temp_sample_counts = temp_df['sample_count'].values;
#                 s_bin_array = temp_df['bin'].values
#                 extra_bins = len(bin_array) - len(s_bin_array)
#                 if extra_bins > 0:
#                     s_bin_array = np.concatenate((s_bin_array, [0 for i in range(0, extra_bins)]))
# 
#                 s = 0;
#                 g = 0;
#                 global_counts = []
#                 sample_counts = []
#                 for i in range(0, len(bin_array)):
#                     if g_bin_array[i-g] != bin_array[i]:
#                         global_counts.append(0);
#                         g += 1;
#                     else:
#                         global_counts.append(temp_global_counts[i-g])
# 
#                     if s_bin_array[i-s] != bin_array[i]:
#                         sample_counts.append(0);
#                         s += 1;
#                     else:
#                         sample_counts.append(temp_sample_counts[i-s])
# 
#                 control_counts = np.array(global_counts)-np.array(sample_counts);
#                 generate_bin_values()
#                 odds_ratios, top_MHs, bottom_MHs, top_test_MHs, bottom_test_MHs, upper_cis, lower_cis = stm.generate_or_mh_ci_stats(sample_counts, global_counts, control_counts);
#                 results_identifier = directory + "_" + column + "_" + str(x);
#    #             plm.plot_odds_ratio(trimmed_bin_array, trimmed_ratios, 0, [], [], trimmed_lower_cis, trimmed_upper_cis, max_population-bin_size*2, results_identifier, "Odds Ratio", new_path);
# 
#                 label = column.title() + " " + str(x);
#                 wrm.write_random_weighting_table(results_file, label, ranges, base_random_values, random_values, bin_array, odds_ratios);
# 
# =============================================================================
            
            test_results_file.close()
            results_file.close()
            

    def generate_confounder_analysis(self, population_data, original_target_list, base_path, directory):


        df_loaded = self.dataframe_loaded;

        # Regular analysis
        if population_data is None:
            print "Population data is none"
        if original_target_list is None:
            print"target list is none"
        if base_path is None:
            print "base path is none"
        if directory is None:
            print "directory is none"
        results=self.generate_results(population_data, original_target_list, base_path, directory+"_base", True)
        if results is None:
            print "No results from generate_results"
        dataframe, globals_dataframe, bin_array, sample_counts, global_counts, control_counts, orig_odds_ratios, top_MHs, bottom_MHs, top_test_MHs, bottom_test_MHs, likelihood_ratios, p_samples, p_globals, p_likelihood_ratios, tests_passed = self.generate_results(population_data, original_target_list, base_path, directory+"_base", True)
        # Date analysis
        minimum_date = self.min_date;
        maximum_date = self.max_date;

        date_bands = dict();
        date_keys = []
        start_date = maximum_date;
        end_date = maximum_date - 10000;
        while(True):
            if end_date < minimum_date:
                end_date = minimum_date;

            new_dir = "date_" + str(start_date) + "-" + str(end_date) + "_" + directory;
            new_target_list, f = tam.filter_targets_for_date(original_target_list, end_date, start_date, "");

            new_df = dataframe[(dataframe.period < start_date) & (dataframe.period >= end_date)]
            new_gdf = globals_dataframe[(globals_dataframe.period < start_date) & (globals_dataframe.period >= end_date)]

            self.dataframe = new_df;
            self.globals_dataframe = new_gdf;
            self.dataframe_loaded = True;

            try:
                df, gdf, bin_array, sample_counts, global_counts, control_counts, odds_ratios, top_MHs, bottom_MHs, top_test_MHs, bottom_test_MHs, likelihood_ratios, p_samples, p_globals, p_likelihood_ratios = self.generate_results(population_data, new_target_list, base_path, new_dir, True)

                label = str(int(start_date)) + " to " + str(int(end_date))
                date_keys.append(label);
                date_bands[label] = [bin_array, sample_counts, global_counts, control_counts, odds_ratios, top_MHs, bottom_MHs, top_test_MHs, bottom_test_MHs]
            except ValueError:
                print("Insufficient data for this band")

            if end_date == minimum_date:
                break;
            start_date = end_date
            end_date -= 10000;
        
        # Latitude analysis

        minimum_lat = self.min_lat;
        maximum_lat = self.max_lat;

        latitude_bands = dict();
        latitude_keys = [];
        start_lat = maximum_lat;
        end_lat = maximum_lat-20;
        while(True):
            if end_lat < minimum_lat:
                end_lat = minimum_lat;
            new_dir = "lat_" + str(int(start_lat)) + "-" + str(int(end_lat)) + "_" + directory;
            new_target_list, f = tam.filter_targets_for_latitude(original_target_list, end_lat, start_lat, "");
            
            new_df = dataframe[(dataframe.latitude < start_lat) & (dataframe.latitude >= end_lat)]
            new_gdf = globals_dataframe[(globals_dataframe.latitude < start_lat) & (globals_dataframe.latitude >= end_lat)]

            self.dataframe = new_df;
            self.globals_dataframe = new_gdf;
            self.dataframe_loaded = True;

            try:
                df, gdf, bin_array, sample_counts, global_counts, control_counts, odds_ratios, top_MHs, bottom_MHs, top_test_MHs, bottom_test_MHs, likelihood_ratios, p_samples, p_globals, p_likelihood_ratios = self.generate_results(population_data, new_target_list, base_path, new_dir, True)

                label = str(int(start_lat)) + " to " + str(int(end_lat))
                latitude_keys.append(label)
                latitude_bands[label] = [bin_array, sample_counts, global_counts, control_counts, odds_ratios, top_MHs, bottom_MHs, top_test_MHs, bottom_test_MHs]
            except ValueError:
                print("No Geographic Points Fall in this area")
                
            if end_lat == minimum_lat:
                break;
            start_lat = end_lat
            end_lat -= 20;
# I have eliminated generation of confounder table - this should allow analysis to go forward

# =============================================================================
#         results_path = os.path.join(base_path, "results")
#         if not os.path.exists(results_path):
#             os.makedirs(results_path)
# 
#         new_path = os.path.join(results_path, directory)
#         if not os.path.exists(new_path):
#             os.makedirs(new_path)
#         
#         results_filename= os.path.join(new_path, directory + "_results.csv") 
# 
#         print("Confounder Analysis Results: " + results_filename)
#         a_file = open(results_filename, 'w');
# 
#         lat_or_MHs, lat_MH_stats, lat_MH_ps = stm.get_confounder_analysis_values(latitude_keys, latitude_bands)
#         date_or_MHs, date_MH_stats, date_MH_ps = stm.get_confounder_analysis_values(date_keys, date_bands)
# 
#         wrm.write_confounder_analysis_table(a_file, "Latitude Confounder Analysis", latitude_bands, latitude_keys, lat_or_MHs, lat_MH_stats, lat_MH_ps);
#         wrm.write_confounder_analysis_table(a_file, "Date Confounder Analysis", date_bands, date_keys, date_or_MHs, date_MH_stats, date_MH_ps);
#         a_file.close();
# 
#         plm.plot_crude_or_vs_mh_or(bin_array, orig_odds_ratios, lat_or_MHs, directory+"_latitude", new_path)
#         plm.plot_crude_or_vs_mh_or(bin_array, orig_odds_ratios, date_or_MHs, directory+"_date", new_path)
# 
#         self.dataframe_loaded = df_loaded
# 
# =============================================================================
        return "Generated Results";
    
   

    def generate_results(self, population_data, original_target_list, base_path, directory, is_confounder_analysis):


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

        header_labels = ['population_data', 'max_for_uninhabited', 'date_window', 'critical_distance', "filters_applied","globals"]
        header_values = [population_data.name, max_for_uninhabited, self.date_window, self.critical_distance, self.filters_applied,self.globals]
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
        clustered_target_list, self.dataframe, self.globals_dataframe = tam.process_targets(self.base_path, population_data, original_target_list, self.dataframe, self.globals_dataframe, self.globals, self.dataframe_loaded, self.clustering_on, self.date_window, self.critical_distance, self.critical_time, directory, self.min_date_window, self.min_lat, self.max_lat, self.min_date, self.max_date)

        if self.dataframe.empty or len(clustered_target_list)<10:
            print "Not enough samples in target area"
            f2.write("Not enough samples in Target Areas")
            f2.close()
            return "Not enough samples in target area"

        #####################
        # Process dataframe #
        #####################
        # - gets statistics (means, growth coefficients) of dataframe
        # - filters target list/dataframe by removing clusters with 0 sample means 
        # - This is ugly - we report medians in write module and we calculate them again here.
        all_sample_mediams, all_global_medians, growth_coefficients, samples_gt_globals, n_targets_gt_0, self.dataframe,growth_samples_gt_globals = stm.process_dataframe(self.dataframe, max_for_uninhabited)
       
        ########################################
        # Write filtered clustered target list #
        ########################################
        wrm.write_label(f2, "Target list after clustering")
        wrm.write_cluster_table(f2, self.dataframe, growth_coefficients)
            ################################
#         # Compute and write statistics #
#         ################################
#         # - binomial test
#         # - wilcoxon
        wrm.write_label(f2, "Statistics")
        p_binomial=binom_test(samples_gt_globals,n_targets_gt_0,0.5)
        stats_header_labels = ["Number of cases median samples>median globals", "Number of targets", "pBinomial"]
        stats_header_values = [samples_gt_globals, n_targets_gt_0, p_binomial]
        wrm.write_information(f2, stats_header_labels, stats_header_values, "   ")
        p_binomial=binom_test(growth_samples_gt_globals,n_targets_gt_0,0.5)
        stats_header_labels = ["Number of cases growth samples>growth globals", "Number of targets", "pBinomial"]
        stats_header_values = [growth_samples_gt_globals, n_targets_gt_0, p_binomial]
        wrm.write_information(f2, stats_header_labels, stats_header_values, "   ")
        # plm.plot_confounder_likelihood_ratio(self.dataframe, original_target_list)
        
        #################
        # Generate bins #
        #################
        # - extracts bin values
        # - write bin values to file
        tam.generate_merged_dataframe(base_path,directory,self.dataframe,self.globals_dataframe)
        stm.write_results(f2,directory,new_path,self.dataframe, self.globals_dataframe, population_data,self.min_globals, self.min_p)
        bin_array, sample_counts, global_counts, control_counts, odds_ratios, lower_cis, upper_cis, top_MHs, bottom_MHs, top_test_MHs, bottom_test_MHs, likelihood_ratios, p_samples, p_globals, p_controls, p_likelihood_ratios,minimum_globals = stm.generate_bin_values(self.dataframe, self.globals_dataframe, population_data)
       # wrm.write_bin_table(f2, bin_array, sample_counts, global_counts, control_counts, odds_ratios, lower_cis, upper_cis, likelihood_ratios, p_samples, p_globals, p_controls, p_likelihood_ratios,minimum_globals)
       
# =============================================================================
#         
#     
#     
 
#         t_wilcox,p_wilcox=wilcoxon(all_sample_means,all_global_means)
#         f2.write( 'Wilcoxon stat for sample vs globals means whole period:'+str(float(t_wilcox))+ '   p='+str(float(p_wilcox))+'\n')

#      
    
#         #Only include data with  globals > minimum_globals in statistical tests#
#         ###################################
        trimmed_bin_array=[]
        trimmed_p_samples=[]
        trimmed_p_globals=[]
        trimmed_sample_counts=[]
        trimmed_global_counts=[]
        trimmed_likelihood_ratios=[]
        for i in range(0, len(bin_array)):
             if global_counts>self.min_globals:
                 trimmed_bin_array.append(bin_array[i])
                 trimmed_p_samples.append(p_samples[i])
                 trimmed_p_globals.append(p_globals[i])
                 trimmed_sample_counts.append(sample_counts[i])
                 trimmed_global_counts.append(global_counts[i])
                 trimmed_likelihood_ratios.append(likelihood_ratios[i])
        if len(trimmed_global_counts)<len(global_counts)/2:
            f2.write('insufficient non-zero bins for analysis')
            return ('insufficient non-zero bins for analysis')
            
            
#         #######################
#         # Test distributions match qualitative predictions from model
#         #######################
        tests_passed=0
#         #######################
#         # Test distributions are different
#         #######################
        t_wilcox,p_wilcox=wilcoxon(trimmed_p_samples,trimmed_p_globals)
        f2.write( 'Wilcoxon stat for samples vs globals :'+str(float(t_wilcox))+ '   p='+str(float(p_wilcox))+'\n')
        ks_d,ks_p=ks_2samp(trimmed_p_samples,trimmed_p_globals)
        f2.write( 'KS test  for samples vs globals :'+str(float(ks_d))+ '   p='+str(float(ks_p))+'\n')
        if ks_p<self.min_p:
             f2.write('The two distribitions are significantly different p<0.001'+'\n')
             tests_passed=tests_passed+1
        
#             
#         ##################################
#         # Detect and display threshold and below curve #
#         ##################################
#         
# =============================================================================
#         threshold,successes,trials,p_threshold, below_curve, p_below_curve=stm.detect_threshold(trimmed_bin_array, trimmed_sample_counts,trimmed_global_counts)
#         wrm.write_label(f2,"Threshold analysis"+'\n')
#         f2.write('Threshold: '+str(threshold)+'\n')
#         f2.write('Successes: '+str(successes)+'\n')
#         f2.write('Trials: '+str(trials)+'\n')
#         f2.write('p: '+str(p_threshold)+'\n')
#         f2.write('below_curve: '+str(below_curve)+'\n')
#         f2.write('p_below_curve: '+str(p_below_curve)+'\n')
#         if p_threshold<self.min_p:
#              f2.write('There is a significant threshold effect'+'\n')
#              tests_passed=tests_passed+1
#         if p_below_curve<self.min_p:
#              f2.write('Samples_curve is significantly below globals curve '+'\n')
#              tests_passed=tests_passed+1
# =============================================================================
#         
#         
#         #################
#         # Fit data to logit curve #
#         #################
# =============================================================================
#         logit_results=stm.fit_to_logit(trimmed_bin_array, trimmed_sample_counts, trimmed_global_counts)
#         logit_predictions=stm.generate_logit_predictions(trimmed_bin_array,logit_results.params)
# #        
# =============================================================================
#        
#         
#         
#         # Generate graphs and statistics #
#         ##################################
#         # - plots likelihood ratios and sites graphs
#         # stm.generate_stats_for_ratios(bin_array, likelihood_ratios, sample_counts, global_counts, p_likelihood_ratios, p_samples, p_globals, f2, "p Likelihood Ratio", new_path, directory)
#          # would like to get rid of this routine
#         #stm.generate_stats_for_ratios(bin_array, likelihood_ratios, sample_counts, global_counts, odds_ratios, p_samples, p_globals, f2, "Odds Ratio", new_path, directory, (population_data.max_population-population_data.bin_size*2), lower_cis=lower_cis, upper_cis=upper_cis)
# 
#         # - plots p_graphs and write statistics (binomial and wilcoxon)
#         #threshold_binomial, threshold, threshold_success_count, threshold_trial_count, threshold_samples, threshold_controls = stm.generate_p_threshold_and_binomial(p_samples, p_controls, bin_array)
#        # logit_results=stm.fitToLogit(bin_array, sample_counts, global_counts)
#         plm.plot_p_graphs(trimmed_bin_array, trimmed_p_samples, trimmed_p_globals, directory, new_path)
#         plm.plot_cumulative_p_graphs(trimmed_bin_array, trimmed_p_samples, trimmed_p_globals, directory, new_path)
#         plm.plot_detection_frequencies (trimmed_bin_array, trimmed_likelihood_ratios, logit_predictions,  population_data.max_population-population_data.bin_size*2, directory, "detection_frequencies", new_path)
# =============================================================================
#         wrm.write_label(f2,"Logistic fit"+'\n')
#         f2.write("Intercept: "+str(logit_results.params[1])+'\n')
#         f2.write("Coefficient: "+str(logit_results.params[0])+'\n')
#         f2.write("AIC: "+str(logit_results.aic)+'\n')
#         f2.write("Pearson Chi2: "+str(logit_results.pearson_chi2)+'\n')
# =============================================================================
# # =============================================================================
# #         t_threshold_wilcoxon, p_threshold_wilcoxon = wilcoxon(threshold_controls, threshold_samples)
# #         wrm.write_label(f2, "Statistics for threshold bins")
# #         f2.write("Threshold: " + str(threshold))
# #         stats_header_labels = ["Number of successes", "Number of targets", "pBinomial"]
# #         stats_header_values = [threshold_success_count, threshold_trial_count, threshold_binomial]
# #         wrm.write_information(f2, stats_header_labels, stats_header_values, "   ")
# #         f2.write( 'Wilcoxon stat for pControls vs pSamples:'+str(float(t_threshold_wilcoxon))+ '   p='+str(float(p_threshold_wilcoxon))+'\n')
# # 
# # =============================================================================
#         # - plots targets and globals on a map
#         plm.plot_targets_on_map(self.dataframe, self.globals_dataframe, new_path, directory)
# 
# 
#         f2.close()
# 
# =============================================================================
        if(is_confounder_analysis):
            print("returning CA")
            return self.dataframe, self.globals_dataframe, bin_array, sample_counts, global_counts, control_counts, odds_ratios, top_MHs, bottom_MHs, top_test_MHs, bottom_test_MHs, likelihood_ratios, p_samples, p_globals, p_likelihood_ratios, tests_passed
        else:
            return "Generated results", tests_passed
# =============================================================================
# =============================================================================

def plot_population_by_time(population_data_name, time):
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

    mp.plot_population_by_time(population_data, time);

def plot_min_densities_in_time_range(population_data_name, time_from, time_to, min_density):
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

    mp.plot_min_densities_in_time_range(population_data, time_from, time_to, min_density);


def run_experiment(results_path, target_list_file, output_directory, population_data_name="Eriksson", the_globals="All", date_window=24, user_max_for_uninhabited=1, clustering_on = False, critical_distance=0, filter_date_before=-1, filter_not_direct=False, filter_not_exact=False, filter_not_figurative=False, filter_not_controversial = False, perform_cross_validation=False, number_of_kfolds = 100,  min_date_window=0, critical_time=10000, filter_min_date=-1, filter_max_date=-1, filter_min_lat=-1, filter_max_lat=-1, processed_targets=False, is_confounder_analysis=False,reweighting=False):
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
            

    mp.set_date_window(date_window)
    mp.set_min_date_window(min_date_window)
    if user_max_for_uninhabited != -1:
        mp.set_user_max_for_uninhabited(user_max_for_uninhabited)
    mp.set_clustering(clustering_on)
    mp.set_critical_distance(critical_distance)
    mp.set_critical_time(critical_time)
    mp.set_globals(the_globals)
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
        mp.set_filters_applied(filters_applied)

    # mp.set_perform_cross_validation(perform_cross_validation)
    # if perform_cross_validation:
    #     mp.set_number_of_kfolds(number_of_kfolds)
    #     mp.set_minimum_likelihood_ratio(minimum_likelihood_ratio)

    if not reweighting:
        if is_confounder_analysis:
            return mp.generate_confounder_analysis(population_data, target_list, results_path, output_directory)
        else:
            return mp.generate_results(population_data, target_list, results_path, output_directory, False)
    if reweighting:
       mp.generate_reweighted_data(population_data, target_list, base_path, output_directory)

