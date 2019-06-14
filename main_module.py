"""Convert specific txt files to netcdf format"""

# pylint: disable=R0914
# THIS SEEMS LIKE UP TO DATE VERSION
# from __future__ import division
import numpy as np
import os, gc;
import sys
import json
import target_module as tam
import population_data_module as pdm
import stats_module as stm
import plot_module as plm
import write_module as wrm
import pandas as pd
from classes_module import PopulationData
from datetime import datetime
pd.options.mode.chained_assignment = None 



#==============================================================================
# 'plt.style.use('ggplot')
#==============================================================================

class MainProgram:


    def __init__(self):
        self.base_path = os.getcwd();
        self.parameters_folder = os.path.join(self.base_path, "experiment_parameters");
        self.targets_folder = os.path.join(self.base_path,"targets");


    ##################
    # Setter Methods #
    ##################


    def set_population_data_active(self, is_active, key):
        self.population_data_sources[key].is_active = is_active;
    

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

    def generate_results(self, parameters_filename="default_experiment_param.txt"):

        parameters_filepath = os.path.join(self.parameters_folder, parameters_filename);
        parameters = json.load(open(parameters_filepath));

        population_data = pdm.load_population_data_source(self.base_path, parameters['population_data']);

        target_filepath = os.path.join(self.targets_folder, parameters['target_file'])
        target_list = tam.read_target_list_from_csv(target_filepath);
        
        #####################################
        # Create directory and results file #
        #####################################
        results_path = os.path.join(self.base_path, "results")
        if not os.path.exists(results_path):
            os.makedirs(results_path)

        directory = parameters["results_directory"];
        new_path = os.path.join(results_path, directory)
        if not os.path.exists(new_path):
            os.makedirs(new_path)
        results_filename= os.path.join(new_path, directory + "_results.csv") 

        f2= open(results_filename, 'w')


        ############################
        # Write header information #
        ############################
        dateTime=str(datetime.now())
        f2.write('Date: '+dateTime)
        f2.write('\n')

        wrm.write_parameters(f2, parameters)

        target_list, targets_dataframe, globals_dataframe = tam.process_targets(self.base_path, population_data, target_list, parameters)

        if targets_dataframe.empty or len(target_list)<10:
            f2.write("Not enough sites in Target Areas")
            f2.close()
            return "Not enough sites in target area"

        print("Processing sites and controls dataframe...")
        targets_dataframe = stm.process_dataframe(targets_dataframe)

        print("Saving merged sites and globals dataframes...")
        merged_dataframe = tam.generate_merged_dataframe(self.base_path, directory, targets_dataframe, globals_dataframe, parameters['save_processed_targets']);
       
        ########################################
        # Write filtered clustered target list #
        ########################################
        print("Writing target list...")
        wrm.write_label(f2, "Filtered Target List")
        wrm.write_target_table(f2, targets_dataframe, population_data.time_window)
        
        #################
        # Generate bins #
        #################
        # - extracts bin values
        # - write bin values to file
        print("Writing bins...")
        bin_values_df = stm.generate_bin_values_dataframe(targets_dataframe, globals_dataframe, parameters['bin_size'], parameters['max_population'], parameters['min_globals'])
        wrm.write_bin_table(f2, bin_values_df, parameters['min_globals'])


        #################
        # Stat Analysis #
        #################
        # - n, median, mean, std of samples and globals
        # - K-S2 test for bin distribution of samples (p_samples) vs. bin distribution of globals (p_globals)
        print("Calculating statistics...")
        stat_dictionary, trimmed_bin_values_df = stm.generate_statistics(targets_dataframe, globals_dataframe, bin_values_df, parameters['min_globals'])
        
        if stat_dictionary is None:
            f2.write('insufficient non-zero bins for analysis');
            return 'insufficient non-zero bins for analysis';
        
        print("Writing analysis...")
        wrm.write_analysis(f2, stat_dictionary, parameters['min_p']);

        
        ###############
        # Plot Graphs #
        ###############

        print("Generating graphs...")
        plm.plot_stat_graphs(stat_dictionary, trimmed_bin_values_df, population_data, parameters['bin_size'], parameters['max_population'], directory, new_path);

        # plots targets and globals on a map
        plm.plot_targets_on_map(targets_dataframe, globals_dataframe, new_path, directory)
        
        ###############
        # Compare likelihoods of epidemiological, linear and constant models 
        ###############
        print("Computing likelihoods Models")
        models=('richard', 'linear','constant')
        max_likelihood=np.zeros(len(models))
        for i in range(0,len(models)):
            print("model= " + models[i])
            max_lambda, max_zetta, max_eps, max_likelihood[i], opt_threshold=stm.compute_likelihood_model(directory, results_path, population_data,merged_dataframe, models[i], parameters)
            write_likelihood_results(f2,max_lambda, max_zetta, max_eps, max_likelihood[i], opt_threshold,models[i] );
            gc.collect();
        epid_over_linear=np.exp(max_likelihood[0]-max_likelihood[1])
        epid_over_constant=np.exp(max_likelihood[0]-max_likelihood[2])
        wrm.write_label(f2,'Bayes factors')
        f2.write( 'Bayes factor epidemiological over linear='+'{:.3g}'.format(epid_over_linear)+'\n')
        f2.write( 'Bayes factor epidemiological over constant='+'{:.3g}'.format(epid_over_constant)+'\n')
        f2.close();
        return max_likelihood

        
      

def write_likelihood_results(aFile,max_lambda, max_zetta, max_eps, max_likelihood, interpolated_lambdas,model ):
        wrm.write_label(aFile, "Results of max likelihood analysis for "+model+" model")
        if model=='epidemiological' or model=='richard':
            aFile.write('Relative lambda 0.025='+ '{:.2f}'.format(interpolated_lambdas[0])+"\n")
            aFile.write('Relative lambda 0.25='+ '{:.2f}'.format(interpolated_lambdas[1])+"\n")
            aFile.write('Relative lambda 0.5='+ '{:.2f}'.format(interpolated_lambdas[2])+"\n")
            aFile.write('Relative lambda 0.75='+ '{:.2f}'.format(interpolated_lambdas[3])+"\n")
            aFile.write('Relative lambda 0.975='+ '{:.2f}'.format(interpolated_lambdas[4])+"\n")
            threshold_low=interpolated_lambdas[0]**2
            threshold_high=interpolated_lambdas[4]**2
            aFile.write('0.025 CI for threshold='+ '{:.2f}'.format(threshold_low)+"\n")
            aFile.write('0.975 CI for threshold='+ '{:.2f}'.format(threshold_high)+"\n")
        if model=='epidemiological' or model=='linear' or model=='richard':
            aFile.write("Max eps="+'{:.5f}'.format(max_eps)+"\n")
#            aFile.write("Max comm="+'{:.2f}'.format(max_comm)+"\n")
        aFile.write("Max zetta="+'{:.7f}'.format(max_zetta)+"\n")
        aFile.write("Max likelihood="+'{:.0f}'.format(max_likelihood)+"\n")
        if model=='epidemiological':
            k=3
        else:
            if model=='richard':
                k=3
            else:
                if model=='linear':
                    k=2
                else:
                    if model=='constant':
                        k=1                               
        aFile.write("AIC="+ '{:.2f}'.format(2*k-2*max_likelihood)+"\n")
    

def run_experiment(parameters_filename=""):

    mp = MainProgram()
    if parameters_filename == "":
        mp.generate_results();
    else:
        mp.generate_results(parameters_filename=parameters_filename)
    gc.collect();

