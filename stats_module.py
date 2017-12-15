import os
import numpy as np
import plot_module as plm
import write_module as wrm
import matplotlib.pyplot as plt
from classes_module import Target, PopulationData
from scipy.optimize import curve_fit
from scipy.stats import linregress, chisquare, binom_test, levene
from sklearn.metrics import r2_score
from sklearn.model_selection import RepeatedKFold


def process_dataframe(data, max_for_uninhabited):
    all_sample_means=[]
    all_control_means=[]
    growth_coefficients = []
    samples_gt_controls = 0
    n_targets_gt_0 = 0

    max_cluster_id = data['cluster_id'].max()
    if np.isnan(max_cluster_id):
        max_cluster_id = 0


    for i in range(0,max_cluster_id+1):

        ######################
        # For each cluster.. #
        ######################
        cluster_df = data[data.cluster_id==i]
        # consider only data where densities > max_for_uninhabited
        cluster_df = cluster_df[cluster_df.density > max_for_uninhabited]

        #########################################
        # Get the mean for samples and controls #
        #########################################
        # Get every sample in the cluster..
        sample_cluster_df = cluster_df[cluster_df.type == 's']
        # take the mean of these samples
        sample_mean = sample_cluster_df['density'].mean()

        #if the sample mean is zero we do not consider the sample or the related control
        if np.isnan(sample_mean): 
            sample_mean=0
            # we delete the cluster from the original dataframe
            data = data[data.cluster_id != i]
            continue
        n_targets_gt_0 += 1

        # Get every control in the cluster..
        control_cluster_df = cluster_df[cluster_df.type == 'c']
        # take the mean of these controls
        control_mean = control_cluster_df['density'].mean()
        if np.isnan(control_mean):
            control_mean=0

        all_sample_means.append(sample_mean)
        all_control_means.append(control_mean)


        ###############################
        # Compute growth coefficients #
        ###############################

        # extract all periods and all population as arrays from the samples dataframe
        sample_times = sample_cluster_df['period'].values
        sample_populations = sample_cluster_df['density'].values

        # extract all periods and all population as arrays from the controls dataframe
        control_times = control_cluster_df['period'].values
        control_populations = control_cluster_df['density'].values

        # compute growth coefficients for samples and controls
        growth_coefficients.append(compute_growth_coefficient(sample_times, sample_populations))
        growth_coefficients.append(compute_growth_coefficient(control_times, control_populations))
            
        if sample_mean>control_mean:
            samples_gt_controls += 1

    return all_sample_means, all_control_means, growth_coefficients, samples_gt_controls, n_targets_gt_0, data

def compute_growth_coefficient(times, populations):
    if len(times)>=2:
        for i in range(0, len(populations)):
            if np.isnan(populations[i]):
                populations[i] = 0
        slope, intercept, r_value, p_value, std_err = linregress(times, populations)
        return slope 
    else:
        return -1.



def generate_bin_values(dataframe, controls_dataframe, population_data, max_for_uninhabited, minimum_controls):
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

    # controls dataframe
    controls_dataframe['bin_index'] = (controls_dataframe.density/bin_size)-bins_to_omit
    controls_dataframe['bin_index'] = controls_dataframe.bin_index.astype(int)
    controls_dataframe['bin'] = controls_dataframe.bin_index*bin_size+minimum_bin



    bin_array = []
    sample_counts = []
    control_counts = []
    likelihood_ratios = []
    p_samples = []
    p_controls = []
    p_likelihood_ratios = []

    ##############
    # Get Totals #
    ##############
    # total samples by summing contributions
    # total controls by counting rows
    total_samples = dataframe[dataframe.type=='s']['contribution'].sum()
    total_controls = controls_dataframe['density'].count()

    #########################
    # Loop through each bin #
    #########################
    current_bin = minimum_bin
    while(current_bin < max_population):
        
        bin_array.append(current_bin)
        # sample count: for all samples in the bin, sum all contributions
        samples_dataframe = dataframe[dataframe.type=='s']
        current_sample_count = samples_dataframe[samples_dataframe.bin == current_bin]['contribution'].sum()
        if np.isnan(current_sample_count):
            current_sample_count = 0;
        sample_counts.append(current_sample_count)
        
        # control count: count all controls dataframe rows in the bin
        current_control_count = controls_dataframe[controls_dataframe.bin == current_bin]['density'].count()
        if np.isnan(current_control_count):
            current_control_count = 0;
        control_counts.append(current_control_count)

        
        # likelihood ratio: sample_count/control_count
        likelihood_ratio = -1
        if(current_control_count != 0):
            likelihood_ratio = float(current_sample_count)/current_control_count
        likelihood_ratios.append(likelihood_ratio)
        
        # p_sample: sample_count/total_samples
        p_sample = -1
        if total_samples > 0:
            p_sample = float(current_sample_count)/total_samples
        p_samples.append(p_sample)
        
        # p_control: control_count/total_controls
        p_control = -1
        if total_controls > 0:
            p_control = float(current_control_count)/total_controls
        p_controls.append(p_control)

        p_likelihood_ratio = -1
        if p_control > 0:
            p_likelihood_ratio = float(p_sample)/p_control
        if current_control_count <= minimum_controls:
            p_likelihood_ratio = "NA" 
        p_likelihood_ratios.append(p_likelihood_ratio)

        current_bin += bin_size

    return bin_array, sample_counts, control_counts, likelihood_ratios, p_samples, p_controls, p_likelihood_ratios


def trim_values(bins, likelihood_ratios, n_samples, n_controls, p_likelihood_ratios, p_samples, p_controls):
    # minimum_value = np.percentile(n_controls, removed_percentile)
    trimmed_bins = []
    trimmed_likelihood_ratios = []
    trimmed_n_samples = []
    trimmed_n_controls = []
    trimmed_p_likelihood_ratios = []
    trimmed_p_samples = []
    trimmed_p_controls = []
    for x in range(0, len(bins)):
        if p_likelihood_ratios[x] != "NA":
            trimmed_bins.append(bins[x])
            trimmed_likelihood_ratios.append(likelihood_ratios[x])
            trimmed_n_samples.append(n_samples[x])
            trimmed_n_controls.append(n_controls[x])
            trimmed_p_likelihood_ratios.append(p_likelihood_ratios[x])
            trimmed_p_samples.append(p_samples[x])
            trimmed_p_controls.append(p_controls[x])

    return trimmed_bins, trimmed_likelihood_ratios, trimmed_n_samples, trimmed_n_controls, trimmed_p_likelihood_ratios, trimmed_p_samples, trimmed_p_controls

def divide_graph(bins, samples, controls, p_likelihood_ratios, divider):
    left_samples_array = []
    left_controls_array = []
    left_p_array = []
    right_samples_array = []
    right_controls_array = []
    right_p_array = []
    for x in range(0, len(bins)):
        if bins[x] < divider:
            left_samples_array.append(samples[x])
            left_controls_array.append(controls[x])
            left_p_array.append(p_likelihood_ratios[x])
        else:
            right_samples_array.append(samples[x])
            right_controls_array.append(controls[x])
            right_p_array.append(p_likelihood_ratios[x])

    return left_samples_array, left_controls_array, left_p_array, right_samples_array, right_controls_array, right_p_array

def compute_divided_graph_variance(left_array, right_array):
    left_variance = np.var(left_array)
    right_variance = np.var(right_array)
    lev_stat, lev_p = levene(left_array, right_array)
    return left_variance, right_variance, lev_stat, lev_p

def compute_divided_graph_ratio(left_samples_array, left_controls_array, right_samples_array, right_controls_array):

    left_samples_sum = float(sum(left_samples_array))
    right_samples_sum = float(sum(right_samples_array))
    
    left_controls_sum = float(sum(left_controls_array))
    right_controls_sum = float(sum(right_controls_array))

    left_p_samples = left_samples_sum/(left_samples_sum+right_samples_sum)
    right_p_samples = right_samples_sum/(left_samples_sum+right_samples_sum)
    left_p_controls = left_controls_sum/(left_controls_sum+right_controls_sum)
    right_p_controls = right_controls_sum/(left_controls_sum+right_controls_sum)

    left_p_likelihood_ratio = -1
    if left_p_controls != 0:
       left_p_likelihood_ratio = left_p_samples/left_p_controls
    right_p_likelihood_ratio = -1
    if right_p_controls != 0:
        right_p_likelihood_ratio = right_p_samples/right_p_controls

    return left_p_likelihood_ratio, right_p_likelihood_ratio

def generate_stats_for_likelihood_ratio(bins,likelihood_ratios,n_samples,n_controls, p_likelihood_ratios, p_samples, p_controls, a_file, file_path, directory): 

    ###############
    # Trim Arrays #
    ###############
    # trimmed arrays: values where controls <= minimum_controls

    trimmed_bins, trimmed_likelihood_ratios, trimmed_n_samples, trimmed_n_controls, trimmed_p_likelihood_ratios, trimmed_p_samples, trimmed_p_controls = trim_values(bins, likelihood_ratios, n_samples, n_controls, p_likelihood_ratios, p_samples, p_controls)

    ################################
    # Fitting and Comparing Models #
    ################################

    # linear
    slope, intercept, r, p, stderr=linregress(trimmed_bins,trimmed_p_likelihood_ratios)
    linear_predicted_p_likelihood_ratios =np.asarray(linear_model2(trimmed_bins,intercept,slope))

    # threshold
    max_lambda_tau = max(trimmed_bins)
    threshold_parameters, threshold_covariance = threshold_fit2(trimmed_bins, trimmed_p_likelihood_ratios, max_lambda_tau)
    threshold_predicted_p_likelihood_ratios = threshold_model2(trimmed_bins, threshold_parameters[0], threshold_parameters[1], threshold_parameters[2])

    r2_linear = r2_score(trimmed_p_likelihood_ratios, linear_predicted_p_likelihood_ratios)
    r2_threshold = r2_score(trimmed_p_likelihood_ratios, threshold_predicted_p_likelihood_ratios)

    chi_linear, p = chisquare(linear_predicted_p_likelihood_ratios, trimmed_p_likelihood_ratios)
    chi_threshold, p = chisquare(threshold_predicted_p_likelihood_ratios, trimmed_p_likelihood_ratios)

    ########################
    # Computing Statistics #
    ########################
    # - Divide graph by lambda_tau
    #   - get value for lambda_tau 
    #   - get arrays left and right of lambda_tau
    # - get variances of arrays
    # - get levene score of arrays


    lambda_tau = threshold_parameters[1]

    left_samples_array, left_controls_array, left_p_array, right_samples_array, right_controls_array, right_p_array = divide_graph(trimmed_bins, trimmed_n_samples, trimmed_n_controls, trimmed_p_likelihood_ratios, lambda_tau)

    left_variance, right_variance, lev_stat, lev_p = compute_divided_graph_variance(left_p_array, right_p_array)


    left_p_likelihood_ratio, right_p_likelihood_ratio = compute_divided_graph_ratio(left_samples_array, left_controls_array, right_samples_array, right_controls_array)


    ###############
    # Plot graphs #
    ###############
    plm.plot_likelihood_ratio(trimmed_bins, trimmed_p_likelihood_ratios, 0, linear_predicted_p_likelihood_ratios, "linear", directory, file_path)

    plm.plot_likelihood_ratio(trimmed_bins, trimmed_p_likelihood_ratios, lambda_tau, threshold_predicted_p_likelihood_ratios, "threshold", directory, file_path)


    ############################
    # Write Statistics to File #
    ############################
    
    a_file.write('*********************' +'\n')
    a_file.write('curve fitting'+'\n')
    a_file.write('*********************' +'\n')
    a_file.write('LMS linear model for p_likelihood_ratio:'+"{:8.7f}".format(lms(trimmed_p_likelihood_ratios,linear_predicted_p_likelihood_ratios))+'\n')
    a_file.write('LMS threshold model for p_likelihood_ratios:'+"{:8.7f}".format(lms(trimmed_p_likelihood_ratios,threshold_predicted_p_likelihood_ratios))+'\n')
    a_file.write('Coefficients for linear fit: slope=' + str(slope) + " intercept=" + str(intercept) + "\n")
    a_file.write('Coefficients for threshold fit: ' + "scaling_factor=" + str(threshold_parameters[0]) + " lambda_tau=" + str(threshold_parameters[1]) + " alpha=" + str(threshold_parameters[2]) + "\n")
    a_file.write("r2_threshold: " + str(r2_threshold) + "  r2_linear: " + str(r2_linear) + "\n")
    a_file.write("chi_threshold: " + str(chi_threshold) + "  chi_linear: " + str(chi_linear) + "\n")

    a_file.write('*************************' +'\n')
    a_file.write('Divided Graph Statistics'+'\n')
    a_file.write('*************************' +'\n')
    a_file.write('Variance (lambda_tau as divider): left=' + str(left_variance) + "  right=" + str(right_variance) + "  levene_stat=" + str(lev_stat) + " levene_p=" + str(lev_p) + "\n")
    a_file.write('p_likelihood_ratio (lambda_tau as divider): left=' + str(left_p_likelihood_ratio) + "  right=" + str(right_p_likelihood_ratio) + "\n")


def generate_p_threshold_and_binomial(p_samples, p_controls, bin_array):
    binomial = 100
    threshold = 0
    threshold_trial_count = 0
    threshold_success_count = 0
    threshold_samples = []
    threshold_controls = []

    #########################
    # For every bin value.. #
    #########################
    for i in range(1, len(bin_array)):
        success = 0
        trials = 0

        ########################
        # Count success/trials #
        ########################
        for x in range(0, i):
            trials += 1
            if p_samples[x] <= p_controls[x]:
                success += 1
    
        ###########################################
        # Compute binomial for success and trials #
        ###########################################
        temp = binom_test(success, trials, 0.5)

        ################################################
        # Save values if binomial is smaller than prev #
        ################################################
        if temp < binomial:
            binomial = temp
            threshold = bin_array[i]
            threshold_trial_count = trials
            threshold_success_count = success
            threshold_samples = p_samples[:i]
            threshold_controls = p_controls[:i]

    return binomial, threshold, threshold_success_count, threshold_trial_count, threshold_samples, threshold_controls


def lms(data_y,predicted_y):
    error=0
    for i in range(0,len(data_y)-1):
        error=error+(data_y[i]-predicted_y[i])**2
    return error 
    
def threshold_fit(data_x,data_y):
    p_opt,p_cov=curve_fit(threshold_model, data_x,data_y)
    return p_opt, p_cov 

def threshold_fit2(data_x,data_y, max_lambda_tau):
    p_opt,p_cov=curve_fit(threshold_model2, data_x,data_y, bounds=([0, 0, 0], [np.inf, max_lambda_tau, 1]))
    return p_opt, p_cov

        
def linear_model(x,beta):
    results=[]
    for an_x in x: 
        result=float(an_x*beta)
        results.append(result)
    return (results)

def linear_model2(x,alpha,beta,):
    results=[]
    for an_x in x: 
        result=float(alpha+an_x*beta)
        results.append(result)
        
    return (results)

def threshold_model(x, n_pop,lambda_tau):
    results=[]
    for an_x in x: 
        result=n_pop*(1-lambda_tau*1/an_x)
        if result<0:
            result=0
        results.append(result)
    return(results)

def threshold_model2(x, scaling_factor, lambda_tau, alpha):
    results=[]
    for an_x in x: 
        if an_x - lambda_tau < 0:
            results.append(alpha)
        else:
            result=alpha+scaling_factor*(an_x - lambda_tau)
            if result<0:
                result=alpha
            results.append(result)
    return(results)

def test_threshold_against_linear(train_p_likelihood_ratios, train_bin_values, test_p_likelihood_ratios, test_bin_values, max_lambda_tau):
    # linear
    train_slope, train_intercept, train_r, train_p, train_stderr=linregress(train_bin_values,train_p_likelihood_ratios)
    test_linear_predicted=np.asarray(linear_model2(test_bin_values,train_intercept,train_slope))

    # threshold
    train_parameters, train_covariance = threshold_fit2(train_bin_values,train_p_likelihood_ratios, max_lambda_tau)
    test_threshold_predicted=np.asarray(threshold_model2(test_bin_values, train_parameters[0], train_parameters[1], train_parameters[2]))
  

    r2_linear = r2_score(test_p_likelihood_ratios, test_linear_predicted)
    r2_threshold = r2_score(test_p_likelihood_ratios, test_threshold_predicted)

    chi_linear, p = chisquare(test_linear_predicted, test_p_likelihood_ratios)
    chi_threshold, p  = chisquare(test_threshold_predicted, test_p_likelihood_ratios)

    return train_parameters, train_slope, train_intercept, r2_threshold, r2_linear, chi_threshold, chi_linear, test_threshold_predicted, test_linear_predicted


def cross_validation(file_path, dataframe, controls_dataframe, target_list, population_data, number_of_kfolds, max_for_uninhabited, minimum_controls):
    r2_linear_array = []
    r2_threshold_array = []
    chi_linear_array = []
    chi_threshold_array = []
    slopes = []
    intercepts = []
    parameter_1 = []
    parameter_2 = []
    parameter_3 = []


    bin_file = open(os.path.join(file_path, "bin_file.csv"), "w")
    cluster_list = []
    for x in range(0, len(target_list)):
        if len(target_list[x]) > 0:
            cluster_list.append(target_list[x][0].cluster_id)
    cluster_list = np.array(cluster_list)

    k_fold = RepeatedKFold(n_splits=2, n_repeats=number_of_kfolds)
    i = 0
    for train_index, test_index in k_fold.split(cluster_list):
        train_dataframe = dataframe[dataframe['cluster_id'].isin(cluster_list[train_index])]
        test_dataframe = dataframe[dataframe['cluster_id'].isin(cluster_list[test_index])]

        train_sample_means, train_control_means, train_growth_coefficients, train_samples_gt_controls, train_n_targets_gt_0, train_dataframe = process_dataframe(train_dataframe, max_for_uninhabited)
        train_bin_values, train_n_samples, train_n_controls, train_likelihood_ratios, train_p_samples, train_p_controls, train_p_likelihood_ratios = generate_bin_values(train_dataframe, controls_dataframe, population_data, max_for_uninhabited, minimum_controls)

        test_sample_means, test_control_means, test_growth_coefficients, test_samples_gt_controls, test_n_targets_gt_0, test_dataframe = process_dataframe(test_dataframe, max_for_uninhabited)
        test_bin_values, test_n_samples, test_n_controls, test_likelihood_ratios, test_p_samples, test_p_controls, test_p_likelihood_ratios = generate_bin_values(test_dataframe, controls_dataframe, population_data, max_for_uninhabited, minimum_controls)

        wrm.write_bin_table(bin_file, test_bin_values, test_n_samples, test_n_controls, test_likelihood_ratios, test_p_samples, test_p_controls, test_p_likelihood_ratios)

        test_bins, test_likelihood_ratios, test_n_samples, test_n_controls, test_p_likelihood_ratios, test_p_samples, test_p_controls = trim_values(test_bin_values, test_likelihood_ratios, test_n_samples, test_n_controls, test_p_likelihood_ratios, test_p_samples, test_p_controls)
        train_bins, train_likelihood_ratios, train_n_samples, train_n_controls, train_p_likelihood_ratios, train_p_samples, train_p_controls = trim_values(train_bin_values, train_likelihood_ratios, train_n_samples, train_n_controls, train_p_likelihood_ratios, train_p_samples, train_p_controls)

        max_lambda_tau = max(train_bins)
        train_parameters, train_slope, train_intercept, r2_threshold, r2_linear, chi_threshold, chi_linear, test_threshold_predicted, test_linear_predicted= test_threshold_against_linear(train_p_likelihood_ratios, train_bins, test_p_likelihood_ratios, test_bins, max_lambda_tau) 

        r2_linear_array.append(r2_linear)
        r2_threshold_array.append(r2_threshold)
        chi_linear_array.append(chi_linear)
        chi_threshold_array.append(chi_threshold)
        parameter_1.append(train_parameters[0])
        parameter_2.append(train_parameters[1])
        parameter_3.append(train_parameters[2])
        slopes.append(train_slope)
        intercepts.append(train_intercept)
        print(i)
        i+=1

    bin_file.close()


    plm.plot_boxplot([np.array(r2_linear_array), np.array(r2_threshold_array)], ["linear", "threshold"], "R2 - Linear vs Threshold", file_path)
    plm.plot_boxplot([np.array(parameter_1), np.array(parameter_2), np.array(parameter_3)], ['Param 1', 'Param 2', 'Param 3'], "Parameters", file_path)

    plm.plot_boxplot([np.array(chi_linear_array), np.array(chi_threshold_array)], ["linear", "threshold"], "Chi - Linear vs Threshold", file_path)

    wrm.write_validation_file(file_path, "r2", r2_linear_array, r2_threshold_array, slopes, intercepts, parameter_1, parameter_2, parameter_3)
    wrm.write_validation_file(file_path, "chi", chi_linear_array, chi_threshold_array, slopes, intercepts, parameter_1, parameter_2, parameter_3)

    return r2_linear_array, r2_threshold_array