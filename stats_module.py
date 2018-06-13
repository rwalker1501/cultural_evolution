import numpy as np
import math
import plot_module as plm
import matplotlib.pyplot as plt
from classes_module import Target, PopulationData
from scipy.optimize import curve_fit
from scipy import stats;
from scipy.stats import linregress, chisquare, binom_test, levene
from sklearn.metrics import r2_score
from sklearn.model_selection import RepeatedKFold


def process_dataframe(data, max_for_uninhabited):

    # print(data)
    all_sample_means=[]
    all_global_means=[]
    growth_coefficients = []
    samples_gt_globals = 0
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
        # Get the mean for samples and globals #
        #########################################
        # Get every sample in the cluster..
        sample_cluster_df = cluster_df[cluster_df.type == 's']
        # take the mean of these samples
        sample_mean = sample_cluster_df['density'].mean()

        #if the sample mean is zero we do not consider the sample or the related global
        if np.isnan(sample_mean): 
            sample_mean=0
            # we delete the cluster from the original dataframe
            data = data[data.cluster_id != i]
            continue
        n_targets_gt_0 += 1

        # Get every global in the cluster..
        global_cluster_df = cluster_df[cluster_df.type == 'c']
        # take the mean of these globals
        global_mean = global_cluster_df['density'].mean()
        if np.isnan(global_mean):
            global_mean=0

        all_sample_means.append(sample_mean)
        all_global_means.append(global_mean)


        ###############################
        # Compute growth coefficients #
        ###############################

        # extract all periods and all population as arrays from the samples dataframe
        sample_times = sample_cluster_df['period'].values
        sample_populations = sample_cluster_df['density'].values

        # extract all periods and all population as arrays from the globals dataframe
        global_times = global_cluster_df['period'].values
        global_populations = global_cluster_df['density'].values

        # compute growth coefficients for samples and globals
        growth_coefficients.append(compute_growth_coefficient(sample_times, sample_populations))
        growth_coefficients.append(compute_growth_coefficient(global_times, global_populations))
            
        if sample_mean>global_mean:
            samples_gt_globals += 1

    # print(data)
    return all_sample_means, all_global_means, growth_coefficients, samples_gt_globals, n_targets_gt_0, data

def compute_growth_coefficient(times, populations):
    if len(times)>=2:
        for i in range(0, len(populations)):
            if np.isnan(populations[i]):
                populations[i] = 0
        slope, intercept, r_value, p_value, std_err = linregress(times, populations)
        return slope 
    else:
        return -1.



def generate_bin_values(dataframe, globals_dataframe, population_data, max_for_uninhabited):
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



    bin_array = []
    sample_counts = []
    global_counts = []
    control_counts = [];
    likelihood_ratios = []
    p_samples = []
    p_globals = []
    p_likelihood_ratios = []

    ##############
    # Get Totals #
    ##############
    # total samples by summing contributions
    # total globals by counting rows
    total_samples = dataframe[dataframe.type=='s']['contribution'].sum()
    total_globals = globals_dataframe['density'].count()


    p_for_filter=float(total_samples/total_globals)
    q_for_filter=1-p_for_filter
    # compute min globals necessary to get confidence of 0.95 on likelihood ratio with range of 0.0025
    top_term=1.96**2*p_for_filter*q_for_filter/0.001**2
    bottom_term=1+((1.96**2)*p_for_filter*q_for_filter)/(0.001**2*total_globals)
    minimum_globals=int(top_term/bottom_term)
    # print(top_term)
    # print(bottom_term)
   # minimum_globals=500   This is a temporary override on computed value - to be removed
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
        
        # global count: count all globals dataframe rows in the bin
        current_global_count = globals_dataframe[globals_dataframe.bin == current_bin]['density'].count()
        if np.isnan(current_global_count):
            current_global_count = 0;
        global_counts.append(current_global_count)

        control_count = current_global_count - current_sample_count;
        control_counts.append(control_count);
        
        # likelihood ratio: sample_count/global_count
        likelihood_ratio = -1
        if(current_global_count != 0):
            likelihood_ratio = float(current_sample_count)/current_global_count
        likelihood_ratios.append(likelihood_ratio)
        
        # p_sample: sample_count/total_samples
        p_sample = -1
        if total_samples > 0:
            p_sample = float(current_sample_count)/total_samples
        p_samples.append(p_sample)
        
        # p_global: global_count/total_globals
        p_global = -1
        if total_globals > 0:
            p_global = float(current_global_count)/total_globals
        p_globals.append(p_global)


        p_likelihood_ratio = -1
        if p_global > 0:
            p_likelihood_ratio = float(p_sample)/p_global

        # print(current_global_count)
        # print(minimum_globals)            
        if current_global_count <= minimum_globals:
            p_likelihood_ratio = "NA"
        p_likelihood_ratios.append(p_likelihood_ratio)

        current_bin += bin_size

    odds_ratios, top_MHs, bottom_MHs, top_test_MHs, bottom_test_MHs, upper_cis, lower_cis = generate_or_mh_ci_stats(sample_counts, global_counts, control_counts);

    return bin_array, sample_counts, global_counts, control_counts, odds_ratios, lower_cis, upper_cis, top_MHs, bottom_MHs, top_test_MHs, bottom_test_MHs, likelihood_ratios, p_samples, p_globals, p_likelihood_ratios,minimum_globals

def generate_or_mh_ci_stats(sample_counts, global_counts, control_counts):
    odds_ratios = [];
    top_MHs = [];
    bottom_MHs = [];
    top_test_MHs = [];
    bottom_test_MHs = [];
    upper_cis = [];
    lower_cis = []

    total_samples = sum(sample_counts);
    total_globals = sum(global_counts);
    total_controls = sum(control_counts);

    for i in range(0, len(sample_counts)):
        current_sample_count = sample_counts[i];
        current_global_count = global_counts[i];
        current_control_count = control_counts[i];
        
        odds_ratio = -1;
        if total_samples > 0 and total_controls > 0:
            odds_ratio = float((current_sample_count/current_control_count)/((total_samples - current_sample_count)/(total_controls - current_control_count)));
        # if current_global_count <= minimum_globals:
        #     odds_ratio = "NA";
        odds_ratios.append(odds_ratio)

        upper_ci = "NA";
        lower_ci = "NA";
        if odds_ratio != "NA" and current_sample_count > 0 and current_control_count > 0:
            ci_a = math.exp(math.log(odds_ratio) + 1.96*math.sqrt(1/current_sample_count + 1/(total_samples-current_sample_count) + 1/current_control_count + 1/(total_controls-current_control_count)));
            ci_b = math.exp(math.log(odds_ratio) - 1.96*math.sqrt(1/current_sample_count + 1/(total_samples-current_sample_count) + 1/current_control_count + 1/(total_controls-current_control_count)));
            if ci_a > ci_b:
                upper_ci = ci_a;
                lower_ci = ci_b;
            else:
                upper_ci = ci_b;
                lower_ci = ci_a;
        upper_cis.append(upper_ci);
        lower_cis.append(lower_ci);

        top_MH = -1;
        bottom_MH = -1;
        top_test_MH = -1;
        bottom_test_MH = -1;
        if total_globals > 0:
            top_MH = (current_sample_count*(total_controls-current_control_count))/total_globals;
            bottom_MH = (current_control_count*(total_samples-current_sample_count))/total_globals;
            # (samples-(globals*totalSamples)/totalGlobals)^2
            top_test_MH = current_sample_count-((current_global_count*total_samples)/total_globals)
            # (totalSamples*totalControls*globals*(totalGlobals-globals)/(totalGlobals^2*(totalGlobals-1)
            a = float(total_samples)/float(total_globals)
            b = float(total_controls)/float(total_globals)
            c = float(current_global_count)
            d = float(total_globals-current_global_count)/float(total_globals-1)
            bottom_test_MH = a*b*c*d;
        top_MHs.append(top_MH);
        bottom_MHs.append(bottom_MH)
        top_test_MHs.append(top_test_MH);
        bottom_test_MHs.append(bottom_test_MH);

    return odds_ratios, top_MHs, bottom_MHs, top_test_MHs, bottom_test_MHs, upper_cis, lower_cis


def get_confounder_analysis_values(keys, values):
    
    sample_bin_array = values[keys[0]][0];

    or_MHs = [];
    MH_statistics = [];
    MH_statistic_p_values = [];

    top_or_MH = [];
    bottom_or_MH = [];
    top_MH_statistics = [];
    bottom_MH_statistics = [];
    for bin in sample_bin_array:
        top_or_MH.append(0);
        bottom_or_MH.append(0);
        top_MH_statistics.append(0);
        bottom_MH_statistics.append(0);
    
    for key in keys:
        # data = [bin_array, sample_counts, global_counts, control_counts, odds_ratios, top_MHs, bottom_MHs, top_test_MHs, bottom_test_MHs]
        data = values[key];
        top_MHs = data[5];
        bottom_MHs = data[6];
        top_test_MHs = data[7];
        bottom_test_MHs = data[8];
        for j in range(0, len(top_MHs)):
            top_or_MH[j] += top_MHs[j];
            bottom_or_MH[j] += bottom_MHs[j]
            top_MH_statistics[j] += top_test_MHs[j];
            bottom_MH_statistics[j] += bottom_test_MHs[j];

    for i in range(0, len(sample_bin_array)):
        or_MH = top_or_MH[i]/bottom_or_MH[i];
        or_MHs.append(or_MH);

        MH_statistic = (top_MH_statistics[i]*top_MH_statistics[i])/bottom_MH_statistics[i]
        MH_statistics.append(MH_statistic);

        MH_statistic_p_values.append(1-stats.chi2.cdf(MH_statistic, 1));

    return or_MHs, MH_statistics, MH_statistic_p_values;






def trim_values(bins, likelihood_ratios, n_samples, n_globals, p_likelihood_ratios, p_samples, p_globals, lower_cis, upper_cis):
    # minimum_value = np.percentile(n_globals, removed_percentile)
    trimmed_bins = []
    trimmed_likelihood_ratios = []
    trimmed_n_samples = []
    trimmed_n_globals = []
    trimmed_p_likelihood_ratios = []
    trimmed_p_samples = []
    trimmed_p_globals = []
    trimmed_lower_cis = [];
    trimmed_upper_cis = [];
    for x in range(0, len(bins)):
        if upper_cis[x] != "NA" and lower_cis != "NA" and float(upper_cis[x]) - float(lower_cis[x]) < 1:
            trimmed_bins.append(bins[x])
            trimmed_likelihood_ratios.append(likelihood_ratios[x])
            trimmed_n_samples.append(n_samples[x])
            trimmed_n_globals.append(n_globals[x])
            trimmed_p_likelihood_ratios.append(p_likelihood_ratios[x])
            trimmed_p_samples.append(p_samples[x])
            trimmed_p_globals.append(p_globals[x])
            if len(lower_cis) > 0:
                trimmed_lower_cis.append(lower_cis[x]);
                trimmed_upper_cis.append(upper_cis[x]);

    return trimmed_bins, trimmed_likelihood_ratios, trimmed_n_samples, trimmed_n_globals, trimmed_p_likelihood_ratios, trimmed_p_samples, trimmed_p_globals, trimmed_lower_cis, trimmed_upper_cis

def divide_graph(bins, samples, the_globals, p_likelihood_ratios, divider):
    left_samples_array = []
    left_globals_array = []
    left_p_array = []
    right_samples_array = []
    right_globals_array = []
    right_p_array = []
    for x in range(0, len(bins)):
        if bins[x] < divider:
            left_samples_array.append(samples[x])
            left_globals_array.append(the_globals[x])
            left_p_array.append(p_likelihood_ratios[x])
        else:
            right_samples_array.append(samples[x])
            right_globals_array.append(the_globals[x])
            right_p_array.append(p_likelihood_ratios[x])

    return left_samples_array, left_globals_array, left_p_array, right_samples_array, right_globals_array, right_p_array

def compute_divided_graph_variance(left_array, right_array):
    left_variance = np.var(left_array)
    right_variance = np.var(right_array)
    lev_stat, lev_p = levene(left_array, right_array)
    return left_variance, right_variance, lev_stat, lev_p

def compute_divided_graph_ratio(left_samples_array, left_globals_array, right_samples_array, right_globals_array):

    left_samples_sum = float(sum(left_samples_array))
    right_samples_sum = float(sum(right_samples_array))
    
    left_globals_sum = float(sum(left_globals_array))
    right_globals_sum = float(sum(right_globals_array))


    left_p_samples = left_samples_sum/(left_samples_sum+right_samples_sum)
    right_p_samples = right_samples_sum/(left_samples_sum+right_samples_sum)
    left_p_globals = left_globals_sum/(left_globals_sum+right_globals_sum)
    right_p_globals = right_globals_sum/(left_globals_sum+right_globals_sum)

    left_p_likelihood_ratio = -1
    if left_p_globals != 0:
       left_p_likelihood_ratio = left_p_samples/left_p_globals
    right_p_likelihood_ratio = -1
    if right_p_globals != 0:
        right_p_likelihood_ratio = right_p_samples/right_p_globals

    return left_p_likelihood_ratio, right_p_likelihood_ratio

def generate_stats_for_ratios(bins,likelihood_ratios,n_samples,n_globals, ratios, p_samples, p_globals, a_file, label, file_path, directory, lower_cis = [], upper_cis = []): 

    ###############
    # Trim Arrays #
    ###############
    # trimmed arrays: values where globals <= minimum_globals

    trimmed_bins, trimmed_likelihood_ratios, trimmed_n_samples, trimmed_n_globals, trimmed_ratios, trimmed_p_samples, trimmed_p_globals, trimmed_lower_cis, trimmed_upper_cis = trim_values(bins, likelihood_ratios, n_samples, n_globals, ratios, p_samples, p_globals, lower_cis, upper_cis)

    ################################
    # Fitting and Comparing Models #
    ################################

    # linear with no intercept
    linear_parameters, p_cov =linear_fit(trimmed_bins,trimmed_ratios)
    linear_slope=linear_parameters[0]
    print "linear_slope=", linear_slope
    linear_predicted_ratios =linear_model(trimmed_bins,linear_slope)

    # threshold
    threshold_parameters,p_cov=threshold_fit(trimmed_bins,trimmed_ratios)
    threshold_threshold=threshold_parameters[0]
    threshold_gamma=threshold_parameters[1]
    threshold_slope=threshold_parameters[2]
    threshold_predicted_ratios = threshold_model(trimmed_bins, threshold_threshold, threshold_gamma,threshold_slope)

    r2_linear = r2_score(trimmed_ratios, linear_predicted_ratios)
    r2_threshold = r2_score(trimmed_ratios, threshold_predicted_ratios)

    chi_linear, p = chisquare(linear_predicted_ratios, trimmed_ratios)
    chi_threshold, p = chisquare(threshold_predicted_ratios, trimmed_ratios)

    ########################
    # Computing Statistics #
    ########################
    # - Divide graph by lambda_tau
    #   - get value for lambda_tau 
    #   - get arrays left and right of lambda_tau
    # - get variances of arrays
    # - get levene score of arrays


    lambda_tau = threshold_threshold

    left_samples_array, left_globals_array, left_p_array, right_samples_array, right_globals_array, right_p_array = divide_graph(trimmed_bins, trimmed_n_samples, trimmed_n_globals, trimmed_ratios, lambda_tau)

    left_variance, right_variance, lev_stat, lev_p = compute_divided_graph_variance(left_p_array, right_p_array)


    left_ratio, right_ratio = compute_divided_graph_ratio(left_samples_array, left_globals_array, right_samples_array, right_globals_array)


    ###############
    # Plot graphs #
    ###############
    # plm.plot_likelihood_ratio(trimmed_bins, trimmed_p_likelihood_ratios, 0, linear_predicted_p_likelihood_ratios, "linear", directory, file_path)

    # plm.plot_likelihood_ratio(trimmed_bins, trimmed_p_likelihood_ratios, threshold_threshold, threshold_predicted_p_likelihood_ratios, "threshold", directory, file_path)
    if label == "Odds Ratio":
        plm.plot_odds_ratio(trimmed_bins, trimmed_ratios, 0, linear_predicted_ratios, threshold_predicted_ratios, trimmed_lower_cis, trimmed_upper_cis, directory, label, file_path)
    else:
        plm.plot_ratio(trimmed_bins, trimmed_ratios, 0, linear_predicted_ratios, threshold_predicted_ratios, directory, label, file_path)




    ############################
    # Write Statistics to File #
    ############################
    
    a_file.write('**************************' +'\n')
    a_file.write('curve fitting: ' + label +'\n')
    a_file.write('**************************' +'\n')
    a_file.write('LMS linear model for multiplier:'+"{:8.7f}".format(lms(trimmed_ratios,linear_predicted_ratios))+'\n')
    a_file.write('LMS threshold model for multiplier:'+"{:8.7f}".format(lms(trimmed_ratios,threshold_predicted_ratios))+'\n')
    a_file.write('Coefficients for linear fit: slope=' + str(linear_slope) + "\n")
    a_file.write('Coefficients for threshold fit: gamma='+str(threshold_gamma)+' threshold=' + str(lambda_tau) + ' slope=' + str(threshold_slope) + '\n')
    a_file.write("r2_threshold: " + str(r2_threshold) + "  r2_linear: " + str(r2_linear) + "\n")
    a_file.write("chi_threshold: " + str(chi_threshold) + "  chi_linear: " + str(chi_linear) + "\n")

    a_file.write('*************************' +'\n')
    a_file.write('Divided Graph Statistics'+'\n')
    a_file.write('*************************' +'\n')
    a_file.write('Variance (lambda_tau as divider): left=' + str(left_variance) + "  right=" + str(right_variance) + "  levene_stat=" + str(lev_stat) + " levene_p=" + str(lev_p) + "\n")
    a_file.write('ratio (lambda_tau as divider): left=' + str(left_ratio) + "  right=" + str(right_ratio) + "\n")


def generate_p_threshold_and_binomial(p_samples, p_globals, bin_array):
    binomial = 100
    threshold = 0
    threshold_trial_count = 0
    threshold_success_count = 0
    threshold_samples = []
    threshold_globals = []

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
            if p_samples[x] <= p_globals[x]:
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
            threshold_globals = p_globals[:i]

    return binomial, threshold, threshold_success_count, threshold_trial_count, threshold_samples, threshold_globals


def lms(data_y,predicted_y):
    error=0
    for i in range(0,len(data_y)-1):
        error=error+(data_y[i]-predicted_y[i])**2
    return error 
    
def linear_fit(data_x,data_y):
    p_opt,p_cov=curve_fit(linear_model, data_x,data_y)
    return p_opt, p_cov   

def threshold_fit(data_x,data_y):
    if len(data_x)>3:
        p_opt,p_cov=curve_fit(threshold_model, data_x,data_y,p0=[1000.0,0.1,0.0002], bounds=([0, 0, 0], [np.inf, 25, 1]))
    else:
        p_opt=[0,0,0]
        p_cov=0
    return p_opt, p_cov 

# =============================================================================
# def threshold_fit2(data_x,data_y, max_lambda_tau):
#     p_opt,p_cov=curve_fit(threshold_model2, data_x,data_y, bounds=([0, 0, 0], [np.inf, 25, 1]))
#     return p_opt, p_cov
# =============================================================================

     

def linear_model(x,beta):
# This function implements a linear model with no intercept
    results=[]
    for an_x in x: 
        result=float(an_x*beta)
        results.append(result)
    return (results)

# =============================================================================
# def threshold_model(x,alpha,beta,):
# # This functions a conventional linear model with an intercept but all values below 0 are set to zero
#     results=[]
#     for an_x in x: 
#         result=float(alpha+an_x*beta)
#         if result<0:
#             result=0
#         results.append(result)
#         
#     return (results)
# =============================================================================

def threshold_model(x,threshold,gamma,beta):
# This instantiates a threshold model in which y grows linearly with slope b for values of x >threshold
    results=[]
    for an_x in x: 
        if an_x>=threshold:
            result=float((an_x-threshold)*beta+gamma)
        else:
            result=gamma
        results.append(result)
        
    return (results)



# =============================================================================
# def threshold_model(x, n_pop,lambda_tau):
#     results=[]
#     for an_x in x: 
#         result=n_pop*(1-lambda_tau*1/an_x)
#         if result<0:
#             result=0
#         results.append(result)
#     return(results)
# 
# def threshold_model2(x, scaling_factor, lambda_tau, alpha):
#     results=[]
#     for an_x in x: 
#         if an_x - lambda_tau < 0:
#             results.append(alpha)
#         else:
#             result=alpha+scaling_factor*(an_x - lambda_tau)
#             if result<0:
#                 result=alpha
#             results.append(result)
#     return(results)
# =============================================================================
