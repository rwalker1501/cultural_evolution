import numpy as np
import math
import copy
import plot_module as plm
import matplotlib.pyplot as plt
import statsmodels.api as sm
import statsmodels.formula.api as smf
import write_module as wrm
from mpmath import mpf, exp
from classes_module import Target, PopulationData
from scipy.optimize import curve_fit
from scipy import stats;
from scipy.optimize import minimize
from scipy.stats import linregress, chisquare, binom_test, levene,wilcoxon,ks_2samp
from sklearn.metrics import r2_score
from sklearn.model_selection import RepeatedKFold
from sklearn.linear_model import LogisticRegression
from operator import itemgetter




def process_dataframe(data, max_for_uninhabited):

    # print(data)
    all_sample_medians=[]
    all_global_medians=[]
    growth_coefficients = []
    samples_gt_globals = 0
    growth_samples_gt_globals=0
    samples_positive_growth=0
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
        # Get the median for samples and globals #
        #########################################
        # Get every sample in the cluster..
        sample_cluster_df = cluster_df[cluster_df.type == 's']
        # take the mean of these samples
        sample_median = sample_cluster_df['density'].median()

        #if the sample mean is zero we do not consider the sample or the related global
        if np.isnan(sample_median): 
            sample_median=0
            # we delete the cluster from the original dataframe
            data = data[data.cluster_id != i]
            continue
        n_targets_gt_0 += 1

        # Get every global in the cluster..
        global_cluster_df = cluster_df[cluster_df.type == 'c']
        # take the mean of these globals
        global_median = global_cluster_df['density'].median()
        if np.isnan(global_median):
            global_median=0
        if sample_median>global_median:
            samples_gt_globals += 1
        all_sample_medians.append(sample_median)
        all_global_medians.append(global_median)


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
        growth_coefficient_samples=compute_growth_coefficient(sample_times, sample_populations)
        growth_coefficient_globals=compute_growth_coefficient(global_times, global_populations)
        if growth_coefficient_samples>growth_coefficient_globals:
           growth_samples_gt_globals+=1 
        if growth_coefficient_samples>0:
            samples_positive_growth+=1
        growth_coefficients.append(growth_coefficient_samples)
        growth_coefficients.append(growth_coefficient_globals)           
    # print(data)
    return all_sample_medians, all_global_medians, growth_coefficients, samples_gt_globals, n_targets_gt_0, data,growth_samples_gt_globals,samples_positive_growth

def compute_growth_coefficient(times, populations):
    if len(times)>=2:
        for i in range(0, len(populations)):
            if np.isnan(populations[i]):
                populations[i] = 0
        slope, intercept, r_value, p_value, std_err = linregress(times, populations)
# times increase towards past - so growing population gives negative slope - we invert slope to give intuitive reading 
        return slope * -1
    else:
        return -1.
    
# =============================================================================
# def detect_threshold(bins,n_samples,n_globals):
#     cum_samples_list=[]
#     cum_globals_list=[]
#     cum_samples=0
#     cum_globals=0
#  #   print 'In detect threshold'
#     for i in range(0,len(bins)):
#         cum_samples=cum_samples+n_samples[i]
#         cum_globals=cum_globals+n_globals[i]
#         cum_samples_list.append(cum_samples)
#         cum_globals_list.append(cum_globals)
#     total_samples=cum_samples_list[len(bins)-1]
#     total_globals=cum_globals_list[len(bins)-1]
#     print 'total_samples: ',str(total_samples)
#     print 'total_globals: ', str(total_globals)
#     p_sample=total_samples/total_globals
#     min_p=1
#     below_curve=0
# #    print 'p_sample: ',str(p_sample)
#     for i in range(0,len(bins)):
#         p=stats.binom_test(cum_samples_list[i],cum_globals_list[i],p_sample)
#  #       print 'Threshold: ',str(bins[i]), 'samples: ',str(cum_samples_list[i]),' globals: ',str(cum_globals_list[i]),'p: ',str(p),'p_samples',str(cum_samples_list[i]/float(total_samples)),'p_globals: ',str(cum_globals_list[i]/float(total_globals))
#         if p<min_p:
#             min_p=p
#             threshold=bins[i]
#             successes=cum_samples_list[i]
#             trials=cum_globals_list[i]
#         if cum_samples_list[i]/float(total_samples)<cum_globals_list[i]/float(total_globals):
#             below_curve=below_curve+1
#     p_below_curve=stats.binom_test(below_curve,len(bins))
# #    print "Threshold: ",str(threshold)
# #    print "Successes: ",str(successes)
# #    print "Trials: ", str(trials)
# #    print "p: ", str(min_p)
# #    print "below_curve: ",str(below_curve)
# #    print "p_below_curve: ",str(p_below_curve)
#     return(threshold,successes,trials,min_p,below_curve,p_below_curve)
#     
# =============================================================================
        
        
            
            

def fit_to_logit(bins, sample_counts, global_counts,bin_size):
# This needs to be checked
    add=bin_size/2
    bins=[x+add for x in bins]
    bins=sm.add_constant(bins, prepend=False)
 #   print "bins"
 #   for i in range(0,5):
 #       aRecord=bins[i]
 #       for j in range(0,len(aRecord)):
 #           print bins[j]," "
 #       print "\nl"
    sample_globals=[]
    for i in range(0,len(sample_counts)):
         # I add 0.000001 to all counts to make fit work - temporary fix for bug in statsModels
        sample_globals.append([sample_counts[i]+0.000001,global_counts[i]+0.000001])
    glm_binom = sm.GLM(sample_globals, bins, family=sm.families.Binomial())
    res = glm_binom.fit()
 #   print(res.summary())
    return(res)
    
def fit_to_linear (bins, sample_counts, global_counts):
 #   bins=sm.add_constant(bins, prepend=False)
    probs=[]
    for  i in range(0,len(sample_counts)):
        aProb=sample_counts[i]/global_counts[i]
        probs.append(aProb)
    res=sm.OLS(probs, bins).fit()
#    print(res.summary())
    return(res)
        
def regress_model_with_innovation(parameters, *args):
#   This function computes the logLikelihood of our main model
#   Save initial parameter guesses
#   A lot of type problems with numpy arrays here - needs cleaning up
    
    samples_in=np.asarray(args[0])
    globals_in=np.asarray (args[1])
    densities_in=np.asarray(args[2])+50
 #   print "densities="
 #   print densities_in
    death=parameters[0]
    mu=parameters[1]
    zetta=parameters[2]
    iStar=np.asarray((mu/(mu+death))*densities_in)
#    print "istar="
#    print iStar
#   I am not sure there shouldn't be a log in second term in following expression
    logLik=-(np.sum(np.multiply(samples_in,np.log(zetta*iStar)))-np.sum(np.multiply(globals_in,iStar)))

    return logLik

def fit_model_with_innovation (samples_in, globals_in, densities_in):
    initial_death=35
    initial_mu=0.24
    initial_zetta=4.6*10**-5
    init_parameters=[initial_death,initial_mu,initial_zetta]
    #  method does not support bounds on parameters which must all be positive - have to change method
    results=minimize(regress_model_with_innovation, init_parameters, args=(samples_in, globals_in, densities_in),method='nelder-mead')
    print '*************'
    print 'x results from model'
    print '*************'
    print results.x
    return results    
    
    
def generate_logit_predictions(bins,params):
    predictions=[]
    for a_Bin in bins:
        a_prediction=(math.exp(params[1]+params[0]*a_Bin))/(1+(math.exp(params[1]+params[0]*a_Bin)))
        predictions.append(a_prediction)
    return (predictions)
    
def generate_linear_predictions(bins,params):
    predictions=[]
    for a_Bin in bins:
 #       a_prediction=(params[1]+params[0]*a_Bin)
        a_prediction=(params[0]*a_Bin)
        predictions.append(a_prediction)
    return(predictions)

def generate_bin_values(dataframe, globals_dataframe, population_data):
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
    #we add bin_size/2 to get midpoint of each bin
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
# =============================================================================
#     p_for_filter=float(total_samples/total_globals)
#     minimum_globals=int(1/p_for_filter)
#     # This is for eriksson - will need to change for Timmermann
#     margin_of_error=0.0003
#     q_for_filter=1-p_for_filter
#     # compute min globals necessary to get confidence of 0.95 on likelihood ratio with MOE of 0.0001
#     x=1.96**2*p_for_filter*q_for_filter/margin_of_error**2
#     minimum_globals=total_globals*x/(x+total_globals-1)
#  We get rid of bins with no globals
    minimum_globals=1
# =============================================================================
# =============================================================================
#     top_term=1.96**2*p_for_filter*q_for_filter/0.001**2
#     bottom_term=1+((1.96**2)*p_for_filter*q_for_filter)/(0.0001**2*total_globals)
# =============================================================================
    
    # Don't allow zero values for minimum globals
    if minimum_globals==0:
        minimum_globals=1
    # print(top_term)
    # print(bottom_term)
   # minimum_globals=500   This is a temporary override on computed value - to be removed
    #########################
    # Loop through each bin . data untrimmed - would be better to filter here#
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

        current_bin += bin_size

    odds_ratios, top_MHs, bottom_MHs, top_test_MHs, bottom_test_MHs, upper_cis, lower_cis = generate_or_mh_ci_stats(sample_counts, global_counts, control_counts);

    p_controls = (np.array(control_counts)/sum(control_counts)).tolist();

    for i in range(0, len(bin_array)):
        p_sample = p_samples[i]
        p_control = p_controls[i]
        p_likelihood_ratio = -1
        if p_control > 0:
            p_likelihood_ratio = float(p_sample)/p_control
        p_likelihood_ratios.append(p_likelihood_ratio)

    return bin_array, sample_counts, global_counts, control_counts, odds_ratios, lower_cis, upper_cis, top_MHs, bottom_MHs, top_test_MHs, bottom_test_MHs, likelihood_ratios, p_samples, p_globals, p_controls, p_likelihood_ratios, minimum_globals

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
        
        odds_ratio = float('nan');
        if total_samples > 0 and total_controls > 0:
            odds_ratio = float((current_sample_count/current_control_count)/((total_samples - current_sample_count)/(total_controls - current_control_count)));
        # if current_global_count <= minimum_globals:
        #     odds_ratio = "NA";
        odds_ratios.append(odds_ratio)

        upper_ci = "NA";
        lower_ci = "NA";
        if not math.isnan(odds_ratio) and current_sample_count > 0 and current_control_count > 0:
            ci_a=0
            ci_b=0
          #   ci_a = math.exp(math.log(odds_ratio) + 1.96*math.sqrt(1/current_sample_count + 1/(total_samples-current_sample_count) + 1/current_control_count + 1/(total_controls-current_control_count)));
          #   ci_b = math.exp(math.log(odds_ratio) - 1.96*math.sqrt(1/current_sample_count + 1/(total_samples-current_sample_count) + 1/current_control_count + 1/(total_controls-current_control_count)));
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
            top_MH=0
            bottom_MH=0
            #top_MH = (current_sample_count*(total_controls-current_control_count))/total_globals;
            #bottom_MH = (current_control_count*(total_samples-current_sample_count))/total_globals;
            # (samples-(globals*totalSamples)/totalGlobals)^2
            #top_test_MH = current_sample_count-((current_global_count*total_samples)/total_globals)
            top_test_MH = 0
            # (totalSamples*totalControls*globals*(totalGlobals-globals)/(totalGlobals^2*(totalGlobals-1)
           # a = float(total_samples)/float(total_globals)
           # b = float(total_controls)/float(total_globals)
           # c = float(current_global_count)
           # d = float(total_globals-current_global_count)/float(total_globals-1)
            bottom_test_MH =0
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






def trim_values(bins, likelihood_ratios, n_samples, n_globals, ratios, p_samples, p_globals, lower_cis, upper_cis):
    
    # trimmed_bins = bins
    # trimmed_likelihood_ratios = likelihood_ratios
    # trimmed_n_samples = n_samples
    # trimmed_n_globals = n_globals
    # trimmed_p_likelihood_ratios = p_likelihood_ratios
    # trimmed_p_samples = p_samples
    # trimmed_p_globals = p_globals
    # trimmed_lower_cis = lower_cis;
    # trimmed_upper_cis = upper_cis;

    trimmed_bins = []
    trimmed_likelihood_ratios = []
    trimmed_n_samples = []
    trimmed_n_globals = []
    trimmed_ratios = []
    trimmed_p_samples = []
    trimmed_p_globals = []
    trimmed_lower_cis = [];
    trimmed_upper_cis = [];
    for x in range(0, len(bins)):
        if not math.isnan(ratios[x]):
            trimmed_bins.append(bins[x])
            trimmed_likelihood_ratios.append(likelihood_ratios[x])
            trimmed_n_samples.append(n_samples[x])
            trimmed_n_globals.append(n_globals[x])
            trimmed_ratios.append(ratios[x])
            trimmed_p_samples.append(p_samples[x])
            trimmed_p_globals.append(p_globals[x])
            if len(lower_cis) > 0:
                trimmed_lower_cis.append(lower_cis[x]);
                trimmed_upper_cis.append(upper_cis[x]);

    return trimmed_bins, trimmed_likelihood_ratios, trimmed_n_samples, trimmed_n_globals, trimmed_ratios, trimmed_p_samples, trimmed_p_globals, trimmed_lower_cis, trimmed_upper_cis

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

# =============================================================================
# def generate_stats_for_ratios(bins,likelihood_ratios,n_samples,n_globals, ratios, p_samples, p_globals, a_file, label, file_path, directory, max_xaxis, lower_cis = [], upper_cis = []): 
#     # There is a lot of code here we don't use
#     ###############
#     # Trim Arrays #
#     ###############
#     # trimmed arrays: values where globals <= minimum_globals
# 
#     trimmed_bins, trimmed_likelihood_ratios, trimmed_n_samples, trimmed_n_globals, trimmed_ratios, trimmed_p_samples, trimmed_p_globals, trimmed_lower_cis, trimmed_upper_cis = trim_values(bins, likelihood_ratios, n_samples, n_globals, ratios, p_samples, p_globals, lower_cis, upper_cis)
# 
# 
#     ################################
#     # Fitting and Comparing Models #
#     ################################
# 
#     # linear with no intercept
#     linear_parameters, p_cov =linear_fit(trimmed_bins,trimmed_ratios)
#     linear_slope=linear_parameters[0]
#     print "linear_slope=", linear_slope
#     linear_predicted_ratios =linear_model(trimmed_bins,linear_slope)
# 
#     # threshold
#     threshold_parameters,p_cov=threshold_fit(trimmed_bins,trimmed_ratios)
#     threshold_threshold=threshold_parameters[0]
#     threshold_gamma=threshold_parameters[1]
#     threshold_slope=threshold_parameters[2]
#     threshold_predicted_ratios = threshold_model(trimmed_bins, threshold_threshold, threshold_gamma,threshold_slope)
# 
#     r2_linear = r2_score(trimmed_ratios, linear_predicted_ratios)
#     r2_threshold = r2_score(trimmed_ratios, threshold_predicted_ratios)
# 
#     chi_linear, p = chisquare(linear_predicted_ratios, trimmed_ratios)
#     chi_threshold, p = chisquare(threshold_predicted_ratios, trimmed_ratios)
# 
#     ########################
#     # Computing Statistics #
#     ########################
#     # - Divide graph by lambda_tau
#     #   - get value for lambda_tau 
#     #   - get arrays left and right of lambda_tau
#     # - get variances of arrays
#     # - get levene score of arrays
# 
# 
#     lambda_tau = threshold_threshold
# 
#     left_samples_array, left_globals_array, left_p_array, right_samples_array, right_globals_array, right_p_array = divide_graph(trimmed_bins, trimmed_n_samples, trimmed_n_globals, trimmed_ratios, lambda_tau)
# 
#     left_variance, right_variance, lev_stat, lev_p = compute_divided_graph_variance(left_p_array, right_p_array)
# 
# 
#     left_ratio, right_ratio = compute_divided_graph_ratio(left_samples_array, left_globals_array, right_samples_array, right_globals_array)
# 
#     ###############
#     # Plot graphs #
#     ###############
#     # plm.plot_likelihood_ratio(trimmed_bins, trimmed_p_likelihood_ratios, 0, linear_predicted_p_likelihood_ratios, "linear", directory, file_path)
# 
#     # plm.plot_likelihood_ratio(trimmed_bins, trimmed_p_likelihood_ratios, threshold_threshold, threshold_predicted_p_likelihood_ratios, "threshold", directory, file_path)
#     logit_results=fit_to_logit(bins, n_samples, n_globals)
#     predictions=generate_logit_predictions(bins,logit_results.params)
#     if label == "Odds Ratio":
#   #      plm.plot_odds_ratio(trimmed_bins, trimmed_ratios, predictions, 0, linear_predicted_ratios, threshold_predicted_ratios, trimmed_lower_cis, trimmed_upper_cis, max_xaxis, directory, label, file_path)
#   #   Still need to trim predictions and fix label
#         plm.plot_detection_frequencies (trimmed_bins, trimmed_likelihood_ratios, predictions,  max_xaxis, directory, label, file_path)
#     else:
#         plm.plot_ratio(trimmed_bins, trimmed_ratios, 0, linear_predicted_ratios, threshold_predicted_ratios,max_xaxis, directory, label, file_path)
# 
# 
# 
# 
#     ############################
#     # Write Statistics to File #
#     ############################
#     
#     a_file.write('**************************' +'\n')
#     a_file.write('curve fitting: ' + label +'\n')
#     a_file.write('**************************' +'\n')
#     a_file.write('LMS linear model for multiplier:'+"{:8.7f}".format(lms(trimmed_ratios,linear_predicted_ratios))+'\n')
#     a_file.write('LMS threshold model for multiplier:'+"{:8.7f}".format(lms(trimmed_ratios,threshold_predicted_ratios))+'\n')
#     a_file.write('Coefficients for linear fit: slope=' + str(linear_slope) + "\n")
#     a_file.write('Coefficients for threshold fit: gamma='+str(threshold_gamma)+' threshold=' + str(lambda_tau) + ' slope=' + str(threshold_slope) + '\n')
#     a_file.write("r2_threshold: " + str(r2_threshold) + "  r2_linear: " + str(r2_linear) + "\n")
#     a_file.write("chi_threshold: " + str(chi_threshold) + "  chi_linear: " + str(chi_linear) + "\n")
# 
#     a_file.write('*************************' +'\n')
#     a_file.write('Divided Graph Statistics'+'\n')
#     a_file.write('*************************' +'\n')
#     a_file.write('Variance (lambda_tau as divider): left=' + str(left_variance) + "  right=" + str(right_variance) + "  levene_stat=" + str(lev_stat) + " levene_p=" + str(lev_p) + "\n")
#     a_file.write('ratio (lambda_tau as divider): left=' + str(left_ratio) + "  right=" + str(right_ratio) + "\n")
# 
# =============================================================================

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

# =============================================================================
# def fit_to_threshold_model(data_x,sample_counts, global_counts):
#     data_y=[]
#     for i in range(0,len(sample_counts)):
#         data_y.append(float(sample_counts[i]/global_counts[i]))
#     if len(data_x)>3:
# #p0=[1000.0,0.1,0.0002],
# # guessed threshold needs to be reset for Timmermann data
#  #       p_opt,p_cov=curve_fit(threshold_model, data_x,data_y, p0=[1400.0,0.001,0.002],bounds=([0, 0, 0]), diag=(1./data_x.mean(),1./data_y.mean()), [np.inf, 25, 1])
#         diag1=1./np.asarray(data_x).mean()
#         diag2=1./(np.asarray(data_y).mean())
#         print "diag1=",str(diag1)
#         print "diag2=",str(diag2)
#         p_opt,p_cov=curve_fit(threshold_model, data_x,data_y,p0=[1400.0,0.0000001])
#     else:
#         p_opt=[0,0]
#         p_cov=0
#     return p_opt, p_cov 
# =============================================================================

def logit_values_using_GLM(sample_tuples, global_tuples):

    unique_densities = [i[0] for i in sample_tuples];

    values=sm.add_constant(unique_densities, prepend=False)

    sample_globals = []
    likelihood_ratios = [];
    for i in range(len(sample_tuples)):
        sample_count = sample_tuples[i][1];
        global_count = global_tuples[i][1];
         # I add 0.000001 to all counts to make fit work - temporary fix for bug in statsModels
        sample_globals.append([sample_count+0.000001,global_count+0.000001])

    glm_binom = sm.GLM(sample_globals, values, family=sm.families.Binomial())
    res = glm_binom.fit()
    print(res.params)

    loss = generate_logit_predictions(unique_densities, res.params);
    return unique_densities, res.params[0], res.params[1], loss; 

def generate_tuples_from_dataframe(merged_dataframe):
    unique_densities = np.sort(merged_dataframe.density.unique());
    zero_array = [0 for i in unique_densities];
    sample_base_dict = dict(zip(unique_densities,zero_array));
    global_base_dict = dict(zip(unique_densities,zero_array));

    merged_dataframe['count'] = 0;
    sample_dataframe = merged_dataframe[merged_dataframe.type == 's'];
    global_dataframe = merged_dataframe[merged_dataframe.type == 'g'];
    
    sample_tuples = generate_sorted_dataframe_tuples(sample_base_dict, sample_dataframe);
    global_tuples = generate_sorted_dataframe_tuples(global_base_dict, global_dataframe);

    return sample_tuples, global_tuples;

def logit_values(sample_tuples, global_tuples):

    x = [];
    y = [];
    for i in range(len(sample_tuples)):
        density = sample_tuples[i][0];
        sample_count = sample_tuples[i][1];
        global_count = global_tuples[i][1];

        # sample > 0 check, all densities
        x.append([density]);
        y.append(1 if sample_count > 0 else 0);

    x = np.array(x);
    y = np.array(y);

    logit = LogisticRegression();
    logit.fit(x,y);
    coef = logit.coef_;
    intercept = logit.intercept_;
    loss = logistic_model_2(x*coef + intercept).ravel();
    return x.ravel(), y, coef, intercept, loss;


    
def logistic_model(x):
    return 1 / (1 + np.exp(-x));

def logistic_model_2(x):
    return (np.exp(x))/(1+(np.exp(x)))

def generate_sorted_dataframe_tuples(base_dict, dataframe):
    print("Generate sorted dataframe tuples")
    dataframe['count'] = 0;
    dataframe = dataframe.groupby(['density']).count().reset_index();
    # print(dataframe);
    dataframe_dict = dict(zip(dataframe['density'], dataframe['count']));
    base_dict.update(dataframe_dict)
    sorted_dataframe_tuples = sorted(base_dict.iteritems(), key=itemgetter(0));

    return sorted_dataframe_tuples;


def generate_cumulated_densities(base_dict, dataframe):
    sorted_dataframe_tuples = generate_sorted_dataframe_tuples(base_dict, dataframe);

    cum_values = [0 for x in base_dict.keys()];
    cum_values[0] = sorted_dataframe_tuples[0][1];
    for i in range(1, len(sorted_dataframe_tuples)):
        cum_values[i] = cum_values[i-1] + sorted_dataframe_tuples[i][1];

    return cum_values;



def generate_cumulated_detection_frequency(merged_dataframe):
    unique_densities = np.sort(merged_dataframe.density.unique());
    zero_array = [0 for x in unique_densities];

    samples_base_dict = dict(zip(unique_densities,zero_array));
    globals_base_dict = dict(zip(unique_densities,zero_array));

    samples = merged_dataframe[merged_dataframe.type == 's']
    cum_samples = generate_cumulated_densities(samples_base_dict, samples)

    the_globals = merged_dataframe[merged_dataframe.type == 'g']
    cum_globals = generate_cumulated_densities(globals_base_dict, the_globals);


    print("unique_densities")
    print(unique_densities)
    print("Cumulative samples")
    print(cum_samples);
    print("Cumulative globals")
    print(cum_globals);
    cum_det_freq = [];
    for i in range(len(cum_globals)):
        cum_det_freq.append(float(cum_samples[i])/float(cum_globals[i]));

    return unique_densities, cum_samples, cum_globals, cum_det_freq;

def optimized_pettitt_test(data_x, data_y):
    u=[]
    max_u=0
    threshold=0
    if len(data_x)>3:
        total_t=len(data_x)
        
        u.append(0);
        for j in range(1, total_t):
            u[0]=u[0]+np.sign(data_y[0]-data_y[j])
        if abs(u[0]) > max_u:
            max_u=abs(u[0])
            threshold=data_x[0];
            
        
        for t in range(1,total_t-1):
            to_remove = 0;
            # print("REMOVE");
            for i in range(0, t):
                # print(str(i) + " " + str(t));
                to_remove += np.sign(data_y[i]-data_y[t])
            
            to_add = 0;
            # print("ADD");
            for j in range(t+1, total_t):
                # print(str(t) + " " + str(j));
                to_add += np.sign(data_y[t]-data_y[j]);
            
            u.append(0);
            u[t] = u[t-1] - to_remove + to_add;
            
            if abs(u[t])>max_u:
                max_u=abs(u[t])
                threshold=data_x[t];
        num_p = float(-6*max_u**2);
        denom_p = float(total_t**3+total_t**2);
        print("num: " + str(num_p))
        print("denom: " + str(denom_p))
        p=2*exp(num_p/denom_p);
        print("p: " + str(p))

        # print(u);
    else:
        threshold=0
        p=1
    return threshold,p 

def pettitt_test(data_x, data_y):
    threshold = 0
    max_u = 0
    u = []
    if len(data_x)>3:
        total_t=len(data_x)
        for t in range(0,total_t):
            u.append(0)
            for i in range(0,t):
                for j in range(t+1,total_t):
                        u[t]=u[t]+np.sign(data_y[i]-data_y[j])
            if abs(u[t])>max_u:
                max_u=abs(u[t])
                threshold=data_x[t]
        p=2*math.exp((-6*max_u**2)/(total_t**3+total_t**2))                      
    else:
        threshold=0
        p=1
    return threshold,p 

def detect_threshold(data_x,sample_counts, global_counts,bin_size):
    
# This finds a threshold using pettitt's test
    add=bin_size/2
    data_x=[x+add for x in data_x]
    u=[]
    data_y=[]
    max_u=0
    threshold=0
    for i in range(0,len(sample_counts)):
        if global_counts[i]==0:
            print "global counts =0", "i=",str(i),'\n'
        data_y.append(float(sample_counts[i]/global_counts[i]))
    if len(data_x)>3:
        total_t=len(sample_counts)
        for t in range(0,total_t):
            u.append(0)
            for i in range(0,t):
                for j in range(t+1,total_t):
                        u[t]=u[t]+np.sign(data_y[i]-data_y[j])
  #          print "candidate threshold=",str(data_x[t]), " u=",str(u[t]),'/n'
            if abs(u[t])>max_u:
                max_u=abs(u[t])
                threshold=data_x[t]
        p=2*math.exp((-6*max_u**2)/(total_t**3+total_t**2))
 #       print "threshold=",str(threshold)," p=", str(float(p)),"/n"                         
    else:
        threshold=0
        p=1
    return threshold,p 



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

# =============================================================================
# def threshold_model(x,threshold,beta, zeta):
# # This instantiates a  model exactly based on the theoretical model
#     results=[]
#     for an_x in x: 
#         if an_x>=threshold:
#             istar=float((an_x-threshold)*beta)
#             result=1-(1-zeta)**istar
#         else:
#             result=0
#         results.append(result)
#         
#     return (results)
# =============================================================================

def threshold_model(x,threshold,beta):
# This instantiates a  linear threshold model appoximating the theoretical model
    results=[]
    for an_x in x: 
        if an_x>=threshold:
            result=(an_x-threshold)*beta
        else:
            result=0
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
    
def write_results(aFile,anIdentifier, aPath,dataframe, globals_dataframe,population_data, min_globals, min_p):
    
    bin_array, sample_counts, global_counts, control_counts, odds_ratios, lower_cis, upper_cis, top_MHs, bottom_MHs, top_test_MHs, bottom_test_MHs, likelihood_ratios, p_samples, p_globals, p_controls, p_likelihood_ratios,minimum_globals = generate_bin_values(dataframe, globals_dataframe, population_data)
    wrm.write_bin_table(aFile, bin_array, sample_counts, global_counts, control_counts, odds_ratios, lower_cis, upper_cis, likelihood_ratios, p_samples, p_globals, p_controls, p_likelihood_ratios,minimum_globals)
    print 'minimum globals=',minimum_globals
    
    ################################
    # Compute and write statistics #
    ################################
    # - binomial test
    # - wilcoxon
    wrm.write_label(aFile, "Statistics")
    #p_binomial=binom_test(samples_gt_globals,n_targets_gt_0,0.5)
    #stats_header_labels = ["Number of successes", "Number of targets", "pBinomial"]
    #stats_header_values = [samples_gt_globals, n_targets_gt_0, p_binomial]
    #wrm.write_information(aFile, stats_header_labels, stats_header_values, "   ")

    #t_wilcox,p_wilcox=wilcoxon(all_sample_means,all_global_means)
    #f2.write( 'Wilcoxon stat for sample vs globals means whole period:'+str(float(t_wilcox))+ '   p='+str(float(p_wilcox))+'\n')


    ##################################
    #Only include data with above minimum count of globals in statistical tests#
    ###################################
    trimmed_bin_array=[]
    trimmed_p_samples=[]
    trimmed_p_globals=[]
    trimmed_sample_counts=[]
    trimmed_global_counts=[]
    trimmed_likelihood_ratios=[]
    for i in range(0, len(bin_array)):
        if global_counts[i]>=minimum_globals:
            trimmed_bin_array.append(bin_array[i])
            trimmed_p_samples.append(p_samples[i])
            trimmed_p_globals.append(p_globals[i])
            trimmed_sample_counts.append(sample_counts[i])
            trimmed_global_counts.append(global_counts[i])
            trimmed_likelihood_ratios.append(likelihood_ratios[i])
    if len(trimmed_global_counts)<len(global_counts)/2:
        aFile.write('insufficient data for analysis')
        return('Insufficient data for analysis')
    #######################
    # Test distributions match qualitative predictions from model
    #######################
    tests_passed=0
    #######################
    # Test distributions are different
    #######################
    t_wilcox,p_wilcox=wilcoxon(trimmed_p_samples,trimmed_p_globals)
    aFile.write( 'Wilcoxon stat for samples vs globals with full controls:'+str(float(t_wilcox))+ '   p='+str(float(p_wilcox))+'\n')
    ks_d,ks_p=ks_2samp(trimmed_p_samples,trimmed_p_globals)
    aFile.write( 'KS test  for samples vs globals with full controls:'+str(float(ks_d))+ '   p='+str(float(ks_p))+'\n')
    if ks_p<min_p:
        aFile.write('The two distribitions are significantly different p<0.001'+'\n')
        tests_passed=tests_passed+1
        
  #  *******************************
  #  fit to curve with innovation - only output is on screen
  #  *******************************
    opt_results=fit_model_with_innovation(sample_counts,global_counts,bin_array)
        
    ##################################
    # Detect and display threshold and below curve #
    ##################################
    
# =============================================================================
#     threshold,successes,trials,p_threshold, below_curve, p_below_curve=detect_threshold(trimmed_bin_array, trimmed_sample_counts,trimmed_global_counts)
#     wrm.write_label(aFile,"Threshold analysis"+'\n')
#     aFile.write('Threshold: '+str(threshold)+'\n')
#     aFile.write('Successes: '+str(successes)+'\n')
#     aFile.write('Trials: '+str(trials)+'\n')
#     aFile.write('p: '+str(p_threshold)+'\n')
#     aFile.write('below_curve: '+str(below_curve)+'\n')
# =============================================================================
    threshold,p_threshold=detect_threshold(trimmed_bin_array,trimmed_sample_counts,trimmed_global_counts,population_data.bin_size)
    aFile.write('Threshold: '+str(threshold)+'\n')
    aFile.write('p: '+str(p_threshold)+'\n')
# =============================================================================
#     aFile.write('p_below_curve: '+str(p_below_curve)+'\n')
#     if p_threshold<min_p:
#         aFile.write('There is a significant threshold effect'+'\n')
#         tests_passed=tests_passed+1
#     if p_below_curve<min_p:
#         aFile.write('Samples_curve is significantly below globals curve '+'\n')
#   #      tests_passed=tests_passed+1
# =============================================================================
    ##################################
    #  comnpute medians
    ##################################
    median_samples=dataframe[dataframe.type=='s']['density'].median()
    median_globals=globals_dataframe ['density'].median()
    if median_samples>median_globals:
        aFile.write('Median density for samples>median density for globals'+'\n')
        tests_passed=tests_passed +1
    #########################
    # Write medians
    #########################
    wrm.write_label(aFile,"Medians"+'\n')
    aFile.write('Median density for sites: '+str(median_samples)+'\n')
    aFile.write('Median density for globals: '+str(median_globals)+'\n')
    
    #################
    # Fit data to logit curve #
    #################
    logit_results=fit_to_logit(trimmed_bin_array, trimmed_sample_counts, trimmed_global_counts,population_data.bin_size)
    logit_predictions=generate_logit_predictions(trimmed_bin_array,logit_results.params)
    if logit_results.params[0]>0:
        aFile.write('Logit curve positive')
        tests_passed=tests_passed+1
        
    
    #################
    # Fit data to linear curve #
    #################
# =============================================================================
#     linear_results=fit_to_linear(trimmed_bin_array, trimmed_sample_counts, trimmed_global_counts)
#     linear_predictions=generate_linear_predictions(trimmed_bin_array,linear_results.params)
# =============================================================================
   
    #################
    # Fit data to threshold model #
    #################
# =============================================================================
#     pOpt,pCov=fit_to_threshold_model(trimmed_bin_array, trimmed_sample_counts, trimmed_global_counts)
#     threshold_predictions=threshold_model(trimmed_bin_array, pOpt[0],pOpt[1])
#  #   print "pOpt0=",str(pOpt[0])+'\nl'
# #    print "pOpt1=",str(pOpt[1])+'\nl'
# #    print "pOpt2=",str(pOpt[2])+'\nl'
#     wrm.write_label(aFile,"Threshold model"+'\n')
#     aFile.write('Estimated threshold from curve fitting'+str(pOpt[0])+'\n')
#     aFile.write('Estimated coefficient for population'+str(pOpt[1])+'\n')
# #    aFile.write('Estimated zeta'+str(pOpt[2])+'\n')
# =============================================================================
    
# =============================================================================
#      #################
#     # comparison LMS for logit and threshold curve
#     #*******************
#     wrm.write_label(aFile,"Comparison of LMS"+'\n')
#     aFile.write('Likelihood ratio for logit fit')
#     aFile.write('LMS Logit:' + str(lms(logit_predictions,likelihood_ratios))+'\n')
#     aFile.write('LMS Threshold:'+ str(lms(threshold_predictions, likelihood_ratios))+'\n')
# =============================================================================
    
    # Generate graphs and statistics #
    ##################################
    # - plots likelihood ratios and sites graphs
    # stm.generate_stats_for_ratios(bin_array, likelihood_ratios, sample_counts, global_counts, p_likelihood_ratios, p_samples, p_globals, f2, "p Likelihood Ratio", new_path, directory)
     # would like to get rid of this routine
    #stm.generate_stats_for_ratios(bin_array, likelihood_ratios, sample_counts, global_counts, odds_ratios, p_samples, p_globals, f2, "Odds Ratio", new_path, directory, (population_data.max_population-population_data.bin_size*2), lower_cis=lower_cis, upper_cis=upper_cis)

    # - plots p_graphs and write statistics (binomial and wilcoxon)
    #threshold_binomial, threshold, threshold_success_count, threshold_trial_count, threshold_samples, threshold_controls = stm.generate_p_threshold_and_binomial(p_samples, p_controls, bin_array)
   # logit_results=stm.fitToLogit(bin_array, sample_counts, global_counts)
    plm.plot_p_graphs(trimmed_bin_array, trimmed_p_samples, trimmed_p_globals,population_data.bin_size, anIdentifier, aPath)
    plm.plot_cumulative_p_graphs(trimmed_bin_array, trimmed_p_samples, trimmed_p_globals, population_data.bin_size,median_samples,median_globals,threshold,anIdentifier, aPath)
    plm.plot_detection_frequencies (trimmed_bin_array, trimmed_likelihood_ratios, logit_predictions, population_data.bin_size, threshold,population_data.max_population-population_data.bin_size*2, anIdentifier, "detection_frequencies", aPath)
    wrm.write_label(aFile,"Logistic fit"+'\n')
    aFile.write("Intercept: "+str(logit_results.params[1])+'\n')
    aFile.write("Coefficient: "+str(logit_results.params[0])+'\n')
    aFile.write("AIC: "+str(logit_results.aic)+'\n')
    aFile.write("Pearson Chi2: "+str(logit_results.pearson_chi2)+'\n')
# =============================================================================
#     wrm.write_label(aFile,"Linear fit"+'\n')
#    # aFile.write("Intercept: "+str(linear_results.params[1])+'\n')
#     aFile.write("Coefficient: "+str(linear_results.params[0])+'\n')
#     aFile.write("AIC: "+str(linear_results.aic)+'\n')
# # =============================================================================
# =============================================================================
#         t_threshold_wilcoxon, p_threshold_wilcoxon = wilcoxon(threshold_controls, threshold_samples)
#         wrm.write_label(f2, "Statistics for threshold bins")
#         f2.write("Threshold: " + str(threshold))
#         stats_header_labels = ["Number of successes", "Number of targets", "pBinomial"]
#         stats_header_values = [threshold_success_count, threshold_trial_count, threshold_binomial]
#         wrm.write_information(f2, stats_header_labels, stats_header_values, "   ")
#         f2.write( 'Wilcoxon stat for pControls vs pSamples:'+str(float(t_threshold_wilcoxon))+ '   p='+str(float(p_threshold_wilcoxon))+'\n')
# 
# =============================================================================
    # - plots targets and globals on a map
    plm.plot_targets_on_map(dataframe, globals_dataframe, aPath, anIdentifier)
    wrm.write_label(aFile,"Test summary"+'\n')
    aFile.write("Tests passed:"+str(tests_passed)+'\n')
    
    


    

   # if(is_confounder_analysis):
   #     print("returning CA")
   #     return self.dataframe, self.globals_dataframe, bin_array, sample_counts, global_counts, control_counts, odds_ratios, top_MHs, bottom_MHs, top_test_MHs, bottom_test_MHs, likelihood_ratios, p_samples, p_globals, p_likelihood_ratios, tests_passed
   # else:
   #     return "Generated results", tests_passed
   
    
    
