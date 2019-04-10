from __future__ import division #This means division always gives a floating result
import numpy as np
import pandas as pd;
from scipy.stats import linregress, ks_2samp,poisson;
from scipy.sparse import spdiags
from math import *
import sys
import cPickle as pkl
import json
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection





def compute_likelihood_model(working_data_filename,samples_counts, controls_counts, globals_counts,samples_counts2, controls_counts2,rho_bins_4_python,bin_boundaries2_4_python,bin_width,rho_bins,work_data_available,low_res=False):
 # fix parameter values for scan
 # when the low_res parameter is set to true, the system produces low_res graphs. Used for system testing and exploratory testing
    if low_res:
        lambda_v=np.arange(25,50,1)
        eps_v=np.linspace(0,0.1,num=11,endpoint=False)
        zetta_v=np.exp(np.linspace(log(1e-5),log(1e-3),num=11,endpoint=False))  
    else:
        lambda_v=np.arange(25,50,0.1)
        eps_v=np.linspace(0,0.1,num=101,endpoint=False)
        zetta_v=np.exp(np.linspace(log(1e-5),log(1e-3),num=101,endpoint=False))      
    print 'length binboundaries2_4_python=',len(bin_boundaries2_4_python)
    rho_bins2=bin_boundaries2_4_python+bin_width/2
    print 'length rho_bins2=', len(rho_bins2)
    if not work_data_available:
        n_controls=np.sum(controls_counts)  
        n_samples=np.sum(samples_counts) 
        n_globals=n_controls+n_samples # not sure these are needed
        rho_bins2=bin_boundaries2_4_python[0:len(bin_boundaries2_4_python)-1]+bin_width/2 #This is EC2 not sure this is correct - depends on usage later on, needs to be adjusted to take account of nature of bins
        l_shift=n_samples*(log(float(n_samples)/float(n_controls))-1)
        print 'l_shift=',l_shift
  #      lambda_v=np.arange(25,50,0.1)
        
  #      lambda_v=np.arange(2,4,0.25) #This is completely different from Eriksson code
  #      eps_v=np.linspace(0,0.1,num=10,endpoint=False)  #shortened for debugging purposes
        
        n_lambda=len(lambda_v)
        n_eps=len(eps_v)
        n_zetta=len(zetta_v)
 #       y_acc=np.linspace(0,10e-2,num=1001) #COLUMN VECTOR
        y_acc=np.linspace(0,1e-4,num=401) #COLUMN VECTOR. In Tindbergen program he had different values that generated a lot of zeros. These in turn created problems when we had to divide by last element in yy
        acc=np.zeros((len(y_acc),len(rho_bins)))
        anElement=np.zeros(len(rho_bins))
        print 'an element completed'
        aRow=[]
        lnL=np.zeros((n_lambda,n_eps,n_zetta))
        print 'lnl completed'
   #     zetta_opt=0*lnL - this is never used - have comented it.
        sqrt_rho_bins=np.sqrt(rho_bins) #These are the values we are computing - rhobins_4_python are intervals for histogram only. In original program were inside loop. Have moved it outside
   #     for i_lambda in range (0,n_lambda-1):
        for i_lambda in range (0,n_lambda):
            my_lambda=float(lambda_v[i_lambda]) 
            bin_zeros=np.zeros(len(sqrt_rho_bins)) #in python I can't compare a vector with a scalar
            pI=np.maximum(bin_zeros,1-my_lambda/sqrt_rho_bins)  #COLUMN VECTOR
    #        for i_zetta in range(0,n_zetta-1):
            for i_zetta in range(0,n_zetta):
                print "i_lambda,i_zetta", i_lambda,",",i_zetta
    #            for i_eps in range (0, n_eps-1):
                for i_eps in range (0, n_eps):
                     pObs=np.zeros(len(pI)).astype(float)                
                     pObs=zetta_v[i_zetta]*((1-eps_v[i_eps])*pI+eps_v[i_eps]) #COLUMN VECTOR (scalars * a column vector)
                     pObs_small=np.zeros(len(pObs))
                     pObs_small.fill(1e-20)
                     pObs=np.maximum(pObs,pObs_small)
                     pObs=pObs.astype(float)
                     log_samples=np.dot(samples_counts,np.log(pObs)) 
                     log_controls=np.dot(controls_counts,np.log(1-pObs)) 
                     LL=log_samples+log_controls
                     L=np.exp(LL-l_shift)
                     lnL[i_lambda,i_eps,i_zetta]=LL # This is different from original datastructure. Will require change of later code. I could also assign using an array op.
                     zz=pObs  #should be able to get rid of this
                     rhs=np.floor(1+zz/y_acc[1]).astype(int) #This is original code - yields a 1-based index
               #       print 'rhs=',rhs
                     len_y_acc=np.array(len(y_acc))
                     len_y_acc.fill(len(y_acc))
                     i_acc=np.minimum(len_y_acc,rhs) #vThis yields column vector of indexes corresponding to different values of pObs. ector length =401. Maximum value of index =400 (zero based vector).I am keeping it 1-based
                     for i in range(0,len(i_acc)):
                         x_coord=i_acc[i]-1
                         y_coord=i
                         acc[x_coord,y_coord]=acc[x_coord,y_coord]+L
# =============================================================================
#                      a_range=np.arange(0,len(rho_bins)) #I have taken off the minus 1 in the original - because ranges in matlab are inclusive. Here they are not. But I am not confident
#                      second_term=a_range* len(acc[:,0])
#                      i_s=((i_acc-1).T+second_term).T +1  #This is original code - but doesn't take account of nature of matlab indices - 
#                      
#                      for i in range (0,len(i_acc)): #I am suspcious about this - not sure x and y are right way round. is[0]=-1. I don't like this either
#                          x_coord=(i_s[i]-1)/29 #all this is correct if i_s starts at zero
#                          x_coord=x_coord.astype(int) #This is temporary - will need to generate it dynamicallyx_coord=i_s[i] % 29  #This is temporary - will need to generate it dynamically. Not sure on division by 28
#                          y_coord=(i_s[i]-1)%29
#  #                        print 'i_s,x,y=(',i_s[i],',',x_coord,',',y_coord,')'
#                          acc[x_coord,y_coord]=acc[x_coord,y_coord]+L
#                      print '.'
# =============================================================================
# loops seem to conclude correctly - but not certain of value of acc  
    max_acc=acc.max(axis=0) *1e-10   #largest accumulated likelihood for a given rho multiplied by a small constant - gives roughly constant result
    acc_plus=(acc+max_acc)
    yy=np.cumsum(acc_plus,axis=0) #This will give me cumulated likelihood for each column as in mathlab but is going to lead to shape problems. May need transpose
    yy_pre_div=yy
    np.seterr(divide='ignore',invalid='ignore')
    yy=np.divide(yy,yy[len(yy)-1:]) # Should give me cumulated likelihood up to 1 - seems to work
    pred_int=np.zeros((len(rho_bins),6)) #Not sure about size of this - in the original looks like an empty matrix - NOW LESS SURE
    for k in range (0,len(rho_bins)):
        first_sub=np.array((0)) # 0;
        second_sub=yy[0:(len(yy)-1),k] #yy(1:end-1,k)] I am suspicious of this. Sometimes suddenly jumps to 1,
        data_x=np.hstack((first_sub,second_sub)) #[0; yy(1:end-1,k)]. This looks OK
        interpolated=np.interp([0.025, 0.25, 0.5, 0.75, 0.975], data_x,y_acc) #Not sure I have interpreted this correctly. 
        #interpolated=np.interp([0.025, 0.25, 0.5, 0.75, 0.975], y_acc,data_y)
        # python interpolation function says data_x must be increasing - but here it is non-monotonic
        first_sub2=y_acc[0:(len(y_acc)-1)]+y_acc[1:len(y_acc)] #*(yacc(1:end-1)+yacc(2:end)))
        second_sub2=acc[0:len(acc)-1,k] #Acc(1:end-1,k))This is mostly zeros in current version
        third_sub2=np.sum(acc[0:(len(acc)-1),k],axis=0) #sum(Acc(1:end-1,k))]; yields a single small float
        term2_1=np.dot(first_sub2,second_sub2) #unsure about this. Is it a dot or a matrix multiplicaiton. It also works as a matrix multiplication but that is not what I am using now
        if(term2_1==0) and (third_sub2==0):
            term2=0
        else:
            term2=((0.5*term2_1/third_sub2))
            #term2=((0.5*test/third_sub2))
        test=np.hstack((interpolated,term2)) #This is one dimensional - correct - but not sure why it is doing this. Adds a 6th element which is not in sequence with the others and which seems to be  never used.
        
        pred_int[k,:]=test 
    scale = (2/sqrt(3))/100 #  convert from hexagon pop size to Timmermann units, inds/100 km^2
    lambda_v = lambda_v*sqrt(scale)
    rho_bins = rho_bins*scale
    rho_bins2=rho_bins2*scale
  #     Figure 1 - in the end we will move this into plot library
    patches=[]
    term1=rho_bins
    term2=np.flip(rho_bins,0)
    h1_x=np.hstack((term1,term2))
    h1_y=np.hstack((pred_int[:,0],np.flip(pred_int[:,4],0)))
    h1_array=np.vstack((h1_x,h1_y)).T
    h1 = Polygon(h1_array,linewidth=1,edgecolor='r',facecolor='none',closed=True,antialiased=True)
    patches.append(h1)
    h2_x = np.hstack((rho_bins,np.flip(rho_bins,0))) 
    h2_y=np.hstack((pred_int[:,1],np.flip(pred_int[:,3],0)))
    h2_array=np.vstack((h2_x,h2_y)).T
    h2 = Polygon(h2_array,linewidth=1,edgecolor='g',facecolor='none',closed=True,antialiased=True)
    patches.append(h2)
    fig1 = plt.figure();
    ax = fig1.add_subplot(111)
 #   p = PatchCollection(patches, alpha=0.4)
#    ax.add_collection(p)
    ax.add_patch(h1)
    ax.add_patch(h2)
    ax.plot(rho_bins,pred_int[:,2],linewidth=2, color='blue',antialiased=True)
    ax.plot(rho_bins2,samples_counts2/(samples_counts2+controls_counts2),color='black', marker='o', markersize=6,linestyle='None',antialiased=True)
    
    plt.show
  #     Figure 2 - in the end we will move this into plot library
    fig2 = plt.figure();
    lnlminusmax=lnL-np.amax(lnL)
    exp_lnlminusmax=np.exp(lnlminusmax)
    dim1=np.mean(exp_lnlminusmax,axis=2)  #up to here - look at definition of dimension
    p_lambda = np.squeeze(np.mean(dim1,axis=1))
    p_lambda=np.true_divide(p_lambda,np.trapz(p_lambda,lambda_v))
    ax2=fig2.add_subplot(111)
    ax2.plot(lambda_v,p_lambda)
    plt.xlabel('lambda')
    plt.show
     #     Figure 3 - in the end we will move this into plot library
    fig3 = plt.figure();
    dim1=np.mean(exp_lnlminusmax,axis=2) 
    print 'shape dim1 fig3=',np.shape(dim1)
    p_eps = np.squeeze(np.mean(dim1,axis=0))
    print 'shape p_eps=',np.shape(p_eps)
    print "p_eps=", p_eps
    p_eps=np.true_divide(p_eps,np.trapz(p_eps,eps_v))
    ax3=fig3.add_subplot(111)
    ax3.plot(eps_v,p_eps)
    plt.xlabel('epsilon')
    plt.show
     #     Figure 4 - in the end we will move this into plot library
    fig4 = plt.figure();
    dim1=np.mean(exp_lnlminusmax,axis=0) 
    print 'shape dim1 fig4=',np.shape(dim1)
    p_zetta = np.squeeze(np.mean(dim1,axis=0))
    print 'shape p_zetta=', np.shape (p_zetta) #p_zetta is identical to p_eps - this can't possibly be right . but dim1 seems to be different
    print "p_zettta=", p_zetta
    trapz=np.trapz(p_zetta,np.log10(zetta_v))
    print "trapz=",trapz
    p_zetta=np.true_divide(p_zetta,trapz)
    x_data=np.log10(zetta_v)
    trapz=np.trapz(p_zetta,np.log10(zetta_v))
    y_data=np.true_divide(p_zetta,trapz)
    print 'x_data=', x_data
    ax4=fig4.add_subplot(111)
    ax4.set_xlabel("log10 zetta")
    print "y_data=",y_data
    ax4.plot(x_data, y_data)
    plt.show
    
    
   
        
        
        
        
        
        
        
    
    
    
    
    
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        
                    
                    
                
            
            
        
        
        
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
     
    
    
    
    

def process_dataframe(dataframe):

    conditions = [];
    samples_growth_coefficients = []
    valid_ids = []

    cluster_ids = dataframe.cluster_id.unique();


    for cluster_id in cluster_ids:
        cluster_df = dataframe[dataframe.cluster_id==cluster_id]
        sample_cluster_df = cluster_df[cluster_df.type == 's']
        if np.isnan(sample_cluster_df['density'].median()):
            print("Removing cluster " + str(cluster_id));
            dataframe = dataframe[dataframe.cluster_id != cluster_id];
            continue;

        # extract all periods and all population as arrays from the samples dataframe
        sample_times = sample_cluster_df['period'].values
        sample_populations = sample_cluster_df['density'].values

        # compute growth coefficients for samples
        growth_coefficient_samples=compute_growth_coefficient(sample_times, sample_populations)

        valid_ids.append(cluster_id);
        samples_growth_coefficients.append(growth_coefficient_samples)

    for cluster_id in valid_ids:
        conditions.append((dataframe['cluster_id'] == cluster_id));

    dataframe['samples_growth_coefficient'] = np.select(conditions, samples_growth_coefficients);
    return dataframe;

def compute_growth_coefficient(times, populations):
    if len(times)>=2:
        for i in range(0, len(populations)):
            if np.isnan(populations[i]):
                populations[i] = 0
        slope, intercept, r_value, p_value, std_err = linregress(times, populations)
        return slope 
    else:
        return -1.


def generate_bin_values_dataframe(dataframe, globals_dataframe, population_data, minimum_globals):
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
    likelihood_ratios = []
    p_samples = []
    p_globals = []

    ##############
    # Get Totals #
    ##############
    # total samples by summing contributions
    # total globals by counting rows
    total_samples = dataframe[dataframe.type=='s']['contribution'].sum()
    total_globals = globals_dataframe['density'].count()
    
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

    df = pd.DataFrame({'bin_array': bin_array, 'sample_counts': sample_counts, 'global_counts': global_counts, 'likelihood_ratios': likelihood_ratios, 'p_samples': p_samples, 'p_globals': p_globals})

    return df;

def generate_statistics(dataframe, globals_dataframe, bin_values_df, minimum_globals):

    trimmed_bin_values_df = bin_values_df[bin_values_df.global_counts > minimum_globals];
    trimmed_bin_values_df['cum_p_samples'] = trimmed_bin_values_df.p_samples.cumsum();
    trimmed_bin_values_df['cum_p_globals'] = trimmed_bin_values_df.p_globals.cumsum();

    
    if len(trimmed_bin_values_df.index) < len(bin_values_df.index)/2:
        return None, None;

    stat_dictionary = {};

    stat_dictionary['trimmed_bin_values_df'] = trimmed_bin_values_df;

    stat_dictionary['total_samples'] = dataframe[dataframe.type=='s']['density'].sum()
    stat_dictionary['total_globals'] = globals_dataframe ['density'].sum()

    stat_dictionary['median_samples'] = dataframe[dataframe.type=='s']['density'].median()
    stat_dictionary['median_globals'] = globals_dataframe ['density'].median()


    stat_dictionary['mean_samples'] = dataframe[dataframe.type=='s']['density'].mean()
    stat_dictionary['mean_globals'] = globals_dataframe ['density'].mean()

    stat_dictionary['std_samples'] = dataframe[dataframe.type=='s']['density'].std()
    stat_dictionary['std_globals'] = globals_dataframe ['density'].std();


    trimmed_p_samples = trimmed_bin_values_df['p_samples'].values;
    trimmed_p_globals = trimmed_bin_values_df['p_globals'].values;


    ks_d,ks_p= ks_2samp(trimmed_p_samples,trimmed_p_globals)

    stat_dictionary['ks_d'] = ks_d;
    stat_dictionary['ks_p'] = ks_p;

    return stat_dictionary, trimmed_bin_values_df;




def write_results(aFile,anIdentifier, aPath,dataframe, globals_dataframe,population_data, min_globals, min_p):
    
    bin_array, sample_counts, global_counts, likelihood_ratios, p_samples, p_globals = generate_bin_values(dataframe, globals_dataframe, population_data, min_globals);
    wrm.write_bin_table(aFile, bin_array, sample_counts, global_counts, likelihood_ratios, p_samples, p_globals, min_globals)
    
    ################################
    # Compute and write statistics #
    ################################
    # - binomial test
    # - wilcoxon
    wrm.write_label(aFile, "Statistics")


    #######################
    # Statistics
    #######################

    wrm.write_label(aFile, "Statistics");

    total_samples=dataframe[dataframe.type=='s']['density'].sum()
    total_globals=globals_dataframe ['density'].sum()
    aFile.write('Total sites: '+str(total_samples)+'\n')
    aFile.write('Total globals: '+str(total_globals)+'\n\n')

    median_samples=dataframe[dataframe.type=='s']['density'].median()
    median_globals=globals_dataframe ['density'].median()
    aFile.write('Median density for sites: '+str(median_samples)+'\n')
    aFile.write('Median density for globals: '+str(median_globals)+'\n\n')


    mean_samples=dataframe[dataframe.type=='s']['density'].mean()
    mean_globals=globals_dataframe ['density'].mean()
    aFile.write('Mean density for sites: '+str(mean_samples)+'\n')
    aFile.write('Mean density for globals: '+str(mean_globals)+'\n\n')

    std_samples=dataframe[dataframe.type=='s']['density'].std()
    std_globals=globals_dataframe ['density'].std()
    aFile.write('Mean density for sites: '+str(median_samples)+'\n')
    aFile.write('Mean density for globals: '+str(median_globals)+'\n')


    #######################
    # Test distributions are different
    #######################

    wrm.write_label(aFile, "K-S2 Test\n")

    ks_d,ks_p=ks_2samp(trimmed_p_samples,trimmed_p_globals)
    aFile.write( 'KS test  for samples vs globals with full controls:'+str(float(ks_d))+ '   p='+str(float(ks_p))+'\n')
    if ks_p<min_p:
        aFile.write('The two distribitions are significantly different p<0.001'+'\n')

    ks_d,ks_p=ks_2samp(trimmed_p_samples,trimmed_p_globals)
    f2.write( 'KS test for p_samples vs p_globals :'+str(float(ks_d))+ '   p='+str(float(ks_p))+'\n')
    if ks_p<self.min_p:
         f2.write('The two distribitions are significantly different p<0.001'+'\n')    
        
    # plot graphs
    plm.plot_p_graphs(trimmed_bin_array, trimmed_p_samples, trimmed_p_globals,population_data.bin_size, anIdentifier, aPath)
    plm.plot_cumulative_p_graphs(trimmed_bin_array, trimmed_p_samples, trimmed_p_globals, population_data.bin_size,median_samples,median_globals, anIdentifier, aPath)
    plm.plot_detection_frequencies (trimmed_bin_array, trimmed_likelihood_ratios, population_data.bin_size, population_data.max_population-population_data.bin_size*2, anIdentifier, "detection_frequencies", aPath)

    # - plots targets and globals on a map
    plm.plot_targets_on_map(dataframe, globals_dataframe, aPath, anIdentifier)
    

    
