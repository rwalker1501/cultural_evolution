
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
import os
from math import *
from matplotlib.pyplot import cm 
from matplotlib.font_manager import FontProperties
from mpl_toolkits.basemap import Basemap
from matplotlib.patches import Polygon
# import seaborn


def plot_stat_graphs(stat_dictionary, bin_values_df, population_data, directory, new_path):

    bin_array = bin_values_df['bin_array'];
    p_samples = bin_values_df['p_samples'];
    p_globals = bin_values_df['p_globals'];
    cum_p_samples = bin_values_df['cum_p_samples'];
    cum_p_globals = bin_values_df['cum_p_globals'];
    likelihood_ratios = bin_values_df['likelihood_ratios'];
    median_samples = stat_dictionary['median_samples'];
    median_globals = stat_dictionary['median_globals'];

    # plot graphs
    plot_p_graphs(bin_array, p_samples, p_globals, population_data.bin_size, directory, new_path)
    plot_cumulative_p_graphs(bin_array, cum_p_samples,cum_p_globals, population_data.bin_size, median_samples, median_globals, directory, new_path)
    plot_detection_frequencies(bin_array, likelihood_ratios, population_data.bin_size, population_data.max_population-population_data.bin_size*2, directory, new_path)



def plot_detection_frequencies(bins, actual_ratios, bin_size, max_xaxis, identifier, file_path):
    fig = plt.figure();
    add=bin_size/2
    bins=[x+add for x in bins]
    label = "detection_frequencies";
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(bins,actual_ratios,'bo',label='Observations')
    plt.ylabel(label)
    plt.xlabel("Population density per cell")
    plt.legend(title= "Legend")
    y_lim=max(actual_ratios)
    plt.gca().set_ylim(0,y_lim)
    fig.savefig(os.path.join(file_path, str(identifier) + "_" + label + ".png"))
    plt.close()

def plot_p_graphs(bins, p_samples, p_globals, bin_size,identifier, file_path):
    add=bin_size/2
    bins=[x+add for x in bins]
    fig = plt.figure()
    ax = fig.add_subplot(111)
    # ax.set_title(title)
    ax.plot(bins,p_samples,'b', label="Samples")
    ax.plot(bins, p_globals,'r',label="Globals")
    plt.gca().set_ylim(0, 0.40)
    
    plt.ylabel("Relative detection frequency")
    plt.xlabel("Population density per cell")
    plt.legend(title= "Legend")
    fig.savefig(os.path.join(file_path, str(identifier) + "_p_graph.png"))
    plt.close()
    
def plot_cumulative_p_graphs(bins, cum_p_samples, cum_p_globals, bin_size, median_samples,median_globals,identifier, file_path):
    add=bin_size/2
    bins=[x+add for x in bins]
    fig = plt.figure()
    ax = fig.add_subplot(111)

    ax.plot(bins,cum_p_samples,'b', label="Samples")
    ax.plot(bins, cum_p_globals,'r',label="Globals")
    plt.gca().set_ylim(0, 1.0);
    ax.axvline(median_samples, color='b', linestyle='--',label="Median samples")
    ax.axvline(median_globals, color='r', linestyle='--', label="Median globals")
    plt.ylabel("Cumulated relative detection frequency")
    plt.xlabel("Population density per cell")
    plt.legend(title= "Legend")
    fig.savefig(os.path.join(file_path, str(identifier) + "_cum_p_graph.png"))
    plt.close()


def plot_targets_on_map(dataframe, controls_dataframe, dir_path, directory):
    my_map = Basemap(projection='robin', lat_0=57, lon_0=-135,
        resolution = 'l', area_thresh = 1000.0,
        llcrnrlon=-136.25, llcrnrlat=56,
        urcrnrlon=-134.25, urcrnrlat=57.75)
     
    my_map.drawcoastlines()
    my_map.drawcountries()
    my_map.fillcontinents(color='lightgray')
    my_map.drawmapboundary()
     
    my_map.drawmeridians(np.arange(0, 360, 30))
    my_map.drawparallels(np.arange(-90, 90, 30))


    grouped_dataframe = dataframe[dataframe.type=='s']
    dir_exact_lats = grouped_dataframe[(grouped_dataframe.is_dir == True) & (grouped_dataframe.is_exact == True)]['latitude'].values
    dir_exact_lons = grouped_dataframe[(grouped_dataframe.is_dir == True) & (grouped_dataframe.is_exact == True)]['longitude'].values
    not_dir_exact_lats = grouped_dataframe[~((grouped_dataframe.is_dir == True) & (grouped_dataframe.is_exact == True))]['latitude'].values
    not_dir_exact_lons = grouped_dataframe[~((grouped_dataframe.is_dir == True) & (grouped_dataframe.is_exact == True))]['longitude'].values



    controls_lats = controls_dataframe['latitude'].values
    controls_lons = controls_dataframe['longitude'].values


    t_x,t_y = my_map(dir_exact_lons, dir_exact_lats)
    f_x,f_y = my_map(not_dir_exact_lons, not_dir_exact_lats)
    c_x, c_y = my_map(controls_lons, controls_lats)
    my_map.plot(c_x, c_y, 'ro', markersize=3, label="globals")
    my_map.plot(t_x, t_y, 'yo', markersize=3, label="dir and exact")
    my_map.plot(f_x, f_y, 'bo', markersize=3, label="not dir and exact")
    
    fontP = FontProperties()
    fontP.set_size('small')
    plt.legend(title="Legend", loc='lower center', prop=fontP, ncol=3)
    plt.savefig(os.path.join(dir_path, directory + "_target_on_map.png"))
    plt.close()

def plot_densities_on_map_by_time_point(population_data, time):
    plt.figure(figsize=(14,8))
    my_map = Basemap(projection='robin', lat_0=57, lon_0=-135,
    resolution = 'l', area_thresh = 1000.0,
    llcrnrlon=-136.25, llcrnrlat=56,
    urcrnrlon=-134.25, urcrnrlat=57.75)
     
    my_map.drawcoastlines()
    my_map.drawcountries()
    my_map.fillcontinents(color='lightgray', zorder=0)
    my_map.drawmapboundary()
     
    my_map.drawmeridians(np.arange(0, 360, 30))
    my_map.drawparallels(np.arange(-90, 90, 30))

    time_index = -1
    for i in range(0, len(population_data.time_array)):
        if time == population_data.time_array[i]*population_data.time_multiplier:
            time_index = i
            break

    lats = []
    lons = []
    dens = []
    for i in range(0, len(population_data.density_array[time_index])):
        density = population_data.density_array[time_index][i]
        if density != 0:
            lats.append(population_data.lat_array[i])
            lons.append(population_data.lon_array[i])
            dens.append(density)

    cmap = cm.get_cmap('jet')

    dens = np.array(dens).astype(float)/3000;
    x,y = my_map(lons, lats)
    rgba = cmap(dens);



    ax = my_map.scatter(x, y, marker='o', c=rgba, s=3, zorder=1)


    sm = cm.ScalarMappable(cmap=cmap)
    sm.set_array([])

    sm.set_clim([0, 3000])
    plt.colorbar(sm, orientation='horizontal', pad=0.03, aspect=50)


    plt.savefig("densitites_on_map_" + population_data.name + "_" + str(time) + "Kya.png");
    
    map_path = get_map_file_path(population_data.name.lower() + "_densitites_on_map_" + str(time) + ".png");
    print("Map filepath: " + map_path);
    plt.savefig(map_path, dpi=500);
    plt.close();


def plot_densities_on_map_by_range(population_data, min_density, max_density, start_time, end_time):
    my_map = Basemap(projection='robin', lat_0=57, lon_0=-135,
    resolution = 'l', area_thresh = 1000.0,
    llcrnrlon=-136.25, llcrnrlat=56,
    urcrnrlon=-134.25, urcrnrlat=57.75)
     
    my_map.drawcoastlines()
    my_map.drawcountries()
    my_map.fillcontinents(color='lightgray')
    my_map.drawmapboundary()
     
    my_map.drawmeridians(np.arange(0, 360, 30))
    my_map.drawparallels(np.arange(-90, 90, 30))

    start_time_index = -1
    end_time_index = -1
    for i in range(0, len(population_data.time_array)):
        if start_time == population_data.time_array[i]*population_data.time_multiplier:
            start_time_index = i
        if end_time == population_data.time_array[i]*population_data.time_multiplier:
            end_time_index = i
        if start_time_index != -1 and end_time_index != -1:
            break

    if start_time_index > end_time_index:
        temp = end_time_index
        end_time_index = start_time_index
        start_time_index = temp

    lats = []
    lons = []
    dens = []
    for i in range(start_time_index, end_time_index + 1):
        dens_array = population_data.density_array[i]
        for j in range(0, len(dens_array)):
            density = dens_array[j]
            if density > min_density and density < max_density:
                # print(density)
                lats.append(population_data.lat_array[j])
                lons.append(population_data.lon_array[j])
                dens.append(density)

    print(start_time_index)
    print(end_time_index)
    x,y = my_map(lons, lats)
    my_map.plot(x, y, 'yo', markersize=1)

    map_path = get_map_file_path(population_data.name.lower() + "_densities_on_den_range_" + str(min_density) + "-" + str(max_density) + "_and_time_range_" + str(start_time) + "-" + str(end_time) + ".png");
    print("Map filepath: " + map_path);
    plt.savefig(map_path, dpi=500)
    plt.close()

def plot_densities_on_map_by_time_range(population_data, start_time, end_time):
    my_map = Basemap(projection='robin', lat_0=57, lon_0=-135,
    resolution = 'l', area_thresh = 1000.0,
    llcrnrlon=-136.25, llcrnrlat=56,
    urcrnrlon=-134.25, urcrnrlat=57.75)
     
    my_map.drawcoastlines()
    my_map.drawcountries()
    my_map.fillcontinents(color='lightgray')
    my_map.drawmapboundary()
     
    my_map.drawmeridians(np.arange(0, 360, 30))
    my_map.drawparallels(np.arange(-90, 90, 30))

    start_time_index = -1
    end_time_index = -1
    for i in range(0, len(population_data.time_array)):
        if start_time == population_data.time_array[i]*population_data.time_multiplier:
            start_time_index = i
        if end_time == population_data.time_array[i]*population_data.time_multiplier:
            end_time_index = i
        if start_time_index != -1 and end_time_index != -1:
            break

    if start_time_index > end_time_index:
        temp = end_time_index
        end_time_index = start_time_index
        start_time_index = temp

    mid_lats = []
    mid_lons = []
    high_lats = []
    high_lons = []
    for i in range(start_time_index, end_time_index + 1):
        dens_array = population_data.density_array[i]
        for j in range(0, len(dens_array)):
            density = dens_array[j]
            if density > 3000 and density < 5400:
                # print(density)
                mid_lats.append(population_data.lat_array[j])
                mid_lons.append(population_data.lon_array[j])
            if density > 5400 and density < 6000:
                # print(density)
                high_lats.append(population_data.lat_array[j])
                high_lons.append(population_data.lon_array[j])

    print(start_time_index)
    print(end_time_index)
    x,y = my_map(mid_lons, mid_lats)
    my_map.plot(x, y, 'yo', markersize=1)
    x,y = my_map(high_lons, high_lats)
    my_map.plot(x, y, 'ro', markersize=1)

    map_path = get_map_file_path(population_data.name.lower() + "_densities_time_range_" + str(start_time) + "-" + str(end_time) + ".png");
    print("Map filepath: " + map_path);
    plt.savefig(map_path, dpi=500)
    plt.close()
    
def plot_maximum_likelihood(acc,rho_bins,rho_bins2,y_acc, lambda_v, opt_threshold, samples_counts2, controls_counts2, scale, directory,file_path):
    max_acc=acc.max(axis=0) *1e-10   #largest accumulated likelihood for a given rho multiplied by a small constant - gives roughly constant result
    acc_plus=(acc+max_acc)
    yy=np.cumsum(acc_plus,axis=0) #This will give me cumulated likelihood for each column as in mathlab but is going to lead to shape problems. May need transpose
    np.seterr(divide='ignore',invalid='ignore') #probably not necessary
    yy=np.divide(yy,yy[len(yy)-1:]) # Should give me cumulated likelihood up to 1 - seems to work
    pred_int=np.zeros((len(rho_bins),6)) #Not sure about size of this - in the original looks like an empty matrix - NOW LESS SURE
    for k in range (0,len(rho_bins)):
# =============================================================================
#         first_sub=np.array((0)) # 0;
#         second_sub=yy[0:(len(yy)-1),k] #yy(1:end-1,k)] I am suspicious of this. Sometimes suddenly jumps to 1,
# =============================================================================
        data_x=np.hstack((np.array((0)),yy[0:(len(yy)-1),k])) #[0; yy(1:end-1,k)]. This looks OK
        interpolated=np.interp([0.025, 0.25, 0.5, 0.75, 0.975], data_x,y_acc) #Not sure I have interpreted this correctly. 
# =============================================================================
#         first_sub2=y_acc[0:(len(y_acc)-1)]+y_acc[1:len(y_acc)] #*(yacc(1:end-1)+yacc(2:end)))
#         second_sub2=acc[0:len(acc)-1,k] #Acc(1:end-1,k))This is mostly zeros in current version
#         third_sub2=np.sum(acc[0:(len(acc)-1),k],axis=0) #sum(Acc(1:end-1,k))]; yields a single small float
# =============================================================================
        term2_1=np.dot(y_acc[0:(len(y_acc)-1)]+y_acc[1:len(y_acc)],acc[0:len(acc)-1,k]) #unsure about this. Is it a dot or a matrix multiplicaiton. It also works as a matrix multiplication but that is not what I am using now
        if(term2_1==0) and (np.sum(acc[0:(len(acc)-1),k],axis=0)==0):
            term2=0
        else:
            term2=((0.5*term2_1/np.sum(acc[0:(len(acc)-1),k],axis=0)))
        pred_int[k,:]=np.hstack((interpolated,term2)) #This is one dimensional - correct - but not sure why it is doing this. Adds a 6th element which is not in sequence with the others and which seems to be  never used.
# =============================================================================
#     lambda_v = lambda_v*sqrt(scale)
#     rho_bins = rho_bins*scale
#     rho_bins2=rho_bins2*scale
# =============================================================================
    patches=[]
    term1=rho_bins
    term2=np.flip(rho_bins,0)
    h1_x=np.hstack((term1,term2))
    h1_y=np.hstack((pred_int[:,0],np.flip(pred_int[:,4],0)))
    h1_array=np.vstack((h1_x,h1_y)).T
    h1 = Polygon(h1_array,linewidth=1,edgecolor='none',facecolor=(1,0.9,0.9),closed=True,antialiased=True)
    patches.append(h1)
    h2_x = np.hstack((rho_bins,np.flip(rho_bins,0))) 
    h2_y=np.hstack((pred_int[:,1],np.flip(pred_int[:,3],0)))
    h2_array=np.vstack((h2_x,h2_y)).T
    h2 = Polygon(h2_array,linewidth=1,edgecolor='none',facecolor=(1,0.7,0.7),closed=True,antialiased=True)
    patches.append(h2)
    fig1 = plt.figure();
    ax = fig1.add_subplot(111)
    plt.gca().set_xlim(0, 30);
 #   p = PatchCollection(patches, alpha=0.4)
#    ax.add_collection(p)
    ax.add_patch(h1)
    ax.add_patch(h2)
    ax.axvline(opt_threshold, color='g', linestyle='--',label="Threshold")
    ax.plot(rho_bins,pred_int[:,2],linewidth=1, color='blue',antialiased=True)
    ax.plot(rho_bins2,samples_counts2/(samples_counts2+controls_counts2),color='black', marker='o', markersize=3,linestyle='None',antialiased=True)
    fig_path=os.path.join(file_path, str(directory)) + "/"+directory+"_fit to model.png"
    fig1.savefig(fig_path)
    plt.close()
    
def plot_parameter_values(lnL,lambda_v, zetta_v, eps_v,directory,file_path):
    #      Figure 2
    fig2 = plt.figure();
    lnlminusmax=lnL-np.amax(lnL)
    exp_lnlminusmax=np.exp(lnlminusmax)
    dim1=np.mean(exp_lnlminusmax,axis=2)  #up to here - look at definition of dimension
    p_lambda = np.squeeze(np.mean(dim1,axis=1))
    p_lambda=np.true_divide(p_lambda,np.trapz(p_lambda,lambda_v))
    ax2=fig2.add_subplot(111)
    ax2.plot(lambda_v,p_lambda)
    plt.xlabel('lambda')
    fig_path=os.path.join(file_path, str(directory)) + "/"+directory+"_lambda.png"
    fig2.savefig(fig_path)
    plt.show
     #     Figure 3 - in the end we will move this into plot library
    fig3 = plt.figure();
    dim1=np.mean(exp_lnlminusmax,axis=2) 
    p_eps = np.squeeze(np.mean(dim1,axis=0))
    p_eps=np.true_divide(p_eps,np.trapz(p_eps,eps_v))
    ax3=fig3.add_subplot(111)
    ax3.plot(eps_v,p_eps)
    plt.xlabel('epsilon')
    fig_path=os.path.join(file_path, str(directory)) + "/"+directory+"_eps.png"
    fig3.savefig(fig_path)
    plt.show
     #     Figure 4 - in the end we will move this into plot library
    fig4 = plt.figure();
    dim1=np.mean(exp_lnlminusmax,axis=0) 
    p_zetta = np.squeeze(np.mean(dim1,axis=0))
    trapz=np.trapz(p_zetta,np.log10(zetta_v))
    p_zetta=np.true_divide(p_zetta,trapz)
    x_data=np.log10(zetta_v)
    trapz=np.trapz(p_zetta,np.log10(zetta_v))
    y_data=np.true_divide(p_zetta,trapz)
    ax4=fig4.add_subplot(111)
    ax4.set_xlabel("log10 zetta")
    ax4.plot(x_data, y_data)
    fig_path=os.path.join(file_path, str(directory)) + "/"+directory+"_zetta.png"
    fig4.savefig(fig_path)
    plt.show

def get_map_file_path(filename):
    base_path = os.getcwd();
    map_path = os.path.join(base_path, "maps")
    if not os.path.exists(map_path):
        os.makedirs(map_path);
    file_path = os.path.join(map_path, filename);
    return file_path

if __name__ == '__main__':
        main()


