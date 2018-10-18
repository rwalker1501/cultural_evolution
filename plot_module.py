
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
import os
from matplotlib.pyplot import cm 
from matplotlib.font_manager import FontProperties
from mpl_toolkits.basemap import Basemap
# import seaborn

    
def plot_crude_or_vs_mh_or(bin_array, crude_or, mh_or, identifier, file_path):
    fig = plt.figure()
    ax = fig.add_subplot(111)

    ax.plot(bin_array, crude_or, color='blue', label="Crude OR");
    ax.plot(bin_array, mh_or, color='coral', label="MH OR");
    plt.ylabel("OR")
    plt.xlabel("Population density")
    plt.gca().set_ylim(0, 7)

    plt.legend(title= "Legend")
    fig.savefig(os.path.join(file_path, str(identifier) + "_crude_vs_mh.png"))
    plt.close()

def plot_time_clusters(time, labels):
    n = np.array(labels).max() + 1
    colors = cm.rainbow(np.linspace(0,1,n))

    fig = plt.figure(figsize=(11,8))
    ax = fig.add_subplot(111)
    for i in range(0, len(time)):
        ax.scatter(time[i], 5, c=colors[labels[i]])

    plt.show()




def plot_histogram(values, titles, title, file_path):
    fig = plt.figure(figsize=(11,8))
    ax = fig.add_subplot(111)
    bins=np.histogram(np.hstack((values[0],values[1])), bins=50)[1]
    colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
    mean_patches = []
    std_patches = []

    for i in range(0, len(values)):
        plt.hist(values[i], bins, alpha=0.5, label=titles[i])
        mean = round(values[i].mean(),2)
        mean_patches.append(mpatches.Patch(label=mean, alpha=0.5, color=colors[i]))
        std = round(values[i].std(),2)
        std_patches.append(mpatches.Patch(label=std, alpha=0.5, color=colors[i]))


    main_legend = plt.legend(title= "Legend")
    mean_legend = plt.legend(loc=(1.01, 0.2), handles=mean_patches, title="Mean")
    std_legend = plt.legend(loc=(1.01, 0), handles=std_patches, title="STD")
    ax.add_artist(main_legend)
    ax.add_artist(mean_legend)
    ax.add_artist(std_legend)
    png_path = os.path.join(file_path , "histogram-" + title + ".png")
    # plt.tight_layout()
    fig.savefig(png_path)
    plt.close()

def plot_boxplot(values, titles, title, file_path):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    # ax.set_title(title)
    ax.boxplot(values)
    ticks = [i+1 for i in range(0, len(titles))]
    plt.xticks(ticks, titles)
    plt.ylabel("value")
    filename = title.lower().replace(" ", "_")
    fig.savefig(os.path.join(file_path, filename+".png"))

def plot_ratio(bins, actual_ratios, lambda_tau, linear_predicted_ratios, threshold_predicted_ratios, identifier, label, file_path):
    fig = plt.figure()
    ax = fig.add_subplot(111)

    plt.scatter(bins, actual_ratios)
    ax.axvline(lambda_tau, color='k', linestyle='--')
    ax.axhline(1, color='k', linestyle='--')
    linear = ax.plot(bins, linear_predicted_ratios, 'b', label="linear");
    threshold = ax.plot(bins, threshold_predicted_ratios, 'g', label="threshold");
    plt.ylabel(label)
    plt.xlabel("Population density")
    plt.gca().set_ylim(bottom=0)

    label = label.lower();
    label = label.replace(" ", "_");
    fig.savefig(os.path.join(file_path, str(identifier) + "_" + label + ".png"))
    plt.close()

def plot_odds_ratio(bins, actual_ratios, predictions,lambda_tau, linear_predicted_ratios, threshold_predicted_ratios, lower_cis, upper_cis, max_xaxis, identifier, label, file_path):
    fig = plt.figure()
    ax = fig.add_subplot(111)

    y_cis_lower = [];
    y_cis_upper = [];
    for i in range(0, len(bins)):
        lower_ci = lower_cis[i];
        upper_ci = upper_cis[i];
        ratio = actual_ratios[i];
        if lower_ci == "NA":
            y_cis_lower.append(0);
            y_cis_upper.append(0);
        else:
            lower = ratio - lower_ci;
            upper = upper_ci - ratio;
            y_cis_lower.append(lower);
            y_cis_upper.append(upper);
    # have commented out error bars
    #ax.errorbar(bins, actual_ratios, yerr=[y_cis_lower, y_cis_upper], fmt="o", capsize=3)
    # Got rid of threshold line
    plt.plot(bins,actual_ratios,'bo')
    # ax.axvline(lambda_tau, color='k', linestyle='--')
    # ax.axhline(1, color='k', linestyle='--')
    # linear = ax.plot(bins, linear_predicted_ratios, 'b', label="linear");
    # threshold = ax.plot(bins, threshold_predicted_ratios, 'g', label="threshold");
    plt.ylabel(label)
    plt.xlabel("Population density")
    plt.gca().set_ylim(0, 7)
    plt.gca().set_xlim(0, max_xaxis)

    label = label.lower();
    label = label.replace(" ", "_");
    fig.savefig(os.path.join(file_path, str(identifier) + "_" + label + ".png"))
    plt.close()


    
def plot_logit(bins, det_freq_values, unique_densities, logit_predictions, label, identifier, file_path):
    fig = plt.figure();
    if max(bins)/2 > len(bins):
        plt.scatter(bins, det_freq_values, color="blue");
    else:
        plt.plot(bins, det_freq_values, color="blue");
    # plt.plot(unique_densities, logit_predictions, color="red");
    plt.xlabel("Densities");
    plt.ylabel("Values");

    fig.tight_layout();
    # x_lim1 = max(bins);
    # x_lim2 = max(unique_densities);
    # x_lim = max(x_lim1, x_lim2);
    # plt.gca().set_xlim(0,max(bins)/2);
    y_lim=max(det_freq_values);
    # y_lim2=max(logit_predictions);
#    y_lim3=max(threshold_predictions)
    # y_lim=max(y_lim1,y_lim2);
    plt.gca().set_ylim(0,y_lim);
    fig.savefig(os.path.join(file_path, identifier + "_" + label + ".png"));
    plt.close();



def plot_detection_frequencies (bins, actual_ratios, logit_predictions, bin_size,predicted_threshold,max_xaxis, identifier, label, file_path):
    fig = plt.figure()

# =============================================================================
#     for i in range(0, len(bins)):
#         ratio = actual_ratios[i];
# =============================================================================
    add=bin_size/2
    bins=[x+add for x in bins]
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(bins,actual_ratios,'bo',label='Observations')
    ax.plot(bins,logit_predictions,'r',label='Logit curve')
#    ax.plot(bins,threshold_predictions,'g',label='Theoretical model')
    ax.axvline(predicted_threshold, color='k', linestyle='--',label='Estimated_threshold')
    #ax.axhline(1, color='k', linestyle='--')
    # threshold = ax.plot(bins, threshold_predictions, 'b', label="linear");
    # threshold = ax.plot(bins, threshold_predicted_ratios, 'g', label="threshold");
    plt.ylabel(label)
    plt.xlabel("Population density per cell")
#    plt.gca().set_xlim(0, max_xaxis)
    plt.legend(title= "Legend")
    #max yaxis needs to be set dynamically
    y_lim1=max(actual_ratios)
    y_lim2=max(logit_predictions)
#    y_lim3=max(threshold_predictions)
    y_lim=max(y_lim1,y_lim2)
    plt.gca().set_ylim(0,y_lim)
    label = label.lower();
    label = label.replace(" ", "_");
    fig.savefig(os.path.join(file_path, str(identifier) + "_" + label + ".png"))
    plt.close()


# def plot_likelihood_ratio(bins, actual_multipliers, lambda_tau, predicted_multipliers, label, identifier, file_path):
#     fig = plt.figure()
#     ax = fig.add_subplot(111)

#     plt.scatter(bins, actual_multipliers, label=label)
#     if label=="threshold":
#         ax.axvline(lambda_tau, color='k', linestyle='--')
#     ax.axhline(1, color='k', linestyle='--')
#     ax.plot(bins, predicted_multipliers, 'g')
#     plt.ylabel("p Likelihood Ratio")
#     plt.xlabel("Population density")
#     plt.gca().set_ylim(bottom=0)
#     fig.savefig(os.path.join(file_path, str(identifier) + "_" + label + "_ratio.png"))
#     plt.close()

def plot_sites_graph(bins, actual_sites, predicted_sites, label, identifier, title, file_path):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    # ax.set_title(title)
    ax.plot(bins,actual_sites,'r', label="actual")
    if label=="threshold":
        ax.plot(bins, predicted_sites,'g',label=label)
    else:
        ax.plot(bins, predicted_sites,'b',label=label)
    plt.ylabel("Count")
    plt.xlabel("Population density")
    plt.legend(title= "Legend")
    fig.savefig(os.path.join(file_path, str(identifier) + "_" + label + "_sites.png"))
    plt.close()

def plot_p_graphs(bins, p_samples, p_globals, bin_size,identifier, file_path):
    add=bin_size/2
    bins=[x+add for x in bins]
    fig = plt.figure()
    ax = fig.add_subplot(111)
    # ax.set_title(title)
    ax.plot(bins,p_samples,'b', label="Samples")
    ax.plot(bins, p_globals,'r',label="Globals")
    # ax.axvline(threshold, color='k', linestyle='--')
    plt.gca().set_ylim(0, 0.40)

    # text = "Threshold: " + str(threshold) + "\nSuccesses: " + str(success) + "\n Trials: " + str(trials) + "\nBinomial: " + str(binomial)
    
    plt.ylabel("Relative detection frequency")
    plt.xlabel("Population density per cell")
    plt.legend(title= "Legend")
    # plt.text(2, 2, text, horizontalalignment='center', verticalalignment='center', transform = ax.transAxes)
    fig.savefig(os.path.join(file_path, str(identifier) + "_p_graph.png"))
    plt.close()
    
def plot_cumulative_p_graphs(bins, p_samples, p_globals, bin_size, median_samples,median_globals,threshold,identifier, file_path):
    add=bin_size/2
    bins=[x+add for x in bins]
    fig = plt.figure()
    ax = fig.add_subplot(111)
    # ax.set_title(title)
    cum_samples=0
    cum_p_samples=[]
    for a_p_Sample in p_samples:
        cum_samples=cum_samples+a_p_Sample
        cum_p_samples.append(cum_samples)
    # normalize values because some p_samples may be missing
# =============================================================================
#     normalizer=1/cum_samples
#     cum_p_samples=(x * normalizer for x in cum_p_samples)
# =============================================================================
    cum_globals=0
    cum_p_globals=[]
    for a_p_globals in p_globals:
        cum_globals=cum_globals+a_p_globals
        cum_p_globals.append(cum_globals)
# =============================================================================
#     normalizer=1/cum_globals
#     cum_p_globals=(x * normalizer for x in cum_p_globals)
# =============================================================================
    
    ax.plot(bins,cum_p_samples,'b', label="Samples")
    ax.plot(bins, cum_p_globals,'r',label="Globals")
    # ax.axvline(threshold, color='k', linestyle='--')
    plt.gca().set_ylim(0, 1.0)

    # text = "Threshold: " + str(threshold) + "\nSuccesses: " + str(success) + "\n Trials: " + str(trials) + "\nBinomial: " + str(binomial)
    ax.axvline(median_samples, color='b', linestyle='--',label="Median samples")
    ax.axvline(median_globals, color='r', linestyle='--', label="Median globals")
    ax.axvline(threshold, color='g', linestyle='--', label="Threshold")
    plt.ylabel("Cumulated relative detection frequency")
    plt.xlabel("Population density per cell")
    plt.legend(title= "Legend")
    # plt.text(2, 2, text, horizontalalignment='center', verticalalignment='center', transform = ax.transAxes)
    fig.savefig(os.path.join(file_path, str(identifier) + "_cum_p_graph.png"))
    plt.close()

def plot_growth_coefficient_boxplot(data, file_path):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    # ax.set_title("Growth Coefficients")
    sample = data[0::2]
    controls = data[1::2]
    agg_data = [sample, controls]
    ax.boxplot(agg_data)
    plt.ylabel("value")
    plt.xticks([1,2],['Samples', 'Controls'])
    fig.savefig(os.path.join(file_path,"gc_boxplot.png"))
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


    # grouped_dataframe = dataframe[dataframe.type=='s'].groupby("cluster_id").first()
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

def plot_densities_on_map_by_time(population_data, time):
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

    #why is this divided by 6000 - don't understand
    dens = np.array(dens).astype(float)/3000;
    x,y = my_map(lons, lats)
    rgba = cmap(dens);

    # white = np.array([[1,1,1, 0] for i in range(0, len(dens))])
    # vector = white - rgba
    # temp_rgba = rgba + vector*0.4; 
    # rgba[:, :-1] = temp_rgba[:, 0:3]

    # print(rgba)


    ax = my_map.scatter(x, y, marker='o', c=rgba, s=3, zorder=1)


    sm = cm.ScalarMappable(cmap=cmap)
    sm.set_array([])
    #This used to be divided by 6000
    sm.set_clim([0, 3000])
    plt.colorbar(sm, orientation='horizontal', pad=0.03, aspect=50)

    # plt.title('Densities on Map: ' + population_data.name + " " + str(time) + "BP")

    plt.savefig("densitites_on_map_" + population_data.name + "_" + str(time) + "Kya.png");
    plt.close();

def plot_min_densities_in_time_range(population_data, time_from, time_to, min_density):
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

    time_indices = []
    for i in range(0, len(population_data.time_array)):
        time = population_data.time_array[i]*population_data.time_multiplier
        if time <= time_from and time >= time_to:
            time_indices.append(i)

    temp_time_densities = []
    for time_index in time_indices:
        temp_time_densities.append(population_data.density_array[time_index]);

    densities = np.mean(np.array(temp_time_densities), axis = 0);

    lons = []
    lats = []

    for i in range(0, len(densities)):
        if densities[i] >= min_density:
            lats.append(population_data.lat_array[i])
            lons.append(population_data.lon_array[i])

    x,y = my_map(lons, lats)
    my_map.plot(x, y, 'ro', markersize=5)

    plt.savefig("high_densitites_on_map_time_range" + str(time_to) + "-" + str(time_from));
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

    plt.savefig("densities_on_range_" + str(min_density) + "-" + str(max_density) + "_and_" + str(start_time) + "-" + str(end_time) + ".png", dpi=500)
    plt.close()

def plot_densities_on_map_snapshot(population_data, start_time, end_time):
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

    plt.savefig("densities_snapshot_" + str(start_time) + "-" + str(end_time) + ".png", dpi=500)
    plt.close()


def plot_boxplot_from_bin_file(bin_file_path):
    bins = []
    samples = []
    controls = []

    bin_file = open(bin_file_path)

    line = "-1"
    while line != '':
        line = bin_file.readline()
        line = bin_file.readline()
        line = bin_file.readline()
        line = bin_file.readline()

        line = bin_file.readline()
        line_values = line.split(";")
        bins.append(int(line_values[0]))
        samples.append(int(line_values[1]))
        controls.append(int(line_values[2]))

    dataframe = pd.DataFrame({'bins': bins, 'samples': samples, 'controls': controls})



if __name__ == '__main__':
        main()
        # df = pd.read_csv("r2_200_cross_validation_results.csv", sep=",")
        # r2_linear = df['r2_linear'].values
        # r2_threshold = df['r2_threshold'].values
        # values = [r2_linear, r2_threshold]
        # titles = ["linear", "threshold"]
        # file_path = os.getcwd()
        # plot_histogram(values, titles, "histo", file_path)


