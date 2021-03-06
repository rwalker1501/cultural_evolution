
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
import os
from matplotlib.pyplot import cm 
from matplotlib.font_manager import FontProperties
from mpl_toolkits.basemap import Basemap
# import seaborn

    
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


def plot_likelihood_ratio(bins, actual_multipliers, lambda_tau, predicted_multipliers, label, identifier, file_path):
    fig = plt.figure()
    ax = fig.add_subplot(111)

    plt.scatter(bins, actual_multipliers, label=label)
    if label=="threshold":
        ax.axvline(lambda_tau, color='k', linestyle='--')
    ax.axhline(1, color='k', linestyle='--')
    ax.plot(bins, predicted_multipliers, 'g')
    plt.ylabel("p Likelihood Ratio")
    plt.xlabel("Population density")
    plt.gca().set_ylim(bottom=0)
    fig.savefig(os.path.join(file_path, str(identifier) + "_" + label + "_ratio.png"))
    plt.close()

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

def plot_p_graphs(bins, p_samples, p_controls, threshold, identifier, file_path):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    # ax.set_title(title)
    ax.plot(bins,p_samples,'b', label="samples")
    ax.plot(bins, p_controls,'r',label="controls")
    ax.axvline(threshold, color='k', linestyle='--')

    # text = "Threshold: " + str(threshold) + "\nSuccesses: " + str(success) + "\n Trials: " + str(trials) + "\nBinomial: " + str(binomial)
    
    plt.ylabel("p")
    plt.xlabel("Population density")
    plt.legend(title= "Legend")
    # plt.text(2, 2, text, horizontalalignment='center', verticalalignment='center', transform = ax.transAxes)
    fig.savefig(os.path.join(file_path, str(identifier) + "_p_graph.png"))
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


    grouped_dataframe = dataframe[dataframe.type=='s'].groupby("cluster_id").first()
    lats = grouped_dataframe['latitude'].values
    lons = grouped_dataframe['longitude'].values


    controls_lats = controls_dataframe['latitude'].values
    controls_lons = controls_dataframe['longitude'].values


    x,y = my_map(lons, lats)
    c_x, c_y = my_map(controls_lons, controls_lats)
    my_map.plot(c_x, c_y, 'ro', markersize=3)
    my_map.plot(x, y, 'go', markersize=3)

    plt.savefig(os.path.join(dir_path, directory + "_target_on_map.png"))
    plt.close()

def plot_densities_on_map_by_time(population_data, time):
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

    x,y = my_map(lons, lats)
    my_map.plot(x, y, 'go', markersize=10)

    plt.show()


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


