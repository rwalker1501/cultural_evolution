
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
import os
from matplotlib.pyplot import cm 
from matplotlib.font_manager import FontProperties
from mpl_toolkits.basemap import Basemap
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

def get_map_file_path(filename):
    base_path = os.getcwd();
    map_path = os.path.join(base_path, "maps")
    if not os.path.exists(map_path):
        os.makedirs(map_path);
    file_path = os.path.join(map_path, filename);
    return file_path

if __name__ == '__main__':
        main()


