import os
import csv
import numpy as np
import pandas as pd
from os import listdir
from copy import deepcopy
from os.path import isfile, join
from clusterer import ClusterAnalysis
from classes_module import Target, PopulationData

def process_targets(base_path, population_data, original_target_list, dataframe, controls_dataframe, controls, dataframe_loaded, clustering_on, date_window, critical_distance, critical_time, directory, min_date_window):
    
    #########################################
    # Cluster targets and group in clusters #
    #########################################

    my_cluster_analyzer=ClusterAnalysis()
    if clustering_on:
        clustered_list = my_cluster_analyzer.cluster_targets_by_dist_and_time(original_target_list, critical_distance, critical_time)
    else:
        clustered_list=original_target_list #uses default clustering used when list read from file

    # transform target list into a 2D array target list that groups clusters
    folded_target_list = fold_target_list(clustered_list)
    # clustered_list, folded_target_list = limit_target_list_to_oldest(folded_target_list)


    ###################################################
    # Extract dataframe and save as processed targets #
    ###################################################
    # - saves extracted dataframe as <directory>_dataframe.csv
    # - saves controls as <directory>_controls_df.csv
    # - saves target list as <directory>_targets.csv
    if dataframe_loaded is False:

        processed_targets_dir = os.path.join(base_path, "processed_targets")

        if not os.path.exists(processed_targets_dir):
            os.makedirs(processed_targets_dir)

        # extract controls dataframe depending on controls parameter
        if controls == "Australia":
            controls_dataframe = load_bin_controls_for_australia(population_data, clustered_list, date_window)
        elif controls == "France and Spain":
            controls_dataframe = load_bin_controls_for_francespain(population_data, clustered_list, date_window)
        elif controls == "Trial Latitudes":
            controls_dataframe = load_bin_controls_for_trial_latitudes(population_data, clustered_list, date_window)
        elif controls == "No Empty Lats":
            controls_dataframe = load_bin_controls_no_empty_lat(population_data, clustered_list, date_window)
        elif controls == "Trial Latitudes 2":
            controls_dataframe = load_bin_controls_for_trial_latitudes2(population_data, clustered_list, date_window)
        elif controls == "All":
            controls_dataframe = load_bin_controls_for_all(population_data, clustered_list, date_window)
        # save controls dataframe in processed_targets folder
        controls_dataframe_filename = os.path.join(processed_targets_dir, directory + "_controls_df.csv") 
        controls_dataframe.to_csv(controls_dataframe_filename, sep=";")

        # extract dataframe
        extracted_dataframe = extract_dataframe(population_data, folded_target_list, date_window, min_date_window)
        
        # save dataframe in processed_targets folder
        dataframe_filename = os.path.join(processed_targets_dir, directory + "_dataframe.csv")
        extracted_dataframe.to_csv(dataframe_filename, sep=";")
        dataframe = extracted_dataframe


        # save targets in processed_targets folder 
        targets_filename = directory + "_targets"
        save_target_list_to_csv(clustered_list, processed_targets_dir, targets_filename)


    return folded_target_list, dataframe, controls_dataframe

def limit_target_list_to_oldest(target_list):
    new_target_list = []
    new_folded_target_list = []

    for cluster in target_list:
        max_location = ""
        max_date = -1
        for target in cluster:
            if target.date_from > max_date:
                max_location = target.location
                max_date = target.date_from

        for target in cluster:
            if target.location == max_location:
                new_target_list.append(target)
                new_folded_target_list.append([target])

    print(len(new_target_list))
    return new_target_list, new_folded_target_list

def load_processed_targets(base_path, filename):
    processed_targets_dir = os.path.join(base_path, "processed_targets")

    chosen_file = filename + '_dataframe.csv'
    dataframe_filename = os.path.join(processed_targets_dir, chosen_file)
    dataframe = pd.read_csv(dataframe_filename, sep=";")

    controls_dataframe_filename = os.path.join(processed_targets_dir, filename + "_controls_df.csv")
    controls_dataframe = pd.read_csv(controls_dataframe_filename, sep=";")

    targets_filename = os.path.join(processed_targets_dir, filename + '_targets')
    target_list = read_target_list_from_csv(targets_filename)

    return target_list, dataframe, controls_dataframe

def clear_processed_targets(base_path):
    processed_targets_dir = os.path.join(base_path, "processed_targets")
    filenames = [f for f in os.listdir(processed_targets_dir) if isfile(join(processed_targets_dir, f))]
    for name in filenames:
        print("Removing: " + name)
        filepath = os.path.join(processed_targets_dir, name)
        os.remove(filepath)
    print("\n")


def extract_dataframe(population_data, target_list, date_window, min_date_window):

        
    lat_np = population_data.lat_array
    lon_np = population_data.lon_array
    time_np = population_data.time_array
    den_np = population_data.density_array
    time_multiplier = population_data.time_multiplier


    ############################
    # Create Time Dictionaries #
    ############################
    # time_dict[some_time] = index_in_time_array
    # next_time_dict[some_time] = next_time
    time_dict, next_time_dict = create_time_dictionaries(time_np, time_multiplier, population_data.ascending_time)
             
    latlon_length = len(lat_np)
    periods = []
    densities = []
    cluster_ids = []
    pseudo_types = []
    types = []
    locations = []
    contributions = []
    latitudes = []
    longitudes = []
    target_lat = []
    target_lon = []

    #######################
    # Loop through latlon #
    #######################
    for latlon_ind in range(0, latlon_length):
        if latlon_ind % 500 == 0:
            print('.   ')        
        
        #   skip columns (latlon) where the sum of column is 0
        #   In other words, skip latlon if sum across time is 0
        densities_in_latlonind = den_np[:, latlon_ind]
        sum_of_densities = np.nansum(densities_in_latlonind)      
        if np.isnan(sum_of_densities) or sum_of_densities <= 0:
            continue

        lat=lat_np[latlon_ind]
        lon=lon_np[latlon_ind]

        #######################################
        # For every target in every cluster.. #
        #######################################
        for cluster in target_list:
            number_of_targets_in_cluster = len(cluster)
            for target in cluster:


                ############################################
                # Loop from date_to to date_from of target #
                ############################################
                date_from = target.date_from
                date_to = target.date_to
                if date_from - date_to < date_window:
                    date_to = target.date_from - date_window;

                time = date_to
                if  date_to % time_multiplier != 0:
                    # get next value divisible by time_multiplier
                    # For example: 
                    #   time_multiplier = 25
                    #   date_to = 20
                    #   time = 20 + (25 - 20) = 25
                    time = date_to + (time_multiplier - date_to % time_multiplier)
                
                while(time <= date_from and time != -1):
                    # small loop is created here so that code is not repeated
                    for is_control in range(0,2):
                        # if time and lat/lon is in target conditions
                        if in_target_location(target, time, lat, lon, is_control == 1, date_window):
                            try:
                                time_index = time_dict[time]
                                density = den_np[time_index][latlon_ind]
                                if np.isnan(density):
                                    density=0
                                # save information only if density > 0
                                if density > 0:
                                    target_lat.append(target.orig_lat)
                                    target_lon.append(target.orig_lon)
                                    latitudes.append(lat)
                                    longitudes.append(lon)
                                    periods.append(time)
                                    densities.append(density)
                                    cluster_ids.append(target.cluster_id)
                                    locations.append(target.location)
                                    contributions.append(float(1)/number_of_targets_in_cluster)

                                    if is_control == 0:
                                        pseudo_types.append('a')
                                        types.append('s')
                                    else:
                                        pseudo_types.append('b')
                                        types.append('c')
                            except KeyError:
                                time = -1
                    try:
                        time = next_time_dict[time]
                    except KeyError:
                        time = -1
        

    new_df = pd.DataFrame({'location': locations, 'density': densities, 'period': periods, 'latitude': latitudes, 'longitude': longitudes, 'target_lat': target_lat, 'target_lon': target_lon, 'cluster_id':cluster_ids, 'pseudo_type':pseudo_types, 'type':types, 'contribution': contributions})
    return new_df



def create_time_dictionaries(time_np, time_multiplier, ascending_time):
    time_dict = dict()
    next_time_dict = dict()
    for i in range(-1, len(time_np)):
        time = -1
        if i > -1:
            time = time_np[i]*time_multiplier
            time_dict[time] = i;

        if ascending_time:
            # from 0 to len(time_np)-2
            if i > -1 and i < len(time_np)-1:
                next_time_dict[time_np[i]*time_multiplier] = time_np[i+1]*time_multiplier
            #the last time in the array (no next time)
            elif i == len(time_np)-1:
                next_time_dict[time_np[i]*time_multiplier] = -1
        else: 
            if i < len(time_np)-1:
                next_time_dict[time_np[i+1]*time_multiplier] = time

    return time_dict, next_time_dict

def load_bin_controls_for_all(population_data, target_list, date_window):

    lat_np = population_data.lat_array
    lon_np = population_data.lon_array
    time_np = population_data.time_array
    den_np = population_data.density_array
    time_multiplier = population_data.time_multiplier

    pseudo_target_list = []

    minimum_latitude = 90
    maximum_latitude = -90
    minimum_date = 150001
    maximum_date = -1
    ######################################
    # Getting the minimum/maximum values #
    ######################################
    for target in target_list:
        if target.orig_lat > maximum_latitude:
            maximum_latitude = target.orig_lat
        if target.orig_lat < minimum_latitude:
            minimum_latitude = target.orig_lat
        if target.date_from > maximum_date:
            maximum_date = target.date_from
        if target.date_from < minimum_date:
            minimum_date = target.date_from
    minimum_latitude = int(round(minimum_latitude))
    maximum_latitude = int(round(maximum_latitude)) 
    maximum_date=maximum_date+date_window
    print('minimum date=', minimum_date)
    print('maximum date=',maximum_date)
    print('minimum latitude=', minimum_latitude)
    print('maximum latitude=',maximum_latitude)
    lat = (minimum_latitude+maximum_latitude)/2
    time = 0
    temp_target = Target(lat,0,maximum_latitude,-1, minimum_latitude,1,"location",maximum_date, time,"country","date_of_reference","is_direct","calibrated","kind","figurative","source","is_controversial",0)

    df = extract_dataframe(population_data, [[temp_target]], maximum_date, minimum_date)
    df = df[df.type=='c']
    del df['location']
    del df['cluster_id']
    del df['pseudo_type']
    del df['type']
    del df['contribution']

    return df

def load_bin_controls_for_australia(population_data, target_list, date_window):
    minimum_date = 150001
    maximum_date = -1
    ######################################
    # Getting the minimum/maximum values #
    ######################################
    for target in target_list:
        if target.date_from > maximum_date:
            maximum_date = target.date_from
        if target.date_from < minimum_date:
            minimum_date = target.date_from
    maximum_date += date_window

    latitude = -25.27
    lat_nw = -11.17
    lat_se = -39.16
    longitude = 133.77
    lon_nw = 112.14
    lon_se = 154.86
    time = 0
    temp_target = Target(latitude,longitude,lat_nw,lon_nw,lat_se,lon_se,"location",maximum_date, time,"country","date_of_reference","is_direct","calibrated","kind","figurative","source","is_controversial",0)

    df = extract_dataframe(population_data, [[temp_target]], maximum_date, minimum_date)
    df = df[df.type=='s']
    del df['location']
    del df['cluster_id']
    del df['pseudo_type']
    del df['type']
    del df['contribution']
    return df


def load_bin_controls_for_francespain(population_data, target_list, date_window):
    minimum_date = 150001
    maximum_date = -1
    ######################################
    # Getting the minimum/maximum values #
    ######################################
    for target in target_list:
        if target.date_from > maximum_date:
            maximum_date = target.date_from
        if target.date_from < minimum_date:
            minimum_date = target.date_from
    maximum_date += date_window

    latitude = 42
    lat_nw = 50.84
    lat_se = 35
    longitude = 0
    lon_nw = 350
    lon_se = 7
    time = 0
    target_east = Target(latitude,longitude,lat_nw,lon_nw,lat_se,359,"location",maximum_date, time,"country","date_of_reference","is_direct","calibrated","kind","figurative","source","is_controversial",0)
    target_west = Target(latitude,longitude,lat_nw,0,lat_se,lon_se,"location",maximum_date, time, "country","date_of_reference","is_direct","calibrated","kind","figurative","source","is_controversial",0)

    df = extract_dataframe(population_data, [[target_east, target_west]], maximum_date, minimum_date)

    df = df[df.type=='s']
    del df['location']
    del df['cluster_id']
    del df['pseudo_type']
    del df['type']
    del df['contribution']
    return df

def load_bin_controls_for_trial_latitudes(population_data, target_list, date_window):
    df = load_bin_controls_for_all(population_data, target_list, date_window)
    new_df = df.loc[~df.latitude.between(-10, 40)]
    return new_df

def load_bin_controls_no_empty_lat(population_data, target_list, date_window):
    df = load_bin_controls_for_all(population_data, target_list, date_window)
    new_df = df.loc[~df.latitude.between(0, 20)]
    return new_df

def load_bin_controls_for_trial_latitudes2(population_data, target_list, date_window):
    df = load_bin_controls_no_empty_lat(population_data, target_list, date_window)
    new_df = df.loc[~df.latitude.between(-30, -10)]
    new_df = df.loc[~df.latitude.between(40, 50)]
    return new_df

def in_target_location(target, time, lat, lon, is_control, date_window):
# Used to deal to test if a particular hexagon with a particular density matches at least one of a list of targets in time and location
#Used to avoid duplicates when dealing with geographical clusters of targets
# Time still has to be checked because some targets could match location but not time
    lat_nw=target.lat_nw
    # if we are testing whether a hexagon is in controls we use all possible longitudes
    #would be better not to control logitudes at all
    if is_control:
        lon_nw=target.lon_se  
    else:
        lon_nw=target.lon_nw
    lat_se=target.lat_se
    lon_se=target.lon_se
    date_from=target.date_from
    date_to=target.date_to
    if date_from - date_to < date_window:
        date_to = target.date_from - date_window;

    if lat<-90 or lat>90:
        print('latitude out of range')
        sys.exit()
    if lon<0 or lon>360:
        print('longitude out of range')
        sys.exit()
    if time<=date_from:
        if time>=date_to:
            if lat<=lat_nw:
                if lat>=lat_se:
                    if lon_in_bin(lon_nw,lon_se,lon):
                        return(True)
    return(False)

def lon_in_bin(target_nw,target_se,lon):
 # calculates if a longtitude lies between lon_nw and lon_se moving eastwards
 # it would be possible to calculate clockArithmetic just once for whole of computeStats
 # but this way it is hidden and doesn't pollute the code
 # print 'In lon_in_bin longNW=',lon_nw,'lon_se=',lon_se,'target=',target
 # check if target lies between lon_nw and lonSW  
    if lon>=target_nw and lon<=target_se:
            normal_in_bin=True
            # print(lon, target_nw, target_se)
    else:
            normal_in_bin=False
# check if we are doing W to E (normal arithmetic) or E to west (abnormal)
    if target_nw<target_se:
        normal_arithmetic=True
    else:
        normal_arithmetic=False
    if normal_arithmetic:
        return(normal_in_bin)
    else:
        return(not(normal_in_bin))

def read_target_list_from_csv(filename):
# reads list of targets from a CSV file - each is assigned a sequential cluster ID which can later be replaced during clustering
    
    target_list=[]

    with open(filename +".csv",'rb') as target_file:
        my_reader=csv.reader(target_file,delimiter=',')
        next(my_reader,None) # skip the headers
        cluster_id=0
        for row in my_reader:
            location=row[0]
            location = location.replace('\n', ' ').replace('\r', '') #avoids misreading by csv file reader
            lat=float(row[1])
            lon=float(row[2])
            if lon<0:
                lon=lon+360
            # date_to=None
            lat_nw=lat+0.5
            lon_nw=lon-1
            lat_se=lat-0.5
            lon_se=lon+1
            date_from=0;
            date_to=1000000;
            dates = row[3].split(";");
            for string_date in dates:
                date = int(string_date);
                if date > date_from:
                    date_from = date;
                if date < date_to:
                    date_to = date;
            country=row[4]
            date_of_reference=row[5]
            is_direct=row[6]
            calibrated=row[7]
            kind=row[8]
            figurative=row[9]
            source=row[10]
            is_controversial=row[11]
            target=Target(lat,lon,lat_nw,lon_nw,lat_se,lon_se,location,date_from, date_to, country,date_of_reference,is_direct,calibrated,kind,figurative,source,is_controversial,cluster_id)
            cluster_id=cluster_id+1
            target_list.append(target)
    target_file.close
    return(target_list)


def save_target_list_to_csv(target_list, directory, filename):
    old_dir = os.getcwd() 
    os.chdir(directory)
    with open(filename +".csv",'wb') as target_file:
        my_writer=csv.writer(target_file,delimiter=',')
        headers='location,latitude,longitude,date,country,date of reference,is_direct,calibrated,kind,figurative,source,is_controversial,cluster_id'
        my_writer.writerow(headers)
        for target in target_list:
            print('saving ',target.location)
            row=[]
            row.append(target.location)
            row.append(target.orig_lat)
            row.append(target.orig_lon)
            row.append(str(target.date_from) + ";" + str(target.date_to));
            row.append(target.country)
            row.append(target.date_of_reference)
            row.append(target.is_direct)
            row.append(target.calibrated)
            row.append(target.kind)
            row.append(target.figurative)
            row.append(target.source)
            row.append(target.is_controversial)
            row.append(target.cluster_id)
            my_writer.writerow(row)
    target_file.close
    os.chdir(old_dir)


def fold_target_list(target_list):
    # This function creates a target list where each item is a list of targets belonging to the same cluster.
    #To be executed only after clustering
    
    max_target=max(target.cluster_id for target in target_list)
    new_list=[]
    for i in range(0,max_target+1):
        new_list.append([])
    for target in target_list:
        new_list[target.cluster_id].append(deepcopy(target))
    return(new_list)

def filter_targets_for_date_before(target_list, value, filters_applied):
    filtered_list=[target for target in target_list if target.date_from>=value]
    filters_applied=filters_applied+" Excluded targets more recent than "+str(int(value))+" BP;"
    return filtered_list, filters_applied
    
def filter_targets_for_not_direct(target_list, filters_applied):
    print("length of original list=",len(target_list))
    filtered_list=[target for target in target_list if target.is_direct=='Yes']
    print("length of filtered list=",len(filtered_list))
    filters_applied=filters_applied+" "+"Excluded targets with indirect measurements;"
    return filtered_list, filters_applied

def filter_targets_for_not_figurative(target_list, filters_applied):
    print("length of original list=",len(target_list))
    filtered_list=[target for target in target_list if target.figurative=='Yes']
    print("length of filtered list=",len(filtered_list))
    filters_applied=filters_applied+" "+"Excluded targets with non-figurative art;"
    return filtered_list, filters_applied

def filter_targets_for_not_controversial(target_list, filters_applied):
    print("length of original list=",len(target_list))
    filtered_list=[target for target in target_list if target.is_controversial=='No']
    print("length of filtered list=",len(filtered_list))
    filters_applied=filters_applied+" "+"Excluded targets with disputed measurements;"
    return filtered_list, filters_applied

def filter_targets_for_date(target_list, minimum_date, maximum_date, filters_applied):
    filtered_targets = []
    for target in target_list:
        if target.date_from >= minimum_date and target.date_from <= maximum_date:
            filtered_targets.append(target)

    filters_applied += " Only targets within" + str(minimum_date) + " - " + str(maximum_date) + " date"

    return filtered_targets, filters_applied

def filter_targets_for_abs_latitude(target_list, minimum_lat, maximum_lat, filters_applied):
    filtered_targets = []
    for target in target_list:
        abs_lat = abs(target.orig_lat)
        if abs_lat >= minimum_lat and abs_lat <= maximum_lat:
            filtered_targets.append(target)

    filters_applied += " Only targets within" + str(minimum_lat) + " - " + str(maximum_lat) + " absolute latitude"

    return filtered_targets, filters_applied

def filter_targets_for_latitude(target_list, minimum_lat, maximum_lat, filters_applied):
    filtered_targets = []
    for target in target_list:
        lat = target.orig_lat
        if lat >= minimum_lat and lat <= maximum_lat:
            filtered_targets.append(target)

    filters_applied += " Only targets within" + str(minimum_lat) + " - " + str(maximum_lat) + " latitude"

    return filtered_targets, filters_applied

def print_target(target, date_window):
    print(target.location,': cluster_id:',target.cluster_id, ' lat_nw:',target.lat_nw,' lon_nw:',target.lon_nw,'lat_se:',target.lat_se,'lon_se:',target.lon_se,'date_from:',target.date_from,'  date_to:',target.date_from+date_window)
    