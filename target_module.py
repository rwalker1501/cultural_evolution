import os
import csv
import re
import numpy as np
import pandas as pd
from os import listdir
from copy import deepcopy
from os.path import isfile, join
from clusterer import ClusterAnalysis
from classes_module import Target, PopulationData

def process_targets(base_path, population_data, original_target_list, dataframe, globals_dataframe, parameters, max_for_uninhabited, directory):
    
    #########################################
    # Cluster targets and group in clusters #
    #########################################

    my_cluster_analyzer=ClusterAnalysis()
    if parameters['clustering_on']:
        clustered_list = my_cluster_analyzer.cluster_targets_by_dist_and_time(original_target_list, parameters['critical_distance'], parameters['critical_time'])
    else:
        clustered_list=original_target_list #uses default clustering used when list read from file

    # transform target list into a 2D array target list that groups clusters
    folded_target_list = fold_target_list(clustered_list)
    # clustered_list, folded_target_list = limit_target_list_to_oldest(folded_target_list)

    ###################################################
    # Extract dataframe and save as processed targets #
    ###################################################
    # - saves extracted dataframe as <directory>_dataframe.csv
    # - saves the_globals as <directory>_globals_df.csv
    # - saves target list as <directory>_targets.csv
    if parameters['dataframe_loaded'] is False:

        processed_targets_dir = os.path.join(base_path, "processed_targets")

        if not os.path.exists(processed_targets_dir):
            os.makedirs(processed_targets_dir)

        date_window = parameters['date_window'];
        date_lag = parameters['date_lag'];
        min_date = parameters['min_date'];
        max_date = parameters['max_date'];
        min_lat = parameters['min_lat'];
        max_lat = parameters['max_lat'];

        # extract the_globals dataframe depending on the_globals parameter
        the_globals = parameters['globals_type']
        brute_all_globals_df = load_all_globals_brute(population_data, min_lat, max_lat, min_date, max_date, max_for_uninhabited)
        if the_globals == "Australia":
            globals_dataframe = load_bin_globals_for_australia(brute_all_globals_df);
        elif the_globals == "France and Spain":
            globals_dataframe = load_bin_globals_for_francespain(brute_all_globals_df);
        elif the_globals == "All Except Equatorials":
            globals_dataframe = load_bin_globals_equatorials_exclusive(brute_all_globals_df);
        elif the_globals == "All":
            globals_dataframe = brute_all_globals_df;

        globals_dataframe_filename = os.path.join(processed_targets_dir, directory + "_globals_df.csv") 
        
        print("Saving globals dataframe...")
        # WRITE PROCESSED GLOBALS
        globals_dataframe.to_csv(globals_dataframe_filename, sep=";")

        # extract dataframe rank
        print("Extracting sites from targets...")
        new_df = extract_dataframe(population_data, folded_target_list, max_for_uninhabited, date_window, date_lag)
        
        # get closest site to target
        new_df['distance'] = (new_df['latitude']-new_df['target_lat'])*(new_df['latitude']-new_df['target_lat']) + (new_df['longitude']-new_df['target_lon'])*(new_df['longitude']-new_df['target_lon'])
        new_df['rank'] = new_df.groupby(['target_lat', 'target_lon', 'period'])['distance'].rank(method="first",ascending=True);
        
        samples_df = new_df[(new_df['type'] == 's') & (new_df['rank'] < 2)];
        samples_df = samples_df.groupby(['density', 'latitude', 'longitude', 'period', 'type']).first().reset_index();
        controls_df = new_df[new_df['type'] == 'c'];

        dataframe = pd.concat([samples_df, controls_df]);

        print("Saving sites dataframe...")
        # save dataframe in processed_targets folder
        dataframe_filename = os.path.join(processed_targets_dir, directory + "_dataframe.csv")
        
        # WRITE PROCESSED TARGET DATEFRAME
        dataframe.to_csv(dataframe_filename, sep=";",quoting=csv.QUOTE_NONNUMERIC) #this is an addition to get file written in good format for excel

        print("Saving target list...")
        # save targets in processed_targets folder 
        targets_filename = directory + "_targets"
        save_target_list_to_csv(clustered_list, processed_targets_dir, targets_filename)
        # exit()


    return folded_target_list, dataframe, globals_dataframe

def reset_cluster_id_values(target_list):
    for i in range(0, len(target_list)):
        target_list[i].cluster_id = i;
    return target_list;

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
    return new_target_list, new_folded_target_list

def load_processed_targets(base_path, filename):
    processed_targets_dir = os.path.join(base_path, "processed_targets")

    chosen_file = filename + '_dataframe.csv'
    dataframe_filename = os.path.join(processed_targets_dir, chosen_file)
    dataframe = pd.read_csv(dataframe_filename, sep=";")

    globals_dataframe_filename = os.path.join(processed_targets_dir, filename + "_globals_df.csv")
    globals_dataframe = pd.read_csv(globals_dataframe_filename, sep=";")

    targets_filename = os.path.join(processed_targets_dir, filename + '_targets')
    target_list = read_target_list_from_csv(targets_filename)

    return target_list, dataframe, globals_dataframe

def clear_processed_targets(base_path):
    processed_targets_dir = os.path.join(base_path, "processed_targets")
    filenames = [f for f in os.listdir(processed_targets_dir) if isfile(join(processed_targets_dir, f))]
    for name in filenames:
        print("Removing: " + name)
        filepath = os.path.join(processed_targets_dir, name)
        os.remove(filepath)
    print("\n")


def extract_dataframe(population_data, target_list, max_for_uninhabited, date_window, date_lag):

        
    lat_np = population_data.lat_array
    lon_np = population_data.lon_array
    time_np = population_data.time_array
    den_np = population_data.density_array
    time_multiplier = population_data.time_multiplier
    density_multiplier=population_data.density_multiplier
    smallest_time = min(time_np)*time_multiplier

    ############################
    # Create Time Dictionaries #
    ############################
    # time_dict[some_time] = index_in_time_array
    # next_time_dict[some_time] = next_time
    time_dict, next_time_dict = create_time_dictionaries(time_np, time_multiplier, population_data.ascending_time)
             
    latlon_length = len(lat_np)
    target_location = [];
    target_date_from = [];
    target_date_to = [];
    target_lats = [];
    target_lons = [];
    periods = []
    densities = []
    pseudo_types = []
    types = []
    contributions = []
    latitudes = []
    longitudes = []
    dirs = []
    exacts = []
    cluster_ids = []

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
            # if cluster[0].cluster_id != 18:
            #     continue;
            number_of_targets_in_cluster = len(cluster)
            for target in cluster:


                ###################################################
                # Loop from date_from to date_from + dw of target #
                ###################################################
                date_from = target.date_from + date_lag + date_window
                # date_to = target.date_to
                date_to = target.date_from + date_lag

                time = date_to
                if smallest_time > time:
                    time = smallest_time
                if  time % time_multiplier != 0:
                    # get next value divisible by time_multiplier
                    # For example: 
                    #   time_multiplier = 25
                    #   date_to = 40
                    #   time = 40 + (25 - 40%25) = 40 + (25 - 15) = 50
                    time = time + (time_multiplier - time % time_multiplier)
                
                # # made exclusive (used to include date if divisible by multiplier)
                # time = time + (time_multiplier - time % time_multiplier)
                # made exclusive (used to be time <= date_from)
                while(time <= date_from and time != -1):
                    # small loop is created here so that code is not repeated
                    for is_global in range(0,2):
                        # if time and lat/lon is in target conditions
                        if in_target_location(target, lat, lon, is_global == 1):
                            try:
                                time_index = time_dict[time]
                                density = den_np[time_index][latlon_ind]*density_multiplier
                                if np.isnan(density):
                                    density=0
                                # save information only if density > max_for_uninhabited
                                if density > max_for_uninhabited:
                                    cluster_ids.append(target.cluster_id);
                                    target_location.append(target.location);
                                    target_date_from.append(target.date_from);
                                    target_date_to.append(target.date_to);
                                    latitudes.append(lat)
                                    longitudes.append(lon)
                                    target_lats.append(target.orig_lat);
                                    target_lons.append(target.orig_lon);

                                    periods.append(time)
                                    densities.append(density)
                                    contributions.append(float(1)/number_of_targets_in_cluster)

                                    if target.is_direct == "Yes":
                                        dirs.append(True)
                                    else:
                                        dirs.append(False)

                                    if target.age_estimation.lower() == "exact age":
                                        exacts.append(True)
                                    else:
                                        exacts.append(False)

                                    if is_global == 0:
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
    
    new_df = pd.DataFrame({'target_location': target_location, 'target_date_from': target_date_from, 'target_date_to': target_date_to, 'cluster_id': cluster_ids, 'density': densities, 'period': periods, 'latitude': latitudes, 'longitude': longitudes, 'target_lat': target_lats, 'target_lon': target_lons, 'pseudo_type':pseudo_types, 'type':types, 'contribution': contributions, 'is_dir': dirs, 'is_exact': exacts})
    

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
            # print(ascending_time)
            if i < len(time_np)-1:
                next_time_dict[time_np[i+1]*time_multiplier] = time

    return time_dict, next_time_dict


def load_all_globals_brute(population_data, min_lat, max_lat, min_date, max_date, max_for_uninhabited):
    lat_np = population_data.lat_array
    lon_np = population_data.lon_array
    time_np = population_data.time_array
    den_np = population_data.density_array
    time_multiplier = population_data.time_multiplier
    density_multiplier = population_data.density_multiplier;


    # MADE EXCLUSIVE
    # time_mask = (time_np*time_multiplier <= max_date) & (time_np*time_multiplier >= min_date);
    time_mask = (time_np*time_multiplier < max_date) & (time_np*time_multiplier > min_date);
    
    print("Min date: " + str(min_date));
    print("Max date: " + str(max_date));
    print("Min lat: " + str(min_lat));
    print("Max lat: " + str(max_lat));


    # MADE EXCLUSIVE
    latlon_mask = (lat_np <= max_lat) & (lat_np >= min_lat);
    # latlon_mask = (lat_np < max_lat) & (lat_np > min_lat);

    mask = latlon_mask[np.newaxis,:] & time_mask[:, np.newaxis];

    latlon_length = len(lat_np);
    time_length = len(time_np);
    indices = np.reshape(range(latlon_length*time_length), [-1, latlon_length])
    valid_ind = indices[mask];


    # print("Latlon: " + str(latlon_length))
    # print("Time: " + str(time_length))
    # print("density: " + str(len(den_np)));
    # print(indices[time_length-1][latlon_length-1]/(latlon_length))
    # print(latlon_length*time_length/(time_length-1))

    print("Masking lat, lon, time..")
    periods = time_np[valid_ind/(latlon_length)]*time_multiplier

    valid_latlon_ind = valid_ind%latlon_length
    latitudes = lat_np[valid_latlon_ind]
    longitudes = lon_np[valid_latlon_ind]

    print("Masking densities");
    densities = den_np[mask]*density_multiplier;

    print("Generating dataframe...")
    new_df = pd.DataFrame({'density': densities, 'period': periods, 'latitude': latitudes, 'longitude': longitudes})
    print("Filtering...")
    new_df = new_df[new_df.density > max_for_uninhabited];
    print("Globals dataframe generated.")

    return new_df


def load_bin_globals_for_australia(brute_all_globals_df):
    
    # Boundaries
    latitude = -25.27
    max_lat = -11.17
    min_lat = -39.16
    longitude = 133.77
    min_lon = 112.14
    max_lon = 154.86

    df = brute_all_globals_df[(brute_all_globals_df.latitude >= min_lat) & (brute_all_globals_df.latitude <= max_lat) & (brute_all_globals_df.longitude >= min_lon) & (brute_all_globals_df.longitude <= max_lon)]
    return df


def load_bin_globals_for_francespain(brute_all_globals_df):

    # Boundaries 
    latitude = 42
    max_lat = 50.84
    min_lat = 35
    longitude = 0
    # 0 in between, special condition
    min_lon = 350
    max_lon = 7

    df = brute_all_globals_df[(brute_all_globals_df.latitude >= min_lat) & (brute_all_globals_df.latitude <= max_lat) & ((brute_all_globals_df.longitude >= min_lon) | (brute_all_globals_df.longitude <= max_lon))]
    return df

def load_bin_globals_equatorials_exclusive(brute_all_globals_df):
    equatorial_north = 20;
    equatorial_south = -10;
    df = brute_all_globals_df[~((brute_all_globals_df.latitude <= equatorial_north) & (brute_all_globals_df.latitude >= equatorial_south))];

    return df;

def in_target_location(target, lat, lon, is_global):
# Used to deal to test if a particular hexagon with a particular density matches at least one of a list of targets in time and location
#Used to avoid duplicates when dealing with geographical clusters of targets
# Time still has to be checked because some targets could match location but not time
    lat_nw=target.lat_nw
    lat_se=target.lat_se
    lon_se=target.lon_se
    lon_nw=target.lon_nw

    if lat<-90 or lat>90:
        print('latitude out of range')
        sys.exit()
    if lon<0 or lon>360:
        print('longitude out of range')
        sys.exit()
    
    if lat<=lat_nw:
        if lat>=lat_se:
            if is_global:
                return True;
            if lon_in_bin(lon_nw,lon_se,lon):
                return(True)
            # if lon >= lon_nw and lon <= lon_se:
            #     return True;
    return False

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

def read_target_list_format1_from_csv(filename):
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
            lat_nw=lat+1
            lon_nw=lon-2
            lat_se=lat-1
            lon_se=lon+2
            date=int(float(row[3]))
            country=row[4]
            date_of_reference=row[5]
            is_direct=row[6]
            calibrated=row[7]
            kind=row[8]
            figurative=row[9]
            source=row[10]
            is_controversial=row[11]
            target=Target(lat,lon,lat_nw,lon_nw,lat_se,lon_se,location,date,date,country,date_of_reference,is_direct,calibrated,kind,figurative,source,is_controversial,cluster_id)
            cluster_id=cluster_id+1
            target_list.append(target)
    target_file.close
    return(target_list)

def read_target_list_from_csv(filename):

    target_list=[]

    with open(filename +".csv",'rb') as target_file:
        my_reader=csv.reader(target_file,delimiter=',')
        next(my_reader,None) # skip the headers
        cluster_id=0
        for row in my_reader:
            location=row[0]
            location = location.replace('\n', ' ').replace('\r', '') 
            lat=float(re.sub("[^0-9.-]", "", row[1]))
            lon=float(re.sub("[^0-9.-]", "", row[2]))
            if lon<0:
                lon=lon+360
            # date_to=None
            lat_nw=lat+1
            lon_nw=lon-2
            lat_se=lat-1
            lon_se=lon+2
            date_from = int(re.sub("[^0-9.-]", "", row[3]));
            date_to = date_from;
            if row[4] == "Modern":
                date_to = 0;
            elif row[4] == "Single sample" or len(row[4]) == 0:
                date_to = date_from;
            else:
                date_to = int(re.sub("[^0-9.-]", "", row[4]))

            country=row[5]
            is_direct="Yes" if row[6]=="Direct" else "No"
            age_est = row[7]
            calibrated="Yes" if row[8]=="Yes" else "No"
            kind=row[9]
            figurative=row[10];
            if row[10] == "Unknown" or row[10] == "-":
                figurative = "No";
            elif len(row[10]) > 2:
                figurative = "Yes" if row[10][0:3] == "Yes" else "No"
            source="Primary"
            is_controversial="No"
            target=Target(lat,lon,lat_nw,lon_nw,lat_se,lon_se,location,date_from, date_to, country,is_direct,calibrated,kind,figurative,source,is_controversial, age_est, cluster_id)
            cluster_id=cluster_id+1
            target_list.append(target)
    target_file.close
    return(target_list)


def save_target_list_to_csv(target_list, directory, filename):
    old_dir = os.getcwd() 
    os.chdir(directory)
    with open(filename +".csv",'wb') as target_file:
        my_writer=csv.writer(target_file,delimiter=',')
        headers='location,latitude,longitude,date_from, date_to, country, is_direct,age_estimation,calibrated,kind, figurative, cluster_id'
        my_writer.writerow(headers)
        for target in target_list:
            row=[]
            row.append(target.location)
            row.append(target.orig_lat)
            row.append(target.orig_lon)
            row.append(target.date_from);
            row.append(target.date_to);
            row.append(target.country)
            row.append(target.is_direct)
            row.append(target.age_estimation);
            row.append(target.calibrated)
            row.append(target.kind)
            row.append(target.figurative)
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

def filter_targets_for_not_exact_age(target_list, filters_applied):
    print("length of original list=",len(target_list))
    filtered_list=[target for target in target_list if target.age_estimation.lower()=='exact age']
    print("length of filtered list=",len(filtered_list))
    filters_applied=filters_applied+" "+"Only include targets with exact age;"
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
        if target.date_to >= minimum_date and target.date_from <= maximum_date:
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

def get_min_max_date_of_targets(target_list):
    minimum_date = 1000000;
    maximum_date = 0;
    for target in target_list:
        if target.date_from > maximum_date:
            maximum_date = target.date_from;
        if target.date_to < minimum_date:
            minimum_date = target.date_to;
    return minimum_date, maximum_date;


def get_min_max_latitude_of_targets(target_list):
    minimum_lat = 1000;
    maximum_lat = -1000;
    for target in target_list:
        if target.orig_lat > maximum_lat:
            maximum_lat = target.orig_lat;
        if target.orig_lat < minimum_lat:
            minimum_lat = target.orig_lat;
    return minimum_lat, maximum_lat;


def create_binned_column(dataframe, new_column_name, base_column_name, interval):

    conditions = [];
    choices = [];
    min_val = int(dataframe[base_column_name].min()/interval);
    max_val = int(dataframe[base_column_name].max()/interval) + 1;

    for i in range(min_val, max_val):
        lower_bound = i*interval;
        upper_bound = i*interval + (interval-1);
        condition = ((dataframe[base_column_name] >= lower_bound) & (dataframe[base_column_name] <= upper_bound));
        choice = str(lower_bound) + "-" + str(upper_bound);
        conditions.append(condition);
        choices.append(choice);

    dataframe[new_column_name] = np.select(conditions, choices);

    return dataframe

def generate_merged_dataframe(base_path,directory, dataframe, globals_dataframe):
    processed_targets_dir = os.path.join(base_path, "processed_targets")
    merged_df_filename = os.path.join(processed_targets_dir, directory + "_merged_df.csv") 
    if not os.path.exists(processed_targets_dir):
        os.makedirs(processed_targets_dir)
    temp_globals_df = globals_dataframe.copy();
    temp_samples_df = dataframe.copy();
    temp_samples_df = temp_samples_df[temp_samples_df.type == 's'];

     # DELETE COLUMNS
    temp_samples_df.drop(["cluster_id","samples_growth_coefficient", "contribution", "distance", "is_dir", "is_exact", "pseudo_type", "rank", "target_date_from", "target_date_to", "target_lat", "target_location", "target_lon","type"], axis=1, inplace = True)

    # is_sample
    temp_globals_df['is_sample'] = 0;    
    temp_samples_df['is_sample'] = 1;

    # MERGE
    to_concat = [temp_globals_df, temp_samples_df]
    merged_df = pd.concat(to_concat, sort=True);
    
    # abs lat
    merged_df['abs_latitude'] = merged_df['latitude'].abs();

    # binned_period
    merged_df = create_binned_column(merged_df, 'binned_period', 'period', 1000);

    # binned_period
    merged_df = create_binned_column(merged_df, 'binned_latitude', 'latitude', 10);


    merged_df.to_csv(merged_df_filename, sep=";")
    print("merged df filename: " + merged_df_filename);
    return merged_df