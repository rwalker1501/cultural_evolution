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

def process_targets(base_path, population_data, target_list, parameters):

    min_date = parameters['min_date'];
    max_date = parameters['max_date'];
    min_lat = parameters['min_lat'];
    max_lat = parameters['max_lat'];
    max_for_uninhabited = parameters['max_for_uninhabited'];
    globals_type = parameters['globals_type'];



    target_list = filter_targets_for_latitude(target_list, min_lat, max_lat);
    target_list = filter_targets_for_date(target_list, min_date, max_date);
    if parameters['remove_not_direct_targets']:
        target_list = filter_targets_for_not_direct(target_list);
    if parameters['remove_not_exact_age_targets']:
        target_list = filter_targets_for_not_exact_age(target_list);
    if parameters['remove_not_figurative_targets']:
        target_list = filter_targets_for_not_figurative(target_list);
    if parameters['remove_not_controversial_targets']:
        target_list = filter_targets_for_not_controversial(target_list);



    globals_dir = os.path.join(base_path, "globals");
    if not os.path.exists(globals_dir):
        os.makedirs(globals_dir);
    filenames_in_globals = [f for f in os.listdir(globals_dir) if isfile(join(globals_dir,f))]

    globals_filename = (globals_type + "_lat_" + str(min_lat) + "-" + str(max_lat) + "_date_" + str(min_date) + "-" + str(max_date) +  "_mfu_" + str(max_for_uninhabited) + ".csv").lower().replace(" ", "_");
    globals_dataframe_path = os.path.join(globals_dir, globals_filename);
    if globals_filename not in filenames_in_globals:
        if globals_type == "Australia":
                globals_dataframe = load_bin_globals_for_australia(population_data, min_date, max_date, max_for_uninhabited)
        elif globals_type=="No equatorials":
            globals_dataframe = load_bin_globals_for_no_equatorials(population_data, min_lat, max_lat,min_date, max_date,max_for_uninhabited)
        elif globals_type == "France and Spain":
            globals_dataframe = load_bin_globals_for_francespain(population_data, min_date, max_date, max_for_uninhabited)
        elif globals_type == "All":
                globals_dataframe = load_all_globals_brute(population_data, min_lat, max_lat, min_date, max_date, max_for_uninhabited)

        print("Generating globals...")
        globals_dataframe.to_csv(globals_dataframe_path, sep=";")
    else:
        print("Reading globals file...")
        globals_dataframe = pd.read_csv(globals_dataframe_path, sep=";");

    ###################################################
    # Extract dataframe and save as processed targets #
    ###################################################
    # - saves extracted dataframe as <directory>_dataframe.csv
    # - saves the_globals as <directory>_globals_df.csv
    # - saves target list as <directory>_targets.csv
    loaded_processed_targets = False
    if parameters['use_processed_targets']:
        loaded_processed_targets, target_list, dataframe, globals_dataframe = load_processed_targets(base_path, directory, globals_type, min_lat, min_date, max_date, max_for_uninhabited)
        if loaded_processed_targets is False:
            print("Load processed targets failed. Processing targets with existing parameters..")


    if loaded_processed_targets is False:

        processed_targets_dir = os.path.join(base_path, "processed_targets")

        if not os.path.exists(processed_targets_dir):
            os.makedirs(processed_targets_dir)

        # extract dataframe rank
        print("Extracting sites from targets...")
        new_df = extract_dataframe(population_data, target_list, max_for_uninhabited)
        
        # get closest site to target
        new_df['distance'] = (new_df['latitude']-new_df['target_lat'])*(new_df['latitude']-new_df['target_lat']) + (new_df['longitude']-new_df['target_lon'])*(new_df['longitude']-new_df['target_lon'])
        new_df['rank'] = new_df.groupby(['target_lat', 'target_lon', 'period'])['distance'].rank(method="first",ascending=True);
        
        samples_df = new_df[(new_df['type'] == 's') & (new_df['rank'] < 2)];
        samples_df = samples_df.groupby(['density', 'latitude', 'longitude', 'period', 'type']).first().reset_index();
        controls_df = new_df[new_df['type'] == 'c'];

        dataframe = pd.concat([samples_df, controls_df]);

        if parameters["save_processed_targets"]:
            print("Saving sites dataframe...")
            # save dataframe in processed_targets folder
            dataframe_path = os.path.join(processed_targets_dir, directory + "_dataframe.csv")
            # WRITE PROCESSED TARGET DATEFRAME
            dataframe.to_csv(dataframe_path, sep=";",quoting=csv.QUOTE_NONNUMERIC) 

            print("Saving target list...")
            # save targets in processed_targets folder 
            targets_filename = directory + "_targets"
            save_target_list_to_csv(target_list, processed_targets_dir, targets_filename)

    return target_list, dataframe, globals_dataframe

def load_processed_targets(base_path, filename, globals_type, min_lat, min_date, max_date, max_for_uninhabited):

    processed_targets_dir = os.path.join(base_path, "processed_targets")
    if not os.path.exists(processed_targets_dir):
        os.makedirs(processed_targets_dir);

    dataframe_filename = filename + '_dataframe.csv'
    filenames_in_processed_targets = [f for f in os.listdir(processed_targets_dir) if isfile(join(processed_targets_dir,f))]

    globals_dir = os.path.join(base_path, "globals");
    if not os.path.exists(globals_dir):
        os.makedirs(globals_dir);

    globals_filename = (globals_type + "_lat_" + str(min_lat) + "-" + str(max_lat) + "_date_" + str(min_date) + "-" + str(max_date) +  "_mfu_" + str(max_for_uninhabited) + ".csv").lower().replace(" ", "_");
    
    filenames_in_globals = [f for f in os.listdir(globals_dir) if isfile(join(globals_dir,f))]
    
    if dataframe_filename not in filenames_in_processed_targets or globals_dataframe_filename not in filenames_in_globals:
        return False, [], [], [];


    dataframe_path = os.path.join(processed_targets_dir, dataframe_filename)
    dataframe = pd.read_csv(dataframe_path, sep=";")

    globals_dataframe_path = os.path.join(globals_dir, globals_dataframe_filename);
    globals_dataframe = pd.read_csv(globals_dataframe_path, sep=";")

    targets_filename = os.path.join(processed_targets_dir, filename + '_targets')
    target_list = read_target_list_from_csv(targets_filename)

    return True, target_list, dataframe, globals_dataframe

def extract_dataframe(population_data, target_list, max_for_uninhabited):

        
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
    
    # print(time_dict[14200])
    latlon_length = len(lat_np)
    target_ids = []
    target_location = [];
    target_date_from = [];
    target_date_to = [];
    target_lats = [];
    target_lons = [];
    periods = []
    densities = []
    pseudo_types = []
    types = []
    latitudes = []
    longitudes = []
    dirs = []
    exacts = []

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
        for key, target in target_list.iteritems():

            ###################################################
            # Loop from date_from to date_from + dw of target #
            ###################################################
            date_from = target.date_from + population_data.time_window
            # date_to = target.date_to
            date_to = target.date_from

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
                                target_ids.append(key);
                                target_location.append(target.location);
                                target_date_from.append(target.date_from);
                                target_date_to.append(target.date_to);
                                latitudes.append(lat)
                                longitudes.append(lon)
                                target_lats.append(target.orig_lat);
                                target_lons.append(target.orig_lon);

                                periods.append(time)
                                densities.append(density)

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

                                has_site = True;
                        except KeyError:
                            print("Error occurred")
                            time = -1
                try:
                    time = next_time_dict[time]
                except KeyError:
                    time = -1
    
    new_df = pd.DataFrame({'target_id': target_ids, 'target_location': target_location, 'target_date_from': target_date_from, 'target_date_to': target_date_to, 'density': densities, 'period': periods, 'latitude': latitudes, 'longitude': longitudes, 'target_lat': target_lats, 'target_lon': target_lons, 'pseudo_type':pseudo_types, 'type':types, 'is_dir': dirs, 'is_exact': exacts})
    

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

    # INCLUSIVE
    time_mask = (time_np*time_multiplier <= max_date) & (time_np*time_multiplier >= min_date);
    
    print("Min date: " + str(min_date));
    print("Max date: " + str(max_date));
    print("Min lat: " + str(min_lat));
    print("Max lat: " + str(max_lat));


    # INCLUSIVE
    latlon_mask = (lat_np <= max_lat) & (lat_np >= min_lat);

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

def load_bin_globals_for_no_equatorials (population_data, min_lat, max_lat, min_date, max_date, max_for_uninhabited):
    df = load_all_globals_brute(population_data, min_lat, max_lat, min_date, max_date, max_for_uninhabited)
    new_df = df.loc[~df.latitude.between(-10, 20)]
    return new_df

def load_globals_for_all(population_data, target_list, date_window, min_lat, max_lat, min_date, max_date):

    lat_np = population_data.lat_array
    lon_np = population_data.lon_array
    time_np = population_data.time_array
    den_np = population_data.density_array
    time_multiplier = population_data.time_multiplier

    pseudo_target_list = []

    minimum_latitude = min_lat
    maximum_latitude = max_lat
    minimum_date = min_date
    maximum_date = max_date

    
    # Using date_from to date_from+dw
    date_window = maximum_date - minimum_date
    maximum_date = minimum_date

    # ######################################
    # # Getting the minimum/maximum values #
    # ######################################
    # for target in target_list:
    #     if target.orig_lat > maximum_latitude:
    #         maximum_latitude = target.orig_lat
    #     if target.orig_lat < minimum_latitude:
    #         minimum_latitude = target.orig_lat
    #     if target.date_from > maximum_date:
    #         maximum_date = target.date_from
    #     if target.date_from < minimum_date:
    #         minimum_date = target.date_from
    # minimum_latitude = int(round(minimum_latitude))
    # maximum_latitude = int(round(maximum_latitude)) 
    # maximum_date=maximum_date+date_window
    print('minimum date=', minimum_date)
    print('maximum date=',maximum_date)
    print('minimum latitude=', minimum_latitude)
    print('maximum latitude=',maximum_latitude)
    lat = (minimum_latitude+maximum_latitude)/2
    time = 0
    temp_target = Target(lat,0,maximum_latitude,-2, minimum_latitude,2,"location",maximum_date, minimum_date,"country","is_direct","calibrated","kind","figurative","source","is_controversial","age_estimation", 0)

    df = extract_dataframe(population_data, [[temp_target]], date_window)
    return df

def load_bin_globals_for_all(population_data, target_list, date_window, min_lat, max_lat, min_date, max_date):

    df = load_globals_for_all(population_data, target_list, date_window, min_lat, max_lat, min_date, max_date)
    df = df[df.type=='c']
    del df['location']
    del df['cluster_id']
    del df['pseudo_type']
    del df['type']
    del df['contribution']
    del df['is_dir'];
    del df['is_exact'];

    return df

def load_bin_globals_for_australia(population_data, min_date, max_date, max_for_uninhabited):
    df = load_all_globals_brute(population_data, -40, -11, min_date, max_date, max_for_uninhabited)
    new_df = df.loc[(df.latitude.between(-39.16,-11.17)) & (df.longitude.between(112.14,154.86))]
    print(new_df)
    return new_df

def load_bin_globals_for_australia_old(population_data, min_date, max_date, max_for_uninhabited):
    # maximum_latitude = -500;
    # minimum_latitude = 500
    minimum_date = min_date
    maximum_date = max_date
    
    # Using date_from to date_from+dw
    date_window = maximum_date - minimum_date
    maximum_date = minimum_date

    ######################################
    # Getting the minimum/maximum values #
    ######################################
    # for target in target_list:
    #     if target.orig_lat > maximum_latitude:
    #         maximum_latitude = target.orig_lat
    #     if target.orig_lat < minimum_latitude:
    #         minimum_latitude = target.orig_lat
    #     if target.date_from > maximum_date:
    #         maximum_date = target.date_from
    #     if target.date_from < minimum_date:
    #         minimum_date = target.date_from

    latitude = -25.27
    lat_nw = -11.17
    lat_se = -39.16
    longitude = 133.77
    lon_nw = 112.14
    lon_se = 154.86
    time = minimum_date
    temp_target = Target(latitude,longitude,lat_nw,lon_nw,lat_se,lon_se,"location",maximum_date, time,"country","is_direct","calibrated","kind","figurative","source","is_controversial","age_estimation", 0)

    df = extract_dataframe(population_data, {'au':temp_target}, max_for_uninhabited)
    df = df[df.type=='s']
    del df['target_location']
    del df['target_date_from']
    del df['target_date_to']
    del df['target_id']
    del df['type']
    del df['is_dir'];
    del df['is_exact'];
    return df


def load_bin_globals_for_francespain(population_data, min_date, max_date, max_for_uninhabited):
    df = load_all_globals_brute(population_data, 34, 60, min_date, max_date, max_for_uninhabited)
    new_df = df.loc[(df.latitude.between(35,50.84)) & ((df.longitude.between(0, 7)) | (df.longitude.between(350, 360)))]
    return new_df

def load_bin_globals_for_francespain_old(population_data, min_date, max_date, max_for_uninhabited):
    maximum_latitude = -500;
    minimum_latitude = 500
    minimum_date = min_date
    maximum_date = max_date


    # Using date_from to date_from+dw
    date_window = maximum_date - minimum_date
    maximum_date = minimum_date

    # ######################################
    # # Getting the minimum/maximum values #
    # ######################################
    # for target in target_list:
    #     if target.orig_lat > maximum_latitude:
    #         maximum_latitude = target.orig_lat
    #     if target.orig_lat < minimum_latitude:
    #         minimum_latitude = target.orig_lat
    #     if target.date_from > maximum_date:
    #         maximum_date = target.date_from
    #     if target.date_from < minimum_date:
    #         minimum_date = target.date_from

    latitude = 42
    lat_nw = 50.84
    lat_se = 35
    longitude = 0
    lon_nw = 350
    lon_se = 7
    time = minimum_date
    target_east = Target(latitude,longitude,lat_nw,lon_nw,lat_se,359,"east_location",maximum_date, time,"country","is_direct","calibrated","kind","figurative","source","is_controversial","age_estimation", 0)
    target_west = Target(latitude,longitude,lat_nw,0,lat_se,lon_se,"west_location",maximum_date, time, "country","is_direct","calibrated","kind","figurative","source","is_controversial","age_estimation", 0)

    df = extract_dataframe(population_data, {'fr1':target_east, 'fr2': target_west}, max_for_uninhabited)

    df = df[df.type=='s']
    del df['target_location']
    del df['target_date_from']
    del df['target_date_to']
    del df['target_id']
    del df['type']
    del df['is_dir'];
    del df['is_exact'];
    return df

def load_bin_globals_for_trial_latitudes(population_data, target_list, date_window, min_lat, max_lat, min_date, max_date):
    df = load_bin_globals_for_all(population_data, target_list, date_window, min_lat, max_lat, min_date, max_date)
    new_df = df.loc[~df.latitude.between(-10, 40)]
    return new_df

def load_bin_globals_no_empty_lat(population_data, target_list, date_window, min_lat, max_lat, min_date, max_date):
    df = load_bin_globals_for_all(population_data, target_list, date_window, min_lat, max_lat, min_date, max_date)
    new_df = df.loc[~df.latitude.between(0, 20)]
    return new_df

def load_bin_globals_for_trial_latitudes2(population_data, target_list, date_window, min_lat, max_lat, min_date, max_date):
    df = load_bin_globals_no_empty_lat(population_data, target_list, date_window, min_lat, max_lat, min_date, max_date)
    new_df = df.loc[~df.latitude.between(-30, -10)]
    new_df = df.loc[~df.latitude.between(40, 50)]
    return new_df

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

    target_list={}

    with open(filename +".csv",'rb') as target_file:
        my_reader=csv.reader(target_file,delimiter=',')
        next(my_reader,None) 
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


            target_id = location + "(lat: " + str(lat) + ", lon: " + str(lon) + ", date_from: " + str(date_from) + ")"

            target=Target(target_id, lat,lon,lat_nw,lon_nw,lat_se,lon_se,location,date_from, date_to, country,is_direct,calibrated,kind,figurative,source,is_controversial, age_est)

            target_list[target_id] = target
    target_file.close
    return target_list


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

    
def filter_targets_for_not_direct(target_list):
    print("length of original list=",len(target_list))

    filtered_list = {};
    for key, target in target_list.iteritems():
        if target.is_direct == 'Yes':
            filtered_list[key] = target;

    print("length of filtered list=",len(filtered_list))
    
    return filtered_list

def filter_targets_for_not_exact_age(target_list):
    print("length of original list=",len(target_list))

    filtered_list = {};
    for key, target in target_list.iteritems():
        if target.age_estimation.lower()=='exact age':
            filtered_list[key] = target;

    print("length of filtered list=",len(filtered_list))

    return filtered_list

def filter_targets_for_not_figurative(target_list):
    print("length of original list=",len(target_list))
    

    filtered_list = {};
    for key, target in target_list.iteritems():
        if target.figurative == 'Yes':
            filtered_list[key] = target;

    print("length of filtered list=",len(filtered_list))

    return filtered_list

def filter_targets_for_not_controversial(target_list):
    print("length of original list=",len(target_list))
    
    filtered_list = {};
    for key, target in target_list.iteritems():
        if target.is_controversial == 'No':
            filtered_list[key] = target;

    print("length of filtered list=",len(filtered_list))

    return filtered_list

def filter_targets_for_date(target_list, minimum_date, maximum_date):
    filtered_list = {};
    for key, target in target_list.iteritems():
        if target.date_from >= minimum_date and target.date_from <= maximum_date:
            filtered_list[key] = target

    return filtered_list

def filter_targets_for_abs_latitude(target_list, minimum_lat, maximum_lat, filters_applied):
    filtered_targets = []
    for target in target_list:
        abs_lat = abs(target.orig_lat)
        if abs_lat >= minimum_lat and abs_lat <= maximum_lat:
            filtered_targets.append(target)

    filters_applied += " Only targets within" + str(minimum_lat) + " - " + str(maximum_lat) + " absolute latitude"

    return filtered_targets, filters_applied

def filter_targets_for_latitude(target_list, minimum_lat, maximum_lat):
    filtered_list = {};
    for key, target in target_list.iteritems():
        lat = target.orig_lat
        if lat >= minimum_lat and lat <= maximum_lat:
            filtered_list[key] = target;
    return filtered_list



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

def generate_merged_dataframe(base_path, directory, dataframe, globals_dataframe, save_processed_targets):

    temp_globals_df = globals_dataframe.copy();
    temp_samples_df = dataframe.copy();
    temp_samples_df = temp_samples_df[temp_samples_df.type == 's'];

     # DELETE COLUMNS
    temp_samples_df.drop(["target_id","samples_growth_coefficient", "distance", "is_dir", "is_exact", "pseudo_type", "rank", "target_date_from", "target_date_to", "target_lat", "target_location", "target_lon","type"], axis=1, inplace = True)

    # is_sample
    temp_globals_df['is_sample'] = 0;    
    temp_samples_df['is_sample'] = 1;

    # MERGE
    to_concat = [temp_globals_df, temp_samples_df]
    merged_df = pd.concat(to_concat);
    
    # abs lat
    merged_df['abs_latitude'] = merged_df['latitude'].abs();

    # binned_period
    merged_df = create_binned_column(merged_df, 'binned_period', 'period', 10000);

    # binned_period
    merged_df = create_binned_column(merged_df, 'binned_latitude', 'latitude', 10);

    if save_processed_targets:
        processed_targets_dir = os.path.join(base_path, "processed_targets")
        merged_df_filename = os.path.join(processed_targets_dir, directory + "_merged_df.csv") 
        if not os.path.exists(processed_targets_dir):
            os.makedirs(processed_targets_dir)
        merged_df.to_csv(merged_df_filename, sep=";")
        print("Saved merged dataframe: " + merged_df_filename);
        
    return merged_df