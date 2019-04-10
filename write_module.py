import numpy as np

def write_label(a_file, label):
    for i in range(0, len(label)+8):
        a_file.write('*')
    a_file.write('\n')
    a_file.write(label + "\n")
    for i in range(0, len(label)+8):
        a_file.write('*')
    a_file.write('\n')

def write_headers(a_file, headers, delimiter):
    for i in range(0, len(headers)):
        a_file.write(str(headers[i]))
        if i != len(headers)-1:
            a_file.write(delimiter)
    a_file.write('\n')
    
def write_confounder_analysis_table(a_file, label, dictionary, keys, or_MHs, MH_stats, MH_stat_ps):

    headers = [""];
    small_headers = ["Bin values"]
    for key in keys:
        headers.append(key);
        headers.append("");
        headers.append("");
        headers.append("");
        headers.append("");
        headers.append("");
        headers.append("");
        headers.append("");
        small_headers.append("Samples");
        small_headers.append("Globals");
        small_headers.append("Controls");
        small_headers.append("Odds Ratio");
        small_headers.append("top ORMH");
        small_headers.append("bottom ORMH");
        small_headers.append("top test MH");
        small_headers.append("bottom test MH");
    headers.append("OR MH");
    headers.append("MH Statistics");
    headers.append("MH Statistics p-value");

    write_label(a_file, label);
    write_headers(a_file, headers, ";");
    write_headers(a_file, small_headers, ";");

    bin_array = dictionary[keys[0]][0];

    for i in range(0, len(bin_array)):
        a_file.write(str(bin_array[i]) + ";");
        for j in range(0, len(keys)):
            key = keys[j];
            values = dictionary[key];
            a_file.write(str(values[1][i]) + ";");
            a_file.write(str(values[2][i]) + ";");
            a_file.write(str(values[3][i]) + ";");
            a_file.write(str(values[4][i]) + ";");
            a_file.write(str(values[5][i]) + ";");
            a_file.write(str(values[6][i]) + ";");
            a_file.write(str(values[7][i]) + ";");
            a_file.write(str(values[8][i]) + ";");
        a_file.write(str(or_MHs[i]) + ";");
        a_file.write(str(MH_stats[i]) + ";");
        a_file.write(str(MH_stat_ps[i]) + "\n");


def write_random_weighting_table(a_file, label, column_keys, column_weights, column_weights_normalized, bin_array, odds_ratios):
    write_label(a_file, "Keys and Weights For " + label);
    a_file.write("Number of Keys: " + str(len(column_keys)) + "\n");

    column_weights_normalized = np.round(column_weights_normalized, 6)

    column_weights = column_weights.tolist();
    column_weights_normalized = column_weights_normalized.tolist();

    column_keys.insert(0, "KEYS");
    column_weights.insert(0, "WEIGHTS");
    column_weights_normalized.insert(0, "NORMALIZED_WEIGHTS");

    write_headers(a_file, column_keys, ";");
    write_headers(a_file, column_weights, ";");
    write_headers(a_file, column_weights_normalized, ";");

    write_label(a_file, "Bins and Odds Ratios For " + label);
    headers = ['bin', 'odds_ratio'];
    write_headers(a_file, headers, ";");

    for i in range(0, len(bin_array)):
        a_file.write(str(bin_array[i]) + ";");
        a_file.write(str(odds_ratios[i]) + "\n");



def write_table(a_file, label, headers, values, delimiter):
    write_label(a_file, label)
    write_headers(a_file, headers, delimiter)
    
    rows = len(values)
    columns = len(values[0])
    for j in range(0, columns):
        for i in range(0, rows):
            a_file.write(str(values[i][j]))
            if i != rows-1:
                a_file.write(delimiter)
        a_file.write('\n')

def write_target(a_file, target, date_window):
    a_file.write('"'+target.location+'"'+';')
    a_file.write(str(int(target.cluster_id))+';')
    a_file.write(str(float(target.orig_lat))+';')
    a_file.write(str(float(target.orig_lon))+';')
    a_file.write(str(int(target.date_from))+';')
    a_file.write(str(int(target.date_to)))
    a_file.write('\n')

def write_target_table(a_file, targetList, date_window):
    a_file.write('Location;ClusterID;Latitude;Longitude;DateFrom;DateTo\n')
    for target in targetList:
        write_target(a_file,target,date_window)
    # a_file.write('\n')
    
def write_information(a_file, labels, values, delimiter):
    information_str=""
    for i in range(0, len(labels)):
        information_str=information_str+str(labels[i]) + "=" + str(values[i])+" "
    a_file.write(information_str)
    a_file.write('\n')

def write_parameters(a_file, parameters):
    for key, value in parameters.iteritems():
        if isinstance(value, basestring):
            a_file.write(key + ": " + value);
        else:
            a_file.write(key + ": " + str(value));
        a_file.write("\n");

def write_cluster_table(a_file, dataframe, parameters):
    date_window = parameters['date_window'];
    date_lag = parameters['date_lag'];
    cluster_headers=["Name of Site", "Latitude", "Longitude", "TargetDateFrom", "TargetDateTo", "AnalysisDateFrom", "AnalysisDateTo", "Direct", "Exact", "Population density", "Controls population density", "Growth Coefficient"]
    write_headers(a_file, cluster_headers, ';')

    ################################
    # Group by cluster id and type #
    # Modified to report medians instead of means
    ################################
    # - aggregate by mean
    #new_data = dataframe.groupby(['cluster_id', 'pseudo_type']).mean().reset_index()
    new_data = dataframe.groupby(['cluster_id', 'pseudo_type']).mean().reset_index()
    temp_types = new_data['pseudo_type'].values
    cluster_ids = new_data.cluster_id.unique();
    # print("ClusterID")
    # print(cluster_ids)

    for cluster_id in cluster_ids:
        cluster_df = dataframe[dataframe.cluster_id == cluster_id] 
        # print(cluster_df)
        location = cluster_df['target_location'].values[0];
        latitude = cluster_df['target_lat'].values[0];
        longitude = cluster_df['target_lon'].values[0];
        date_from = cluster_df['target_date_from'].values[0];
        date_to = cluster_df['target_date_to'].values[0];
        date_from_analysis = date_from + date_lag;
        date_to_analysis = date_from + date_lag + date_window;


        direct = 'not direct';
        if cluster_df['is_dir'].values[0]:
            direct = 'direct'

        exact = 'not exact';
        if cluster_df['is_exact'].values[0]:
            exact = 'exact'

        sample_mean = cluster_df[cluster_df.type == 's']['density'].values[0];
        controls_mean = cluster_df[cluster_df.type == 'c']['density'].values[0];

        sample_growth_coefficient = cluster_df['samples_growth_coefficient'].values[0]

        a_file.write("\"" + str(location) + "\";")
        # a_file.write(str(cluster_ids[i]) + ";")
        a_file.write(str(latitude) + ";")
        a_file.write(str(longitude) + ";")
        a_file.write(str(date_from) + ";")
        a_file.write(str(date_to) + ";")
        a_file.write(str(date_from_analysis) + ";")
        a_file.write(str(date_to_analysis) + ";")
        a_file.write(str(direct) + ";")
        a_file.write(str(exact) + ";")
        a_file.write(str(sample_mean) + ";")
        a_file.write(str(controls_mean) + ";")
        a_file.write(str(sample_growth_coefficient))
        a_file.write("\n")
        
def write_bin_table(a_file, bin_values_df, minimum_globals):

    write_label(a_file, "Distribution of values for samples and globals - Minimum_globals="+str(minimum_globals))

    columns = ['Bin value', 'Samples', 'Globals', 'Detection Frequency', 'Relative Frequency of Sites', 'Relative Frequency of Globals']
    write_headers(a_file,columns,";")
    bin_array = bin_values_df['bin_array'].values
    sample_counts = bin_values_df['sample_counts'].values
    global_counts = bin_values_df['global_counts'].values
    likelihood_ratios = bin_values_df['likelihood_ratios'].values
    p_samples = bin_values_df['p_samples'].values
    p_globals = bin_values_df['p_globals'].values

    for i in range(0, len(bin_array)):
        a_file.write(str(bin_array[i]) + ';')
        a_file.write(str(sample_counts[i]) + ';')
        a_file.write(str(global_counts[i]) + ';')
        a_file.write(str(likelihood_ratios[i]) + ';')
        a_file.write(str(p_samples[i]) + ';')
        a_file.write(str(p_globals[i]) + ";")
        a_file.write("\n")

def write_analysis(f2, stat_dictionary, min_p):

    write_label(f2, "Statistics")

    f2.write('Total sites: '+str(stat_dictionary['total_samples'])+'\n')
    f2.write('Total globals: '+str(stat_dictionary['total_globals'])+'\n\n')

    f2.write('Median density for sites: '+str(stat_dictionary['median_samples'])+'\n')
    f2.write('Median density for globals: '+str(stat_dictionary['median_globals'])+'\n\n')

    f2.write('Mean density for sites: '+str(stat_dictionary['mean_samples'])+'\n')
    f2.write('Mean density for globals: '+str(stat_dictionary['mean_globals'])+'\n\n')

    f2.write('Standard deviation of density for sites: '+str(stat_dictionary['std_samples'])+'\n')
    f2.write('Standard deviation of density for globals: '+str(stat_dictionary['std_globals'])+'\n\n')

    ks_p = stat_dictionary['ks_p']
    ks_d = stat_dictionary['ks_d']

    f2.write('K S2 test for p_samples vs p_globals:\n');
    f2.write('    ks_d: ' + str(float(ks_d)) + '\n')
    f2.write('    ks_p: ' + str(float(ks_p)) + '\n')
    if ks_p < min_p:
        f2.write('    The two distribitions are significantly different p<0.001'+'\n')

def write_lat_bin_table(a_file, bin_array, table, label):
    write_label(a_file, "Distribution of values for " + label + " per latitude range")

    # Write bins
    a_file.write(";;")
    for i in range(0, len(bin_array)):
        a_file.write(str(bin_array[i]))
        if i != len(bin_array) - 1:
            a_file.write(";")
    a_file.write("\n")

    min_latitudes = [x for x in range(-90, 90, 10)]
    max_latitudes = [x for x in range(-80, 91, 10)]
    lat_arr_length = len(min_latitudes)

    for i in range(0, lat_arr_length):
        a_file.write(str(min_latitudes[i]) + ";" + str(max_latitudes[i]) + ";")
        values = table[i]
        for j in range(0, len(values)):
            a_file.write(str(values[j]))
            if j != len(values) - 1:
                a_file.write(";")
        a_file.write("\n")

def write_validation_file(base_path, title, score_linear_array, score_threshold_array, slopes, intercepts, parameter_1, parameter_2, parameter_3):
    valid_file = open(base_path+ "/" + title+".csv", "w")
    print(base_path)
    headers = [title+"_score_linear_array", title+"_score_threshold_array", "slopes", "intercepts", "parameter_1", "parameter_2", "parameter_3"]
    write_headers(valid_file, headers, ";")

    for i in range(0, len(score_linear_array)):
        valid_file.write(str(score_linear_array[i]) + ";")
        valid_file.write(str(score_threshold_array[i]) + ";")
        valid_file.write(str(slopes[i]) + ";")
        valid_file.write(str(intercepts[i]) + ";")
        valid_file.write(str(parameter_1[i]) + ";")
        valid_file.write(str(parameter_2[i]) + ";")
        valid_file.write(str(parameter_3[i]) + "\n")

    valid_file.close()

