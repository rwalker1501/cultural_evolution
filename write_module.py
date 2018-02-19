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
        a_file.write(headers[i])
        if i != len(headers)-1:
            a_file.write(delimiter)
    a_file.write('\n')
    

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
    a_file.write(str(int(target.date_from+date_window)))
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


def write_cluster_table(a_file, dataframe, growth_coefficients):
    cluster_headers=["First Location in cluster", "Cluster", "Latitude", "Longitude", "Sample/Control", "Mean Population Density over Period", "N samples over period", "growthCoefficient"]
    write_headers(a_file, cluster_headers, ';')

    ################################
    # Group by cluster id and type #
    ################################
    # - aggregate by mean
    new_data = dataframe.groupby(['cluster_id', 'pseudo_type']).mean().reset_index()

    temp_types = new_data['pseudo_type'].values
    cluster_ids = new_data['cluster_id'].values
    means = new_data['density'].values
    # aggregate by sum
    counts = dataframe.groupby(['cluster_id', 'pseudo_type']).sum()['contribution'].values

    for i in range(0, len(cluster_ids)):
        # get first target
        location = dataframe[dataframe.cluster_id == cluster_ids[i]].head(1)['location'].values[0]
        longitude = dataframe[dataframe.cluster_id == cluster_ids[i]].head(1)['target_lon'].values[0]
        latitude = dataframe[dataframe.cluster_id == cluster_ids[i]].head(1)['target_lat'].values[0]
        a_file.write("\"" + str(location) + "\";")
        a_file.write(str(cluster_ids[i]) + ";")
        a_file.write(str(latitude) + ";")
        a_file.write(str(longitude) + ";")

        if temp_types[i]=='a':
            a_file.write("Sample;")
        else:
            a_file.write('Control;')
        a_file.write(str(means[i]) + ";")
        a_file.write(str(counts[i]) + ";")
        a_file.write(str(growth_coefficients[i]))
        a_file.write("\n")
# =============================================================================
#         if i != 0 and i%2==1:
#             a_file.write("\n")   #I was getting an unexplained CR on every two lines - think this is responsible.
# =============================================================================

def write_bin_table(a_file, bin_array, sample_counts, control_counts, likelihood_ratios, p_samples, p_controls, p_likelihood_ratios,minimum_controls):

    print 'minimum controls in write_bin_table=',str(minimum_controls)
    write_label(a_file, "Distribution of values for samples and controls - Minimum-controls="+str(minimum_controls))

    columns = ['Bin value', 'Samples', 'Controls', 'Likelihood Ratios', 'pSamples', 'pControls', 'p Likelihood Ratios']
    write_headers(a_file,columns,";")

    for i in range(0, len(bin_array)):
        a_file.write(str(bin_array[i]) + ';')
        a_file.write(str(sample_counts[i]) + ';')
        a_file.write(str(control_counts[i]) + ';')
        a_file.write(str(likelihood_ratios[i]) + ';')
        a_file.write(str(p_samples[i]) + ';')
        a_file.write(str(p_controls[i]) + ";")
        a_file.write(str(p_likelihood_ratios[i]))
        a_file.write("\n")


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

