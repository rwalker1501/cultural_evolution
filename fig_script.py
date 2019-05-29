#!/usr/bin/env python
# -*- coding: utf-8 -*-  
import os
import main_module as mm

#--------------------------#
# run_experiment parameters
#--------------------------#
# results_path: directory where results folder will be created
# target_list_file: filename of target list in tests folder (no extension)
# output_directory: directory of results (will be in results folder)
# population_data_name="Eriksson"
# globals="All"
# date_window=10000
# user_max_for_uninhabited=-1
# clustering_on = False
# critical_distance=1
# filter_date_before=-1
# filter_not_direct=False
# filter_not_exact=False;
# filter_not_figurative=False
# filter_not_controversial = False
# perform_cross_validation=False
# number_of_kfolds = 100
# minimum_likelihood_ratio = 0.0000001
# min_date_window=0
# critical_time=10000
# filter_min_date=-1
# filter_max_date=-1
# filter_min_lat=-1
# filter_max_lat=-1
# processed_targets = False
# is_confounder_analysis = False

#-----------------------------------#
# plot_population_by_time parameters
#-----------------------------------#
# population_data_name: name of population data
# time: time in BP


#--------#
# SCRIPT #
#--------#
base_path = os.getcwd()


#mm.run_experiment(base_path, "rpv12 no equatorials","rpv12 no equatorials NEWHR",globals_type="no equatorials",date_window=24,low_res=False);
#mm.run_experiment(base_path, "rpv12 no equatorials exact direct", "rpv12 no equatorials exact direct NEW",globals_type="No equatorials",date_window=24,low_res=False);
#mm.run_experiment(base_path, "rpv12 no equatorials", "rpv12 no equatorials",globals_type="No equatorials",date_window=24,low_res=False);
#mm.run_experiment(base_path, "rpv12 no equatorials", "rpv12 no equatorials NEWHR",globals_type="No equatorials",date_window=24,low_res=False);
#mm.run_experiment(base_path, "rpv12 no equatorials", "rpv12_timmermann", population_data_name="Timmermann",globals_type="No equatorials",date_window=999,low_res=True);
# Figure 1A: map of all sites and controls. Use a different colour (green) for sites where we have "direct" and "exact age" compared to all other sites
# Fig 2A: pGraph for world, Eriksson
# Fig 2B: OR graph for world, Eriksson
#
#mm.run_experiment(base_path, "rpv12", "eriksson full",date_window=24,low_res=False);
# Here we need to add experiment where we use only consensus densities

# # # Fig 1B: map of population densities with Eriksson model at 25000 BP
#mm.plot_population_by_time("Eriksson", 25000);


# # FIG2E: OR graph for Australia, Eriksson
#mm.run_experiment(base_path, "rpv12_fr_sp", "rpv12_fr_sp", globals_type="France and Spain", date_window=24,low_res=False)
#mm.run_experiment(base_path, "rpv12_au_exact_direct", "rpv12_au_exact_direct", globals_type="Australia", date_window=24,low_res=False)
#mm.run_experiment(base_path, "rpv12 no equatorials exact direct", "rpv12_lat60-40", filter_min_lat=40, filter_max_lat=59.9999,date_window=24,low_res=False)
#mm.run_experiment(base_path, "rpv12 no equatorials exact direct", "rpv12_lat20-40", filter_min_lat=20, filter_max_lat=39.9999,date_window=24,low_res=False)
#mm.run_experiment(base_path, "rpv12 no equatorials exact direct", "rpv12_lat-10--30", filter_min_lat=-29.9999, filter_max_lat=-10,date_window=24,low_res=False)
#mm.run_experiment(base_path, "rpv12 no equatorials exact direct", "rpv12_lat-30--50", filter_min_lat=-49.9999, filter_max_lat=-30,date_window=24,low_res=False)
mm.run_experiment(base_path, "rpv12 no equatorials", "rpv12_dat 10k-45k", globals_type="No equatorials",filter_min_date=10000, filter_max_date=45000,date_window=24,low_res=False)
mm.run_experiment(base_path, "rpv12 no equatorials", "rpv12_dat<10k", filter_min_date=0, filter_max_date=9999,globals_type="No equatorials", date_window=24,low_res=False)



# FIG2G: pGraph for world, Timmermann (as in Fig 2A)
# FIG2H: OR graph for world, Timmermann (as in Fig2B)
#mm.run_experiment(base_path, "rpv12", "rpv12_timmermann", population_data_name="Timmermann",date_window=999);
# mm.run_experiment(base_path, "rpv12", "rpv12_timmermann",  processed_targets = True,  population_data_name="Timmermann",date_window=999);
# mm.run_experiment(base_path, "rpv12_test", "rpv12_tim_test", population_data_name="Timmermann",date_window=999);

# SM1 Map showing location of regions with high mean population density in period 5000-50000BP
#mm.plot_min_densities_in_time_range("Eriksson", 50000, 5000, 5000)
# SM2: OR graph for each date band - include empty lats, error bars, no trend lines.
# SM3: OR graph for each latitude band -  include empty lats, error bars, no trend lines.
#mm.run_experiment(base_path, "rpv12", "rpv12_ca", is_confounder_analysis=True, date_window=24,date_window=24,low_res=False);

# SM4: OR graphs for France/Spain with different clustering parameters (distance = 100, 500, 1000 km, time = 1000, 5000, 10000 years)
#mm.run_experiment(base_path, "rpv10_fr_sp", "rpv10_fr_sp_cd100_ct1000", the_globals="France and Spain", clustering_on=True, critical_distance=100, critical_time=1000)
#mm.run_experiment(base_path, "rpv10_fr_sp", "rpv10_fr_sp_cd100_ct5000", the_globals="France and Spain", clustering_on=True, critical_distance=100, critical_time=5000)
#mm.run_experiment(base_path, "rpv10_fr_sp", "rpv10_fr_sp_cd100_ct10000", the_globals="France and Spain", clustering_on=True, critical_distance=100, critical_time=10000)

#mm.run_experiment(base_path, "rpv10_fr_sp", "rpv10_fr_sp_cd500_ct1000", the_globals="France and Spain", clustering_on=True, critical_distance=500, critical_time=1000)
#mm.run_experiment(base_path, "rpv10_fr_sp", "rpv10_fr_sp_cd500_ct5000", the_globals="France and Spain", clustering_on=True, critical_distance=500, critical_time=5000)
#mm.run_experiment(base_path, "rpv10_fr_sp", "rpv10_fr_sp_cd500_ct10000", the_globals="France and Spain", clustering_on=True, critical_distance=500, critical_time=10000)

#mm.run_experiment(base_path, "rpv10_fr_sp", "rpv10_fr_sp_cd1000_ct1000", the_globals="France and Spain", clustering_on=True, critical_distance=1000, critical_time=1000)
#mm.run_experiment(base_path, "rpv10_fr_sp", "rpv10_fr_sp_cd1000_ct5000", the_globals="France and Spain", clustering_on=True, critical_distance=1000, critical_time=5000)
#mm.run_experiment(base_path, "rpv10_fr_sp", "rpv10_fr_sp_cd1000_ct10000", the_globals="France and Spain", clustering_on=True, critical_distance=1000, critical_time=10000)

# OR graphs for Australia with different clustering parameters (distance = 100, 500, 1000 km, time = 1000, 5000, 10000 years)
#mm.run_experiment(base_path, "rpv10_au", "rpv10_au_cd100_ct1000", the_globals="Australia", clustering_on=True, critical_distance=100, critical_time=1000)
#mm.run_experiment(base_path, "rpv10_au", "rpv10_au_cd100_ct5000", the_globals="Australia", clustering_on=True, critical_distance=100, critical_time=5000)
#mm.run_experiment(base_path, "rpv10_au", "rpv10_au_cd100_ct10000", the_globals="Australia", clustering_on=True, critical_distance=100, critical_time=10000)

#mm.run_experiment(base_path, "rpv10_au", "rpv10_au_cd500_ct1000", the_globals="Australia", clustering_on=True, critical_distance=500, critical_time=1000)
#mm.run_experiment(base_path, "rpv10_au", "rpv10_au_cd500_ct5000", the_globals="Australia", clustering_on=True, critical_distance=500, critical_time=5000)
#mm.run_experiment(base_path, "rpv10_au", "rpv10_au_cd500_ct10000", the_globals="Australia", clustering_on=True, critical_distance=500, critical_time=10000)

#mm.run_experiment(base_path, "rpv10_au", "rpv10_au_cd1000_ct1000", the_globals="Australia", clustering_on=True, critical_distance=1000, critical_time=1000)
#mm.run_experiment(base_path, "rpv10_au", "rpv10_au_cd1000_ct5000", the_globals="Australia", clustering_on=True, critical_distance=1000, critical_time=5000)
#mm.run_experiment(base_path, "rpv10_au", "rpv10_au_cd1000_ct10000", the_globals="Australia", clustering_on=True, critical_distance=1000, critical_time=10000)

#mm.run_experiment(base_path, "rpv11_no_big_gaps", "rpv11_no_big_gaps")
#mm.run_experiment(base_path, "rpv11_no_low_pop_australia", "rpv11_no_low_pop_australia")

#mm.run_experiment(base_path, "rpv11_no_big_gaps", "rpv11_no_big_gaps_timmermann", population_data_name="Timmermann")
#mm.run_experiment(base_path, "rpv11_no_low_pop_australia", "rpv11_no_low_pop_australia_timmermann", population_data_name="Timmermann")