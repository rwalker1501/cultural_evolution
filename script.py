import os
import main_module as mm

#--------------------------#
# run_experiment parameters
#--------------------------#
# results_path: directory where results folder will be created
# target_list_file: filename of target list in tests folder (no extension)
# output_directory: directory of results (will be in results folder)
# population_data_name="Eriksson"
# controls="All"
# date_window=10000
# user_max_for_uninhabited=-1
# clustering_on = False
# critical_distance=1
# filter_date_before=-1
# filter_not_direct=False
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
desktop_path = "/Users/richard/Dropbox (HappyFamily)/Richard global sync/EPFL documents/Documents 2017/Richard articles and papers/Basic informational constraint/Empirical studies" #saves results at this address
#base_path = os.getcwd()
#base_path="/Users/rwalker/Dropbox (HappyFamily)/Richard global sync/EPFL documents/Documents 2017/Richard articles and papers/Basic informational constraint/Empirical studies" #saves results at this address
base_path="/Users/richard/Dropbox (HappyFamily)/Richard global sync/EPFL documents/Documents 2017/Richard articles and papers/Basic informational constraint/Empirical studies" #saves results at this address
############
# ERIKSSON #
############

# =============================================================================
#mm.run_experiment(base_path, "rockpaintings v8a", "full_a_e e001", controls="All")
# =============================================================================
# =============================================================================
# mm.run_experiment(base_path, "trial_latitudes", "trial_latitudes_e", controls="Trial Latitudes")
# mm.run_experiment(base_path, "trial_latitudes2", "trial_latitudes2_e", controls="Trial Latitudes 2")
# 
# 
# mm.run_experiment(base_path, "rockpaintings v8", "full_a_e_d0-10k", controls="No Empty Lats", filter_min_date = 0, filter_max_date = 10000)
# mm.run_experiment(base_path, "rockpaintings v8", "full_a_e_d10k-20k", controls="No Empty Lats", filter_min_date = 10000, filter_max_date = 20000)
# mm.run_experiment(base_path, "rockpaintings v8", "full_a_e_d20k-30k", controls="No Empty Lats", filter_min_date = 20000, filter_max_date = 30000)
# mm.run_experiment(base_path, "rockpaintings v8", "full_a_e_d30k-40k", controls="No Empty Lats", filter_min_date = 30000, filter_max_date = 40000)
# mm.run_experiment(base_path, "rockpaintings v8", "full_a_e_d40k-50k", controls="No Empty Lats", filter_min_date = 40000, filter_max_date = 50000)
# mm.run_experiment(base_path, "rockpaintings v8", "full_a_e_d10k-30k", controls="No Empty Lats", filter_min_date = 10000, filter_max_date = 30000)
# 
# 
# mm.run_experiment(base_path, "rockpaintings v8a", "full_a_e_dir", controls="No Empty Lats", filter_not_direct=True)
# mm.run_experiment(base_path, "rockpaintings v8a", "full_a_e_fig", controls="No Empty Lats", filter_not_figurative=True)
# 
# mm.run_experiment(base_path, "rockpaintings v8a", "full_a_e_cd500_ct10000", controls="No Empty Lats", clustering_on=True, critical_distance=500, critical_time=10000)
# 
# 
# ##############
# # Timmermann #
# ##############
# mm.run_experiment(base_path, "rockpaintings v8a", "full_a_t", population_data_name="Timmermann", controls="No Empty Lats")
# mm.run_experiment(base_path, "trial_latitudes", "trial_latitudes_t", population_data_name="Timmermann", controls="Trial Latitudes")
# mm.run_experiment(base_path, "trial_latitudes2", "trial_latitudes2_t", population_data_name="Timmermann", controls="Trial Latitudes 2")
# 
# 
# mm.run_experiment(base_path, "rockpaintings v8a", "full_a_t_d0-10k", population_data_name="Timmermann", controls="No Empty Lats", filter_min_date = 0, filter_max_date = 10000)
# mm.run_experiment(base_path, "rockpaintings v8a", "full_a_t_d10k-20k", population_data_name="Timmermann", controls="No Empty Lats", filter_min_date = 10000, filter_max_date = 20000)
# mm.run_experiment(base_path, "rockpaintings v8a", "full_a_t_d20k-30k", population_data_name="Timmermann", controls="No Empty Lats", filter_min_date = 20000, filter_max_date = 30000)
# mm.run_experiment(base_path, "rockpaintings v8a", "full_a_t_d30k-40k", population_data_name="Timmermann", controls="No Empty Lats", filter_min_date = 30000, filter_max_date = 40000)
# mm.run_experiment(base_path, "rockpaintings v8a", "full_a_t_d40k-50k", population_data_name="Timmermann", controls="No Empty Lats", filter_min_date = 40000, filter_max_date = 50000)
# mm.run_experiment(base_path, "rockpaintings v8a", "full_a_t_d10k-30k", population_data_name="Timmermann", controls="No Empty Lats", filter_min_date = 10000, filter_max_date = 30000)
# 
# 
# mm.run_experiment(base_path, "rockpaintings v8a", "full_a_t_dir", population_data_name="Timmermann", controls="No Empty Lats", filter_not_direct=True)
# mm.run_experiment(base_path, "rockpaintings v8a", "full_a_t_fig", population_data_name="Timmermann", controls="No Empty Lats", filter_not_figurative=True)
# 
# mm.run_experiment(base_path, "rockpaintings v8a", "full_a_t_cd500_ct10000", population_data_name="Timmermann", controls="No Empty Lats", clustering_on=True, critical_distance=500, critical_time=10000)
# =============================================================================



#mm.run_experiment(base_path, "quicktest", "test2", processed_targets=False)
#mm.run_experiment(base_path, "rockpaintings v8", "full_a_e script", processed_targets=False)

# =============================================================================
# =============================================================================
#mm.run_experiment(base_path, "rockpaintings v8a", "full_a_e_lat-40--30", processed_targets=False, filter_min_lat = -40, filter_max_lat = -30)
#mm.run_experiment(base_path, "rockpaintings v8a", "full_a_e_lat-30--20", filter_min_lat = -30, filter_max_lat = -20)
#mm.run_experiment(base_path, "rockpaintings v8a", "full_a_e_lat-20--10", filter_min_lat = -20, filter_max_lat = -10)
#mm.run_experiment(base_path, "rockpaintings v8a", "full_a_e_lat-10-0", processed_targets=False, filter_min_lat = -10, filter_max_lat = 0)
#mm.run_experiment(base_path, "rockpaintings v8a", "full_a_e_lat0-10", filter_min_lat = 0, filter_max_lat = 10)
#mm.run_experiment(base_path, "rockpaintings v8a", "full_a_e_lat10-20", filter_min_lat = 10, filter_max_lat = 20)
#mm.run_experiment(base_path, "rockpaintings v8a", "full_a_e_lat20-30", processed_targets=False, filter_min_lat = 20, filter_max_lat = 30)
#mm.run_experiment(base_path, "rockpaintings v8a", "full_a_e_lat30-40", processed_targets=False, filter_min_lat = 30, filter_max_lat = 40)
#mm.run_experiment(base_path, "rockpaintings v8a", "full_a_e_lat40-50", processed_targets=False, filter_min_lat = 40, filter_max_lat = 50)
#mm.run_experiment(base_path, "rockpaintings v8a", "full_a_e_lat50-60", processed_targets=False, filter_min_lat = 50, filter_max_lat = 60)
# # 
#mm.run_experiment(base_path, "rockpaintings v8", "full_a_e_d0-10k", filter_min_date = 0, filter_max_date = 10000)
#mm.run_experiment(base_path, "rockpaintings v8", "full_a_e_d10k-20k", processed_targets=False, filter_min_date = 10000, filter_max_date = 20000)
#mm.run_experiment(base_path, "rockpaintings v8", "full_a_e_d20k-30k", processed_targets=False, filter_min_date = 20000, filter_max_date = 30000)
#mm.run_experiment(base_path, "rockpaintings v8", "full_a_e_d30k-40k", processed_targets=False, filter_min_date = 30000, filter_max_date = 40000)
#mm.run_experiment(base_path, "rockpaintings v8", "full_a_e_d40k-50k", processed_targets=False, filter_min_date = 40000, filter_max_date = 50000)
#mm.run_experiment(base_path, "rockpaintings v8", "full_a_e_d10k-30k", processed_targets=False, filter_min_date = 10000, filter_max_date = 30000)
# # 
# # 
# =============================================================================
mm.run_experiment(base_path, "france_spain", "fs_e", controls="France and Spain", processed_targets=False)
# =============================================================================
# # # mm.run_experiment(base_path, "france_spain", "fs_e_d0-10k", controls="France and Spain", filter_min_date = 0, filter_max_date = 10000)
# # # mm.run_experiment(base_path, "france_spain", "fs_e_d10k-20k", controls="France and Spain", filter_min_date = 10000, filter_max_date = 20000)
# # # mm.run_experiment(base_path, "france_spain", "fs_e_d20k-30k", controls="France and Spain", filter_min_date = 20000, filter_max_date = 30000)
# # # mm.run_experiment(base_path, "france_spain", "fs_e_d30k-80k", controls="France and Spain", filter_min_date = 30000, filter_max_date = 80000)

mm.run_experiment(base_path, "france_spain", "fs_e_cd100_ct10000", controls="France and Spain", clustering_on=True, critical_distance=100, critical_time=10000)#mm.run_experiment(base_path, "france_spain", "fs_e_cd250_ct10000", controls="France and Spain", clustering_on=True, critical_distance=250, critical_time=10000)
mm.run_experiment(base_path, "france_spain", "fs_e_cd500_ct10000", controls="France and Spain", clustering_on=True, critical_distance=500, critical_time=10000)

mm.run_experiment(base_path, "france_spain", "fs_e_cd100_ct5000", controls="France and Spain", clustering_on=True, critical_distance=100, critical_time=5000)
mm.run_experiment(base_path, "france_spain", "fs_e_cd250_ct5000", controls="France and Spain", clustering_on=True, critical_distance=250, critical_time=5000)
mm.run_experiment(base_path, "france_spain", "fs_e_cd500_ct5000", controls="France and Spain", clustering_on=True, critical_distance=500, critical_time=5000)

mm.run_experiment(base_path, "france_spain", "fs_e_cd100_ct1000", controls="France and Spain", clustering_on=True, critical_distance=100, critical_time=1000)
mm.run_experiment(base_path, "france_spain", "fs_e_cd250_ct1000", controls="France and Spain", clustering_on=True, critical_distance=250, critical_time=1000)
mm.run_experiment(base_path, "france_spain", "fs_e_cd500_ct1000", controls="France and Spain", clustering_on=True, critical_distance=500, critical_time=1000)

mm.run_experiment(base_path, "australia", "au_e_with_gamma", controls="Australia", processed_targets=False)
# # # mm.run_experiment(base_path, "australia", "au_e_d0-10k", controls="Australia", filter_min_date = 0, filter_max_date = 10000)
# # # mm.run_experiment(base_path, "australia", "au_e_d10k-20k", controls="Australia", filter_min_date = 10000, filter_max_date = 20000)
# # # mm.run_experiment(base_path, "australia", "au_e_d20k-30k", controls="Australia", filter_min_date = 20000, filter_max_date = 30000)
# # # mm.run_experiment(base_path, "australia", "au_e_d30k-80k", controls="Australia", filter_min_date = 30000, filter_max_date = 80000)

mm.run_experiment(base_path, "australia", "au_e_cd100_ct10000", controls="Australia", clustering_on=True, critical_distance=100, critical_time=10000)
mm.run_experiment(base_path, "australia", "au_e_cd250_ct10000", controls="Australia", clustering_on=True, critical_distance=250, critical_time=10000)
mm.run_experiment(base_path, "australia", "au_e_cd500_ct10000_with_gamma_error0.001", controls="Australia", clustering_on=True, critical_distance=500, critical_time=10000)

mm.run_experiment(base_path, "australia", "au_e_cd100_ct5000", controls="Australia", clustering_on=True, critical_distance=100, critical_time=5000)
mm.run_experiment(base_path, "australia", "au_e_cd250_ct5000", controls="Australia", clustering_on=True, critical_distance=250, critical_time=5000)
mm.run_experiment(base_path, "australia", "au_e_cd500_ct5000", controls="Australia", clustering_on=True, critical_distance=500, critical_time=5000)

mm.run_experiment(base_path, "australia", "au_e_cd100_ct1000", controls="Australia", clustering_on=True, critical_distance=100, critical_time=1000)
mm.run_experiment(base_path, "australia", "au_e_cd250_ct1000", controls="Australia", clustering_on=True, critical_distance=250, critical_time=1000)
mm.run_experiment(base_path, "australia", "au_e_cd500_ct1000", controls="Australia", clustering_on=True, critical_distance=500, critical_time=1000)


# # # mm.run_experiment(base_path, "no_australia_a", "no_australia_e", processed_targets=True)
# # # mm.run_experiment(base_path, "no_frsp_a", "no_frsp_e")
# # mm.run_experiment(base_path, "no_au_frsp_a", "no_au_frsp_e", processed_targets=True)


# # # # ##############
# # # # # TIMMERMANN #
# # # # ##############

# # mm.run_experiment(base_path, "rockpaintings v8a", "full_a_t", population_data_name="Timmermann", )

# # # mm.run_experiment(base_path, "rockpaintings v8a", "full_a_t_dir", population_data_name="Timmermann", filter_not_direct=True)
# # # mm.run_experiment(base_path, "rockpaintings v8a", "full_a_t_fig", population_data_name="Timmermann", filter_not_figurative=True)

# # mm.run_experiment(base_path, "rockpaintings v8a", "full_a_t_lat-40--30", population_data_name="Timmermann", filter_min_lat = -40, filter_max_lat = -30)
# # # mm.run_experiment(base_path, "rockpaintings v8a", "full_a_t_lat-30--20", population_data_name="Timmermann", filter_min_lat = -30, filter_max_lat = -20)
# # # mm.run_experiment(base_path, "rockpaintings v8a", "full_a_t_lat-20--10", population_data_name="Timmermann", filter_min_lat = -20, filter_max_lat = -10)
# # mm.run_experiment(base_path, "rockpaintings v8a", "full_a_t_lat-10-0", population_data_name="Timmermann", filter_min_lat = -10, filter_max_lat = 0)
# # # mm.run_experiment(base_path, "rockpaintings v8a", "full_a_t_lat0-10", population_data_name="Timmermann", filter_min_lat = 0, filter_max_lat = 10)
# # # mm.run_experiment(base_path, "rockpaintings v8a", "full_a_t_lat10-20", population_data_name="Timmermann", filter_min_lat = 10, filter_max_lat = 20)
# # mm.run_experiment(base_path, "rockpaintings v8a", "full_a_t_lat20-30", population_data_name="Timmermann", filter_min_lat = 20, filter_max_lat = 30)
# # mm.run_experiment(base_path, "rockpaintings v8a", "full_a_t_lat30-40", population_data_name="Timmermann", filter_min_lat = 30, filter_max_lat = 40)
# # # mm.run_experiment(base_path, "rockpaintings v8a", "full_a_t_lat40-50", population_data_name="Timmermann", filter_min_lat = 40, filter_max_lat = 50)
# # mm.run_experiment(base_path, "rockpaintings v8a", "full_a_t_lat50-60", population_data_name="Timmermann", filter_min_lat = 50, filter_max_lat = 60)

# # # mm.run_experiment(base_path, "rockpaintings v8a", "full_a_t_d0-10k", population_data_name="Timmermann", filter_min_date = 0, filter_max_date = 10000)
# # mm.run_experiment(base_path, "rockpaintings v8a", "full_a_t_d10k-20k", population_data_name="Timmermann", filter_min_date = 10000, filter_max_date = 20000)
# # mm.run_experiment(base_path, "rockpaintings v8a", "full_a_t_d20k-30k", population_data_name="Timmermann", filter_min_date = 20000, filter_max_date = 30000)
# # mm.run_experiment(base_path, "rockpaintings v8a", "full_a_t_d30k-40k", population_data_name="Timmermann", filter_min_date = 30000, filter_max_date = 40000)
# # mm.run_experiment(base_path, "rockpaintings v8a", "full_a_t_d40k-50k", population_data_name="Timmermann", filter_min_date = 40000, filter_max_date = 50000)
# # mm.run_experiment(base_path, "rockpaintings v8a", "full_a_t_d10k-30k", population_data_name="Timmermann", filter_min_date = 10000, filter_max_date = 30000)


# # mm.run_experiment(base_path, "france_spain", "fs_t", population_data_name="Timmermann", controls="France and Spain")
# # # mm.run_experiment(base_path, "france_spain", "fs_t_d0-10k", population_data_name="Timmermann", controls="France and Spain", filter_min_date = 0, filter_max_date = 10000)
# # # mm.run_experiment(base_path, "france_spain", "fs_t_d10k-20k", population_data_name="Timmermann", controls="France and Spain", filter_min_date = 10000, filter_max_date = 20000)
# # # mm.run_experiment(base_path, "france_spain", "fs_t_d20k-30k", population_data_name="Timmermann", controls="France and Spain", filter_min_date = 20000, filter_max_date = 30000)
# # # mm.run_experiment(base_path, "france_spain", "fs_t_d30k-80k", population_data_name="Timmermann", controls="France and Spain", filter_min_date = 30000, filter_max_date = 80000)

# # # mm.run_experiment(base_path, "france_spain", "fs_t_cd100_ct10000", population_data_name="Timmermann", controls="France and Spain", clustering_on=True, critical_distance=100, critical_time=10000)
# # # mm.run_experiment(base_path, "france_spain", "fs_t_cd250_ct10000", population_data_name="Timmermann", controls="France and Spain", clustering_on=True, critical_distance=250, critical_time=10000)
# # # mm.run_experiment(base_path, "france_spain", "fs_t_cd500_ct10000", population_data_name="Timmermann", controls="France and Spain", clustering_on=True, critical_distance=500, critical_time=10000)

# # # mm.run_experiment(base_path, "france_spain", "fs_t_cd100_ct5000", population_data_name="Timmermann", controls="France and Spain", clustering_on=True, critical_distance=100, critical_time=5000)
# # # mm.run_experiment(base_path, "france_spain", "fs_t_cd250_ct5000", population_data_name="Timmermann", controls="France and Spain", clustering_on=True, critical_distance=250, critical_time=5000)
# # # mm.run_experiment(base_path, "france_spain", "fs_t_cd500_ct5000", population_data_name="Timmermann", controls="France and Spain", clustering_on=True, critical_distance=500, critical_time=5000)

# # # mm.run_experiment(base_path, "france_spain", "fs_t_cd100_ct1000", population_data_name="Timmermann", controls="France and Spain", clustering_on=True, critical_distance=100, critical_time=1000)
# # # mm.run_experiment(base_path, "france_spain", "fs_t_cd250_ct1000", population_data_name="Timmermann", controls="France and Spain", clustering_on=True, critical_distance=250, critical_time=1000)
# # # mm.run_experiment(base_path, "france_spain", "fs_t_cd500_ct1000", population_data_name="Timmermann", controls="France and Spain", clustering_on=True, critical_distance=500, critical_time=1000)

# # mm.run_experiment(base_path, "australia", "au_t", population_data_name="Timmermann", controls="Australia")
# # # mm.run_experiment(base_path, "australia", "au_t_d0-10k", population_data_name="Timmermann", controls="Australia", filter_min_date = 0, filter_max_date = 10000)
# # # mm.run_experiment(base_path, "australia", "au_t_d10k-20k", population_data_name="Timmermann", controls="Australia", filter_min_date = 10000, filter_max_date = 20000)
# # # mm.run_experiment(base_path, "australia", "au_t_d20k-30k", population_data_name="Timmermann", controls="Australia", filter_min_date = 20000, filter_max_date = 30000)
# # # mm.run_experiment(base_path, "australia", "au_t_d30k-80k", population_data_name="Timmermann", controls="Australia", filter_min_date = 30000, filter_max_date = 80000)

# # # mm.run_experiment(base_path, "australia", "au_t_cd100_ct10000", population_data_name="Timmermann", controls="Australia", clustering_on=True, critical_distance=100, critical_time=10000)
# # # mm.run_experiment(base_path, "australia", "au_t_cd250_ct10000", population_data_name="Timmermann", controls="Australia", clustering_on=True, critical_distance=250, critical_time=10000)
# # # mm.run_experiment(base_path, "australia", "au_t_cd500_ct10000", population_data_name="Timmermann", controls="Australia", clustering_on=True, critical_distance=500, critical_time=10000)

# # # mm.run_experiment(base_path, "australia", "au_t_cd100_ct5000", population_data_name="Timmermann", controls="Australia", clustering_on=True, critical_distance=100, critical_time=5000)
# # # mm.run_experiment(base_path, "australia", "au_t_cd250_ct5000", population_data_name="Timmermann", controls="Australia", clustering_on=True, critical_distance=250, critical_time=5000)
# # # mm.run_experiment(base_path, "australia", "au_t_cd500_ct5000", population_data_name="Timmermann", controls="Australia", clustering_on=True, critical_distance=500, critical_time=5000)

# # # mm.run_experiment(base_path, "australia", "au_t_cd100_ct1000", population_data_name="Timmermann", controls="Australia", clustering_on=True, critical_distance=100, critical_time=1000)
# # # mm.run_experiment(base_path, "australia", "au_t_cd250_ct1000", population_data_name="Timmermann", controls="Australia", clustering_on=True, critical_distance=250, critical_time=1000)
# # # mm.run_experiment(base_path, "australia", "au_t_cd500_ct1000", population_data_name="Timmermann", controls="Australia", clustering_on=True, critical_distance=500, critical_time=1000)


# # # mm.run_experiment(base_path, "no_australia_a", "no_australia_t", population_data_name="Timmermann")
# # # mm.run_experiment(base_path, "no_frsp_a", "no_frsp_t", population_data_name="Timmermann")
# # mm.run_experiment(base_path, "no_au_frsp_a", "no_au_frsp_t", population_data_name="Timmermann")


# ####################
# # Clustering Tests #
# ####################

# mm.run_experiment(base_path, "rockpaintings v8a", "full_a_e_cd100_ct10000", processed_targets=True, clustering_on=True, critical_distance=100, critical_time=10000)
# mm.run_experiment(base_path, "rockpaintings v8a", "full_a_e_cd250_ct10000", processed_targets=True, clustering_on=True, critical_distance=250, critical_time=10000)
# mm.run_experiment(base_path, "rockpaintings v8a", "full_a_e_cd500_ct10000", processed_targets=True, clustering_on=True, critical_distance=500, critical_time=10000)

# mm.run_experiment(base_path, "rockpaintings v8a", "full_a_e_cd100_ct5000", processed_targets=True, clustering_on=True, critical_distance=100, critical_time=5000)
# mm.run_experiment(base_path, "rockpaintings v8a", "full_a_e_cd250_ct5000", processed_targets=True, clustering_on=True, critical_distance=250, critical_time=5000)
# mm.run_experiment(base_path, "rockpaintings v8a", "full_a_e_cd500_ct5000", processed_targets=True, clustering_on=True, critical_distance=500, critical_time=5000)

# mm.run_experiment(base_path, "rockpaintings v8a", "full_a_e_cd100_ct1000", processed_targets=True, clustering_on=True, critical_distance=100, critical_time=1000)
# mm.run_experiment(base_path, "rockpaintings v8a", "full_a_e_cd250_ct1000", processed_targets=True, clustering_on=True, critical_distance=250, critical_time=1000)
# mm.run_experiment(base_path, "rockpaintings v8a", "full_a_e_cd500_ct1000", processed_targets=True, clustering_on=True, critical_distance=500, critical_time=1000)


# # mm.run_experiment(base_path, "rockpaintings v8a", "full_a_t_cd100_ct10000",  population_data_name="Timmermann", clustering_on=True, critical_distance=100, critical_time=10000)
# # mm.run_experiment(base_path, "rockpaintings v8a", "full_a_t_cd250_ct10000",  population_data_name="Timmermann", clustering_on=True, critical_distance=250, critical_time=10000)
# # mm.run_experiment(base_path, "rockpaintings v8a", "full_a_t_cd500_ct10000", population_data_name="Timmermann", clustering_on=True, critical_distance=500, critical_time=10000)

# # mm.run_experiment(base_path, "rockpaintings v8a", "full_a_t_cd100_ct5000",  population_data_name="Timmermann", clustering_on=True, critical_distance=100, critical_time=5000)
# # mm.run_experiment(base_path, "rockpaintings v8a", "full_a_t_cd250_ct5000",  population_data_name="Timmermann", clustering_on=True, critical_distance=250, critical_time=5000)
# # mm.run_experiment(base_path, "rockpaintings v8a", "full_a_t_cd500_ct5000",  population_data_name="Timmermann", clustering_on=True, critical_distance=500, critical_time=5000)

# # mm.run_experiment(base_path, "rockpaintings v8a", "full_a_t_cd100_ct1000",  population_data_name="Timmermann", clustering_on=True, critical_distance=100, critical_time=1000)
# # mm.run_experiment(base_path, "rockpaintings v8a", "full_a_t_cd250_ct1000", population_data_name="Timmermann", clustering_on=True, critical_distance=250, critical_time=1000)
# # mm.run_experiment(base_path, "rockpaintings v8a", "full_a_t_cd500_ct1000",  population_data_name="Timmermann", clustering_on=True, critical_distance=500, critical_time=1000)
