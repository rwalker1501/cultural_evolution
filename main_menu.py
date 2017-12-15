import os
import sys
import random
import target_module as tam
from copy import deepcopy
from os.path import isfile, join
from main_module import MainProgram
from classes_module import Target, PopulationData


class Driver:
	def __init__(self):
		self.main_program = MainProgram()

	def run(self):
		print("Reading files..")
		self.main_program.load_population_data();
		base_path = self.main_program.get_base_path()

		print('')
		self.print_label('Population analysis tool v. 0.1 by Richard Walker')
		user_option = 0
		while user_option != '9':
			target_list = self.main_program.get_current_target_list()
			population_data = self.main_program.get_population_data()
			dataframe_loaded = self.main_program.get_dataframe_loaded()
			filters_applied = self.main_program.get_filters_applied()

			print("\nActive Population Data")
			for pop_data in population_data:
				if pop_data.is_active:
					print("   " + pop_data.name)

			self.print_target_list(target_list)

			print("\nProcessed targets loaded: " + str(dataframe_loaded))
			print("Filters applied: " + str(filters_applied))
			print('-----------------------')
			print('1)    Define list of targets')
			print('2)    Filter target list')
			print('3)    Load processed targets')
			print('4)    Set parameters')
			print('5)    Plot population by time')
			print('6)    Clear processed targets')
			print('7)    Generate results')
			print('8)    Add new population data source')
			print('9)    Exit')
			random.seed(2002)
			user_option = raw_input("Choose an option: ")
			if user_option=='1':
				self.define_target_list(target_list)
			elif user_option=='2':
				self.filter_target_list(target_list)
			elif user_option=='3':
				self.load_processed_targets()
			elif user_option=='4':
				self.set_parameters()
			elif user_option=='5':
				time = self.get_number("Input time: ", 0, 150000)
				self.main_program.plot_population(time)
			elif user_option=='6':
				user_option = ""
				while user_option != 'y' and user_option != 'n':
					user_option = raw_input("Are you sure you want to delete all processed targets? [y/n]: ")
				if user_option == "y":
					tam.clear_processed_targets(base_path)
			elif user_option=='7':
				for population_data_source in population_data:
					if population_data_source.is_active:
						directory=raw_input("Insert directory name: ")
						report = self.main_program.generate_results(population_data_source, target_list, base_path, directory)
						print(report)
			elif user_option=='8':
				self.add_population_data()
			elif user_option!='9':
				print("Invalid option. Try again.")

	def define_target_list(self, some_target_list):
		self.print_label("Define target list")
		print()

		new_list = deepcopy(some_target_list)
		user_option = ""
		while user_option != '3' and user_option != '4':
			self.print_target_list(new_list)

			print("\n------")
			print("1) Read target list")
			print("2) Save target list")
			print("3) Cancel")
			print("4) Save and exit")
			user_option = raw_input("Choose an option: ")
			if user_option=='1':
				base_path = self.main_program.get_base_path()
				filename=self.ask_filename(os.path.join(base_path,"targets"))
				new_list = self.main_program.read_target_list(filename)
			elif user_option=='2':
				filename=raw_input('Insert file name - no extension:  ')
				self.main_program.save_target_list(filename, new_list)
			elif user_option=='3':
				self.main_program.set_target_list(some_target_list)
			elif user_option=='4':
				self.main_program.set_target_list(new_list)


	def filter_target_list(self, some_target_list):
		self.print_label("Filter target list")

		clustering = self.main_program.get_clustering()
		critical_distance = self.main_program.get_critical_distance()
		critical_time = self.main_program.get_critical_time()

		new_list = deepcopy(some_target_list)
		user_option = ""
		while user_option != '9' and user_option != '10':
			filters_applied = self.main_program.get_filters_applied()

			self.print_target_list(new_list)

			print("\nClustering: " + str(clustering))
			print("Critical distance: " + str(critical_distance))
			print("Critical time: " + str(critical_time))
			print("\n------")
			print("1) Exclude targets more recent than...")
			print("2) Exclude targets with indirect measurements")
			print("3) Exclude non-figurative targets")
			print("4) Exclude targets with controversial measurements")
			print("5) Filter targets by date range")
			print("6) Filter targets by latitude range")
			print("7) Turn clustering ON")
			print("8) Turn clustering OFF")
			print("9) Cancel")
			print("10) Save and exit")

			user_option = raw_input("Choose an option: ")
			if user_option=='1':
				value = self.get_number("Exclude all targets more recent than: ", 0, 150000)
				new_list, filters_applied = tam.filter_targets_for_date_before(new_list, value, filters_applied)
			elif user_option=='2':
				new_list, filters_applied = tam.filter_targets_for_not_direct(new_list, filters_applied)
			elif user_option=='3':
				new_list, filters_applied = tam.filter_targets_for_not_figurative(new_list, filters_applied)
			elif user_option=='4':
				new_list, filters_applied = tam.filter_targets_for_not_controversial(new_list, filters_applied)
			elif user_option=='5':
				minimum_date = self.get_number("Insert minimum date: ", 0, 150000)
				maximum_date = self.get_number("Insert maximum date: ", 0, 150000)
				new_list, filters_applied = tam.filter_targets_for_date(new_list, minimum_date, maximum_date, filters_applied)
			elif user_option=='6':
				minimum_latitude = self.get_number("Insert minimum latitude: ", -90, 90)
				maximum_latitude = self.get_number("Insert maximum latitude: ", -90, 90)
				new_list, filters_applied = tam.filter_targets_for_latitude(new_list, minimum_latitude, maximum_latitude, filters_applied)
			elif user_option=='7':
				critical_distance = self.get_number("Insert critical distance: ", 1, 1000)
				critical_time = self.get_number("Insert critical time: ", 1, 150000)
				clustering = True
			elif user_option=='8':
				clustering = False
				critical_distance = 0
			elif user_option=='10':
				self.main_program.set_clustering(clustering)
				self.main_program.set_critical_distance(critical_distance)
				self.main_program.set_critical_time(critical_time)
				self.main_program.set_filters_applied(filters_applied)
				self.main_program.set_target_list(new_list)



	def load_processed_targets(self):
		base_path = self.main_program.get_base_path()
		processed_targets_dir = os.path.join(base_path, "processed_targets")
		filenames = [f for f in os.listdir(processed_targets_dir) if isfile(join(processed_targets_dir, f))]
		
		if len(filenames) == 0:
			print("")
			self.print_label("There are no processed targets.")
		else:
			for name in filenames:
				if '_dataframe.csv' in name:
					name = name.replace('_dataframe.csv', '')
					print(name) 
			print("\n")

			chosen_file = ""
			while not (os.path.isfile(os.path.join(processed_targets_dir, chosen_file))):
				filename=raw_input('Pick a file to be loaded: ')
				chosen_file = filename + '_dataframe.csv'
				print(chosen_file)

			target_list, dataframe, controls_dataframe = tam.load_processed_targets(base_path, filename)
			self.main_program.set_target_list(target_list)
			self.main_program.set_dataframe(dataframe, controls_dataframe)

	def set_parameters(self):
		user_option = ""
		while user_option != '9':
			self.print_label("Set Parameters")
			date_window = self.main_program.get_date_window()
			mfu = self.main_program.get_user_max_for_uninhabited()
			default_mfu = self.main_program.get_default_mfu()
			min_controls = self.main_program.get_minimum_controls()
			perform_cv = self.main_program.get_perform_cross_validation()
			num_kfolds = self.main_program.get_number_of_kfolds()
			controls = self.main_program.get_controls()
			print("")
			print("1) Set active population data")
			print("2) Define date window for target: " + str(date_window))
			print("3) Toggle default max population for areas considered as uninhabited: " + str(default_mfu) )
			print("4) Define max population for areas considered as uninhabited (sets default=False): " + str(mfu))
			print("5) Set minimum controls: " + str(min_controls))
			print("6) Toggle perform cross validation: " + str(perform_cv))
			print("7) Set folds for cross validation: " + str(num_kfolds))
			print("8) Define controls range: " + str(controls))
			print("9) Save and exit")
			user_option = raw_input("Choose an option: ")
			if user_option=='1':
				self.toggle_population_data()
			elif user_option=='2':
				self.main_program.set_date_window(self.get_number('Insert date window: ',0,150000))
			elif user_option=='3':
				self.main_program.set_default_mfu(not default_mfu)
			elif user_option=='4':
				self.main_program.set_user_max_for_uninhabited(self.get_number('Insert max density for uninhabited - for no max write -1: ',-1,5000))
			elif user_option=='5':
				self.main_program.set_minimum_controls(self.get_number('Insert minimum controls: ', 0, 500000))
			elif user_option=='6':
				self.main_program.set_perform_cross_validation(not perform_cv)
			elif user_option=='7':
				self.main_program.set_number_of_kfolds(self.get_number("Insert number of kfolds: ", 1, 200))
			elif user_option=='8':
				controls_option = 1
				self.print_label("Define Controls Range")
				print("1) All")
				print("2) Australia")
				print("3) France and Spain")
				print("4) Trial Latitudes")
				print("5) Trial Latitudes 2")
				print("6) No Empty Lats (New Controls)")
				controls_option = self.get_number("Insert option: ", 1, 6)
				if controls_option == 1:
					self.main_program.set_controls("All")
				elif controls_option == 2:
					self.main_program.set_controls("Australia")
				elif controls_option == 3:
					self.main_program.set_controls("France and Spain")
				elif controls_option == 4:
					self.main_program.set_controls("Trial Latitudes")
				elif controls_option == 5:
					self.main_program.set_controls("Trial Latitudes 2")
				elif controls_option == 6:
					self.main_program.set_controls("No Empty Lats")
			elif user_option!='9':
				print("Invalid option. Try again.")


	def toggle_population_data(self):
		valid_number=False
		user_option=''
		while valid_number is False:
			population_data_sources = self.main_program.get_population_data()
			self.print_label("Toggle Active Population Data")
			valid_number = True

			for i in range(1, len(population_data_sources)):
				population_data = population_data_sources[i]
				to_print = str(i) + ") " + population_data.name + " - " + str(population_data.is_active)
				print(to_print)
			print(str(len(population_data_sources)) + ") Exit")

			user_option = raw_input("Choose an option: ")
			try:
				user_int=int(user_option)

				if user_int > len(population_data_sources):
					print('not a valid number')
					valid_number=False

				if valid_number:
					if user_int == len(population_data_sources):
						return
					user_active = population_data_sources[user_int].is_active
					self.main_program.set_population_data_active(not user_active, user_int)
					valid_number = False

			except ValueError:
				print('not a valid number')
				valid_number=False

	def add_population_data(self):
		base_path = self.main_program.get_base_path()

		pop_data_name = raw_input("Name: ")
		
		pop_data_path = os.path.join(base_path, "population_data")
		pop_data_binary_path = "-"
		while(not os.path.isfile(pop_data_binary_path)):
			filename = raw_input("Insert valid binary data filename: ")
			pop_data_binary_path = os.path.join(pop_data_path, filename)

		pop_data_info_path = "-"
		while(not os.path.isfile(pop_data_info_path)):
			filename = raw_input("Insert valid info data filename: ")
			pop_data_info_path = os.path.join(pop_data_path, filename)

		self.main_program.add_population_data(pop_data_name, pop_data_binary_path, pop_data_info_path)
			

	def print_label(self, label):
		asterisks = ''
		for i in range(0, len(label)+2):
			asterisks += "*"
		print(asterisks)
		print(label)
		print(asterisks)

	def print_target_list(self, target_list):
		if len(target_list) == 0:
			print("<Empty Target List>")
			print("")
		date_window = self.main_program.get_date_window()
		for i in range(0, len(target_list)):
			target = target_list[i]
			print(str(i+1) + ")")
			tam.print_target(target, date_window)


	def ask_filename(self, file_path):
		filename=""
		print('path=',file_path)
		os.chdir(file_path)
		filenames = [f for f in os.listdir(file_path) if isfile(join(file_path,f))]
		print("***********************")
		print("*** AVAILABLE FILES ***")
		print("***********************")
		for name in filenames:
			print(name)
		print("***********************")
		while not (os.path.isfile(filename+'.csv')):
			filename=raw_input('Insert file name - no extension:  ')
			print("filename is:", filename+".csv")
		return filename 
        

	def get_number(self, label, min, max):
		validated = False
		while not validated:
			valid_number = True;
			try:
				value = float(raw_input(label))
			except ValueError:
				print("Not a valid number")
				valid_number = False

			if valid_number and value >= min and value <= max:
				validated = True
				return value
			elif valid_number:
				print("Value out of range")

if __name__ == '__main__':
	dr = Driver()
	dr.run();