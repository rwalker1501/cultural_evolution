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
		self.print_label('Population analysis tool v. 1 by Richard Walker')
		user_option = 0
		while user_option != '9':
			target_list = self.main_program.get_current_target_list()
			population_data = self.main_program.get_population_data()
			parameters = self.main_program.get_parameters();

			print("Active Population Data")
			for key, val in population_data.iteritems():
				if population_data[key].is_active:
					print("   " + key)

			self.print_target_list(target_list)

			print("\nProcessed targets loaded: " + str(parameters['processed_targets']))
			filters_applied = parameters['filters_applied']
			if len(filters_applied) > 0:
				print("Filters applied: " + str(filters_applied));
			else:
				print("Filters applied: None");
			
			print('-----------------------')
			print('1)    Define list of targets')
			print('2)    Filter target list')
			print('3)    Load processed targets')
			print('4)    Set parameters')
			print('5)    Create density map plots')
			print('6)    Clear processed targets')
			print('7)    Generate results')
			print('8)    Add new population data source')
			print('9)   Exit')
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
				self.create_map_plots();
			elif user_option=='6':
				user_option = ""
				while user_option != 'y' and user_option != 'n':
					user_option = raw_input("Are you sure you want to delete all processed targets? [y/n]: ")
				if user_option == "y":
					tam.clear_processed_targets(base_path)
			elif user_option=='7':
				for pop_data_name, pop_data in population_data.iteritems():
					if pop_data.is_active:
						directory=raw_input("Insert directory name: ")
						report = self.main_program.generate_results(pop_data, target_list, base_path, directory)
						print(report);
			elif user_option=='8':
				self.add_population_data()
			elif user_option!='9':
				print("Invalid option. Try again.")

	def create_map_plots(self):
		self.print_label("Create Density Map Plots")
		print()


		population_data = self.main_program.get_population_data()
		print("Active Population Data")
		for pop_data in population_data:
			if pop_data.is_active:
				print("   " + pop_data.name)

		print("\n------")
		print("1) Plot densities on map by time point")
		print("2) Plot densities on map by time range")
		print("3) Plot densities on map by density and time range")
		print("4) Cancel")

		plot_option = self.get_number("Insert option: ", 1, 4)
		if plot_option == 1:
			time = self.get_number("Input time point: ", 0, 150000)
			for population_data_source in population_data:
				if population_data_source.is_active:
					self.main_program.plot_population_by_time(population_data_source, time);
		elif plot_option == 2:
			start_time = self.get_number("Input start time: ", 0, 150000)
			end_time = self.get_number("Input end time: ", 0, 150000)
			for population_data_source in population_data:
				if population_data_source.is_active:
					self.main_program.plot_population_by_time_range(population_data_source, start_time, end_time);
		elif plot_option == 3:
			min_density = self.get_number("Input min density: ", 0, 150000)
			max_density = self.get_number("Input max density: ", 0, 150000)
			start_time = self.get_number("Input start time: ", 0, 150000)
			end_time = self.get_number("Input end time: ", 0, 150000)
			for population_data_source in population_data:
				if population_data_source.is_active:
					self.main_program.plot_population_by_range(population_data_source, min_density, max_density, start_time, end_time);

		if plot_option != 4:
			print("Create another map plot?")
			print("\n------")
			print("1) Yes");
			print("2) No");
			plot_option = self.get_number("Insert option: ", 1, 2);
			if plot_option == 1:
				self.create_map_plots();
			

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

		parameters = self.main_program.get_parameters();
		filters_applied = parameters['filters_applied']

		new_list = deepcopy(some_target_list)
		user_option = ""
		while user_option != '10' and user_option != '11':

			self.print_target_list(new_list)
			if len(filters_applied) > 0:
				print("\nFilters applied: " + str(filters_applied));
			else:
				print("\nFilters applied: None");
			print("\n------")
			print("1) Exclude targets more recent than...")
			print("2) Exclude targets with indirect measurements")
			print("3) Only include targets with exact age")
			print("4) Exclude non-figurative targets")
			print("5) Exclude targets with controversial measurements")
			print("6) Filter targets by date range")
			print("7) Filter targets by latitude range")
			print("8) Cancel")
			print("9) Save and exit")

			user_option = raw_input("Choose an option: ")
			if user_option=='1':
				value = self.get_number("Exclude all targets more recent than: ", 0, 150000)
				new_list, filters_applied = tam.filter_targets_for_date_before(new_list, value, filters_applied)
			elif user_option=='2':
				new_list, filters_applied = tam.filter_targets_for_not_direct(new_list, filters_applied)
			elif user_option=='3':
				new_list, filters_applied = tam.filter_targets_for_not_exact_age(new_list, filters_applied)
			elif user_option=='4':
				new_list, filters_applied = tam.filter_targets_for_not_figurative(new_list, filters_applied)
			elif user_option=='5':
				new_list, filters_applied = tam.filter_targets_for_not_controversial(new_list, filters_applied)
			elif user_option=='6':
				minimum_date = self.get_number("Insert minimum date: ", 0, 150000)
				maximum_date = self.get_number("Insert maximum date: ", 0, 150000)
				new_list, filters_applied = tam.filter_targets_for_date(new_list, minimum_date, maximum_date, filters_applied)
			elif user_option=='7':
				minimum_latitude = self.get_number("Insert minimum latitude: ", -90, 90)
				maximum_latitude = self.get_number("Insert maximum latitude: ", -90, 90)
				new_list, filters_applied = tam.filter_targets_for_latitude(new_list, minimum_latitude, maximum_latitude, filters_applied)
			elif user_option=='9':
				new_list = tam.reset_cluster_id_values(new_list);
				self.main_program.set_parameter('filters_applied', filters_applied);
				self.main_program.set_target_list(new_list)

			new_list = tam.reset_cluster_id_values(new_list);



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

			success, target_list, dataframe, globals_dataframe = tam.load_processed_targets(base_path, filename)
			if not success:
				print("LOADING PROCESSED TARGETS FAILED: Missing files")
				return False
			else:
				self.main_program.set_target_list(target_list)
				self.main_program.set_dataframe(dataframe, globals_dataframe)
				return True

	def set_parameters(self):
		user_option = ""
		while user_option != '10':
			self.print_label("Set Parameters")
			parameters = self.main_program.get_parameters();

			print("")
			print("1) Set active population data")
			print("2) Define date window for target: " + str(parameters['date_window']))
			print("3) Toggle default max population for areas considered as uninhabited: " + str(parameters['default_max_for_uninhabited']) )
			print("4) Define max population for areas considered as uninhabited (sets default=False): " + str(parameters['user_max_for_uninhabited']))
			print("5) Define globals range: " + str(parameters['globals_type']))
			print("6) Set minimum latitude for globals: " + str(parameters['min_lat']))
			print("7) Set maximum latitude for globals: " + str(parameters['max_lat']) )
			print("8) Set minimum date for globals: " + str(parameters['min_date']))
			print("9) Set maximum date for globals: " + str(parameters['max_date']))
			print("10) Save and exit")
			user_option = raw_input("Choose an option: ")
			if user_option=='1':
				self.toggle_population_data()
			elif user_option=='2':
				self.main_program.set_parameter('date_window', self.get_number('Insert date window: ',0,150000))
			elif user_option=='3':
				self.main_program.set_parameter('default_max_for_uninhabited', not parameters['default_max_for_uninhabited'])
			elif user_option=='4':
				self.main_program.set_parameter('user_max_for_uninhabited', self.get_number('Insert max density for uninhabited - for no max write -1: ',-1,5000))
			elif user_option=='5':
				globals_option = 1
				self.print_label("Define globals Range")
				print("1) All")
				print("2) No Equatorials")
				print("3) Australia")
				print("4) France and Spain")
				globals_option = self.get_number("Insert option: ", 1, 6)
				if globals_option == 1:
					self.main_program.set_parameter('globals_type', "All")
				elif globals_option == 2:
					self.main_program.set_parameter('globals_type', "No equatorials")
				elif globals_option == 3:
					self.main_program.set_parameter('globals_type', "Australia")
				elif globals_option == 4:
					self.main_program.set_parameter('globals_type', "France and Spain")
			elif user_option=='6':
				self.main_program.set_parameter('min_lat', self.get_number("Insert minimum latitude: ", -90, 90))
			elif user_option=='7':
				self.main_program.set_parameter('max_lat', self.get_number("Insert maximum latitude: ", -90, 90))
			elif user_option=='8':
				self.main_program.set_parameter('min_date', self.get_number("Insert minimum date: ", 0, 1000000))
			elif user_option=='9':
				self.main_program.set_parameter('max_date', self.get_number("Insert maximum date: ", 0, 1000000))
			elif user_option!='10':
				print("Invalid option. Try again.")


	def toggle_population_data(self):
		valid_number=False
		user_option=''
		while valid_number is False:
			population_data_sources = self.main_program.get_population_data()
			self.print_label("Toggle Active Population Data")
			valid_number = True
			pop_data_names = [];

			i = 0;
			for key, val in population_data_sources.iteritems():
				to_print = str(i) + ") " + key + " - " + str(population_data_sources[key].is_active);
				pop_data_names.append(key)
				print(to_print);
				i+=1;
			print(str(i) + ") Exit")

			user_option = raw_input("Choose an option: ")
			try:
				user_int=int(user_option)

				if user_int > len(pop_data_names):
					print('not a valid number')
					valid_number=False

				if valid_number:
					if user_int == len(pop_data_names):
						return
					pop_data_name = pop_data_names[user_int]
					pop_data_active_status = population_data_sources[pop_data_name].is_active
					self.main_program.set_population_data_active(not pop_data_active_status, pop_data_name)
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
			return
		for i in range(0, len(target_list)):
			target = target_list[i]
			print(str(i+1) + ") ")
			print(target);
		print("\nNumber of Targets: " + str(len(target_list)));



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