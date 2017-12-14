#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 20 11:18:02 2017

@author: rwalker
"""

#==============================================================================
# class dataStore():
#     def __init__(self):
#         self.theData=
#==============================================================================
         
class Target():
    def __init__(self,orig_lat, orig_lon,lat_nw,lon_nw,lat_se,lon_se,location,date_from,country,date_of_reference,is_direct,calibrated,kind,figurative,source,is_controversial,cluster_id):
        # dateTo is legacy - in current versions of code dateTo is calculated as date_from + dateWindow
        self.location=location
        self.lat_nw=lat_nw
        self.lon_nw=lon_nw
        self.lat_se=lat_se
        self.lon_se=lon_se
        self.orig_lat=orig_lat
        self.orig_lon=orig_lon
        self.date_from=date_from
        # self.dateTo=dateTo
        self.country=country
        self.date_of_reference=date_of_reference
        self.is_direct=is_direct
        self.calibrated=calibrated
        self.kind=kind
        self.figurative=figurative
        self.source=source
        self.is_controversial=is_controversial
        self.cluster_id=cluster_id

        
class PopulationData():

    def __init__(self, name, is_active, lat_array, lon_array, time_array, density_array, time_multiplier, bin_size, max_population, max_for_uninhabited, ascending_time):
        
        self.name = name
        self.lat_array = lat_array
        self.lon_array = lon_array
        self.time_array = time_array
        self.density_array = density_array
        self.time_multiplier = time_multiplier
        self.bin_size = bin_size
        self.max_population = max_population
        self.is_active = is_active
        self.max_for_uninhabited = max_for_uninhabited
        self.ascending_time = ascending_time
        
        
        
        
        
        
        
        
        

        