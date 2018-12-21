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
    def __init__(self,orig_lat, orig_lon,lat_nw,lon_nw,lat_se,lon_se,location,date_from, date_to, country,is_direct,calibrated,kind,figurative,source,is_controversial, age_estimation, cluster_id):
        self.location=location
        self.lat_nw=lat_nw
        self.lon_nw=lon_nw
        self.lat_se=lat_se
        self.lon_se=lon_se
        self.orig_lat=orig_lat
        self.orig_lon=orig_lon
        self.date_from=date_from
        self.date_to=date_to
        self.country=country
        self.is_direct=is_direct
        self.calibrated=calibrated
        self.kind=kind
        self.figurative=figurative
        self.source=source
        self.is_controversial=is_controversial
        self.age_estimation = age_estimation;
        self.cluster_id=cluster_id

    def __str__(self):
        ret = self.location + "\n";
        # ret += "    Cluster ID: " + str(self.cluster_id) + "\n"
        ret += "    Actual Latitude: " + str(self.orig_lat) + "\n"
        ret += "    Actual Longitude: " + str(self.orig_lon) + "\n"
        ret += "    Latitude (NW): " + str(self.lat_nw) + "\n"
        ret += "    Latitude (SE): " + str(self.lat_se) + "\n"
        ret += "    Longitude (NW): " + str(self.lon_nw) + "\n"
        ret += "    Longitude (SE): " + str(self.lon_se) + "\n"
        ret += "    Date From: " + str(self.date_from) + "\n"
        ret += "    Date To: " + str(self.date_to)

        return ret;

class PopulationData():

    def __init__(self, name, is_active, lat_array, lon_array, time_array, density_array, time_multiplier, density_multiplier,bin_size, max_population, max_for_uninhabited, ascending_time):
        
        self.name = name
        self.lat_array = lat_array
        self.lon_array = lon_array
        self.time_array = time_array
        self.density_array = density_array
        self.time_multiplier = time_multiplier
        self.density_multiplier=density_multiplier
        self.bin_size = bin_size
        self.max_population = max_population
        self.is_active = is_active
        self.max_for_uninhabited = max_for_uninhabited
        self.ascending_time = ascending_time
        
        
        
        
        
        
        
        
        

        