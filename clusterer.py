import matplotlib.pyplot as plt
import numpy as np
import plot_module as plm
from sklearn import cluster, mixture
from itertools import cycle, islice
from mpl_toolkits.mplot3d import Axes3D  # needed for 3D plots!
from copy import deepcopy
import sys

class ClusterAnalysis:
    
    
    def __init__(self):
        x=10
        self.EARTH_RADIUS = 6371.0  # in km
        
    
    def cluster_targets_by_dist_and_time(self, target_list, critical_dist, critical_time):
        #######################
        # Cluster by distance #
        #######################
        distance_clustered_targets, num_clusters = self.cluster_targets_by_distance(target_list, critical_dist)
        
        # group targets by cluster ids through a 2D array
        # [[targets cluster id 0], [targets with cluster id 1],(and so on)]
        target_groups = [[] for i in range(0, num_clusters)]
        for target in target_list:
            cluster_id = target.cluster_id
            target_groups[cluster_id].append(target)

        ##############################
        # Cluster each group by time #
        ##############################
        # each distance clustered groups would be
        # clustered by time. We then assign new cluster ids
        # for every target in the group
        prev_cluster_id = 0
        new_target_list = []
        for target_group in target_groups:
            # get the time for each target in order of targets
            time = []
            for target in target_group:
                time.append(target.date_from)

            # cluster time array using DBSCAN
            time = np.array(time).reshape(-1, 1)
            dbscan = cluster.DBSCAN(eps=critical_time, min_samples=1)
            dbscan.fit(time)

            # get the labels from clustering by time (starts from 0)
            # note that this is the label only for the single group
            labels = dbscan.labels_.astype(np.int)

            # we assign cluster ids based on the labels
            # and the last cluster id of the previous group
            # For example:
            #   prev_labels = [0, 1, 2, 3, 4, 5]
            #   prev_cluster_id = 5 + 1
            #   new_labels = [0, 1, 2, 3]
            #   new cluster ids: 6 + label -> 6, 7, 8, 9
            #   new prev_label = 9 + 1
            for i in range(0, len(target_group)):
                target = target_group[i]
                target.cluster_id = prev_cluster_id + labels[i]
                new_target_list.append(target)
                print(target.location,"cluster_id:",target.cluster_id)

            prev_cluster_id = prev_cluster_id + max(labels) + 1

        # sort targets by cluster id
        sorted_target_list=sorted(new_target_list,key=lambda an_id: an_id.cluster_id) 
        for target in sorted_target_list:
            print("cluster_id:",target.cluster_id, target.location)

        return sorted_target_list


    def cluster_targets_by_distance(self,target_list,critical_dist):

    
        # Get 3D cartesian coordinates of target list
        xx, yy, zz, lon, lat = self.read_target_list(target_list)
        
        dbscan = cluster.DBSCAN(eps=critical_dist/self.EARTH_RADIUS,  min_samples = 1, algorithm='ball_tree', metric='haversine')
        X = np.stack((np.deg2rad(lon), np.deg2rad(lat)), axis=1)
        dbscan.fit(X)
        y_pred = dbscan.labels_.astype(np.int)
        
        max_cluster=max(y_pred)
        cluster_index=max_cluster +1 # this is for all targets with no cluster (-1 cluster)
        
        new_labels=[]
        for label in y_pred:
            if label==-1:
                new_labels.append(cluster_index)
                cluster_index=cluster_index+1
            else:
                new_labels.append(label)

        i=0
        for target in target_list:
            target.cluster_id=new_labels[i]
            i=i+1
        sorted_target_list=sorted(deepcopy(target_list),key=lambda an_id: an_id.cluster_id)
        return sorted_target_list, cluster_index # cluster index also indicates number of clusters
    

    def read_target_list(self,target_list):
        """
        Reads a target_list
        :return: xx, yy, zz - 3D coordinates (in km); lon, lat - spherical coordinates (in deg)
        """
       
        lat = []
        lon = []
        if len(target_list)==0:
            print("No targets to cluster")
            sys.exit
        for target in target_list:
                lat.append(target.orig_lat)
                lon.append(target.orig_lon)

        lon = np.array(lon)
        lat = np.array(lat)
        X = np.array([list(self.spherical_to_cartesian(r=self.EARTH_RADIUS, lat=np.deg2rad(la_), lon=np.deg2rad(lo_)))
                      for la_, lo_ in zip(lat, lon)])
        xx = X[:, 0]
        yy = X[:, 1]
        zz = X[:, 2]
    
        return xx, yy, zz, lon, lat
    
    
    def spherical_to_cartesian(self,r, lat, lon):
        """
        Transforms spherical coordinates into 3D cartesian ones
        :param r: radius
        :param lat: latitude in radians
        :param lon: longitude in radians
        :return: (x, y, z) - a tuple with the 3D coordinates
        """
        x = r * np.cos(lon) * np.cos(lat)
        y = r * np.sin(lon) * np.cos(lat)
        z = r * np.sin(lat)
        return x, y, z

