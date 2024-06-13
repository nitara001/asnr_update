#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 15 11:30:58 2024

@author: melissacollier
"""
import networkx as nx
import pandas as pd
import os


# set your working directory to a folder that has three subfolders- one called "edgelists", 
# one called "matrice"s and one called "graphmls"

directory = "C:\\Users\\s2607536\\OneDrive - University of Edinburgh\\ASNR\\New_networks"

#create sep directories for matrices and edgelists
edgelist_directory = directory+"/edgelists"
matrix_directory = directory+"/matrices"

#first go through the edgelists
#the following code will work only for edgelists set up where column 1 = source node, column 2 = target node, column 3 = weight
for filename in os.listdir(edgelist_directory):
    if filename == ".DS_Store":
        print("Mac error, skip") # this is a weird thing with hidden files on Macs, ignore
    else:
        print(filename)
        f = os.path.join(edgelist_directory, filename)

        if os.path.isfile(f): #namke sure item is a file not a folder (if you have additional folders)
            df = pd.read_csv(f, index_col=0) # read in the edgelist
            col_names = df.columns

            G = nx.from_pandas_edgelist(df, col_names[0], col_names[1], col_names[2]) #create the network
            filename = filename.replace('.csv', '') #remove the .csv extension for the graphml name
            nx.write_graphml(G, directory+"/graphmls/"+filename+".graphml") # write the network as a graphml to folder called "graphmls"


#now for adjancenct matrices

for filename in os.listdir(matrix_directory):
    
    if filename == ".DS_Store":
        print("Mac error, skip") # this is a weird thing with hidden files on Macs, ignore
    else:
        print(filename)
        f = os.path.join(matrix_directory, filename)

        if os.path.isfile(f):
            df = pd.read_csv(f, index_col=0)
            col_names = df.columns
            adf = df.set_index([col_names]) #make sure the index of the matrix is the name as the column names

            A = nx.from_pandas_adjacency(adf) #create network from ad matrix
            filename = filename.replace('.csv', '')
            nx.write_graphml(A, directory+"/graphmls/"+filename+".graphml")

