#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 15 11:30:58 2024

@author: melissacollier
"""
import networkx as nx
import pandas as pd
import os

# Set your working directory to a folder that has three subfolders: "edgelists", "matrices", and "graphmls"
directory = "C:\\Users\\s2607536\\OneDrive - University of Edinburgh\\ASNR\\New_networks"

# Create separate directories for matrices and edgelists
edgelist_directory = os.path.join(directory, "edgelists")
matrix_directory = os.path.join(directory, "matrices")
graphml_directory = os.path.join(directory, "graphmls")

# Ensure the graphml directory exists
os.makedirs(graphml_directory, exist_ok=True)

# Process edgelists
for filename in os.listdir(edgelist_directory):
    if filename == ".DS_Store":
        print("Mac error, skip")  # Ignore hidden files on Macs
        continue

    print(f"Processing edgelist: {filename}")
    f = os.path.join(edgelist_directory, filename)

    if os.path.isfile(f):
        df = pd.read_csv(f, index_col=0)  # Read in the edgelist
        col_names = df.columns

        # Check if there are at least 3 columns
        if len(col_names) >= 3:
            # Create the network
            G = nx.from_pandas_edgelist(df, col_names[0], col_names[1], col_names[2])
        elif len(col_names) == 2:
            # If there are only 2 columns (unweighted), create the network without weight
            G = nx.from_pandas_edgelist(df, col_names[0], col_names[1])
        else:
            print(f"Invalid format in file {filename}, skipping...")
            continue

        # Remove the .csv extension for the graphml name
        graphml_filename = filename.replace('.csv', '') + ".graphml"
        nx.write_graphml(G, os.path.join(graphml_directory, graphml_filename))

# Process adjacency matrices
for filename in os.listdir(matrix_directory):
    if filename == ".DS_Store":
        print("Mac error, skip")  # Ignore hidden files on Macs
        continue

    print(f"Processing adjacency matrix: {filename}")
    f = os.path.join(matrix_directory, filename)

    if os.path.isfile(f):
        df = pd.read_csv(f, index_col=0)
        col_names = df.columns

        # Ensure the index of the matrix is the same as the column names
        df.index = df.index.astype(str)
        df = df.loc[df.index.intersection(col_names), col_names]

        # Create network from adjacency matrix
        A = nx.from_pandas_adjacency(df)
        graphml_filename = filename.replace('.csv', '') + ".graphml"
        nx.write_graphml(A, os.path.join(graphml_directory, graphml_filename))
