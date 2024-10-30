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

# Define the subdirectories for edgelists, matrices, and graphml outputs
edgelist_directory = os.path.join(directory, "edgelists")
matrix_directory = os.path.join(directory, "matrices")
graphml_directory = os.path.join(directory, "graphmls")

# Ensure the graphml directory exists
os.makedirs(graphml_directory, exist_ok=True)

# Process edgelists
for filename in os.listdir(edgelist_directory):
    if filename.endswith('.csv'):  # Process only CSV files
        print(f"Processing edgelist: {filename}")
        filepath = os.path.join(edgelist_directory, filename)

        if os.path.isfile(filepath):
            try:
                # Read the CSV file into a DataFrame
                df = pd.read_csv(filepath)

                # Ensure at least 2 columns for source and target
                if len(df.columns) >= 2:
                    # If a 'weight' column exists, use it, else just the source-target
                    if 'weight' in df.columns:
                        G = nx.from_pandas_edgelist(df, source=df.columns[0], target=df.columns[1], edge_attr=True)
                    else:
                        G = nx.from_pandas_edgelist(df, source=df.columns[0], target=df.columns[1], edge_attr=df.columns[2:])

                    # Remove the .csv extension for the GraphML name
                    graphml_filename = filename.replace('.csv', '') + ".graphml"

                    # Save the graph to GraphML
                    nx.write_graphml(G, os.path.join(graphml_directory, graphml_filename))
                    print(f"GraphML saved as {graphml_filename}")

                else:
                    print(f"File {filename} does not have enough columns, skipping...")

            except Exception as e:
                print(f"Error processing file {filename}: {e}")

# Process adjacency matrices
for filename in os.listdir(matrix_directory):
    if filename.endswith(".csv") and filename != ".DS_Store":  # Avoid hidden files like .DS_Store
        print(f"Processing adjacency matrix: {filename}")
        filepath = os.path.join(matrix_directory, filename)

        if os.path.isfile(filepath):
            try:
                # Read the CSV as an adjacency matrix
                df = pd.read_csv(filepath, index_col=0)

                # Ensure the index matches the columns for adjacency matrix
                df.index = df.index.astype(str)
                df = df.loc[df.index.intersection(df.columns), df.columns]

                # Create the graph from the adjacency matrix
                A = nx.from_pandas_adjacency(df)

                # Remove the .csv extension for the GraphML name
                graphml_filename = filename.replace('.csv', '') + ".graphml"

                # Save the graph to GraphML
                nx.write_graphml(A, os.path.join(graphml_directory, graphml_filename))
                print(f"GraphML saved as {graphml_filename}")

            except Exception as e:
                print(f"Error processing file {filename}: {e}")
