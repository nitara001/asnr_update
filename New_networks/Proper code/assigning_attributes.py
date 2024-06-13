import pandas as pd
import networkx as nx
import os
attributes_directory = r"C:\Users\s2607536\OneDrive - University of Edinburgh\ASNR\New_networks\Attributes"
graphmls_directory = r"C:\Users\s2607536\OneDrive - University of Edinburgh\ASNR\New_networks\graphmls"
metadata_df = pd.read_csv("C:\\Users\\s2607536\\OneDrive - University of Edinburgh\\ASNR\\New_networks\\metadata.csv")

#function to convert attributes csv to a dictionary
# Function to convert attributes CSV to a dictionary
def csv_to_dict(csv_file):
    df = pd.read_csv(csv_file)
    df['id'] = df['id'].astype(str)  # Convert 'id' column to string
    df.set_index('id', inplace=True)  # Make 'id' the index
    network_attrs = df.to_dict(orient='index')
    return network_attrs


#function to set graph attributes
def add_attributes_to_graph(graph, attributes):
    nx.set_node_attributes(graph, attributes)


# Go through metadata for matching (metadata should contain name of network file (column 1) and name of attribute file (column 2))
for i, row in metadata_df.iterrows():
    graph_file = row['network']
    attribute_file = row['attributes']
    
    # Specify file paths of where the attribute and graph files are
    attribute_path = os.path.join(attributes_directory, f"{attribute_file}.csv")
    graphml_path = os.path.join(graphmls_directory, f"{graph_file}.graphml")
    
    # Read CSV file and convert to dictionary
    attribute_attrs = csv_to_dict(attribute_path)
    
    # Load GraphML file
    G = nx.read_graphml(graphml_path)
    
    # Add attributes to the graph
    add_attributes_to_graph(G, attribute_attrs)
    print(f"Attributes for {graph_file}:")
    for node, attr in G.nodes(data=True):
     print(f"Node {node}: {attr}")
    # Write the graph with attributes to a new GraphML file
    new_graphml_file = os.path.join("C:\\Users\\s2607536\\OneDrive - University of Edinburgh\\ASNR\\New_networks\\graphmls\\graphmls_with_attributes", f"{graph_file}_with_attributes.graphml")

    nx.write_graphml(G, new_graphml_file)


    

