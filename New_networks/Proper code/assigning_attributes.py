import pandas as pd
import networkx as nx
import os

# Directories for attributes and graphml files
attributes_directory = r"C:\Users\s2607536\OneDrive - University of Edinburgh\ASNR\New_networks\Attributes"
graphmls_directory = r"C:\Users\s2607536\OneDrive - University of Edinburgh\ASNR\New_networks\graphmls"
new_graphmls_directory = r"C:\Users\s2607536\OneDrive - University of Edinburgh\ASNR\New_networks\new_graphmls"
metadata_df = pd.read_csv(r"C:\Users\s2607536\OneDrive - University of Edinburgh\ASNR\New_networks\metadata.csv")

os.makedirs(new_graphmls_directory, exist_ok=True)

# Function to convert attributes CSV to a dictionary
def csv_to_dict(csv_file):
    df = pd.read_csv(csv_file)
    id_column = next((col for col in ['id', 'ID', 'node_id', 'Node', 'node'] if col in df.columns), None)
    if id_column is None:
        raise ValueError(f"No valid ID column found in {csv_file}.")
    df[id_column] = df[id_column].astype(str)
    df.set_index(id_column, inplace=True)
    return df.to_dict(orient='index')

# Process all GraphML files in the graphmls directory
for graph_file in os.listdir(graphmls_directory):
    if graph_file.endswith('.graphml'):
        graphml_path = os.path.join(graphmls_directory, graph_file)
        G = nx.read_graphml(graphml_path)

        # Check for matching attribute file
        attribute_file = metadata_df.loc[metadata_df['network'] == graph_file.replace('.graphml', ''), 'attributes']
        if not attribute_file.empty:
            attribute_path = os.path.join(attributes_directory, f"{attribute_file.values[0]}.csv")
            if os.path.exists(attribute_path):
                # Load and add attributes
                attribute_attrs = csv_to_dict(attribute_path)
                nx.set_node_attributes(G, attribute_attrs)
                print(f"Attributes added for {graph_file}.")
                
                # Save the graph with attributes
                new_graphml_file = os.path.join(new_graphmls_directory, graph_file)
                nx.write_graphml(G, new_graphml_file)
                print(f"New graph with attributes saved to {new_graphml_file}.")
            else:
                print(f"Attribute file {attribute_file.values[0]} does not exist.")
        else:
            print(f"No attributes listed for {graph_file}.")
