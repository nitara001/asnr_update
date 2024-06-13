import networkx as nx
import community
import fnmatch
import os
import numpy as np
import csv
import zipfile
import pandas as pd
#####################################################################333
os.chdir('/Users/melissacollier/Library/CloudStorage/GoogleDrive-mac532@georgetown.edu/My Drive/Melissa_Shweta_Files/ASNR/Network_Cleaning_Code')
#########################################################################

qd = pd.read_csv("qualitative_data.csv")
directory = "/Users/melissacollier/Library/CloudStorage/GoogleDrive-mac532@georgetown.edu/My Drive/Melissa_Shweta_Files/ASNR/Network_Cleaning_Code"
g_directory = "/Users/melissacollier/Library/CloudStorage/GoogleDrive-mac532@georgetown.edu/My Drive/Melissa_Shweta_Files/ASNR/Network_Cleaning_Code/new_graphmls"

#first go through the edgelists
#the following code will work only for edgelists set up where column 1 = source node, column 2 = target node, column 3 = weight


#######################################
def calculate_avg_wd(G, partition, n_nodes):

	wdlist = []
	for node1 in G.nodes():
		nbrs = G.neighbors(node1)
		mod1 = partition[node1]
		mod_nbrs = [node2 for node2 in nbrs if partition[node2]==mod1]
		wd = len(mod_nbrs)
		wdlist.append(wd) 
		
	return sum(wdlist)/(1.*n_nodes)
	
#######################################	
#######################################
def calculate_Qmax(G, mod_nodes):


	Lt= sum([G.degree(node) for node in G.nodes()])
	total  =0
	
	for mod in mod_nodes.keys():
		Lk = sum([G.degree(node) for node in mod_nodes[mod]])
		total+= (1.0*Lk/Lt) - (1.0*Lk/Lt)**2 
		
		

	return total

########################################################################
def calculate_jaccard(g1, g2):

	edges1 = g1.edges()
	edges2 = g2.edges()
	w11 = len(list(set(edges1) & set(edges2)))
	w10 = len(list(set(edges1) - set(edges2)))
	w01 =  len(list(set(edges2) - set(edges1)))

	ratio = w11/ (1.*(w11+w10+w01))

	return ratio
#################################################################################
def normalize_edge_weight(Graw):

	totwt = sum(nx.get_edge_attributes(Graw, 'weight').values())
	print ("tot wt = "), totwt
	if totwt!=1:
		G = nx.Graph()
		G.add_nodes_from(Graw.nodes())
	
		for (n1, n2), wt in nx.get_edge_attributes(Graw, 'weight').items():
			G.add_edge(n1,n2)
			G[n1][n2]['weight'] = wt
	
	else: G = Graw.copy()

	return G
	
################################################################################
rows = []
num = 1
for filename in os.listdir(g_directory):

    if filename == ".DS_Store":
        print("Mac error, skip") # this is a weird thing with hidden files on Macs, ignore
    else:
        if filename.endswith(".graphml"):
            
            ### get qualitative data
            qdata= qd.loc[qd['filename']== filename]
            dirname = qdata['dirname'].iloc[0]
            clas = qdata['class'].iloc[0]
            genus = qdata['genus'].iloc[0]
            species = qdata['species'].iloc[0]
            interaction_type = qdata['interaction_type'].iloc[0]
            definition_of_interaction = qdata['definition_of_interaction'].iloc[0]
            edge_wt_type = qdata['edge_wt_type'].iloc[0]
            geographical_location = qdata['geographical_location'].iloc[0]
            data_record_technique = qdata['data_record_technique'].iloc[0]
            time_span = qdata['time_span'].iloc[0]
            resolution = qdata['resolution'].iloc[0]
            time_span_per_day = qdata['time_span_per_day'].iloc[0]
            attributes_available = qdata['attributes_available'].iloc[0]
            population_type = qdata['population_type'].iloc[0]
            citation= qdata['Citation'].iloc[0]
            
            
            f = os.path.join(g_directory, filename)                                #namke sure item is a file not a folder (if you have additional folders)    
	
            Graw= nx.read_graphml(f)
            Graw_new = nx.read_graphml(f)
            Graw.remove_edges_from(nx.selfloop_edges(Graw, data= True))
            n_nodes = len(Graw.nodes())
            n_edges = len(Graw.edges())
            
            
        	#####
        	## if network does not have weights, add a weight of one to all edges
            if len(nx.get_edge_attributes(Graw,'weight'))==0:
                for (n1,n2) in Graw.edges(): Graw[n1][n2]['weight']=1
        	####
        
        	####################################################
        	#if no edges then return NAs
            if n_edges==0:
                elements = [num, filename, "NA", n_nodes, n_edges] + ["NA"]*18
                writer.writerow(elements)
                continue
        	########################################################
        	## remove edges with weight zero
            for (n1, n2) in Graw.edges():
                if Graw[n1][n2]['weight']==0: 
                    print("Removing NULL edge!!!"), (n1, n2), len(Graw.edges()),
                    Graw_new.remove_edge(n1,n2)
                    print(len(Graw_new.edges()))
        	##########################################################
        	##normalize edge weight	
        
            G = normalize_edge_weight(Graw_new)	
           
           ##########
           ##Running this bc no weights and I just need the outputs for these graphs
            #G = Graw_new
            #print(G.degree())
            
        	#########################################################
            is_connect = nx.is_connected(G) 
            num_comp = nx.number_connected_components(G)
            giant_cc = max(list((G.subgraph(c).copy() for c in nx.connected_components(G))))
            #giant_cc = max(nx.connected_component_subgraphs(G), key=len)
            largest_cc = len(giant_cc)
            num_isolates = len(list(nx.isolates(G)))
            density = round(nx.density(G),3)
            deg_list = [G.degree(node) for node in G.nodes()]
            avg_deg = round(np.mean(deg_list),3)
            std_deg = round(np.std(deg_list),3)
        	
            if avg_deg>0:cv_deg = round(float(std_deg)/avg_deg,3)
            else: cv_deg = "NA"
            if std_deg>0:
                first_part = float(n_nodes)/((n_nodes-1)*(n_nodes-2))
                sec_part = sum([(float(deg - avg_deg)/std_deg)**3 for deg in deg_list])
                skew = round(first_part*sec_part,3)
            else: skew="NA"
            wt_deg_list = [G.degree(node, weight="weight") for node in G.nodes()]
            wt_edge_list = [G[n1][n2]['weight'] for (n1,n2) in G.edges()]
            node_stg = round(np.mean(wt_deg_list),3)
            edge_stg = round(np.mean(wt_edge_list),3)
            std_node_stg = round(np.std(wt_deg_list),3)
            std_edge_stg = round(np.std(wt_edge_list),3)
            CV_node_stg = round(float(std_node_stg)/node_stg,3)
            CV_edge_stg = round(float(std_edge_stg)/edge_stg,3)
            high_node_stg = round(max(wt_deg_list),3)
            asrt = round(nx.degree_assortativity_coefficient(G),3)
        	#asrt_wt = round(nx.degree_assortativity_coefficient(G, weight="weight"),3)
            betw_list  = list(nx.betweenness_centrality(G).values())
            betw_wt_list  = list(nx.betweenness_centrality(G, weight="weight").values())
            avg_betw =  round(np.mean(betw_list) ,3)
            std_betw =  round(np.std(betw_list) ,3)
            betw_wt =  round(np.mean(betw_wt_list) ,3)
            std_betw_wt =  round(np.std(betw_wt_list) ,3)
            highest_betw = max(betw_list)
            highest_betw_wt = max(betw_wt_list)
            clstr = round(nx.average_clustering(G),3)
            clstr_wt = round(nx.average_clustering(G, weight="weight"),3)
            trans = round(nx.transitivity(G),3)
        	#avg_jaccard = calculate_avg_jaccard(num, G)
        	#for the rest of the computations, network is required to be connected
            if not nx.is_connected(G): G = max(nx.connected_component_subgraphs(G), key=len)
            G1 = nx.Graph()
            G1.add_nodes_from(G.nodes())
            G1.add_edges_from(G.edges())
        		
            partition = community.best_partition(G1)
            Q = round(community.modularity(partition, G1),3)
            modules = list(set(partition.values()))
            mod_nodes= {}
            for mod in modules: mod_nodes[mod] = [node for node in G1.nodes() if partition[node]==mod]
            Qmax = round(calculate_Qmax(G1, mod_nodes),3)
            coh = calculate_avg_wd(G1, partition, n_nodes)/(1.*avg_deg)
            diam = nx.diameter(G)
            avg_modsize = float(len(G1.nodes()))/len(modules)
            if Qmax>0:Qrel = round(float(Q)/Qmax,3)
            else: Qrel="NA"
        	
                #This part will write a new graphml file that has removed self loops and NULL edges
            
            nx.write_graphml(G, directory+'/clean_graphmls/'+filename)
        	
            
            #This writes the row in your new csv file
            elements = [num, dirname, filename, clas, genus, species, interaction_type, definition_of_interaction, edge_wt_type, geographical_location, population_type, data_record_technique, 
                        time_span, resolution, time_span_per_day, attributes_available, is_connect,num_comp, n_nodes, n_edges, density, max(deg_list), avg_deg, std_deg, cv_deg, skew, node_stg, std_node_stg, CV_node_stg, edge_stg, std_edge_stg, CV_edge_stg, high_node_stg, asrt,  avg_betw, std_betw, highest_betw, betw_wt, std_betw_wt, highest_betw_wt, clstr, clstr_wt, trans,
                        Q, Qmax, Qrel, coh , len(modules), avg_modsize, diam, citation]
            rows.append(elements)
            print(elements)
            num = num+1
            
df = pd.DataFrame(rows, columns = ["Graph#", "dirname", "filename", "class", "genus", "species", "interaction_type", "definition_of_interaction", "edge_wt_type", "geographical_location", "population_type", "data_record_technique", "time_span", "resolution", "time_span_per_day", "attributes_available",
                 "is_connected", "num_components", "nodes", "edges", "network.density", "highest.degree", "avg.degree","std.degree.node",  "CV.degree", "skewness", "avg.node.strength", "std.node.strength",  "CV.node.strength", "avg.edge.strength", "std.edge.strength",  "CV.edge.strength",  "highest.node.strength", "deg.assort",  
                 "avg.betw.centrality", "std.betw.centrality", "highest.betw" ,  "avg.betw.centrality.weighted", "std.betw.centrality.weighted", "highest.betw.weighted" , "clustering", "clustering.weighted", "transitivity", "Q", "Qmax", "Qrel", "cohesion", "num.modules", "avg.modsize", "diameter", "Citation"]) 

df.to_csv("Network_Summary.csv")