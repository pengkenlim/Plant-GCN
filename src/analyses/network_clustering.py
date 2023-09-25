import networkx as nx
from networkx.algorithms.community import label_propagation_communities


# Create an empty graph
G = nx.Graph()

# Add weighted edges
G.add_edge('A', 'B', weight=0.9)
G.add_edge('A', 'D', weight=0.9)
G.add_edge('D', 'B', weight=0.9)
G.add_edge('B', 'C', weight=0.9)
G.add_edge('A', 'C', weight=-1)
G.add_edge('C', 'D', weight=-1)




communities = label_propagation_communities(G)
for i, community in enumerate(communities):
    print(f"Community {i}: {community}")


from sklearn.cluster import FuzzyCMeans