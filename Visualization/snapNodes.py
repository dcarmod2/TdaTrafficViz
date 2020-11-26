import json
import networkx as nx
import numpy as np

bigG = nx.read_gpickle('subG.p')

with open("CycleNodeData_level_-3/CycleNodeData_frame_855_gen_5.json","r") as f:
    contents = json.load(f)

coords = [feature['geometry']['coordinates'] for feature in contents['features']]
#print(coords)
nodeL = list(bigG.nodes)
snapped_nodes = [nodeL[np.argmin([np.linalg.norm(np.array([bigG.nodes[x]['y'],bigG.nodes[x]['x']])-np.array([coord[1],coord[0]])) for x in nodeL])] for coord in coords]

snapped_coords = [[bigG.nodes[node]['y'],bigG.nodes[node]['x']] for node in snapped_nodes]
#print(snapped_coords)

nodefeatureList = [{'type':'Feature', 'properties': {'color':'green'},'geometry':{'type':'Point', 'coordinates':(x,y)}} for y,x in snapped_coords]

nodeFeatureCol = {'type':'FeatureCollection','features':nodefeatureList}

with open('BigCycleCoords.json','w+') as f:
    json.dump(nodeFeatureCol,f)

