import matplotlib
matplotlib.use('Agg')
import networkx as nx
import numpy as np
import os
import matplotlib.pyplot as plt
import pandas as pd
import itertools
import subprocess
import time as timemod
import re
import pickle


def condense_graph(subG,clusterlabelfile, num_levels, level,nodeL):
    with open(clusterlabelfile,"r+") as f:
        lines = f.readlines()
        line = lines[level]
        num_communs = int(re.findall(r'(\d+) nodes',line)[0])

    if level < 0:
        hier = num_levels + level
    else:
        hier = level
    with open("clusters.csv","w+") as f:
        subprocess.check_call([r"./DirectedLouvain/bin/hierarchy","subG.tree","-l",str(hier)], stdout=f)


    cluster_df = pd.read_csv('clusters.csv',index_col = False, names = ['node','community'],sep=' ')

    def convert_inds(ind):
        return nodeL[ind]

    cluster_df['node'] = cluster_df['node'].map(convert_inds)

    for i, tup in cluster_df.iterrows():
        node, comm = tup
        subG.nodes[node]['community'] = comm
        

    print("Condensing graph...")
    tic = timemod.clock()
    condensed_G = nx.DiGraph()
    for i in range(0,num_communs):
        condensed_G.add_node(i)
        nodes = list(cluster_df[cluster_df['community'] == i]['node'])
        condensed_G.nodes[i]['squished_nodes'] = nodes
        temp = subG.subgraph(nodes)
        distances = nx.floyd_warshall(temp,weight='pace')
        avg_pace = np.nanmean([t for src,tardict in distances.items() for dest,t in tardict.items() if src in nodes and dest in nodes and np.isfinite(t)])
        if np.isfinite(avg_pace):
            condensed_G.add_edge(i,i,pace=avg_pace)


    for start,end in itertools.permutations(condensed_G.nodes,2):
        clus1 = cluster_df[cluster_df['community']==start]
        clus2 = cluster_df[cluster_df['community']==end]
        temp = subG.subgraph(list(clus1['node'])+list(clus2['node']))
        distances = nx.floyd_warshall(temp,weight='pace')
        edges_12_pace = np.nanmean([t for src,tardict in distances.items() for dest,t in tardict.items() 
                                if src in set(clus1['node']) and dest in set(clus2['node']) and np.isfinite(t)])
        if np.isfinite(edges_12_pace):
            condensed_G.add_edge(start,end,pace=edges_12_pace) 
    toc = timemod.clock()
    print("Graph condensed in " + str(toc-tic) + " seconds.")
    return cluster_df, condensed_G

def condense_graph_spaths(subG,clusterlabelfile, num_levels, level, nodeL, savefile=None):
    with open(clusterlabelfile,"r+") as f:
        lines = f.readlines()
        line = lines[level]
        num_communs = int(re.findall(r'(\d+) nodes',line)[0])

    if level < 0:
        hier = num_levels + level
    else:
        hier = level
    with open("clusters.csv","w+") as f:
        subprocess.check_call([r"./DirectedLouvain/bin/hierarchy","subG.tree","-l",str(hier)], stdout=f)


    cluster_df = pd.read_csv('clusters.csv',index_col = False, names = ['node','community'],sep=' ')

    def convert_inds(ind):
        return nodeL[ind]

    cluster_df['node'] = cluster_df['node'].map(convert_inds)

    for i, tup in cluster_df.iterrows():
        node, comm = tup
        subG.nodes[node]['community'] = comm
        

    print("Condensing graph...")
    tic = timemod.clock()
    condensed_G = nx.DiGraph()
    for i in range(0,num_communs):
        condensed_G.add_node(i)
        nodes = list(cluster_df[cluster_df['community'] == i]['node'])
        condensed_G.nodes[i]['squished_nodes'] = nodes
        temp = subG.subgraph(nodes)
        distances = nx.floyd_warshall(temp,weight='pace')
        paces = [(t,src,dest) for src,tardict in distances.items() for dest,t in tardict.items() if src in nodes and dest in nodes and np.isfinite(t)]
        avg_pace = np.nanmean([x[0] for x in paces])
        
        if np.isfinite(avg_pace):
            condensed_G.add_edge(i,i,pace=avg_pace)
            quickest_pace_ind = np.argmin([x[0] for x in paces])
            fast_src,fast_dest = paces[quickest_pace_ind][1],paces[quickest_pace_ind][2]
            condensed_G.edges[(i,i)]['spath'] = nx.dijkstra_path(subG,fast_src,fast_dest,weight='pace')

    #set positions
    pos = {node:(np.mean([subG.nodes[cnode]['y'] for cnode in cluster_df[cluster_df['community']==node]['node']]),np.mean([subG.nodes[cnode]['x'] for cnode in cluster_df[cluster_df['community']==node]['node']])) for node in condensed_G.nodes}
    nx.set_node_attributes(condensed_G,pos,'pos')

    for start,end in itertools.permutations(condensed_G.nodes,2):
        clus1 = cluster_df[cluster_df['community']==start]
        clus2 = cluster_df[cluster_df['community']==end]
        temp = subG.subgraph(list(clus1['node'])+list(clus2['node']))
        distances = nx.floyd_warshall(temp,weight='pace')
        edges_12 =[(t,src,dest) for src,tardict in distances.items() for dest,t in tardict.items() if src in set(clus1['node']) and dest in set(clus2['node']) and np.isfinite(t)]
        edges_12_pace = np.nanmean([x[0] for x in edges_12])
        
        if np.isfinite(edges_12_pace):
            condensed_G.add_edge(start,end,pace=edges_12_pace)
            #get closest actual nodes to centroid
            sind = np.argmin([np.linalg.norm(np.array([subG.nodes[x[1]]['y'],subG.nodes[x[1]]['x']])-np.array(condensed_G.nodes[start]['pos'])) for x in edges_12])
            qsrc = edges_12[sind][1]
            eind = np.argmin([np.linalg.norm(np.array([subG.nodes[x[2]]['y'],subG.nodes[x[2]]['x']])-np.array(condensed_G.nodes[end]['pos'])) for x in edges_12])
            qsrc = edges_12[eind][2]
            qsrc,qdest = edges_12[sind][1],edges_12[sind][2]
            condensed_G.edges[(start,end)]['spath'] = nx.dijkstra_path(subG,qsrc,qdest,weight = 'pace')

    
    if savefile:
        nx.write_gpickle(condensed_G,savefile)
    toc = timemod.clock()
    print("Graph condensed in " + str(toc-tic) + " seconds.")
    return cluster_df, condensed_G

def get_largest_comp(filename):
    G = nx.read_gpickle(filename)
    bad_edges = [edge for edge in G.edges if np.isinf(G.edges[edge]['time'])]
    G.remove_edges_from(bad_edges)
    biggest_comp = max(nx.strongly_connected_components(G),key=len)
    
    subgraph = G.subgraph(biggest_comp)
    return subgraph


def get_cost(G,path,attr):
    cost = 0
    for link in path:
        cost += G.edges[link][attr]
    return cost


def safemax(listlike):
    try:
        return max(listlike)
    except ValueError:
        return -1

class chain_toggle_reduce:
    def __init__(self,simplexCoeff,p=2,filtration=None,mark=0,ind=None,chain_ptr=None):
        self.p = p
        self.coeffAugList = simplexCoeff
        self.max_dim = safemax(len(x[1]) for x in self.coeffAugList)-1
        self.mark = mark
        self.filtration = filtration
        self.ind = ind
        self.chain_ptr = chain_ptr
        self.operated_upon = simplexCoeff.copy()
        self.label = repr(self.coeffAugList)
        return
    
    def is_empty(self):
        return self.coeffAugList == []
    
    def get_coeff_inv(self,simpl):
        coeff = [simp.coeffAugList[0][0] for simp in self if simp == simpl][0]
        return (coeff**(self.p-2)) % self.p
        
    def __bdry__(self):
        def simp_bdry(coeffAug):
            coeff,simp = coeffAug
            return [((-1)**i * coeff,simp[:i]+simp[i+1:]) for i in range(0,len(simp))]
            
        return chain_toggle_reduce(sum([simp_bdry(x) for x in self.coeffAugList],[]),
                         p=self.p,filtration=self.filtration)
            
    def reduce(self):
        seenDict = dict()
        for coeff,simp in self.coeffAugList:
            try:
                seenDict[tuple(simp)] += coeff
            except Exception as e:
                seenDict[tuple(simp)] = coeff
        almost_reduced = [(val % self.p,list(key)) for key,val in seenDict.items() if val % self.p != 0]
        reduced = [(coeff,simp) for coeff,simp in almost_reduced if simp != []]
        return chain_toggle_reduce(reduced,self.p,self.filtration,self.mark,self.ind,self.chain_ptr)
            
    def __add__(self,other):
        return chain_toggle_reduce(self.coeffAugList + other.coeffAugList,self.p,self.filtration)
    
    def __rmul__(self,const):
        self.coeffAugList = [((const*coeff) % self.p,simp) for coeff,simp in self.coeffAugList]
        return self
    
    def __lt__(self,other):
        return (len(self.coeffAugList),self.filtration) <= (len(other.coeffAugList),other.filtration)
    
    def __eq__(self,other):
        
        if len(self.coeffAugList) != len(other.coeffAugList):
            return False
        
        
        for i,coeffAug in enumerate(self.coeffAugList):
            if coeffAug[1] != other.coeffAugList[i][1]:
                return False
            
        return True
    
    def __getitem__(self,index):
        return chain_toggle_reduce([self.coeffAugList[index]],self.p,self.filtration)


def get_simplicesGen_min(dist_mat,filts,dim,p):
    simplices = []
    for comb in itertools.combinations(range(0,len(dist_mat)),dim+1):
        relevant_dist = dist_mat[tuple(np.meshgrid(comb,comb))]
        lower = np.tril(relevant_dist)
        upper = np.triu(relevant_dist)
        min_arr = np.minimum(lower,upper.T)
        max_val = np.max(min_arr.flatten())
        for filt in filts:
            if filt > max_val:
                simplices.append(chain_toggle_reduce([(1,list(comb))],p,filt))
                break
        
    return simplices

def get_simplicesGen(dist_mat,filts,dim,p):
    simplices = []
    for comb in itertools.combinations(range(0,len(dist_mat)),dim+1):
        max_val = max(dist_mat[tuple(np.meshgrid(comb,comb))].flatten())
        for filt in filts:
            if filt > max_val:
                simplices.append(chain_toggle_reduce([(1,list(comb))],p,filt))
                break
        
    return simplices

def removePivotRowsGen(simps,simp):
    #cycle = [simp.coeffAugList]
    bdry = simp.__bdry__()
    if bdry.reduce().is_empty():
        return bdry.reduce()
    
    inds = [simps.index(simp) for simp in bdry]
    marked_bdry = chain_toggle_reduce([simps[ind].coeffAugList[0] for ind in inds if simps[ind].mark == 1],2)
    
    while not marked_bdry.reduce().is_empty():
        marked_inds = [simps.index(simp) for simp in marked_bdry.reduce()]
        max_ind = max(marked_inds)
        
        if not simps[max_ind].ind:
            break
        else:
            simp.operated_upon = (chain_toggle_reduce(simp.operated_upon) + (-1)*simps[max_ind].chain_ptr.get_coeff_inv(simps[max_ind]) * chain_toggle_reduce(simps[simps[max_ind].ind].operated_upon)).coeffAugList
            #cycle.append(simps[simps[max_ind].ind].coeffAugList)
            marked_bdry = marked_bdry + (-1)*simps[max_ind].chain_ptr.get_coeff_inv(simps[max_ind])*simps[max_ind].chain_ptr
        
            
    return marked_bdry.reduce()

def computeIntervalsGen(simps):
    dim = len(simps[-1].coeffAugList[0][1])
    interDict = {i : dict() for i in range(0,dim)}
    genDict = {i: dict() for i in range(0,dim)}
    for j in range(0,len(simps)):
        red_bdry = removePivotRowsGen(simps,simps[j])
        if red_bdry.is_empty():
            simps[j].mark = 1
            k = len(simps[j].coeffAugList[0][1])-1
            genDict[k].update({simps[j].label: (simps[j].filtration, chain_toggle_reduce(simps[j].operated_upon.copy()).reduce().coeffAugList)}) 
        else:
            max_ind = max(simps.index(simp) for simp in red_bdry)
            k = len(simps[max_ind].coeffAugList[0][1])-1
            simps[max_ind].ind = j
            simps[max_ind].chain_ptr = red_bdry
            #if simps[max_ind].filtration != simps[j].filtration:
            interDict[k].update({simps[max_ind].label : (simps[max_ind].filtration,simps[j].filtration)})
            
                
    for j in range(0,len(simps)):
        if simps[j].mark and not simps[j].ind:
            k = len(simps[j].coeffAugList[0][1])-1
            interDict[k].update({simps[j].label : (simps[j].filtration,np.inf)})
    return interDict,genDict
