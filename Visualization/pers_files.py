import matplotlib
matplotlib.use('Agg')
import networkx as nx
import numpy as np
import os
import matplotlib.pyplot as plt
import pandas as pd
import itertools
from matplotlib import animation
import time as timemod
import re
import pickle
import json
import shapely.geometry



def num_display(num):
    return "{0:.0f}".format(num)

def convert_val(val):
    if np.isinf(val):
        return 'inf'
    else:
        return np.round(val)

def __min_birth_max_death(persistence, band=0.0):
    """This function returns (min_birth, max_death) from the persistence.

    :param persistence: The persistence to plot.
    :type persistence: list of tuples(dimension, tuple(birth, death)).
    :param band: band
    :type band: float.
    :returns: (float, float) -- (min_birth, max_death).
    """
    # Look for minimum birth date and maximum death date for plot optimisation
    max_death = 0
    min_birth = persistence[0][1][0]
    for interval in reversed(persistence):
        if float(interval[1][1]) != float("inf"):
            if float(interval[1][1]) > max_death:
                max_death = float(interval[1][1])
        if float(interval[1][0]) > max_death:
            max_death = float(interval[1][0])
        if float(interval[1][0]) < min_birth:
            min_birth = float(interval[1][0])
    if band > 0.0:
        max_death += band
    return (min_birth, max_death)


"""
Only 13 colors for the palette
"""
palette = [
    "#ff0000",
    "#00ff00",
    "#0000ff",
    "#00ffff",
    "#ff00ff",
    "#ffff00",
    "#000000",
    "#880000",
    "#008800",
    "#000088",
    "#888800",
    "#880088",
    "#008888",
]


def plot_persistence_barcode_dan(
    persistence=[],
    ax = None,
    persistence_file="",
    alpha=0.6,
    max_intervals=1000,
    max_barcodes=1000,
    inf_delta=0.1,
    inter_ind=None,
    legend=False,
    minbirth_maxdeath = None,
    color_pal = None,
    **kwargs
):
    """This function plots the persistence bar code from persistence values list
    or from a :doc:`persistence file <fileformats>`.

    :param persistence: Persistence intervals values list grouped by dimension.
    :type persistence: list of tuples(dimension, tuple(birth, death)).
    :param ax: Axis on which to plot barcode.
    :type ax: pyplot axes instance.
    :param persistence_file: A :doc:`persistence file <fileformats>` style name
        (reset persistence if both are set).
    :type persistence_file: string
    :param alpha: barcode transparency value (0.0 transparent through 1.0
        opaque - default is 0.6).
    :type alpha: float.
    :param max_intervals: maximal number of intervals to display.
        Selected intervals are those with the longest life time. Set it
        to 0 to see all. Default value is 1000.
    :type max_intervals: int.
    :param inf_delta: Infinity is placed at :code:`((max_death - min_birth) x
        inf_delta)` above :code:`max_death` value. A reasonable value is
        between 0.05 and 0.5 - default is 0.1.
    :type inf_delta: float.
    :param legend: Display the dimension color legend (default is False).
    :type legend: boolean.
    :returns: A matplotlib object containing horizontal bar plot of persistence
        (launch `show()` method on it to display it).
        
    kwargs is passed to ax.set
    """
    try:
        import matplotlib.pyplot as plt
        import matplotlib.patches as mpatches

        if persistence_file is not "":
            if path.isfile(persistence_file):
                # Reset persistence
                persistence = []
                diag = read_persistence_intervals_grouped_by_dimension(
                    persistence_file=persistence_file
                )
                for key in diag.keys():
                    for persistence_interval in diag[key]:
                        persistence.append((key, persistence_interval))
            else:
                print("file " + persistence_file + " not found.")
                return None

        if max_barcodes is not 1000:
            #print("Deprecated parameter. It has been replaced by max_intervals")
            max_intervals = max_barcodes

        if max_intervals > 0 and max_intervals < len(persistence):
            # Sort by life time, then takes only the max_intervals elements
            persistence = sorted(
                persistence,
                key=lambda life_time: life_time[1][1] - life_time[1][0],
                reverse=True,
            )[:max_intervals]

        persistence = sorted(persistence, key=lambda birth: birth[1][0])

        (min_birth, max_death) = minbirth_maxdeath or __min_birth_max_death(persistence)
        ind = 0
        delta = (max_death - min_birth) * inf_delta
        # Replace infinity values with max_death + delta for bar code to be more
        # readable
        infinity = max_death + delta
        axis_start = min_birth - delta
        # Draw horizontal bars in loop
        for inner_ind,interval in enumerate(reversed(persistence)):
            if color_pal:
                color = color_pal[inner_ind]
            else:
                if not inter_ind:
                    color = palette[interval[0]]
                elif inter_ind == inner_ind or inner_ind == len(persistence)+inter_ind:
                    color = palette[-1]
                else:
                    color = palette[interval[0]]
                
            if float(interval[1][1]) != float("inf"):
                # Finite death case
                if ax:
                    ax.barh(
                    ind,
                    (interval[1][1] - interval[1][0]),
                    height=0.8,
                    left=interval[1][0],
                    alpha=alpha,
                    color=color,
                    linewidth=0,
                )
                else:
                    plt.barh(
                        ind,
                        (interval[1][1] - interval[1][0]),
                        height=0.8,
                        left=interval[1][0],
                        alpha=alpha,
                        color=color,
                        linewidth=0,
                    )
            else:
                # Infinite death case for diagram to be nicer
                if ax:
                    ax.barh(
                        ind,
                        (infinity - interval[1][0]),
                        height=0.8,
                        left=interval[1][0],
                        alpha=alpha,
                        color=color,
                        linewidth=0,
                    )
                else:
                    plt.barh(
                        ind,
                        (infinity - interval[1][0]),
                        height=0.8,
                        left=interval[1][0],
                        alpha=alpha,
                        color=color,
                        linewidth=0,
                    )
                
            ind = ind + 1

        if legend:
            dimensions = list(set(item[0] for item in persistence))
            plt.legend(
                handles=[
                    mpatches.Patch(color=palette[dim], label=str(dim))
                    for dim in dimensions
                ],
                loc="lower right",
            )
        ax.set(**kwargs)
        # Ends plot on infinity value and starts a little bit before min_birth
        if ax:
            ax.axis([axis_start, infinity, 0, ind])
        else:
            plt.axis([axis_start, infinity, 0, ind])
            
        return plt

    except ImportError:
        print("This function is not available, you may be missing matplotlib.")

def convert_to_gd_format(pers_list,gen_list,max_dimension=1):
    combined = [[(i, inter),gen_list[i][label]] for i in pers_list for label,inter in pers_list[i].items() if i <= max_dimension]
    return [x[0] for x in combined], [x[1] for x in combined]

def nodify(genList,ref_graph):
    nodeL = list(ref_graph.nodes)
    return [(filt,[tuple([nodeL[ind] for ind in coeffAug[1]]) for coeffAug in chain]) for filt,chain in genList]

def make_graph(simps,filtration,ref_graph):
    G = nx.Graph()
    nodes = [get_corresponding_nodes_from_simp(x,ref_graph)[0] for x in simps if x.filtration <= filtration and x.max_dim == 0]
    edges = [get_corresponding_nodes_from_simp(x,ref_graph) for x in simps if x.filtration <= filtration and x.max_dim == 1]
    
    G.add_nodes_from(nodes)
    G.add_edges_from(edges)
    return G

def get_corresponding_node_from_gen(genList,graph,dim,max_filt):
    nodeL = list(graph.nodes)
    return [[tuple([nodeL[y] for y in x[1]]) for x in dat[1]]  for label,dat in genList[dim].items() if dat[0] < max_filt]

def get_corresponding_nodes_from_simp(simp,graph):
    nodeL = list(graph.nodes)
    return [nodeL[x] for x in simp.coeffAugList[0][1]]


def ani_frame_files(filename,inter_pruned,gen_pruned,G,pos,filtrations,bigG,level,simps,dpi,plot_bars=False,bar_file=None,plot_cycles=False,cycle_file=None,cycle_type='pers',fontsize=30,figsize=(12,12)):
    plt.rcParams.update({'font.size': fontsize})
    num_rows = 3
    fig = plt.figure(figsize=(20,20), constrained_layout=False)
    outer_grid = fig.add_gridspec(num_rows,2)
    ax00 = fig.add_subplot(outer_grid[0,0])
    ax01 = fig.add_subplot(outer_grid[0,1])
    ax10 = fig.add_subplot(outer_grid[1,0])
    ax11 = fig.add_subplot(outer_grid[1,1])
    ax21 = fig.add_subplot(outer_grid[2,1])

    axs = fig.get_axes()
    for ax in axs:
            #ax.set_aspect('equal')
            #ax.get_xaxis().set_visible(False)
        ax.get_yaxis().set_visible(False)
        if ax.is_first_col():
            ax.get_xaxis().set_visible(False)


    #fig.set_size_inches([12,12])


    fig.set_tight_layout(True)
    pers, geners = convert_to_gd_format(inter_pruned,gen_pruned)
    generators = nodify(geners,G)
    data = [(filtration , gen, ind) for filtration in filtrations 
           for ind,gen in 
               enumerate(get_corresponding_node_from_gen(gen_pruned,G,1,max_filt=filtration))
           ]
    minb,maxd = __min_birth_max_death(pers)
    
    def threshmax(val,maxval):
        if val <= maxval:
            return val
        else:
            return maxval
    
    #Use for displaying multiple H_1 gens
    num_gens = 0

    def update_img(mfilt_tup):
        filt_ind, mfilt = mfilt_tup
        relevant_inds = [i for i,gen in enumerate(generators) if gen[0] <= mfilt]
        modified_pers = [(dim_minmax[0],(dim_minmax[1][0],threshmax(dim_minmax[1][1],mfilt))) for ind,dim_minmax in enumerate(pers)
                        if ind in relevant_inds]
        relevant_gens = [gen for filt,gen in generators if filt <= mfilt]

        filteredG = make_graph(simps,mfilt,G)
        nx.set_node_attributes(filteredG,pos,'pos')
        spaths = {edge: G.edges[edge]['spath'] for edge in G.edges}
        for edge in filteredG.edges:
            spath_in_condensed = nx.dijkstra_path(G,edge[0],edge[1],weight='pace')
            actual_spath = []
            for i,src in enumerate(spath_in_condensed[:-2]):
                actual_spath += spaths[(src,spath_in_condensed[i+1])]
                actual_spath += nx.dijkstra_path(bigG,actual_spath[-1],spaths[(spath_in_condensed[i+1],spath_in_condensed[i+2])][0])
            actual_spath += spaths[(spath_in_condensed[-2],spath_in_condensed[-1])]
            filteredG.edges[edge]['spath'] = actual_spath

        # create rips graph features

        for node in filteredG.nodes:
            filteredG.nodes[node]['x'] = G.nodes[node]['pos'][1]
            filteredG.nodes[node]['y'] = G.nodes[node]['pos'][0]

        coord_list = [(filteredG.nodes[node]['y'], filteredG.nodes[node]['x']) for node in filteredG.nodes]
        nodefeatureList = [{'type':'Feature', 'properties':{},'geometry':{'type':'Point','coordinates': (x,y)}} for y,x in coord_list]

        nodeFeatureCol = {'type':'FeatureCollection','features': nodefeatureList}
        with open('RipsGraphNodes_level_'+str(level)+'/RipsGraphNodes_frame_'+str(filt_ind)+'.json','w+') as f:
            json.dump(nodeFeatureCol,f)

        edgefeatureList = []


        for edge in filteredG.edges:
            try:
                edgefeatureList.append({'type': 'Feature', 'properties': {}, 'geometry': shapely.geometry.mapping(filteredG.edges[edge]['geometry'])})
            except KeyError:
                source,tar = edge
                lats,lons = filteredG.nodes[source]['y'], filteredG.nodes[source]['x']
                latt,lont = filteredG.nodes[tar]['y'], filteredG.nodes[tar]['x']
                edgefeatureList.append({'type': 'Feature', 'properties': {}, 'geometry': {'type':'LineString', 'coordinates': ((lons,lats),(lont,latt))}})

        edgeFeatureCol = {'type':'FeatureCollection','features': edgefeatureList}
        with open('RipsGraphEdges_level_'+str(level)+'/RipsGraphEdges_frame_'+str(filt_ind)+'.json','w+') as f:
            json.dump(edgeFeatureCol,f)
       
        ax00.clear()
        nx.draw_networkx(filteredG,pos,node_size=20,with_labels=False,ax=ax00)

        ax01.clear()
        plot_persistence_barcode_dan(modified_pers,ax01,title='Total Barcode',xlabel='Pace Filtration',ylabel='Generators',
                            inter_ind=None,minbirth_maxdeath = (minb,maxd))

        # save a separate copy of total barcode
        tempfig,tempax = plt.subplots(figsize=figsize)
        tempax.set_yticks([])
        plot_persistence_barcode_dan(modified_pers,tempax,title='Total Barcode',xlabel='Pace Filtration',ylabel='Generators',
                            inter_ind=None,minbirth_maxdeath = (minb,maxd))

        plt.savefig('TotPers_level_'+str(level)+'/TotPers_frame_'+str(filt_ind)+'.png',bbox_inches='tight')
        plt.close(tempfig)


        zero_pers = [x for x in modified_pers if x[0] == 0]
        comps = [list(x) for x in nx.connected_components(filteredG)]
        relevant_gens_zero = [gen for gen in relevant_gens if len(gen[0]) == 1]
        color_pal = list(reversed(['#%02x%02x%02x'%(20*i % 255,30**i % 255, 40*i % 255) for i,_ in enumerate(zero_pers)]))


        ax10.clear()
        for i,gen in enumerate(relevant_gens_zero):
            comp = [x for x in comps if gen[0][0] in x][0]
            subgraph = filteredG.subgraph(comp)
            temp = subgraph.copy()
            temp.remove_edges_from(subgraph.edges)
            nx.draw_networkx(temp,pos,node_size=200,with_labels=False,ax=ax10,node_color='#%02x%02x%02x'%(20*i % 255,30**i % 255, 40*i % 255))
            coord_list = [(temp.nodes[node]['y'], temp.nodes[node]['x']) for node in temp.nodes]
            nodefeatureList = [{'type':'Feature', 'properties':{'color':'#%02x%02x%02x'%(20*i % 255,30**i % 255, 40*i % 255)},'geometry':{'type':'Point','coordinates': (x,y)}} for y,x in coord_list]

            nodeFeatureCol = {'type':'FeatureCollection','features': nodefeatureList}
            with open('H0gens_level_'+str(level)+'/H0gens_frame_'+str(filt_ind)+'_gen_'+str(i)+'.json','w+') as f:
                json.dump(nodeFeatureCol,f)

        
            

        # 
        ax11.clear()
        plot_persistence_barcode_dan(zero_pers,ax11,title='$H_0$ Barcode',xlabel='Pace Filtration',ylabel='Generators',inter_ind=None,minbirth_maxdeath = (minb,maxd),color_pal=color_pal)

        tempfig,tempax = plt.subplots(figsize=figsize)
        tempax.set_yticks([])
        plot_persistence_barcode_dan(zero_pers,tempax,title='$H_0$ Barcode',xlabel='Pace Filtration',ylabel='Generators',inter_ind=None,minbirth_maxdeath = (minb,maxd),color_pal=color_pal)
    
        plt.savefig('H0Plots_level_'+str(level)+'/H0Plots_frame_'+str(filt_ind)+'.png',bbox_inches='tight')
        plt.close(tempfig)


        #
        axs = fig.get_axes()
        for ax in axs[5:]:
            fig.delaxes(ax)
        one_pers = [(x[0],(x[1][0],x[1][1])) for x in modified_pers if x[0] == 1]
        one_impact = [(x[0],(np.log(x[1][0]),np.log(x[1][1]))) for x in modified_pers if x[0] == 1]
        relevant_gens_one = [gen for gen in relevant_gens if len(gen[0]) == 2]
        color_pal = list(reversed(['#%02x%02x%02x'%(20*i % 255,30**i % 255, 40*i % 255) for i,_ in enumerate(one_pers)]))

        num_gens = len(relevant_gens_one)
        inner_grid = outer_grid[2,0].subgridspec(int(np.ceil(np.sqrt(num_gens))),int(np.ceil(np.sqrt(num_gens))))

        for i,gen in enumerate(relevant_gens_one):
            ax = fig.add_subplot(inner_grid[i])
            ax.set_xticks([])
            ax.set_yticks([])
            fig.add_subplot(ax)
            subgraph = filteredG.edge_subgraph(gen)
            temp = subgraph.copy()
            nx.draw_networkx(temp,pos,node_size=20,with_labels=False,ax=ax,edge_color='#%02x%02x%02x'%(20*i % 255,30**i % 255, 40*i % 255))

            coord_list = [(temp.nodes[node]['y'], temp.nodes[node]['x']) for node in temp.nodes]
            nodefeatureList = [{'type':'Feature', 'properties':{'color':'#%02x%02x%02x'%(20*i % 255,30**i % 255, 40*i % 255)},'geometry':{'type':'Point','coordinates': (x,y)}} for y,x in coord_list]

            nodeFeatureCol = {'type':'FeatureCollection','features': nodefeatureList}
            with open('CycleNodeData_level_'+str(level)+'/CycleNodeData_frame_'+str(filt_ind)+'_gen_'+str(i)+'.json','w+') as f:
                json.dump(nodeFeatureCol,f)

            edgefeatureList = []
            spathList = []
            for edge in temp.edges:
                spathList.append(filteredG.edges[edge]['spath'])
                try:
                    edgefeatureList.append({'type': 'Feature', 'properties': {}, 'geometry': shapely.geometry.mapping(temp.edges[edge]['geometry'])})
                except KeyError:
                    source,tar = edge
                    lats,lons = temp.nodes[source]['y'], temp.nodes[source]['x']
                    latt,lont = temp.nodes[tar]['y'], temp.nodes[tar]['x']
                    edgefeatureList.append({'type': 'Feature', 'properties': {'color':'#%02x%02x%02x'%(20*i % 255,30**i % 255, 40*i % 255)}, 'geometry': {'type':'LineString', 'coordinates': ((lons,lats),(lont,latt))}})

            edgeFeatureCol = {'type':'FeatureCollection','features': edgefeatureList}
            with open('CycleEdgeData_level_'+str(level)+'/CycleEdgeData_frame_'+str(filt_ind)+'_gen_'+str(i)+'.json','w+') as f:
                json.dump(edgeFeatureCol,f)

            joinedPathList = []
            for path1,path2 in zip(spathList[:-1],spathList[1:]):
                to_add = nx.dijkstra_path(bigG,path1[-1],path2[0])
                joinedPathList.append(path1)
                joinedPathList.append(to_add)
                
            joinedPathList.append(spathList[-1])
            joinedPathList.append(nx.dijkstra_path(bigG,spathList[-1][-1],spathList[0][0]))

            concatPath = sum(joinedPathList,[])
            pathfeatureList = []
            for src,tar in zip(concatPath[:-1],concatPath[1:]):
                try:
                    pathfeatureList.append({'type': 'Feature', 'properties': {}, 'geometry': shapely.geometry.mapping(bigG.edges[(src,tar)]['geometry'])})
                except KeyError:
                    lats,lons = bigG.nodes[src]['y'], bigG.nodes[src]['x']
                    latt,lont = bigG.nodes[tar]['y'], bigG.nodes[tar]['x']
                    pathfeatureList.append({'type': 'Feature', 'properties': {'color':'#%02x%02x%02x'%(20*i % 255,30**i % 255, 40*i % 255)}, 'geometry': {'type':'LineString', 'coordinates': ((lons,lats),(lont,latt))}})

            pathFeatureCol = {'type':'FeatureCollection','features': pathfeatureList}
      
            with open('CyclePathData_level_'+str(level)+'/CyclePathData_frame_'+str(filt_ind)+'_gen_'+str(i)+'.json','w+') as f:
                json.dump(pathFeatureCol,f)
        
        #
        ax21.clear()
        tempfig,tempax = plt.subplots(figsize=figsize)
        tempax.set_yticks([])
        plot_persistence_barcode_dan(one_pers,tempax,title='$H_1$ Barcode',xlabel='Pace Filtration',ylabel='Generators',inter_ind=None,minbirth_maxdeath = (minb,maxd),color_pal=color_pal)
        plt.savefig('H1Pers_level_'+str(level)+'/H1Pers_frame_'+str(filt_ind)+'.png',bbox_inches='tight')
        plt.close(tempfig)

        if cycle_type == 'pers':
            plot_persistence_barcode_dan(one_pers,ax21,title='$H_1$ Barcode',xlabel='Pace Filtration',ylabel='Generators',inter_ind=None,minbirth_maxdeath = (minb,maxd),color_pal=color_pal)
        elif cycle_type == 'impact':
            plot_persistence_barcode_dan(one_impact,ax21,title='$H_1$ Impact',xlabel='Log Pace Filtration',ylabel='Generators',inter_ind=None,minbirth_maxdeath = (0,np.log(maxd)),color_pal=color_pal)
        return fig.get_axes()

    #legend(loc=0)
    ani = animation.FuncAnimation(fig,update_img,[(i,x) for i,x in enumerate(filtrations) if x<= maxd],interval=200)
    writer = animation.writers['ffmpeg'](fps=5)

    ani.save(filename,writer=writer,dpi=dpi)
    plt.close()

    if plot_bars:
        fig,ax = plt.subplots(figsize=(12,12))
        ax.set_yticks([])
        plot_persistence_barcode_dan(pers,ax,title='$H_0$ and $H_1$ Barcodes',xlabel='Pace Filtration',ylabel='Generators',inter_ind=None,minbirth_maxdeath = (minb,maxd))
        plt.savefig(bar_file, bbox_inches='tight')
        plt.close()

    if plot_cycles:
        fig = plt.figure(figsize=(12,12), constrained_layout=False)
        filteredG = make_graph(simps,maxd,G)
        one_pers = [(x[0],(x[1][0],x[1][1])) for x in pers if x[0] == 1]
        one_impact = [(x[0],(np.log(x[1][0]),np.log(x[1][1]))) for x in pers if x[0] == 1]
        relevant_gens_one = [gen for filt,gen in generators if len(gen[0]) == 2]
        color_pal = list(reversed(['#%02x%02x%02x'%(20*i % 255,30**i % 255, 40*i % 255) for i,_ in enumerate(one_pers)]))

        num_gens = len(relevant_gens_one)
        inner_grid = fig.add_gridspec(int(np.ceil(np.sqrt(num_gens))),int(np.ceil(np.sqrt(num_gens))))

        for i,gen in enumerate(relevant_gens_one):
            ax = fig.add_subplot(inner_grid[i])
            ax.set_xticks([])
            ax.set_yticks([])
            fig.add_subplot(ax)
            subgraph = filteredG.edge_subgraph(gen)
            temp = subgraph.copy()
            nx.draw_networkx(temp,pos,node_size=20,with_labels=False,ax=ax,edge_color='#%02x%02x%02x'%(20*i % 255,30**i % 255, 40*i % 255))

        plt.savefig(cycle_file,bbox_inches='tight')
        plt.close()

        fig,ax = plt.subplots(figsize=(12,12))
        ax.set_yticks([])
        plot_persistence_barcode_dan(one_pers,ax,title='$H_1$ Barcode',xlabel='Pace Filtration',ylabel='Generators',inter_ind=None,minbirth_maxdeath = (minb,maxd),color_pal=color_pal)
        plt.savefig(cycle_file.replace('.png','_bar.png'),bbox_inches='tight')
        plt.close()

        fig,ax = plt.subplots(figsize=(12,12))
        ax.set_yticks([])
        plot_persistence_barcode_dan(one_impact,ax,title='$H_1$ Impact',xlabel='Log Pace Filtration',ylabel='Generators',inter_ind=None,minbirth_maxdeath = (0,np.log(maxd)),color_pal=color_pal)
        plt.savefig(cycle_file.replace('.png','_bar_impact.png'),bbox_inches='tight')
        plt.close()

    return 

def produce_files(inter_pruned,gen_pruned,G,pos,filtrations,bigG,level,simps,dpi):

    pers, geners = convert_to_gd_format(inter_pruned,gen_pruned)
    generators = nodify(geners,G)
    data = [(filtration , gen, ind) for filtration in filtrations 
           for ind,gen in 
               enumerate(get_corresponding_node_from_gen(gen_pruned,G,1,max_filt=filtration))
           ]
    minb,maxd = __min_birth_max_death(pers)
    
    def threshmax(val,maxval):
        if val <= maxval:
            return val
        else:
            return maxval
    
    #Use for displaying multiple H_1 gens
    num_gens = 0

    for filt_ind,mfilt in enumerate( [x for x in filtrations if x<= maxd]):
        relevant_inds = [i for i,gen in enumerate(generators) if gen[0] <= mfilt]
        modified_pers = [(dim_minmax[0],(dim_minmax[1][0],threshmax(dim_minmax[1][1],mfilt))) for ind,dim_minmax in enumerate(pers)
                        if ind in relevant_inds]
        relevant_gens = [gen for filt,gen in generators if filt <= mfilt]
       

        filteredG = make_graph(simps,mfilt,G)
        nx.set_node_attributes(filteredG,pos,'pos')
        spaths = {edge: G.edges[edge]['spath'] for edge in G.edges}
        nx.set_edge_attributes(filteredG,spaths,'spath')

        # Rips complex json data #
        for node in filteredG.nodes:
            filteredG.nodes[node]['x'] = G.nodes[node]['pos'][1]
            filteredG.nodes[node]['y'] = G.nodes[node]['pos'][0]

        coord_list = [(filteredG.nodes[node]['y'], filteredG.nodes[node]['x']) for node in filteredG.nodes]
        nodefeatureList = [{'type':'Feature', 'properties':{},'geometry':{'type':'Point','coordinates': (x,y)}} for y,x in coord_list]

        nodeFeatureCol = {'type':'FeatureCollection','features': nodefeatureList}
        with open('RipsGraphNodes_level_'+str(level)+'/RipsGraphNodes_frame_'+str(filt_ind)+'.json','w+') as f:
            json.dump(nodeFeatureCol,f)

        edgefeatureList = []


        for edge in filteredG.edges:
            try:
                edgefeatureList.append({'type': 'Feature', 'properties': {}, 'geometry': shapely.geometry.mapping(filteredG.edges[edge]['geometry'])})
            except KeyError:
                source,tar = edge
                lats,lons = filteredG.nodes[source]['y'], filteredG.nodes[source]['x']
                latt,lont = filteredG.nodes[tar]['y'], filteredG.nodes[tar]['x']
                edgefeatureList.append({'type': 'Feature', 'properties': {}, 'geometry': {'type':'LineString', 'coordinates': ((lons,lats),(lont,latt))}})

        edgeFeatureCol = {'type':'FeatureCollection','features': edgefeatureList}
        with open('RipsGraphEdges_level_'+str(level)+'/RipsGraphEdges_frame_'+str(filt_ind)+'.json','w+') as f:
            json.dump(edgeFeatureCol,f)

        # Total persistence plot #

        fig,ax = plt.subplots()
        plot_persistence_barcode_dan(modified_pers,ax,title='Total Barcode',xlabel='Pace Filtration',ylabel='Generators',
                            inter_ind=None,minbirth_maxdeath = (minb,maxd))
        plt.savefig('TotPers_level_'+str(level)+'/TotPers_frame_'+str(filt_ind)+'.png')
        plt.close()

        zero_pers = [x for x in modified_pers if x[0] == 0]
        comps = [list(x) for x in nx.connected_components(filteredG)]
        relevant_gens_zero = [gen for gen in relevant_gens if len(gen[0]) == 1]
        color_pal = list(reversed(['#%02x%02x%02x'%(20*i % 255,30**i % 255, 40*i % 255) for i,_ in enumerate(zero_pers)]))


        for i,gen in enumerate(relevant_gens_zero):
            comp = [x for x in comps if gen[0][0] in x][0]
            subgraph = filteredG.subgraph(comp)
            temp = subgraph.copy()
            coord_list = [(temp.nodes[node]['y'], temp.nodes[node]['x']) for node in temp.nodes]
            nodefeatureList = [{'type':'Feature', 'properties':{'color':'#%02x%02x%02x'%(20*i % 255,30**i % 255, 40*i % 255)},'geometry':{'type':'Point','coordinates': (x,y)}} for y,x in coord_list]

            nodeFeatureCol = {'type':'FeatureCollection','features': nodefeatureList}
            with open('H0gens_level_'+str(level)+'/H0gens_frame_'+str(filt_ind)+'_gen_'+str(i)+'.json','w+') as f:
                json.dump(nodeFeatureCol,f)

 #           temp.remove_edges_from(subgraph.edges)
#            nx.draw_networkx(temp,pos,node_size=200,with_labels=False,ax=ax10,node_color='#%02x%02x%02x'%(20*i % 255,30**i % 255, 40*i % 255))

            
        fig,ax = plt.subplots()
        plot_persistence_barcode_dan(zero_pers,ax,title='$H_0$ Barcode',xlabel='Pace Filtration',ylabel='Generators',inter_ind=None,minbirth_maxdeath = (minb,maxd),color_pal=color_pal)
        plt.savefig('H0Plots_level_'+str(level)+'/H0Plots_frame_'+str(filt_ind)+'.png')
        plt.close()

        one_pers = [x for x in modified_pers if x[0] == 1]
        relevant_gens_one = [gen for gen in relevant_gens if len(gen[0]) == 2]
        color_pal = list(reversed(['#%02x%02x%02x'%(20*i % 255,30**i % 255, 40*i % 255) for i,_ in enumerate(one_pers)]))

        for i,gen in enumerate(relevant_gens_one):
            subgraph = filteredG.edge_subgraph(gen)
            temp = subgraph.copy()
            
            coord_list = [(temp.nodes[node]['y'], temp.nodes[node]['x']) for node in temp.nodes]
            nodefeatureList = [{'type':'Feature', 'properties':{'color':'#%02x%02x%02x'%(20*i % 255,30**i % 255, 40*i % 255)},'geometry':{'type':'Point','coordinates': (x,y)}} for y,x in coord_list]

            nodeFeatureCol = {'type':'FeatureCollection','features': nodefeatureList}
            with open('CycleNodeData_level_'+str(level)+'/CycleNodeData_frame_'+str(filt_ind)+'_gen_'+str(i)+'.json','w+') as f:
                json.dump(nodeFeatureCol,f)

            edgefeatureList = []
            spathList = []
            for edge in temp.edges:
                spathList.append(G.edges[edge]['spath'])
                try:
                    edgefeatureList.append({'type': 'Feature', 'properties': {}, 'geometry': shapely.geometry.mapping(temp.edges[edge]['geometry'])})
                except KeyError:
                    source,tar = edge
                    lats,lons = temp.nodes[source]['y'], temp.nodes[source]['x']
                    latt,lont = temp.nodes[tar]['y'], temp.nodes[tar]['x']
                    edgefeatureList.append({'type': 'Feature', 'properties': {'color':'#%02x%02x%02x'%(20*i % 255,30**i % 255, 40*i % 255)}, 'geometry': {'type':'LineString', 'coordinates': ((lons,lats),(lont,latt))}})

            edgeFeatureCol = {'type':'FeatureCollection','features': edgefeatureList}
            with open('CycleEdgeData_level_'+str(level)+'/CycleEdgeData_frame_'+str(filt_ind)+'_gen_'+str(i)+'.json','w+') as f:
                json.dump(edgeFeatureCol,f)

            joinedPathList = []
            for path1,path2 in zip(spathList[:-1],spathList[1:]):
                to_add = nx.dijkstra_path(bigG,path1[-1],path2[0])
                joinedPathList.append(path1)
                joinedPathList.append(to_add)
                
            joinedPathList.append(spathList[-1])
            joinedPathList.append(nx.dijkstra_path(bigG,spathList[-1][-1],spathList[0][0]))

            concatPath = sum(joinedPathList,[])
            pathfeatureList = []
            for src,tar in zip(concatPath[:-1],concatPath[1:]):
                try:
                    pathfeatureList.append({'type': 'Feature', 'properties': {}, 'geometry': shapely.geometry.mapping(bigG.edges[(src,tar)]['geometry'])})
                except KeyError:
                    lats,lons = bigG.nodes[src]['y'], bigG.nodes[src]['x']
                    latt,lont = bigG.nodes[tar]['y'], bigG.nodes[tar]['x']
                    pathfeatureList.append({'type': 'Feature', 'properties': {'color':'#%02x%02x%02x'%(20*i % 255,30**i % 255, 40*i % 255)}, 'geometry': {'type':'LineString', 'coordinates': ((lons,lats),(lont,latt))}})

            pathFeatureCol = {'type':'FeatureCollection','features': pathfeatureList}
      
            with open('CyclePathData_level_'+str(level)+'/CyclePathData_frame_'+str(filt_ind)+'_gen_'+str(i)+'.json','w+') as f:
                json.dump(pathFeatureCol,f)

                   
        fig,ax = plt.subplots()
        plot_persistence_barcode_dan(one_pers,ax,title='$H_1$ Barcode',xlabel='Pace Filtration',ylabel='Generators',inter_ind=None,minbirth_maxdeath = (minb,maxd),color_pal=color_pal)
        plt.savefig('H1Pers_level_'+str(level)+'/H1Pers_frame_'+str(filt_ind)+'.png')
        plt.close()
    
    return
