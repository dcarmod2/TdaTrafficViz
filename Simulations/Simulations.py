import matplotlib
matplotlib.use('Agg')
import re
import pickle
import yaml
from pers_engine import *
from pers_plotter import *


if __name__ == "__main__":

    #Read data from sim_config.yml
    with open('sim_config.yml','r') as f:
        yml = yaml.load(f,Loader=yaml.FullLoader)

    config = yml['Config']
    fontsize = config['plotFontsize']
    genfile = config['genfile']
    logfile = config['logfile']
    movieFilePrefix = config['movieFilePrefix']
    ripsType = config['ripsType']
    cyclePlotType = config['cycleMovieType']
    subG = get_largest_comp(config['graphFile'])
    plt.rcParams.update({'font.size': 22})

    # Assign speed and pace to edges
    for edge in subG.edges:
        subG.edges[edge]['speed'] = subG.edges[edge]['length']/subG.edges[edge]['time']
        subG.edges[edge]['pace'] = subG.edges[edge]['time']/subG.edges[edge]['length']

    # Setup for louvain
    nodeL = list(subG.nodes)
    graphData = [(nodeL.index(edge[0]),nodeL.index(edge[1]),subG.edges[edge]['time'],subG.edges[edge]['pace'],subG.edges[edge]['speed']) for edge in subG.edges]
    datalen = len(graphData)

    with open('subG.txt','w+') as f:
        for i,trip in enumerate(graphData):
            src,dest,time,pace,speed = trip
            if i < datalen - 1:
                f.write("{0:d} {1:d} {2:.9f}\n".format(src,dest,speed))
            else:
                f.write("{0:d} {1:d} {2:.9f}".format(src,dest,speed))

    subprocess.check_call([r"./DirectedLouvain/bin/convert", "-i", "subG.txt", "-o", "subG.bin", "-w", "subG.weights"])

    # run Louvain
    with open("subG.tree","w+") as f:
        subprocess.check_call([r"./DirectedLouvain/bin/community", "subG.bin", "-l", "-1", "-w", "subG.weights"],stdout = f)

    with open("temp.txt","w+") as f:
        subprocess.check_call(["./DirectedLouvain/bin/hierarchy","subG.tree"],stdout=f)

    # get louvain output
    with open("temp.txt","r+") as f:
        lines = f.readlines()
        first = lines[0]
        num_levels = int(re.findall(r'\d+',first)[0]) 

    filtrations_max = []

    #remove old log file
    if os.path.exists(logfile):
        os.remove(logfile)


    ## Start persistent homology calculation ##

    # Get shortest path distances
    print("Calculating distance matrix...")
    tic = timemod.clock()
    dist = np.array(nx.floyd_warshall_numpy(subG,weight='pace'))
    toc = timemod.clock()
    print("Distance matrix calculated in " + str(toc-tic) + " seconds.")
    with open(logfile,"a+") as f:
        f.write("Dist mat: " + str(toc-tic) + " seconds.\n")
    
    filtrations = np.unique(sorted(dist.flatten()))
    if (len(filtrations_max) < len(filtrations)):
        filtrations_max = filtrations.copy()
    
    # form Rips complex
    print("Forming Rips complex...")
    tic = timemod.clock()
    if ripsType == "max":
        rips_comp = [sorted(get_simplicesGen(dist,filtrations_max,k,2))  for k in range(0,3)]
    elif ripsType == "min":
        rips_comp = [sorted(get_simplicesGen_min(dist,filtrations_max,k,2))  for k in range(0,3)]
    else:
        raise Exception('Error: unknown rips configuration options')
    simps = sum(rips_comp,[])
    toc = timemod.clock()
    print("Rips complex formed in " + str(toc-tic) + " seconds.")
    with open(logfile,"a") as f:
        f.write("Rips complex: " + str(toc-tic) + " seconds.\n")
    
    # Calculate persistent homology
    print("Calculating persistent homology...")
    tic = timemod.clock()
    inter,gen = computeIntervalsGen(simps)
    toc = timemod.clock()
    print("Homology and generators calculated in " + str(toc-tic) + " seconds.")
    with open(logfile,"a") as f:
        f.write("Pers Hom: " + str(toc-tic) + " seconds.\n")
    
        
    inter_pruned = {i: {x:y for x,y in inter[i].items() if y[0]!=y[1]} for i in inter.keys()}
    gen_pruned = {i: {x:y for x,y in gen[i].items() if x in inter_pruned[i]} for i in gen.keys()}
    
    #Set up for making movie

    pos = {node: (subG.nodes[node]['y'],subG.nodes[node]['x']) for node in subG.nodes}
    if config['plotGraph']:
        edge_labels = {edge: num_display(subG.edges[edge]['pace']) for edge in subG.edges}
        plt.figure(figsize=(12,12))
        nx.draw_networkx(subG,pos,node_size=20,with_labels=False)
        nx.draw_networkx_edge_labels(subG,pos,edge_labels)
        plt.savefig(config['plotGraphFile'])
        plt.close()
    dpi = 100
    print("Animating...")
    tic = timemod.clock()
    ani_frame(movieFilePrefix + "_level_0.mp4",inter_pruned,gen_pruned,subG,pos,filtrations_max,simps,dpi,config['plotBarcode'],config['plotBarcodeFileBase']+'_level_0.png',plot_cycles = config['plotCycles'], cycle_file = config['plotCycleFileBase'] + '_level_0.png',cycle_type = cyclePlotType,fontsize=fontsize)
    toc = timemod.clock()
    with open(logfile,"a") as f:
        f.write("Animation: " + str(toc-tic) + " seconds.")
    with open(genfile,"wb") as f:
        pickle.dump(gen_pruned,f)

    # Now repeat the above steps with graphs successively condensed by Louvain

    for level in range(-1*num_levels+1,0):
        
        cluster_df,condensed_G = condense_graph(subG,"temp.txt",num_levels,level,nodeL)
        
        print("Plotting...")
        if config['plotCondensed']:
            fig,ax = plt.subplots(figsize=(40,30))
            #fix pos to be centroid
            pos = {node:(np.mean([subG.nodes[cnode]['y'] for cnode in cluster_df[cluster_df['community']==node]['node']]),np.mean([subG.nodes[cnode]['x'] for cnode in cluster_df[cluster_df['community']==node]['node']])) for node in condensed_G.nodes}
            labels = nx.get_edge_attributes(condensed_G,'pace')
            labels_rounded = {key: convert_val(val) for key,val in labels.items()}
            nx.draw_networkx(condensed_G,pos,node_size=20,with_labels=False,ax=ax)
            temp = nx.draw_networkx_edge_labels(condensed_G,pos,edge_labels=labels_rounded)
            plt.savefig(config['plotCondensedFile'].replace('.png','_level_'+str(level)+'.png'),bbox_inches='tight')
            plt.close()
        #in case we want to use speed
        for edge in condensed_G.edges:
            condensed_G.edges[edge]['speed'] = 1/condensed_G.edges[edge]['pace']

        print("Calculating distance matrix...")
        tic = timemod.clock()
        dist = np.array(nx.floyd_warshall_numpy(condensed_G,weight='pace'))
        dist_symm = np.maximum(dist,dist.T)
        toc = timemod.clock()
        print("Distance matrix calculated in " + str(toc-tic) + " seconds.")
        with open(logfile,"a+") as f:
            f.write("Dist mat: " + str(toc-tic) + " seconds.\n")
    
        filtrations = np.unique(sorted(dist.flatten()))
        if (len(filtrations_max) < len(filtrations)):
            filtrations_max = filtrations.copy()
    
        print("Forming Rips complex...")
        tic = timemod.clock()
        simps = sum([sorted(get_simplicesGen(dist,filtrations_max,k,2))  for k in range(0,3)],[])
        toc = timemod.clock()
        print("Rips complex formed in " + str(toc-tic) + " seconds.")
        with open(logfile,"a") as f:
            f.write("Rips complex: " + str(toc-tic) + " seconds.\n")
    
    
        print("Calculating persistent homology...")
        tic = timemod.clock()
        inter,gen = computeIntervalsGen(simps)
        toc = timemod.clock()
        print("Homology and generators calculated in " + str(toc-tic) + " seconds.")
        with open(logfile,"a") as f:
            f.write("Pers Hom: " + str(toc-tic) + " seconds.\n")
    
        
        inter_pruned = {i: {x:y for x,y in inter[i].items() if y[0]!=y[1]} for i in inter.keys()}
        gen_pruned = {i: {x:y for x,y in gen[i].items() if x in inter_pruned[i]} for i in gen.keys()}
    
    #Set up for making movie
        dpi = 100
        print("Animating...")
        tic = timemod.clock()
        ani_frame(movieFilePrefix + "_level_" + str(level) + ".mp4",inter_pruned,gen_pruned,condensed_G,pos,filtrations_max,simps,dpi,config['plotBarcode'],config['plotBarcodeFileBase']+'_level_'+str(level)+'.png',plot_cycles = config['plotCycles'], cycle_file = config['plotCycleFileBase']+'_level_'+str(level)+'.png',cycle_type = cyclePlotType,fontsize=fontsize)
        toc = timemod.clock()
        with open(logfile,"a") as f:
            f.write("Animation: " + str(toc-tic) + " seconds.")
        with open(genfile,"wb") as f:
            pickle.dump(gen_pruned,f)
    
