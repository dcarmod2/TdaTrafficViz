#!/bin/bash

mkdir GraphImages
mkdir MovieFiles
mkdir Barcodes
mkdir GraphNodes_level_0
mkdir GraphEdges_level_0

for i in {1..3}
do
    mkdir "CycleEdgeData_level_-$i"
    mkdir "CycleNodeData_level_-$i"
    mkdir "CyclePathData_level_-$i"
    mkdir "CyclePathDataRev_level_-$i"
    mkdir "H0gens_level_-$i"
    mkdir "H0gens_fullGraph_level_-$i"
    mkdir "H0Plots_level_-$i"
    mkdir "H1Pers_level_-$i"
    mkdir "RipsGraphEdges_level_-$i"
    mkdir "RipsGraphNodes_level_-$i"
    mkdir "GraphEdges_level_-$i"
    mkdir "GraphNodes_level_-$i"
    mkdir "TotPers_level_-$i"
done
