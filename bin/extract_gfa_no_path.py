#!/usr/bin/env python3

import argparse
import sys
from collections import defaultdict

parser = argparse.ArgumentParser()
parser.add_argument("-c", '--cat_bat', type=str,help="contig2classification file")
parser.add_argument("-g", '--gfa', type=str, help="assembly graph file")
parser.add_argument("-t", '--taxa', type=str, help="taxids to keep (coma separated list)")
parser.add_argument("-s", "--size", type=str, help="Filter out contigs smaller than this size")
parser.add_argument("-n", "--nodes", action="store_true", help="Output node ids instead of contigs")
parser.add_argument("--subgraph", type=str, help="File where to store the subgraph (optional)")

args = parser.parse_args()

includes = set()
for taxid in args.taxa.split(","):
    includes.add(taxid)
print(includes)
ctg_list = set()

for line in open(args.cat_bat, "r"):
    tokens = line.split("\t")
    if tokens[1].startswith("no"):
        continue

    for lineage in tokens[3].split(";"):
        if lineage in includes:
            ctg_list.add(tokens[0])

node_to_contig = {}
graph = defaultdict(list)

to_pass = [str(i) for i in ctg_list]
kept = set()

for line in open(args.gfa, "r"):
    if line[0]=="L":
        tokens = line.split('\t')
        print(tokens)
        from_node = str(tokens[1])
        to_node = str(tokens[3])
        if from_node == "62262" or to_node == "62262":
            print("from_node", from_node)
            print("to_node", to_node)
        graph[from_node].append(to_node)
        graph[to_node].append(from_node)


print(graph["62262"])
print("len to_pass", len(to_pass))
while not len(to_pass)==0:
    last = to_pass.pop()
    next_nodes = graph[last]
    #print("next_nodes",last, type(last), next_nodes)
    kept.add(last)

    for node in next_nodes:
        if node in kept:
            continue
        to_pass.append(node)

if args.nodes:
    print("\n".join(str(i) for i in kept))
else:
    selected_contigs = {node_to_contig[i] for i in kept if i in node_to_contig}
    print("\n".join(selected_contigs))

print("n kept", len(kept))

if not args.subgraph is None:
    new_gfa = open(args.subgraph, "w")
    for line in open(args.gfa, "r"):
        if line[0]=="S":
            tokens = line.split('\t')
            if str(tokens[1]) in kept:
                print(line, file=new_gfa)
        elif line[0]=="L":
            tokens = line.split('\t')
            if str(tokens[1]) in kept:
                print(line, file=new_gfa)
    new_gfa.close()
