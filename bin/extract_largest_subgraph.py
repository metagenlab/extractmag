import statistics
import networkx as nx
from Bio import SeqIO
import os
import re
import pandas 
from Bio.SeqUtils import GC

class Assembly():
    def __init__(self, 
                 assembly_path, 
                 scaffolds_or_contigs_fasta=False,
                 scaffolds_or_contigs_paths=False,
                 format='fastg'
                 ):
        
        if format == 'fastg':
            import pyfastg
            self.assembly_nx = pyfastg.parse_fastg(assembly_path)
        elif format == 'GFA':
            self.assembly_nx = self.parse_gfa(assembly_path)
            
        #self.node2contig = self.node_id2contig_name(scaffolds_or_contigs_paths)
        #self.contig2seq = SeqIO.to_dict(SeqIO.parse(scaffolds_or_contigs_fasta, "fasta"))
        
        self.cols = ["#de2d12", "#12de1e", "#2b12de", "#12dbde", "#de2d12", "#12de1e", "#2b12de", "#12dbde"]
        self.subgraphs = {}
        self.subgraphs_index = 0
    
    
    def parse_gfa(self, gfa_file_path):

        '''
        length
        gc
        cov
        
        LN 	i 	Segment length
        RC 	i 	Read count
        FC 	i 	Fragment count
        KC 	i 	k-mer count
        SH 	H 	SHA-256 checksum of the sequence
        UR 	Z 	URI or local file-system path of the sequence. If it does not start with a standard protocol (e.g. ftp), it is assumed to be a local path.
        DP
        GC
        
        '''

        if gfa_file_path.suffix == ".gz":
            import gzip
            gfa_file_obj = gzip.open(gfa_file_path, "rt")
        else:
            gfa_file_obj = open(gfa_file_path, 'r')

        gg_network = nx.DiGraph()

        for line in gfa_file_obj:

            if line[0] == 'H':
                # Header line
                header_list = line.split('\t')
                format_version = header_list[1].strip()

            if line[0] == 'S':
                # Node line
                node_list = line.split('\t')
                node_id = node_list[1]
                node_seq = node_list[2].strip()
                metadata = {}
                for opt_field in node_list[3:]:
                    #print("OPTIONAL", opt_field)
                    label, data_type, value = opt_field.split(":")
                    metadata[label.upper()] = float(value)
                if 'DP' in metadata:
                    cov = metadata['DP']
                else:
                    cov=None
                if 'LN' in metadata:
                    length = metadata['LN']
                else:
                    length = len(node_seq)
                gc = GC(node_seq)
                gg_network.add_node(node_id + '+', sequence=node_seq, gc=gc, length=length, cov=cov)
                gg_network.add_node(node_id + '-', sequence=node_seq, gc=gc, length=length, cov=cov)

        if gfa_file_path.suffix == ".gz":
            import gzip
            gfa_file_obj = gzip.open(gfa_file_path, "rt")
        else:
            gfa_file_obj = open(gfa_file_path, 'r')
            
        for line in gfa_file_obj:

            if line[0] == 'L':
                # link / edge node
                edge_list = line.split('\t')
                edge_source = edge_list[1] + edge_list[2]
                edge_target = edge_list[3] + edge_list[4] 


                if edge_list[2] == '+':
                    edge_source_rev = edge_list[1] + '-'
                else:
                    edge_source_rev = edge_list[1] + '+'
                if edge_list[4] == '+':
                    edge_target_rev = edge_list[3] + '-'
                else:
                    edge_target_rev = edge_list[3] + '+'
                             
                #if edge_list[2] == '+':
                gg_network.add_edge(edge_source, edge_target)
                #gg_network.add_edge(edge_target, edge_source, )
                #else:
                gg_network.add_edge( edge_target_rev, edge_source_rev)
                #gg_network.add_edge(edge_target_rev, edge_source_rev, )

            if line[0] == 'P':
                # Add path info
                path_info = line.split('\t')
                seq_info = path_info[2]
                path_nodes = path_info[1]

        return gg_network


    def get_path_stats(self, subgraph, node_len_cutoff=0):
        '''
        get path gc, median GC and depth
        only keep contigs longer than <node_len_cutoff> to calculate GC and depth
        '''
        cumul_size = sum([ self.assembly_nx.nodes[node]["length"] for node in subgraph])
        gc_median = statistics.median([self.assembly_nx.nodes[node]["gc"] for node in subgraph if self.assembly_nx.nodes[node]["length"] > node_len_cutoff])
        depth_median = statistics.median([self.assembly_nx.nodes[node]["cov"] for node in subgraph if self.assembly_nx.nodes[node]["length"] > node_len_cutoff])
        
        return cumul_size, gc_median, depth_median

    def filter_based_on_depth(self, subgraph, min_depth):
        
        keep = []
        for node in subgraph:
            if self.assembly_nx.nodes[node]["cov"] > min_depth:
                keep.append(node)
        return keep


    def extract_largest_components(self, 
                                   cumul_size_cutoff=10000, 
                                   stats_cutoff=0, 
                                   filter_depth_min_factor=0.1):
        
        components = list(nx.weakly_connected_components( self.assembly_nx))
        
        large_count = 0
        for n,comp in enumerate(components):

            cumul_size, gc_median, depth_median = self.get_path_stats(comp)
            #print(f"One comp {n}:", len(comp), cumul_size, gc_median, depth_median)
            
            if cumul_size > cumul_size_cutoff:

                large_count += 1
                cumul_size, gc_median, depth_median = self.get_path_stats(comp, stats_cutoff)
                # filter low depth contigs
                if filter_depth_min_factor:
                    comp = self.filter_based_on_depth(comp, 
                                                      depth_median*filter_depth_min_factor)

                self.subgraphs[f"{self.subgraphs_index}"] = {"graph": comp,
                                                            "GC_median": gc_median,
                                                            "cumul_size": cumul_size,
                                                            "depth_median": depth_median}

                
    def write_largest(self, out_prefix):
        largest = 0
        g = None
        for sub in self.subgraphs:
            graph_data = self.subgraphs[sub]
            if graph_data["cumul_size"] > largest:
                largest = graph_data["cumul_size"]
                g = graph_data
        with open(f"{out_prefix}_largest_target_subgraph.tsv", "w") as f:
            for node in graph_data["graph"].nodes:
                f.write(node + "\n")
        with open(f"{out_prefix}_largest_target_subgraph_stats.tsv", "w") as f:
                f.write(f'{out_prefix}\t{graph_data["cumul_size"]}\t{graph_data["depth_median"]}\t{graph_data["GC_median"]}\n')
    
    

if __name__ == "__main__":
    
    import argparse
    from pathlib import Path
    
    """Define and immediately parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Validate and transform a tabular samplesheet.",
        epilog="Example: python check_samplesheet.py samplesheet.csv samplesheet.valid.csv",
    )
    parser.add_argument(
        "-g", '--assembly_graph',
        type=Path,
        help="Samtools depth file (TSV format).",
    )
    parser.add_argument(
        "-p", '--prefix',
        type=str,
        default='out',
        help="Output prefix",
    )

    args = parser.parse_args()

    a = Assembly(args.assembly_graph, format=args.format)
    a.extract_largest_components(cumul_size_cutoff=10000, # only keep subgraphs largers than <cumul_size_cutoff>
                                 stats_cutoff=5000) # only use contigs largers than <stats_cutoff> to caldulate median depth
    