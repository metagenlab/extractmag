#!/opt/conda/bin/python

from Bio import SeqIO
from Bio.SeqUtils import GC
import networkx as nx
import statistics
import re 
import pandas
import numpy as np
from Bio import SeqRecord
from Bio import Seq


def calculate_N50(list_of_lengths):
    
    """Calculate N50 for a sequence of numbers.
 
    Args:
        list_of_lengths (list): List of numbers.
 
    Returns:
        float: N50 value.
 
    """
    if len(list_of_lengths) == 0:
        return None
    if len(list_of_lengths) == 1:
        return list_of_lengths[0]
    tmp = []
    for tmp_number in set(list_of_lengths):
            tmp += [tmp_number] * list_of_lengths.count(tmp_number) * tmp_number
    tmp.sort()
 
    if (len(tmp) % 2) == 0:
        median = (tmp[int(len(tmp) / 2) - 1] + tmp[int(len(tmp) / 2)]) / 2
    else:
        median = tmp[int(len(tmp) / 2)]
 
    return median


class Assembly:
    def __init__(self, 
                 gfa,
                 coverm,
                 contig2taxonomy,
                 fasta,
                 summary_rank='p__',
                 output_prefix='out',
                 target_taxon=["Chlamydiota"]):
        
        self.assembly_nx = self.parse_gfa(gfa)
        
        self.subgraphs = {}
        self.subgraphs_index = 0
        self.coverm = pandas.read_csv(coverm, sep="\t")
        self.coverm["Contig"] = self.coverm["Contig"].astype(str)
        self.coverm = self.coverm.set_index("Contig")
        self.ranks = ["p__", "o__","f__","g__", "s__"]
        
        self.taxon2unconnected_contigs = {}
        
        # for spades: mapping with contigs
        self.contig2rec = False

        for node in self.assembly_nx.nodes:
            self.assembly_nx.nodes[node]["cov"] = self.coverm.loc[re.sub("\-|\+", "",node)]["Mean"]
        
        
        # root;d__Bacteria;p__Actinobacteriota;c__Actinomycetia;o__Mycobacteriales;f__Mycobacteriaceae;g__Rhodococcus_B;s__Rhodococcus_B
        self.cat_taxnonomy = pandas.read_csv(contig2taxonomy, sep="\t")
        self.cat_taxnonomy["# contig"] = self.cat_taxnonomy["# contig"].astype(str)
        self.cat_taxnonomy = self.cat_taxnonomy.set_index("# contig")

        #self.cat_taxnonomy["lineage"] = self.cat_taxnonomy.apply(lambda row: self._format_gtdb_lineage(row, target_rank=summary_rank), axis=1)
        for rank in self.ranks:
            self.cat_taxnonomy[rank] = self.cat_taxnonomy.apply(lambda row: self._format_gtdb_lineage(row, target_rank=rank), axis=1)
        #self.cat_taxnonomy.apply(lambda row: row["lineage"].split(";") if not pandas.isna(row["lineage"]) else '', axis=1)
        # parse fasta
        self.fasta = [i for i in SeqIO.parse(fasta, "fasta")]
        self.gc_data = pandas.DataFrame([(record.name, round(GC(record.seq), 2), len(record.seq)) for record in self.fasta], columns=["contig", "GC", "length"]).set_index("contig")
        for node in self.assembly_nx.nodes:
            self.assembly_nx.nodes[node]["gc"] = self.gc_data.loc[re.sub("\-|\+", "",node)]["GC"]
            self.assembly_nx.nodes[node]["length"] = self.gc_data.loc[re.sub("\-|\+", "",node)]["length"]
            #print(type(self.gc_data.loc[re.sub("\-|\+", "",node)]["length"]))
        self.contig_summary = self.coverm.join(self.cat_taxnonomy).join(self.gc_data)
        self.contig_summary_filtered_size = self.contig_summary[self.contig_summary["Length"] > 1000]
        
        self.phylum2cumul_size = self.contig_summary_filtered_size[["p__", "Length"]].groupby("p__").sum().to_dict()["Length"]
        
        self.contig_summary.to_csv(f"{output_prefix}_contigs.tsv", sep="\t")
        

    def _format_gtdb_lineage(self, row, target_rank="p__"):
        if row.classification == 'no taxid assigned':
            return None
        else:
            try:
                classif = row.lineage.split(target_rank)[1].split(";")[0]
            except:
                classif = "failed"
            return classif
    
    def _extract_taxon(self, string, rank='p__'):
        import re 
        string = str(string)

        # root;d__Bacteria;p__Actinobacteriota;c__Actinomycetia;o__Mycobacteriales;f__Mycobacteriaceae;g__Rhodococcus_B;s__Rhodococcus_B
        result = re.search(f"{rank}([^;]+)", string)
        if result:
            return result.group(1)
        else:
            return None
        
    def get_consensus_taxonomy(self, subgraph):
        '''
        Return dict with the most frequent taxon and the percentage of bases supporting the classification
        {"p__": ["Chlamydiae", 99],
        "o__": ["Chlamydiales", 72],
        "f__": ["Chlamydiaceae", 33],}
        '''
        
        cumul_size = sum([ self.assembly_nx.nodes[node]["length"] for node in subgraph])
        tax2size_bp = {}
        for rank in self.ranks:
            tax2size_bp[rank] = {}
        for node in subgraph:
            l = self.assembly_nx.nodes[node]["length"]
            tax = self.cat_taxnonomy.loc[re.sub("\-|\+", "",node)]["lineage"]
            tax_format = [self._extract_taxon(tax, rank) for rank in self.ranks]
            for rank, taxon in zip(self.ranks, tax_format):
                
                # save taxonomy to node
                self.assembly_nx.nodes[node][rank] = taxon
                
                if taxon not in tax2size_bp[rank]:
                    tax2size_bp[rank][taxon] = l 
                else:
                    tax2size_bp[rank][taxon] += l 
        rank2best = {}
        for rank in self.ranks:
            taxo_summary = pandas.DataFrame([(key,value) for key,value in tax2size_bp[rank].items()], columns=["taxonomy", "bp"]).sort_values("bp", ascending=False)
            taxo_summary["percent"] = [(tax_size/cumul_size)*100 for tax_size in taxo_summary["bp"]]
            # exclude missing classification
            if len(taxo_summary) == 1 and taxo_summary.iloc[0].taxonomy == None:
                rank2best[rank] = [taxo_summary.iloc[0].taxonomy, taxo_summary.iloc[0].percent]
            else:
                taxo_summary = taxo_summary.dropna(subset=["taxonomy"])
                rank2best[rank] = [taxo_summary.iloc[0].taxonomy, taxo_summary.iloc[0].percent]
            
        return  rank2best
    
    def extract_subgraph(self, taxon_list, out_prefix):
        
        for taxon in taxon_list:
            df_target = self.contig_summary.dropna(subset=['lineage'])
            df_target = df_target[self.contig_summary['lineage'].astype(str).str.contains(taxon)]
            
            current_nodes = set(df_target.index.tolist())  # our list of known nodes
            checked_nodes = set()  # our list of checked nodes
            node_list = set(df_target.index.tolist()) 
            while len(current_nodes - checked_nodes) > 0:
                # first, randomly select a new node that we know about AND haven't already checked:
                current_node = (current_nodes - checked_nodes).pop()
                # get the neighbors of the node as unique elements:
                neighbors = set(nx.neighbors(self.assembly_nx, current_node))
                # add any results we don't already know about to our list of known nodes:
                current_nodes |= neighbors  
                # save our results for this node:
                node_list |= neighbors
                # finally, mark that we checked the node so we don't check it again:
                checked_nodes.add(current_node)
            subgraph = self.assembly_nx.subgraph(node_list)
            self.to_gfa(subgraph, f"{out_prefix}_{taxon}.gfa")
            return node_list
        
    def to_gfa(self, graph, output_name):
        with open(output_name, "w") as f:
            for node in graph.nodes:
                f.write(graph.nodes[node]["gfa_S"])
            for edge in graph.edges:
                f.write(graph.edges[edge]["gfa_L"])
        
    def get_contigs_from_nodes(self, node_list):
        '''
        get SeqIO record of spades contigs from node list
        '''
        # it can happen that a node is not in the paths file
        contig_set = set([re.sub("'", "", self.node2contig_name[node]) for node in node_list if node in self.node2contig_name])
        rec_set = [self.contig2rec[contig] for contig in contig_set]
        return rec_set
    
    def parse_spades_paths(self, paths_file, fasta_file):
        '''
        parse spades contig file 
        => self.contig2rec
        get correspondance between spades assembly graph nodes and contig names
        => self.node2contig_name
        
        '''
        self.contig2rec = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))
        self.node2contig_name = {}
        with open(paths_file, "r") as f:
            for n,line in enumerate(f):
                if  line.startswith("NODE"):
                    contig_name = line.strip()
                    
                else:
                    node_list = [re.sub("[\+\-;]", "", node) for node in line.strip().split(",")]
                    for node in node_list:
                        self.node2contig_name[node] = contig_name
        
    
    
    def extract_largest_components(self, 
                                   cumul_size_cutoff=10000, 
                                   stats_cutoff=0, 
                                   filter_depth_min_factor=0.1, 
                                   min_contig_size=1000):
        
        components = list(nx.connected_components(self.assembly_nx))
        
        large_count = 0
        for n,comp in enumerate(components):

            cumul_size, cumul_size_filt, gc_median, depth_mean = self.get_path_stats(comp)
            #print(f"One comp {n}:", len(comp), cumul_size, gc_median, depth_median)
            
            if cumul_size > cumul_size_cutoff:
                large_count += 1
                cumul_size, cumul_size_filt, gc_median, depth_mean = self.get_path_stats(comp, stats_cutoff)

                self.save_subgraph(comp, node_len_cutoff=min_contig_size)
    
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

        if gfa_file_path.endswith(".gz"):
            import gzip
            gfa_file_obj = gzip.open(gfa_file_path, "rt")
        else:
            gfa_file_obj = open(gfa_file_path, 'r')

        gg_network = nx.Graph()

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
                gc = round(GC(node_seq), 2)
                #gg_network.add_node(node_id + '+', sequence=node_seq, gc=gc, length=length, cov=cov)
                #gg_network.add_node(node_id + '-', sequence=node_seq, gc=gc, length=length, cov=cov)
                gg_network.add_node(node_id, sequence=node_seq, gc=gc, length=length, cov=cov, gfa_S=line,start=0, end=0)

        if gfa_file_path.endswith(".gz"):
            import gzip
            gfa_file_obj = gzip.open(gfa_file_path, "rt")
        else:
            gfa_file_obj = open(gfa_file_path, 'r')
            
        for line in gfa_file_obj:

            if line[0] == 'L':
                # link / edge node
                edge_list = line.split('\t')
                node_source = edge_list[1] #+ edge_list[2]
                node_target = edge_list[3] #+ edge_list[4] 

                edge_source_oriant = edge_list[2]
                edge_target_oriant = edge_list[4]
                
                # count number of connection ad each ends
                # if oriantation == + normal
                # if oriantation == - reverse complement ==> start is end and end is start
                if edge_source_oriant == '+':
                    gg_network.nodes[node_source]["start"] += 1
                else:
                    gg_network.nodes[node_source]["end"] += 1

                if edge_target_oriant == '+':
                    gg_network.nodes[node_target]["end"] += 1
                else:
                    gg_network.nodes[node_target]["start"] += 1
                       
                #if edge_list[2] == '+':
                gg_network.add_edge(node_source, node_target, gfa_L=line)
                #gg_network.add_edge(edge_target, edge_source, )
                #else:
                #gg_network.add_edge( edge_target, edge_source, gfa_L=line)
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
        cumul_size_filt = sum([ self.assembly_nx.nodes[node]["length"] for node in subgraph if self.assembly_nx.nodes[node]["length"] > node_len_cutoff])
        gc_data = [self.assembly_nx.nodes[node]["gc"] for node in subgraph if self.assembly_nx.nodes[node]["length"] > node_len_cutoff]
        if len(gc_data) > 0:
            gc_median = statistics.median(gc_data)
        else:
            gc_median = None
        depth_data = [self.assembly_nx.nodes[node]["cov"] for node in subgraph if self.assembly_nx.nodes[node]["length"] > node_len_cutoff]
        length_data = [self.assembly_nx.nodes[node]["length"] for node in subgraph if self.assembly_nx.nodes[node]["length"] > node_len_cutoff]
        if len(depth_data) > 0:
            depth_mean = np.average(depth_data, weights=length_data)
        else:
            depth_mean = None
        
        return cumul_size, cumul_size_filt, gc_median, depth_mean
    
    def count_dead_ends(self, network):
        print(type(network), network)
        de = 0
        for node in network.nodes:
            if network.nodes[node]["start"] == 0:
                de += 1
            if network.nodes[node]["end"] == 0:
                de += 1
        return de
    
    
    def save_subgraph(self, subgraph, node_len_cutoff):
        nodes = [re.match(r'\d+', node).group(0) for node in subgraph]
         
        self.subgraphs[self.subgraphs_index] = {}
        
        self.subgraphs[self.subgraphs_index]["nodes"] = nodes
        cumul_size, cumul_size_filt, gc_median, depth_mean = self.get_path_stats(subgraph, node_len_cutoff=node_len_cutoff)
        
        
        # retrieve consensus taxonomy
        rank2top_taxon = self.get_consensus_taxonomy(subgraph) 
        for rank in rank2top_taxon:   
            self.subgraphs[self.subgraphs_index][f"{rank}"] = rank2top_taxon[rank][0]
            self.subgraphs[self.subgraphs_index][f"{rank}_percent_bp"] = rank2top_taxon[rank][1]
        
        self.subgraphs[self.subgraphs_index]["network"] = self.assembly_nx.subgraph(subgraph).copy()
        
        self.subgraphs[self.subgraphs_index]["cumul_size"] = cumul_size
        self.subgraphs[self.subgraphs_index][f"cumul_size_{node_len_cutoff}"] = cumul_size_filt
        self.subgraphs[self.subgraphs_index]["gc_median"] = gc_median
        self.subgraphs[self.subgraphs_index]["depth_mean"] = depth_mean
        
        self.subgraphs_index += 1


    def _get_spaced_colors(self, n):
        
        '''
        Return n spaces colors color
        '''

        max_value = 16581375 #255**3
        interval = int(max_value / n)
        colors = [hex(I)[2:].zfill(6) for I in range(0, max_value, interval)]

        return ['#%02x%02x%02x' % (int(i[:2], 16), int(i[2:4], 16), int(i[4:], 16)) for i in colors]

    def print_color_file(self, filename):
        with open(filename, "w") as f:
            # bandage header
            f.write(f'node,color\n') 
            cols = self._get_spaced_colors(len(self.subgraphs))
            for n,subgraph in enumerate(self.subgraphs):
                graph = self.subgraphs[subgraph]["network"]               
                col = cols[n]
                # save color
                self.subgraphs[subgraph]["col"] = col
                for node in graph.nodes:
                    f.write(f'{node},"{col}"\n')
    
    def _match_taxon(self, row, taxon_list):
        for taxon in taxon_list:
            if taxon in row.lineage:
                return taxon
    
    def print_color_taget_taxon(self, 
                                filename, 
                                target_taxon="Chlamydiota",
                                target_rank='p__',
                                node_list=False):
    
        # filter target 
        if node_list:
            df_target = self.contig_summary.loc[node_list].dropna(subset=['lineage'])
        else:
            df_target = self.contig_summary.dropna(subset=['lineage'])
        
        df_match = df_target[self.contig_summary[target_rank] == target_taxon]
        df_match_other_taxa = df_target[self.contig_summary[target_rank] != target_taxon]
        
        #filter = '|'.join(target_taxons)
        #df_target = df_target[self.contig_summary['lineage'].astype(str).str.contains(filter)]
        #df_target["match"] = df_target.apply(lambda row: self._match_taxon(row, taxon_list=target_taxons), axis=1)
        
        with open(filename, "w") as f:
            # bandage header
            f.write(f'node,color,taxon\n') 
            for n,row in df_match.iterrows():
                f.write(f'{n},"red",{row[target_rank]}\n')  
            for n,row in df_match_other_taxa.iterrows():
                if row[target_rank] == 'failed':
                    f.write(f'{n},"blue",{row[target_rank]}\n')  
                else:
                    f.write(f'{n},"green",{row[target_rank]}\n')  
    
    
    def write_fasta_subgraphs(self, prefix, taxon_list, taxon_rank="p__", tar=False):
        '''
        Write entire subgraph or only circular
        '''
        import tarfile
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        
        if tar:
            tar_arch = tarfile.open(f'{prefix}_replicons.fasta.tar.gz', 'w:gz')
        
        for n, graph_index in enumerate(self.subgraphs):
            if self.subgraphs[graph_index][f"{taxon_rank}"] not in taxon_list:
                continue
            graph = self.subgraphs[graph_index]["network"]
            #contigs = set([self.node2contig[node.strip("-").strip("+")] for node in graph if node.strip("-").strip("+") in self.node2contig ])
            #seqs = [self.contig2seq[contig.strip("'")] for contig in contigs]
            if not self.contig2rec:
                seqs = [SeqRecord(Seq(self.assembly_nx.nodes[node]["sequence"]), id=f"{node}", name="", description="") for node in graph]
            else:
                # if spades contig file
                # convert to contigs
                seqs = self.get_contigs_from_nodes([re.sub("[\+\-]", "", node) for node in self.subgraphs[graph_index]["nodes"]]) 
            SeqIO.write(seqs, f'{prefix}_{self.subgraphs[graph_index][f"{taxon_rank}"]}_{graph_index + 1}_chromosome.fasta', "fasta")
            if tar:
                tar_arch.add(f"{prefix}_{graph_index + 1}.fasta")
        if tar:
            tar_arch.close()

    def extract_unconnected_contigs(self, 
                                    taxon_list, 
                                    taxon_rank="p__",
                                    min_size=1000,
                                    prefix="out"):

        
        for taxon in taxon_list:
            nodes_in_subgraphs = set()
            for n, graph_index in enumerate(self.subgraphs):
                if self.subgraphs[graph_index][f"{taxon_rank}"] != taxon:
                    continue
                graph_nodes = set(self.subgraphs[graph_index]["nodes"])
                nodes_in_subgraphs |= graph_nodes
            
            contig_summary = self.contig_summary[self.contig_summary[taxon_rank].isin(taxon_list)]
            
            large_unconnected_contigs = contig_summary.loc[contig_summary.index.difference(nodes_in_subgraphs)]
            
            large_unconnected_contigs = large_unconnected_contigs[large_unconnected_contigs["Length"] > min_size]
            
            self.taxon2unconnected_contigs[taxon] = large_unconnected_contigs
                
            large_unconnected_contigs.to_csv(f"{prefix}_{taxon}_unconnected_contigs.tsv")
            
            if not self.contig2rec:
                seqs = [SeqRecord(Seq(self.assembly_nx.nodes[node]["sequence"]), id=f"{node}", name="", description="") for node in large_unconnected_contigs.index]
            else:
                # if spades contig file
                # convert to contigs
                seqs = self.get_contigs_from_nodes([re.sub("[\+\-]", "", node) for node in large_unconnected_contigs.index]) 
            
            SeqIO.write(seqs, f"{prefix}_{taxon}_unconnected_contigs.fasta", "fasta")
            
        


    def summary_subgraphs(self, filename, min_node_len=1000, sample_name="sample"):
        with open(filename, "w") as f:
            # bandage header
            taxo = [val for pair in zip(self.ranks, [f"{rank}_percent_bp" for rank in self.ranks]) for val in pair]
            
            h = ["sample", "graph_id", "n_nodes", f"n_nodes_{min_node_len}", "n_contigs",
                 f"n_contigs_{min_node_len}", 
                 "cumul_size", 
                 f"cumul_size_{min_node_len}", 
                 f"N50_{min_node_len}",
                 "cumul_size_spades",
                 f"cumul_size_spades_{min_node_len}", 
                 f"N50_spades_{min_node_len}",
                 f"cumul_size_same_phylum",
                 f"percent_same_phylum",
                 "gc_median", 
                 "depth_mean",
                 "unconnected_contigs_same_phylum", 
                 "unconnected_contigs_same_phylum_length", 
                 "unconnected_contigs_same_phylum_mean_depth",
                 "n_dead_ends",
                 "n_contigs_conflicting_taxonomy",
                 "n_contigs_no_taxonomy",
                 "col"] + taxo
            f.write('\t'.join(h) + '\n') 
            for subgraph in self.subgraphs:
                
                dead_ends = self.count_dead_ends(self.subgraphs[subgraph]["network"])
                
                if self.contig2rec:
                    contigs = self.get_contigs_from_nodes([re.sub("[\+\-]", "", node) for node in self.subgraphs[subgraph]["nodes"]])
                    contigs_filter = [i for i in contigs if len(i) > min_node_len]
                    spades_contigs_len = [len(i) for i in contigs_filter] # only contigs larger than `min_node_len`
                    spades_cumul_size =  sum(spades_contigs_len)
                    spades_N50 = calculate_N50(spades_contigs_len)
                    spades_cumul_size_filt =  sum([len(i) for i in contigs_filter])
                    n_contigs_spades = len(contigs)
                    n_contigs_conflicting_taxonomy = sum([1 for node in self.subgraphs[subgraph]["nodes"] if self.assembly_nx.nodes[node]["p__"] != self.subgraphs[subgraph]["p__"] if self.assembly_nx.nodes[node]["p__"] is not None])
                    n_contigs_no_taxonomy = sum([1 for node in self.subgraphs[subgraph]["nodes"] if self.assembly_nx.nodes[node]['p__'] is None])
                    
                    if self.subgraphs[subgraph]["p__"] is not None:
                        print(self.subgraphs[subgraph]["p__"], type(self.subgraphs[subgraph]["p__"]))
                        cumul_size_same_rank = self.phylum2cumul_size[self.subgraphs[subgraph]["p__"]]
                        fraction_cumul_size = round(self.subgraphs[subgraph][f"cumul_size_{min_node_len}"] / cumul_size_same_rank * 100, 2)
                    else:
                        cumul_size_same_rank = None
                        fraction_cumul_size = None
                    
                else:
                    n_contigs_spades = None
                    contigs_filter = []
                    spades_cumul_size = None
                    spades_cumul_size_filt = None
                    spades_N50 = None
                    
                if self.subgraphs[subgraph]["p__"] in self.taxon2unconnected_contigs:
                    unconnected_same_taxon = self.taxon2unconnected_contigs[self.subgraphs[subgraph]["p__"]]
                    n_unconnected_same_taxon = len(unconnected_same_taxon)
                    cumul_length_unconnected_same_taxon = unconnected_same_taxon["Length"].sum()
                    if len(unconnected_same_taxon["Length"]) > 0:
                        mean_depth_unconnected_same_taxon = np.average(unconnected_same_taxon["Mean"], weights=unconnected_same_taxon["Length"]) 
                    else:
                        mean_depth_unconnected_same_taxon = None
                else:
                    n_unconnected_same_taxon = None
                    cumul_length_unconnected_same_taxon = None
                    mean_depth_unconnected_same_taxon = None
                
                node_length_list = [ int(self.assembly_nx.nodes[node]["length"]) for node in self.subgraphs[subgraph]["nodes"] if self.assembly_nx.nodes[node]["length"] > min_node_len]
                print("node_length_list", node_length_list)
                N50_nodes = calculate_N50(node_length_list)
                
                
                row = [sample_name,
                       subgraph,
                       len(self.subgraphs[subgraph]["nodes"]),
                       len([node  for node in self.subgraphs[subgraph]["nodes"] if self.assembly_nx.nodes[node]["length"] > min_node_len]),
                       n_contigs_spades,
                       len(contigs_filter),
                       self.subgraphs[subgraph]["cumul_size"] ,
                       self.subgraphs[subgraph][f"cumul_size_{min_node_len}"] ,
                       N50_nodes,
                       spades_cumul_size,
                       spades_cumul_size_filt,
                       spades_N50,
                       cumul_size_same_rank,
                       fraction_cumul_size,
                       self.subgraphs[subgraph]["gc_median"] ,
                       self.subgraphs[subgraph]["depth_mean"],
                       n_unconnected_same_taxon,
                       cumul_length_unconnected_same_taxon,
                       mean_depth_unconnected_same_taxon,
                       dead_ends,
                       n_contigs_conflicting_taxonomy,
                       n_contigs_no_taxonomy]
                if 'col' in self.subgraphs[subgraph]:
                    row.append(self.subgraphs[subgraph]["col"])
                else:
                    row.append(None)
                for rank in self.ranks:
                    row += [self.subgraphs[subgraph][rank], 
                            self.subgraphs[subgraph][f"{rank}_percent_bp"]]
                f.write('\t'.join([str(i) for i in row]) + '\n')

    def save_target_taxon_subgraph_nodes(self, 
                                        filename, 
                                        target_rank="p__", 
                                        target_taxon="Chlamydiota"):
        
        with open(filename, "w") as f:
            for subgraph in self.subgraphs:
                if self.subgraphs[subgraph][target_rank] == target_taxon:
                    for node in self.subgraphs[subgraph]["nodes"]:
                        f.write(f'{subgraph}\t{node}\n')
           

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(
        description="Parse ssearch results",
    )
    parser.add_argument(
        "-f", '--fasta',
        type=str,
        help="fasta file",
    )
    parser.add_argument(
        "-t", '--taxonomy',
        type=str,
        help="contigs2taxnonomy CAT BAT file",
    )
    parser.add_argument(
        "-c", '--coverm',
        type=str,
        help="coverm table",
    )
    parser.add_argument(
        "-g", '--gfa',
        type=str,
        help="gfa file",
    )
    parser.add_argument(
        "-p", '--paths',
        type=str,
        help="Spades contigs or scaffolds paths file",
    ) 
    parser.add_argument(
        "-s", '--spades',
        type=str,
        help="Spades contig or scaffold fasta",
    ) 
    parser.add_argument(
        "-o", '--output',
        type=str,
        help="output_name",
    )
        
    args = parser.parse_args()
    
    #summary(args.fasta, args.output, args.identity_cutoff)
    
    A = Assembly(args.gfa,
                 args.coverm,
                 args.taxonomy,
                 args.fasta,
                 "p__",
                 args.output)
    
    # search subgraphs of min 100kb
    A.extract_largest_components(cumul_size_cutoff=100000, 
                                 min_contig_size=1000)
    # print color file with all large subgraphs (for plotting with bandage or other)
    A.print_color_file(f"{args.output}_subgraphs_cols.csv")
    # if spades graph, get mapping with contigs 
    # will output fasta with contigs rather than nodes
    if args.paths and args.spades:
        A.parse_spades_paths(args.paths, args.spades)
        
    A.extract_unconnected_contigs(taxon_list=["Chlamydiota"], 
                                  taxon_rank="p__",
                                  min_size=1000,
                                  prefix=args.output)
    
    # summary table with statistics for subgraphs of min 100kb (all taxons)
    A.summary_subgraphs(f"{args.output}_subgraphs_summary.tsv", 
                        min_node_len=1000,
                        sample_name = args.output)
    
    # color file for bandage with all subgraphs >100kb
    A.save_target_taxon_subgraph_nodes(f"{args.output}_taxon_subgraphs_nodes.tsv")
    

    
    node_list = A.extract_subgraph(["Chlamydiota"], 
                       f"{args.output}")

    # all nodes with the taregt taxonomy
    A.print_color_taget_taxon(f"{args.output}_target_cols.csv", 
                              target_taxon="Chlamydiota",
                              target_rank='p__',
                              node_list = node_list)

    # save large subgraphs with target taxonomy
    # WARNING: will miss plasmids! (too small, not part of the main graph)
    # todo export large contigs with the expected taxonomy (could also search for circular paths)
    A.write_fasta_subgraphs(prefix=args.output, 
                            taxon_list=["Chlamydiota"], 
                            taxon_rank="p__")
