#!/opt/conda/bin/python

from Bio import SeqIO
from Bio.SeqUtils import GC
import networkx as nx
import statistics
import re 
import pandas


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
        self.ranks = ["p__", "o__","f__","g__"]


        for node in self.assembly_nx.nodes:
            self.assembly_nx.nodes[node]["cov"] = self.coverm.loc[re.sub("\-|\+", "",node)]["Mean"]
        
        
        # root;d__Bacteria;p__Actinobacteriota;c__Actinomycetia;o__Mycobacteriales;f__Mycobacteriaceae;g__Rhodococcus_B;s__Rhodococcus_B
        self.cat_taxnonomy = pandas.read_csv(contig2taxonomy, sep="\t")
        self.cat_taxnonomy["# contig"] = self.cat_taxnonomy["# contig"].astype(str)
        self.cat_taxnonomy = self.cat_taxnonomy.set_index("# contig")

        #self.cat_taxnonomy["lineage"] = self.cat_taxnonomy.apply(lambda row: self._format_gtdb_lineage(row, target_rank=summary_rank), axis=1)
        self.cat_taxnonomy["lineage_format"] = self.cat_taxnonomy.apply(lambda row: self._format_gtdb_lineage(row, target_rank=summary_rank), axis=1)
        #self.cat_taxnonomy.apply(lambda row: row["lineage"].split(";") if not pandas.isna(row["lineage"]) else '', axis=1)
        
        # parse fasta
        self.fasta = [i for i in SeqIO.parse(fasta, "fasta")]
        self.gc_data = pandas.DataFrame([(record.name, round(GC(record.seq), 2), len(record.seq)) for record in self.fasta], columns=["contig", "GC", "length"]).set_index("contig")
        for node in self.assembly_nx.nodes:
            self.assembly_nx.nodes[node]["gc"] = self.gc_data.loc[re.sub("\-|\+", "",node)]["GC"]
            self.assembly_nx.nodes[node]["length"] = self.gc_data.loc[re.sub("\-|\+", "",node)]["length"]
            #print(type(self.gc_data.loc[re.sub("\-|\+", "",node)]["length"]))
        self.contig_summary = self.coverm.join(self.cat_taxnonomy).join(self.gc_data)
        self.contig_summary.to_csv(f"{output_prefix}_contigs.tsv", sep="\t")
        

    def _format_gtdb_lineage(self, row, target_rank="p__"):
        if row.classification == 'no taxid assigned':
            return None
        else:
            try:
                classif = row.lineage.split(target_rank)[1].split(";")[0]
            except:
                classif = row.lineage.split(";")[-1].split("__")[-1]
            return classif
    
    def _extract_taxon(self, string, rank='p__'):
        import re 
        string = str(string)

        # root;d__Bacteria;p__Actinobacteriota;c__Actinomycetia;o__Mycobacteriales;f__Mycobacteriaceae;g__Rhodococcus_B;s__Rhodococcus_B
        result = re.search(f"{rank}([A-Za-z]+)", string)
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
                if current_node == "100121":
                    print("n neigbors 100121", len(neighbors), neighbors)
                # add any results we don't already know about to our list of known nodes:
                current_nodes |= neighbors  
                # save our results for this node:
                node_list |= neighbors
                # finally, mark that we checked the node so we don't check it again:
                checked_nodes.add(current_node)
            subgraph = self.assembly_nx.subgraph(node_list)
            self.to_gfa(subgraph, f"{out_prefix}_{taxon}.gfa")
        
    def to_gfa(self, graph, output_name):
        with open(output_name, "w") as f:
            for node in graph.nodes:
                f.write(graph.nodes[node]["gfa_S"])
            for edge in graph.edges:
                f.write(graph.edges[edge]["gfa_L"])
        
        
    
    
    def extract_largest_components(self, 
                                   cumul_size_cutoff=10000, 
                                   stats_cutoff=0, 
                                   filter_depth_min_factor=0.1, 
                                   max_depth=2000000000000000000000):
        
        components = list(nx.connected_components( self.assembly_nx))
        
        large_count = 0
        for n,comp in enumerate(components):

            cumul_size, gc_median, depth_median = self.get_path_stats(comp)
            #print(f"One comp {n}:", len(comp), cumul_size, gc_median, depth_median)
            
            if cumul_size > cumul_size_cutoff:
                large_count += 1
                cumul_size, gc_median, depth_median = self.get_path_stats(comp, stats_cutoff)

                self.save_subgraph(comp)
    
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
                gg_network.add_node(node_id, sequence=node_seq, gc=gc, length=length, cov=cov, gfa_S=line)

        if gfa_file_path.endswith(".gz"):
            import gzip
            gfa_file_obj = gzip.open(gfa_file_path, "rt")
        else:
            gfa_file_obj = open(gfa_file_path, 'r')
            
        for line in gfa_file_obj:

            if line[0] == 'L':
                # link / edge node
                edge_list = line.split('\t')
                edge_source = edge_list[1] #+ edge_list[2]
                edge_target = edge_list[3] #+ edge_list[4] 


                '''
                if edge_list[2] == '+':
                    edge_source_rev = edge_list[1] + '-'
                else:
                    edge_source_rev = edge_list[1] + '+'
                if edge_list[4] == '+':
                    edge_target_rev = edge_list[3] + '-'
                else:
                    edge_target_rev = edge_list[3] + '+'
                '''          
                #if edge_list[2] == '+':
                gg_network.add_edge(edge_source, edge_target, gfa_L=line)
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
        gc_median = statistics.median([self.assembly_nx.nodes[node]["gc"] for node in subgraph if self.assembly_nx.nodes[node]["length"] > node_len_cutoff])
        depth_median = statistics.median([self.assembly_nx.nodes[node]["cov"] for node in subgraph if self.assembly_nx.nodes[node]["length"] > node_len_cutoff])
        
        return cumul_size, gc_median, depth_median
    
    def save_subgraph(self, subgraph):
        nodes = [re.match(r'\d+', node).group(0) for node in subgraph]       
        self.subgraphs[self.subgraphs_index] = {}
        self.subgraphs[self.subgraphs_index]["network"] = subgraph
        self.subgraphs[self.subgraphs_index]["nodes"] = nodes
        cumul_size, gc_median, depth_median = self.get_path_stats(subgraph)
        rank2top_taxon = self.get_consensus_taxonomy(subgraph)
        # retrieve consensus taxonomy
        for rank in rank2top_taxon:   
            self.subgraphs[self.subgraphs_index][f"{rank}"] = rank2top_taxon[rank][0]
            self.subgraphs[self.subgraphs_index][f"{rank}_percent_bp"] = rank2top_taxon[rank][1]
        self.subgraphs[self.subgraphs_index]["cumul_size"] = cumul_size
        self.subgraphs[self.subgraphs_index]["gc_median"] = gc_median
        self.subgraphs[self.subgraphs_index]["depth_median"] = depth_median
        
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
                for node in graph:
                    f.write(f'{node},"{col}"\n')
    
    def _match_taxon(self, row, taxon_list):
        for taxon in taxon_list:
            if taxon in row.lineage:
                return taxon
    
    
    def print_color_taget_taxon(self, filename, target_taxons=["Chlamydiota", "Brevundimonas"]):
    
        filter = '|'.join(target_taxons)
        df_target = self.contig_summary.dropna(subset=['lineage'])
        df_target = df_target[self.contig_summary['lineage'].astype(str).str.contains(filter)]
        df_target["match"] = df_target.apply(lambda row: self._match_taxon(row, taxon_list=target_taxons), axis=1)
        
        nr_taxons =df_target["match"].drop_duplicates().to_list()
        if len(nr_taxons) == 1:
            cols = ["red"]
        elif len(nr_taxons) == 2:
            cols = ["red", "blue"]
        elif len(nr_taxons) == 3:
            cols = ["red", "blue", "green"]
        else:  
            cols = self._get_spaced_colors(len(nr_taxons))
        tax2col = {tax:col for tax,col in zip(nr_taxons, cols)}
        with open(filename, "w") as f:
            # bandage header
            f.write(f'node,color,taxon\n') 
            for n,row in df_target.iterrows():
                f.write(f'{n},"{tax2col[row.match]}",{row.match}\n')     
    
    def write_fasta_subgraphs(self, prefix, tar=False):
        '''
        Write entire subgraph or only circular
        '''
        import tarfile
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        
        if tar:
            tar_arch = tarfile.open(f'{prefix}_replicons.fasta.tar.gz', 'w:gz')
        
        for n, graph_index in enumerate(self.subgraphs):
            graph = self.subgraphs[graph_index]["network"]
            #contigs = set([self.node2contig[node.strip("-").strip("+")] for node in graph if node.strip("-").strip("+") in self.node2contig ])
            #seqs = [self.contig2seq[contig.strip("'")] for contig in contigs]
            
            seqs = [SeqRecord(Seq(self.assembly_nx.nodes[node]["sequence"]), id=f"{node}", name="", description="") for node in graph]
            SeqIO.write(seqs, f"{prefix}_{graph_index + 1}.fasta", "fasta")
            if tar:
                tar_arch.add(f"{prefix}_{graph_index + 1}.fasta")
        tar_arch.close()

    def summary_subgraphs(self, filename):
        with open(filename, "w") as f:
            # bandage header
            taxo = [val for pair in zip(self.ranks, [f"{rank}_percent_bp" for rank in self.ranks]) for val in pair]
            
            h = ["graph_id", "n_nodes", "cumul_size", "gc_median", "depth_median", "col"] + taxo
            f.write('\t'.join(h) + '\n') 
            for subgraph in self.subgraphs:
                row = [subgraph,
                       len(self.subgraphs[subgraph]["nodes"]),
                       self.subgraphs[subgraph]["cumul_size"] ,
                       self.subgraphs[subgraph]["gc_median"] ,
                       self.subgraphs[subgraph]["depth_median"] ]
                if 'col' in self.subgraphs[subgraph]:
                    row.append(self.subgraphs[subgraph]["col"])
                else:
                    row.append(None)
                for rank in self.ranks:
                    row += [self.subgraphs[subgraph][rank], 
                            self.subgraphs[subgraph][f"{rank}_percent_bp"]]
                f.write('\t'.join([str(i) for i in row]) + '\n')
                
    

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
    
    # extract subgraphs of min 100kb
    A.extract_largest_components(cumul_size_cutoff=100000)
    A.print_color_file(f"{args.output}_subgraphs_cols.csv")
    A.summary_subgraphs(f"{args.output}_subgraphs_summary.tsv")
    A.print_color_taget_taxon(f"{args.output}_target_cols.csv")
    A.extract_subgraph(["Chlamydiota", "Brevundimonas"], f"{args.output}")
