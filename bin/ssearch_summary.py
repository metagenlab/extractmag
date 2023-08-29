#!/opt/conda/bin/python

from Bio import SearchIO
import re 
from pathlib import Path






def parse_ssearch(ssearch, output_name, identity_cutoff=50):
    o = open(output_name, "w")
    o.write("sample\tquery\thit\talignment_length\tpercent_identity\tevalue\tbitscore\n")
    sample = ssearch.split("_ssearch.txt")[0]
    records = [i for i in SearchIO.parse(ssearch, 'fasta-m10')]
    for i,record in enumerate(records):
        for n,hit in enumerate(record.hits):
            
            description = hit.hsps[0].query_description
            query_start = hit.hsps[0].query_start
            query_end = hit.hsps[0].query_end
            hit_accession = hit.hsps[0].hit.id
            query_accession = hit.hsps[0].query.id
            alignment_length = hit.hsps[0].aln_span
            
            percent_identity = float(hit.hsps[0].ident_pct)
            evalue = hit.hsps[0].evalue
            bitscore = hit.hsps[0].bitscore
            if percent_identity < identity_cutoff:
                continue
            # 16S_rRNA::NODE_52_length_5592_cov_1.520403:3723-5244(-)
            # 95266;tax=d:Bacteria,p:Bacteroidetes,c:Bacteroidia,o:Bacteroidales,f:Porphyromonadaceae,g:Parabacteroides,s:Parabacteroides
            #depth = re.search("(\d+.\d+)", query_accession.split("_cov_")[1]).group(1)
            o.write(f'{sample}\t{query_accession.split(":")[2]}\t{hit_accession.split("tax=")[1]}\t{alignment_length}\t{percent_identity}\t{evalue}\t{bitscore}\n')

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(
        description="Parse ssearch results",
    )
    parser.add_argument(
        "-i", '--input',
        type=str,
        help="Input ssearch m10",
    )
    parser.add_argument(
        "-o", '--output',
        type=str,
        help="output_name",
    )
    parser.add_argument(
        "-c", '--identity_cutoff',
        type=float,
        help="Identity cutoff for filtering",
        default=50
    )
        
    args = parser.parse_args()
    
    parse_ssearch(args.input, args.output, args.identity_cutoff)