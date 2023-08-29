#!/usr/bin/env python

from Bio import SeqIO
from Bio import SeqUtils

import matplotlib.pyplot as plt
import numpy
import sqlite3


# taxids to keep:    PVC, Chlamydiae]
# to_keep = [1, 131567, 2, 1783257, 204428]
to_keep = ["root", "d__Bacteria", "p__Chlamydiota"]


def extract_lineage(cat_bat_output, output_file):
    cat_bat = open(cat_bat_output, "r")
    output = open(output_file, "w")
    for line_no, line in enumerate(cat_bat):
        if line_no==0:
            # skip header
            continue

        tokens = line.split("\t")
        contig = tokens[0]
        if tokens[1]=="no taxid assigned":
            output.write(f"{tokens[0]}\tunclassified\n")
            continue

        # taxids = (int(i) for i in tokens[2].split(";"))
        taxids = tokens[3].split(";")
        print(taxids)
        i, keep = 0, True
        for taxid, to_check in zip(taxids, to_keep):
            i += 1
            if taxid!=to_check:
                keep = False
                break

        if keep and i==len(to_keep):
            output.write(f"{tokens[0]}\tkeep\n")
        elif keep:
            output.write(f"{tokens[0]}\tcompatible\n")
        else:
            output.write(f"{tokens[0]}\tother\n")
    output.close()


# returns a list of tuples:
# (contig length, average coverage)
def compute_cov_to_length(sam_depths_file, fasta_file, output_file=None):

    # maps name to (GC content/length)
    hsh_contig_data = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        hsh_contig_data[record.id] = (len(record.seq), SeqUtils.GC(record.seq))

    f = open(sam_depths_file, "r")
    hsh_contig_to_coverage = {}
    for line in f:
        contig_name, pos, coverage = line.split("\t")
        if contig_name in hsh_contig_to_coverage:
            hsh_contig_to_coverage[contig_name] += int(coverage)
        else:
            hsh_contig_to_coverage[contig_name] = int(coverage)

    tuples = []

    to_print = []
    for contig_name, ttl_coverage in hsh_contig_to_coverage.items():
        length, gc_content = hsh_contig_data[contig_name]
        coverage = int(ttl_coverage/length)
        tuples.append((length, int(coverage), gc_content))
        to_print.append((contig_name, str(length), str(coverage), str(gc_content)))

    if output_file!=None:
        output = open(output_file, "w")
        output.write("\t".join(["Contig name", "Length", "Coverage", "GC content"]) + "\n")
        output.write("\n".join(["\t".join(i) for i in to_print]))
    return tuples


def merge_contigs_infos(coverage, classifications, output):
    hsh_contig_name_to_info = {}
    to_print = []
    f = open(coverage, "r")
    f.readline()

    line = f.readline()
    while line:
        contig_name, length, coverage, gc = [i.rstrip() for i in line.split("\t")]
        hsh_contig_name_to_info[contig_name] = (length, coverage, gc)
        line = f.readline()

    cat_bat_file = open(classifications, "r")
    line = cat_bat_file.readline()
    while line:
        contig, classification = [i.rstrip() for i in line.split("\t")]
        length, coverage, gc = hsh_contig_name_to_info[contig]
        to_print.append((contig, coverage, length, gc, classification))
        line = cat_bat_file.readline()

    output_file = open(output, "w")
    output_file.write("\t".join(["contig name", "coverage", "length", "gc", "classification"]) + "\n")
    output_file.write("\n".join(["\t".join(i) for i in to_print]))


# returns a dict (String->String[]) mapping contig name
# to node name (useful to interact with bandage)
# NOTE:
#  in the GFA file, contig names have a _1/2 suffix coding for
#  sense/antisense strand (there may be reads from both, that need
#  to be included in the same contig)
def parse_gfa_file(gfa_file, strip_node_sense=False):
    f = open(gfa_file, "r")
    hsh_mapping = {}

    line = f.readline()
    while line and line[0]!="P":
        line = f.readline()

    # first line with path
    while line and line[0]=="P":
        tokens = line.split("\t")
    
        # ignore the _1/2 suffix (see above)
        contig = tokens[1].strip()[:-2]
        nodes = tokens[2].split(",")

        if strip_node_sense:
            stripped_nodes = [node[:-1] for node in nodes]
            nodes = stripped_nodes

        if contig in hsh_mapping:
            hsh_mapping[contig].extend(nodes)
        else:
            hsh_mapping[contig] = nodes
        line = f.readline()

    return hsh_mapping


def get_fasta_from_graph(gfa_filename, node_list_filename, draft_genome_filename, output_filename):
    hsh_contig_to_node = parse_gfa_file(gfa_filename, True)
    hsh_nodes_to_keep = {}
    lst_contigs = []

    with open(node_list_filename, "r") as node_list_file:
        # usually, only one line
        for line in node_list_file:
            nodes_list = line.split(",")
            for node in nodes_list:
                hsh_nodes_to_keep[node.strip()] = True

    for record in SeqIO.parse(draft_genome_filename, "fasta"):
        nodes = hsh_contig_to_node.get(record.id, None)
        if nodes==None:
            continue

        # keep the contig if at least one node is in node_list
        for node in nodes:
            if node in hsh_nodes_to_keep:
                lst_contigs.append(record)
                break

    SeqIO.write(lst_contigs, output_filename, "fasta")


taxo_to_colour = {
    "keep" : "red",
    "unclassified" : "blue",
    "other": "purple",
    "compatible": "green"}


def gen_bandage_csv_file(contig_to_nodes, merged_info_file, csv_file):
    info_file = open(merged_info_file, "r")
    output = open(csv_file, "w")

    output.write(",".join(["node", "color"])+"\n")

    info_file.readline()
    line = info_file.readline()
    while line:
        contig, cov, length, gc, taxo = line.split("\t")

        # hack for now: the contigs name from the gfa are before repeat resolution 
        # by spades: the contig names may thus have changed. For now, ignore those
        # contigs
        nodes = contig_to_nodes.get(contig, [])
        for node in nodes:
            colour = taxo_to_colour.get(taxo.rstrip(), "yellow")
            output.write(f"{node},{colour}\n")
        line = info_file.readline()

    info_file.close()
    output.close()

def plot_with_taxonomy(file_info):
    f = open(file_info, "r")
    f.readline()

    lengths = []
    gcs = []
    colours = []
    coverages = []
    
    line = f.readline()
    while line:
        contig, cov, length, gc, taxo = line.split("\t")
        colour = taxo_to_colour.get(taxo.rstrip(), "yellow")

        lengths.append(int(length))
        gcs.append(int(float(gc)))
        coverages.append(int(cov))
        colours.append(colour)
        line = f.readline()

    plt.figure()
    plt.scatter(lengths, coverages, c=colours)
    plt.xlabel("Length")
    plt.ylabel("Coverage")
    plt.savefig("Length_vs_coverage.png")

    plt.figure()
    plt.scatter(coverages, gcs, c=colours)
    plt.xlabel("Coverage")
    plt.ylabel("GC")
    plt.savefig("Coverage_vs_GC.png")

# blob-ish code: could easily be refactored and simplified
def gen_coverage_from_taxonomy(coverage_fn, fasta, merged_info_fn, name, gen_from_tree=False):
    contigs = {}

    merged_info = open(merged_info_fn, "r")
    merged_info.readline() # skip header
    line = merged_info.readline()
    gcs = []
    cov_to_gc_ratio = []
    long_not_kept = []
    while line:
        contig, cov, length, gc, classification = [i.rstrip() for i in line.split("\t")]
        if int(length) >= 1000 and (gen_from_tree or classification == "keep"):
            gcs.append(float(gc))
            contigs[contig] = [int(cov), int(length), float(gc), classification]
        elif int(length) >= 1000:
            long_not_kept.append(contig)

        line = merged_info.readline()

    statistics_file = f"{name}_statistics.txt"
    length = sum([ctg_length for ctg, (cov, ctg_length, gc, taxo) in contigs.items()])
    mean_coverage = sum([cov*ctg_length for ctg, (cov, ctg_length, gc, taxo) in contigs.items()])/length
    mean_gc = sum([gc*ctg_length for ctg, (cov, ctg_length, gc, taxo) in contigs.items()])/length
    with open(statistics_file, "w") as f:
        f.write("name,coverage,gc,length\n")
        f.write(f"{name},{mean_coverage},{mean_gc},{length}")

    np_gcs = numpy.array(gcs)
    gc_median = numpy.median(np_gcs)
    gc_std = numpy.std(np_gcs)

    cov_values = []
    coverage = open(coverage_fn, "r")
    line = coverage.readline()
    hsh_contig_to_cov_values = {}
    while line:
        contig, position, cov = [i.rstrip() for i in line.split("\t")]
        if contig in contigs:
            cov_values.append(int(cov))
            if contig in hsh_contig_to_cov_values:
                hsh_contig_to_cov_values[contig].append(int(cov))
            else:
                hsh_contig_to_cov_values[contig] = [int(cov)]
        line = coverage.readline()
    out = open("values.txt", "w")
    out.write("\n".join(map(str, cov_values)))

    np_array = numpy.array(cov_values)
    first_quartile = numpy.percentile(cov_values, 25)
    third_quartile = numpy.percentile(cov_values, 75)
    iqr = third_quartile - first_quartile
    no_outliers = np_array[(np_array > first_quartile - 1.5*iqr) & (np_array < third_quartile + 1.5*iqr)]

    mean = numpy.mean(no_outliers)
    median = numpy.median(no_outliers)
    std = numpy.std(no_outliers)

    to_keep = []
    high_coverage = []
    high_gc = []

    high_coverage_threshold = median+2*std
    high_gc_threshold = gc_median + 2*gc_std
    filtered_merged_info_file = open(f"{name}_merged_filtered_infos", "w")
    to_plot_cov = []
    to_plot_len = []
    to_plot_gc = []
    to_plot_annotate = []
    filtered_merged_info_file.write("Contig\tCoverage\tLength\tGC\tclassification\n")
    contigs_colour = []

    filtered_out = []
    for record in SeqIO.parse(fasta, "fasta"):

        if record.id in long_not_kept:
            filtered_out.append(record)

        # keep only Chlamydia contigs longer than 1kpb
        if record.id not in contigs:
            continue

        cov, length, gc, classification = contigs[record.id]
        to_keep.append(record)
        to_print = [record.id, cov, length, gc, classification]
        filtered_merged_info_file.write("\t".join(map(str, to_print)) + "\n")

        if gc>=high_gc_threshold:
            high_gc.append(record)
        if cov>=high_coverage_threshold:
            high_coverage.append(record)
            to_plot_annotate.append([record.id, length, cov])
        to_plot_gc.append(gc)
        to_plot_cov.append(cov)
        to_plot_len.append(length)

    plt.figure()
    plt.scatter(to_plot_len, to_plot_cov)
    plt.xlabel("Length")
    plt.ylabel("Coverage")
    plt.axhline(y=median-2*std)
    # plt.axhline(y=median+2*std)
    # plt.axvline(x=1000)

    # do not annotate for now: makes some of the graphs unreadable
    # for ctg_name, x, y in to_plot_annotate:
    #    display_name = " ".join(ctg_name.split("_")[0:2])
    #    plt.annotate(display_name, (x, y))
    plt.savefig(f"{name}_length_vs_coverage.png")

    plt.figure()
    plt.scatter(to_plot_gc, to_plot_cov)
    plt.xlabel("GC")
    plt.ylabel("Coverage")
    # plt.axvline(x=high_gc_threshold)
    plt.savefig(f"{name}_gc_vs_coverage.png")

    for ctg_name, x, y in to_plot_annotate:
        lst_values = hsh_contig_to_cov_values[ctg_name]
        plt.figure()
        plt.plot(range(0,len(lst_values)),lst_values, markersize=1)
        plt.savefig(f"{name}_{ctg_name}_cov_per_nucleotide.png")

    filtered_contigs_file = open(f"{name}_scaffold_filtered.fasta", "w")
    high_coverage_file = open(f"{name}_high_coverage.fasta", "w")
    high_gc_file = open(f"{name}_high_gc.fasta", "w")
    SeqIO.write(to_keep, filtered_contigs_file, "fasta")
    SeqIO.write(filtered_out, f"{name}_filtered_out.fasta", "fasta")
    SeqIO.write(high_coverage, high_coverage_file, "fasta")
    SeqIO.write(high_gc, high_gc_file, "fasta")

