#! /usr/bin/env python


import os, stat, os.path, time, datetime, subprocess
import re, sys, glob, argparse, csv, shutil, fileinput
from collections import defaultdict
from pprint import pprint
from itertools import islice
import ConfigParser

print("this is zyxue's edits")

def __main__():
    print "Quality and somatic filtering script starts at: %s\n" % datetime.datetime.now()     
    parser = argparse.ArgumentParser(description='Filter variants based on qulaity and somatic filters')
    parser.add_argument('-i','--input_file', help='specify input file', required=True)
    args = vars(parser.parse_args())

    # variant_input_file = "bam_vcf_cnv_path.txt" 
    variant_summary = args['input_file']
    filtered_summary = ".".join([variant_summary, "filtered" ])
    filter_variants(variant_summary, filtered_summary)

    # recalculate occurrence for filtered summary
    out = make_occurrence_dict(filtered_summary)
    gene_patients_occur = out[0]
    var_patients_occur = out[1]
    gene_variants_occur = out[2]
    calculate_occurrence(filtered_summary, gene_patients_occur, var_patients_occur, gene_variants_occur)


def filter_variants(variant_summary, filtered_summary):
    """ d = gene_variant_patients dictionary """
    d = dict()
    with open (filtered_summary,  'wb') as fh:
        writer = csv.writer( fh, delimiter='\t' )
        print "The variant summary file is: %s.\n" % variant_summary
        with open (variant_summary, 'r') as handle:
             records = csv.DictReader(handle,  delimiter='\t')
             headers = records.fieldnames
             writer.writerow(headers)
             for line in records:
                 #print line
                 gene = line['gene']
                 chr = line["chromosome"]
                 pos = line["position"]
                 ref = line["ref_base"]
                 alt = line["alt_base"]
                 in_strelka = line["in_strelka"]
                 content = [line[i] for i in headers]
                 print in_strelka
                 # strelka makes high confidence calls
                 if (in_strelka == 'in_strelka'):
                     print "in strelka, so keep this variant!"
                     writer.writerow(content)
    return [filtered_summary]

def make_occurrence_dict(infile):
    gene_patients = {}
    var_patients = {}
    gene_variants = {}
    with open (infile) as fh:
        records = csv.DictReader(fh,  delimiter='\t')
        for line in records:
            gene = line['gene']
            chr = line["chromosome"]
            pos = line["position"]
            ref = line["ref_base"]
            alt = line["alt_base"]
            patient = line["patient_ID"].split('_')[0]
            variant = '_'.join([gene, chr, pos, ref, alt])
            try:
                gene_patients[gene].append(patient)
            except KeyError:
                gene_patients[gene] = [patient]

            try:
                var_patients[variant].append(patient)
            except KeyError:
                var_patients[variant] = [patient]

            try:
                gene_variants[gene].append(variant)
            except KeyError:
                gene_variants[gene] = [variant]
    # remove duplicate items
    for gene in gene_patients:
        gene_patients[gene] = len(list(set(gene_patients[gene])))  
    #print gene_patients

    for variant in var_patients:
        var_patients[variant] = len(list(set(var_patients[variant])))
    #pprint(var_patients)

    for gene in gene_variants:
        gene_variants[gene] = len(list(set(gene_variants[gene])))
    #print gene_variants
    return [gene_patients, var_patients, gene_variants]


def calculate_occurrence(infile, gene_patients, var_patients, gene_variants):
    outfile = '.'.join([infile, "final", "txt"])
    with open (outfile,  'wb') as fh:
        writer = csv.writer( fh, delimiter='\t' )
        with open (infile, 'r') as handle:
             records = csv.DictReader(handle,  delimiter='\t')
             headers = records.fieldnames
             new_headers = [headers[0]] + ['filtered_num_patients_gene_level', 'filtered_num_INDELs_gene_level'] + headers[1:7] + ['filtered_num_patients_indel_level'] + headers[7:] 
             writer.writerow(new_headers)
             for line in records:
                 gene = line['gene']
                 chr = line["chromosome"]
                 pos = line["position"]
                 ref = line["ref_base"]
                 alt = line["alt_base"]
                 variant = '_'.join([gene, chr, pos, ref, alt])
                 num_gene_patients = gene_patients[gene]
                 num_var_patients = var_patients[variant]
                 num_gene_variants = gene_variants[gene]
                 content = [line[i] for i in headers]
                 final_content = [content[0]] + [num_gene_patients, num_gene_variants] + content[1:7] + [num_var_patients] + content[7:]
                 writer.writerow(final_content)
if __name__ == '__main__':
    __main__()

