#! /usr/bin/env python

import os
import os.path
import datetime
import argparse
import csv
import ConfigParser
import logging
import colorlog

import datetime
print("now is {0}".format(datetime.datetime.now()))

logger = colorlog.getLogger()
logger.setLevel(logging.DEBUG)
handler = logging.StreamHandler()
handler.setLevel(logging.DEBUG)
handler.setFormatter(
    colorlog.ColoredFormatter('%(log_color)s%(levelname)s:%(name)s:%(message)s'))
logger.addHandler(handler)

def __main__():
    logger.info("Quality and somatic filtering script starts at: %s\n" % datetime.datetime.now())
    parser = argparse.ArgumentParser(description='Filter variants based on qulaity and somatic filters')
    parser.add_argument('-i', '--input_file', help='specify input file', required=True)
    parser.add_argument('-f', '--quality_somatic_filters', help='specify quality and somatic filters', required=True)
    parser.add_argument('-p', '--pairing', help='specify if sample paired with matched normal: paired or unpaired', required=True)
    args = vars(parser.parse_args())
    logger.debug('some debug message')
    logger.info('some info message')
    logger.warning('some warning message')
    logger.error('some error message')
    logger.critical('some critical message')

    # variant_input_file = "bam_vcf_cnv_path.txt" 
    variant_summary = args['input_file']
    filters = args['quality_somatic_filters']
    pairing = args['pairing']
    filtered_summary = ".".join([variant_summary, "filtered.tmp"]).replace('.txt', '')
    somatic_summary = ".".join([filtered_summary, "somatic.tmp"]).replace('.txt', '')
    filter_variants(variant_summary, filters, pairing,
                    filtered_summary, somatic_summary)

    # recalculate occurrence for filtered summary
    out = make_occurrence_dict(filtered_summary)
    gene_patients_occur = out[0]
    var_patients_occur = out[1]
    gene_variants_occur = out[2]
    calculate_occurrence(filtered_summary, gene_patients_occur,
                         var_patients_occur, gene_variants_occur)

    # recalculate occurrence for somatic summary
    if (os.path.isfile(somatic_summary)):
        out = make_occurrence_dict(somatic_summary)
        gene_patients_occur = out[0]
        var_patients_occur = out[1]
        gene_variants_occur = out[2]
        calculate_occurrence(somatic_summary, gene_patients_occur,
                             var_patients_occur, gene_variants_occur)
        os.reomve(somatic_summary)

    # remove intermediate files
    os.remove(filtered_summary)
def filter_variants(variant_summary, filters, pairing,
                    filtered_summary, somatic_summary):

    # get filter values
    config = ConfigParser.SafeConfigParser()
    config.read(filters)

    # default quality filtering: 
    DNA_t_cov_thres = float(config.get('quality_filters', 'DNA_t_cov'))
    DNA_t_altC_thres = float(config.get('quality_filters', 'DNA_t_altC'))
    DNA_t_af_thres = float(config.get('quality_filters', 'DNA_t_af'))

    # OR 
    RNA_t_cov_thres = float(config.get('quality_filters', 'RNA_t_cov'))
    RNA_t_altC_thres = float(config.get('quality_filters', 'RNA_t_altC'))
    RNA_t_af_thres = float(config.get('quality_filters', 'RNA_t_af'))
    #misalignment at exon junctions
    RNA_t_percent_thres = float(config.get(
        'quality_filters', 'RNA_t_altref_total_percent'))
    # third allele
    DNA_t_percent_thres = float(config.get(
        'quality_filters', 'DNA_t_altref_total_percent'))


    # default somatic filters
    DNA_n_af_thres = float(config.get('somatic_filters','DNA_n_af'))
    DNA_n_altC_thres = float(config.get('somatic_filters', 'DNA_n_altC'))
    RNA_n_af_thres = float(config.get('somatic_filters', 'RNA_n_af'))
    RNA_n_altC_thres = float(config.get('somatic_filters', 'RNA_n_altC'))

    """ d = gene_variant_patients dictionary """
    d = dict()
    with open(filtered_summary,  'wb') as fh:
        # with open (somatic_summary,  'wb') as fh2:
        writer = csv.writer(fh, delimiter='\t')
        # writer2 = csv.writer( fh2, delimiter='\t' )
        logger.info("The variant summary file is: %s.\n" % variant_summary)
        with open(variant_summary, 'r') as handle:
            records = csv.DictReader(handle,  delimiter='\t')
            headers = records.fieldnames
            writer.writerow(headers)
            # writer2.writerow(headers)
            if (pairing == 'paired'):
                logger.info('Filtering variant with matching normal!')
                for line in records:
                    for item in headers:
                        if(line[item].isdigit()):
                            item = float(line[item])
                        else:
                            item = line[item]
                    """
                    gene = line['gene']
                    chr = line["chromosome"]
                    pos = line["position"]
                    ref = line["ref_base"]
                    alt = line["alt_base"]
                    in_strelka = line["in_strelka"]
                    try:
                        DNA_t_cov = int(line["t_DNA_cov"])
                    except:
                        DNA_t_cov = line["t_DNA_cov"]
                    try:
                        DNA_t_refC = int(line["t_DNA_RefC"])
                    except:
                        DNA_t_refC = line["t_DNA_RefC"]
                    try:
                        DNA_t_altC = int(line["t_DNA_AltC"])
                    except:
                        DNA_t_altC = line["t_DNA_AltC"]
                    try:
                        DNA_t_af = float(line["t_DNA_AF"])
                    except:
                        DNA_t_af = line["t_DNA_AF"]
                    try:
                        RNA_t_cov = int(line["t_RNA_cov"])
                    except:
                        RNA_t_cov = line["t_RNA_cov"]
                    try:
                        RNA_t_refC = int(line["t_RNA_RefC"])
                    except:
                        RNA_t_refC = line["t_RNA_RefC"]
                    try:
                        RNA_t_altC = int(line["t_RNA_AltC"])
                    except:
                        RNA_t_altC = line["t_RNA_AltC"]
                    try:
                        RNA_t_af = float(line["t_RNA_AF"])
                    except:
                        RNA_t_af = line["t_RNA_AF"]
                    RNA_t_altref_total = RNA_t_refC + RNA_t_altC
                    DNA_t_altref_total = DNA_t_refC + DNA_t_altC
                    """
                    content = [line[i] for i in headers]
                    # if (pairing == 'paired'):
                    with open (somatic_summary,  'wb') as fh2:
                        writer2 = csv.writer(fh2, delimiter='\t')
                        writer2.writerow(headers)
                        """
                        try:
                            DNA_n_cov = int(line["n_DNA_cov"])
                            DNA_n_refC = int(line["n_DNA_RefC"])
                            DNA_n_altC = int(line["n_DNA_AltC"])
                            DNA_n_af = float(line["n_DNA_AF"])
                        except KeyError:
                            DNA_n_cov = "na"
                            # DNA_n_refC = "na"
                            DNA_n_altC = "na"
                            DNA_n_af = "na"
                        try:
                            RNA_n_cov = int(line["n_RNA_cov"])
                            # RNA_n_refC = int(line["n_RNA_RefC"])
                            RNA_n_altC = int(line["n_RNA_AltC"])
                            RNA_n_af = float(line["n_RNA_AF"])
                        except KeyError:
                            RNA_n_cov = line["n_RNA_cov"]
                            # RNA_n_refC = line["n_RNA_RefC"]
                            RNA_n_altC = line["n_RNA_AltC"]
                            RNA_n_af = line["n_RNA_AF"]
                            # RNA_n_cov = "na"
                            # RNA_n_refC = "na"
                            # RNA_n_altC = "na"
                            # RNA_n_af = "na"
                        # content = [line[i] for i in headers]
                        # print type(DNA_t_cov_thres), DNA_t_cov_thres
                        # quality filtering
                        # strelka makes high confidence calls
                        """
                        if(in_strelka == 'in_strelka'):
                            # print "in strelka, so keep this variant!"
                            writer.writerow(content)
                            writer2.writerow(content)
                        else:
                            if (DNA_t_cov == 'na'):
                               logger.info("Only transcriptome is sequenced for this tumor! ")
                               if (RNA_t_cov >= RNA_t_cov_thres and
                                   RNA_t_altC >= RNA_t_altC_thres and
                                   RNA_t_af >= RNA_t_af_thres and
                                   RNA_t_altref_total >= RNA_t_percent_thres*RNA_t_cov):
                                    writer.writerow(content)
                                    if ((RNA_n_af <= RNA_n_af_thres) or (RNA_n_altC <= RNA_n_altC_thres)):
                                        writer2.writerow(content)
                            elif (RNA_t_cov == 'na'):
                                logger.info("Only genome is sequenced for this tumor! ")
                                if (DNA_t_cov >= DNA_t_cov_thres and DNA_t_altC >= DNA_t_altC_thres and DNA_t_af >= DNA_t_af_thres and DNA_t_altref_total >= DNA_t_percent_thres*DNA_t_cov):
                                    writer.writerow(content)
                                    if ((DNA_n_af <= DNA_n_af_thres) or (DNA_n_altC <= DNA_n_altC_thres)):
                                        writer2.writerow(content)
        
                            else:
                                # print "Both genome and transcriptome are sequenced for this tumor! "
                                if ((DNA_t_cov >= DNA_t_cov_thres and DNA_t_altC >= DNA_t_altC_thres and DNA_t_af >= DNA_t_af_thres and DNA_t_altref_total >= DNA_t_percent_thres*DNA_t_cov) or
                                    (RNA_t_cov >= RNA_t_cov_thres and RNA_t_altC >= RNA_t_altC_thres and RNA_t_af >= RNA_t_af_thres and RNA_t_altref_total >= RNA_t_percent_thres*RNA_t_cov)):
                                    writer.writerow(content)
                                    
                                    # somatic filters
                                    if (DNA_n_cov == 'na'):
                                        if ((RNA_n_af <= RNA_n_af_thres) or (RNA_n_altC <= RNA_n_altC_thres)):
                                            writer2.writerow(content)
                                    elif (RNA_n_cov == 'na'):
                                        if ((DNA_n_af <= DNA_n_af_thres) or (DNA_n_altC <= DNA_n_altC_thres)):
                                            writer2.writerow(content)
                                    else:
                                        if ((DNA_n_af <= DNA_n_af_thres) or (DNA_n_altC <= DNA_n_altC_thres)) and ((RNA_n_af <= RNA_n_af_thres) or (RNA_n_altC <= RNA_n_altC_thres)):
                                            writer2.writerow(content)
            elif (pairing == 'unpaired'):
                logger.info('Filtering variant without matching normal!')
                for line in records:
                    for item in headers:
                        if(line[item].isdigit()):
                            print('item is:%s,%s' % (line, item))
                            item = float(line[item])
                        else:
                            item = line[item]
                    """
                    gene = line['gene']
                    chr = line["chromosome"]
                    pos = line["position"]
                    ref = line["ref_base"]
                    alt = line["alt_base"]
                    in_strelka = line["in_strelka"]
                    try:
                        DNA_t_cov = int(line["t_DNA_cov"])
                    except:
                        DNA_t_cov = line["t_DNA_cov"]
                    try:
                        DNA_t_refC = int(line["t_DNA_RefC"])
                    except:
                        DNA_t_refC = line["t_DNA_RefC"]
                    try:
                        DNA_t_altC = int(line["t_DNA_AltC"])
                    except:
                        DNA_t_altC = line["t_DNA_AltC"]
                    try:
                        DNA_t_af = float(line["t_DNA_AF"])
                    except:
                        DNA_t_af = line["t_DNA_AF"]
                    try:
                        RNA_t_cov = int(line["t_RNA_cov"])
                    except:
                        RNA_t_cov = line["t_RNA_cov"]
                    try:
                        RNA_t_refC = int(line["t_RNA_RefC"])
                    except:
                        RNA_t_refC = line["t_RNA_RefC"]
                    try:
                        RNA_t_altC = int(line["t_RNA_AltC"])
                    except:
                        RNA_t_altC = line["t_RNA_AltC"]
                    try:
                        RNA_t_af = float(line["t_RNA_AF"])
                    except:
                        RNA_t_af = line["t_RNA_AF"]
                    """
                    RNA_t_altref_total = RNA_t_refC + RNA_t_altC
                    DNA_t_altref_total = DNA_t_refC + DNA_t_altC
                    content = [line[i] for i in headers]
                    if (in_strelka == 'in_strelka'):
                        # print "in strelka, so keep this variant!"
                        writer.writerow(content)
                    else:
                        if (DNA_t_cov == 'na'):
                           logger.info("Only transcriptome is sequenced for this tumor! ")
                           if (RNA_t_cov >= RNA_t_cov_thres and RNA_t_altC >= RNA_t_altC_thres and RNA_t_af >= RNA_t_af_thres and RNA_t_altref_total >= RNA_t_percent_thres*RNA_t_cov):
                                writer.writerow(content)
                        elif (RNA_t_cov == 'na'):
                            logger.info("Only genome is sequenced for this tumor! ")
                            if (DNA_t_cov >= DNA_t_cov_thres and DNA_t_altC >= DNA_t_altC_thres and DNA_t_af >= DNA_t_af_thres and DNA_t_altref_total >= DNA_t_percent_thres*DNA_t_cov):
                                writer.writerow(content)
                        else:
                            # logger.info("Both genome and transcriptome are sequenced for this tumor! ")
                            if ((DNA_t_cov >= DNA_t_cov_thres and DNA_t_altC >= DNA_t_altC_thres and DNA_t_af >= DNA_t_af_thres and DNA_t_altref_total >= DNA_t_percent_thres*DNA_t_cov) or
                                (RNA_t_cov >= RNA_t_cov_thres and RNA_t_altC >= RNA_t_altC_thres and RNA_t_af >= RNA_t_af_thres and RNA_t_altref_total >= RNA_t_percent_thres*RNA_t_cov)):
                                writer.writerow(content)
                         
    return [filtered_summary, somatic_summary]

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
    # outfile = '.'.join([infile, "final", "txt"])
    outfile = infile.replace('.tmp', '.txt')
    with open (outfile,  'wb') as fh:
        writer = csv.writer( fh, delimiter='\t' )
        with open (infile, 'r') as handle:
             records = csv.DictReader(handle,  delimiter='\t')
             headers = records.fieldnames
             new_headers = [headers[0]] + ['filtered_num_patients_gene_level', 'filtered_num_SNVs_gene_level'] + headers[1:7] + ['filtered_num_patients_SNV_level'] + headers[7:] 
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

