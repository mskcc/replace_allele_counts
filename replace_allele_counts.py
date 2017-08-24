#!/usr/bin/python
'''
@Description : This tool helps to replace the allele counts from the caller with the allele counts of GetBaseCountMultiSample
@Created :  05/10/2017
@Updated : 06/01/2017
@author : Ronak H Shah

'''
from __future__ import division
import argparse
import sys
import os
import time
import logging
import multiprocessing

#try:
#    import coloredlogs
#    coloredlogs.install(level='DEBUG')
#except ImportError:
#    logger.warning("replace_allele_counts: coloredlogs is not installed, please install it if you wish to see color in logs on standard out.")
#    pass
try:
    import numpy as np
    import pandas as pd
except ImportError:
    logger.fatal("replace_allele_counts: pandas is not installed, please install pandas as it is required to run the removing process.")
    sys.exit(1)

def main():
   parser = argparse.ArgumentParser(prog='replace_allele_counts.py', description='This tool helps to replace the allele counts from the caller with the allele counts of GetBaseCountMultiSample', usage='%(prog)s [options]')
   parser.add_argument("-v", "--verbose", action="store_true", dest="verbose", help="make lots of noise")
   parser.add_argument("-imaf", "--input-maf", action="store", dest="inputMaf", required=True, type=str, metavar='SomeID.maf', help="Input maf file which needs to be fixed")
   parser.add_argument("-ifill","--fillout", action="store", dest="fillout", required=True, type=str, metavar='SomeID.fillout.maf', help="Input fillout file created by GetBaseCountMultiSample using the input maf")
   parser.add_argument("-omaf","--output-maf", action="store", dest="outputMaf", required=True, type=str, metavar='SomeID.maf', help="Output maf file name")
   parser.add_argument("-o", "--outDir", action="store", dest="outdir", required=False, type=str, metavar='/somepath/output', help="Full Path to the output dir.")
   
   args = parser.parse_args()
   logger.info("Reading Maf...")
   (cvDF,mafDF) =read_maf(args)
   logger.info("Reading Fillout...")
   maf_samples = np.concatenate((mafDF['Tumor_Sample_Barcode'].unique(), mafDF['Matched_Norm_Sample_Barcode'].unique()))
   (filloutDF) =read_fillout(args, maf_samples)
   logger.info("Rewriting fillout lines into maf...")
   group_size = len(mafDF.index) / 10 
   dataframes = []
   pool = multiprocessing.Pool(processes=5)
   m = multiprocessing.Manager()
   q = m.Queue()
   for g, df in mafDF.groupby(np.arange(len(mafDF)) // group_size):
       (cleanDF) = pool.apply_async(replace_allele_count, args=(args, df, filloutDF, q, ))
#       (cleanDF) = replace_allele_count(args, mafDF, filloutDF, q)
   pool.close()
   pool.join()
   while not q.empty():
       dataframes.append(q.get())
   write_output(args,dataframes)
   if(args.verbose):
       logger.info("replace_allele_counts: Finished the run for replacing allele counts.")

def read_maf(args):
    dataDF = pd.read_table(args.inputMaf, comment="#", low_memory=False)
    complex_variant_dataDF = dataDF.loc[dataDF['TYPE'] == "Complex"]
    return(complex_variant_dataDF,dataDF)

def read_fillout(args, maf_samples):
    dataDF = pd.read_table(args.fillout, comment="#", low_memory=False)
    logger.info("sorting by chrom")
    data_by_chrom_and_sample = {}
    for chrom in dataDF['Chromosome'].unique():
        if chrom not in data_by_chrom_and_sample:
            data_by_chrom_and_sample[chrom]={}
        for sample in dataDF['Tumor_Sample_Barcode'].unique():
            #lets do this rough match once instead of every line
            adjusted_sample_name = None
            for mafsample in maf_samples:
                if sample.find(mafsample) > -1:
                    adjusted_sample_name = mafsample
            if not adjusted_sample_name:
                continue
                #fillout can have more samples than maf, it turns out 
                logger.critical("Unable to match fillout sample names with actual maf.")
                logger.critical("Maf samples: %s" % " ".join(maf_samples))
                logger.critical("Unmatched fillout sample id: %s" % sample)
                sys.exit(1)
            #data within the df is untouched, but it is organized in a double level dict by chr, sample name
            #second level sample name matches maf sample name for O(1) into sample specific fillout list 
            data_by_chrom_and_sample[chrom][adjusted_sample_name]=dataDF[(dataDF['Chromosome']== chrom) &
                (dataDF['Tumor_Sample_Barcode']== sample)]

    return (data_by_chrom_and_sample)

def get_index_for_fillout(data_chrom_sample,m_chr,m_start,m_end,m_ref,m_alt, m_sample_name):
    index = None
    filloutDF = data_chrom_sample[m_chr][m_sample_name]
    index = filloutDF[(filloutDF['Start_Position'] == m_start) & 
                           (filloutDF['End_Position'] == m_end) & 
                           (filloutDF['Reference_Allele'] == m_ref) & 
                           (filloutDF['Tumor_Seq_Allele1'] == m_alt)].index.tolist()[0]
    return(index)
    

def replace_allele_count(args,mafDF,filloutDF,q):
    mafDF_copy = mafDF.copy()
    for i_index, i_row in mafDF.iterrows():
        m_chr = i_row.loc['Chromosome']
        m_start = i_row.loc['Start_Position']
        m_end = i_row.loc['End_Position']
        m_type = i_row.loc['TYPE']
        m_vt = i_row.loc['Variant_Type']
        m_ref = (str(i_row.loc['Reference_Allele'])).rstrip()
        m_alt = (str(i_row.loc['Tumor_Seq_Allele2'])).rstrip()
        m_tumor_sample_name = i_row.loc['Tumor_Sample_Barcode']
        m_normal_sample_name = i_row.loc['Matched_Norm_Sample_Barcode']
        
        #if(len(m_ref) >= 2 and len(m_alt) >= 2 and (len(m_ref) != len(m_alt))):
        if(m_type == "Complex"):
            logger.warn("replace_allele_counts: Skipping as complex event: %s, %d, %d, %s, %s", m_chr, m_start, m_end, m_ref, m_alt)
            continue
        else:
            filloutIndexTumor = get_index_for_fillout(filloutDF, m_chr, m_start, m_end, m_ref, m_alt, m_tumor_sample_name)
            filloutIndexNormal = get_index_for_fillout(filloutDF, m_chr, m_start, m_end, m_ref, m_alt, m_normal_sample_name)
            if(filloutIndexTumor != None and filloutIndexNormal != None):
                #tumor_recordToReplace = filloutDF.iloc[[filloutIndexTumor]]
                #normal_recordToReplace = filloutDF.iloc[[filloutIndexNormal]]
                mafDF_copy.set_value(i_index,"t_depth",filloutDF[m_chr][m_tumor_sample_name].get_value(filloutIndexTumor,"t_total_count"))
                mafDF_copy.set_value(i_index,"t_ref_count",filloutDF[m_chr][m_tumor_sample_name].get_value(filloutIndexTumor,"t_ref_count"))
                mafDF_copy.set_value(i_index,"t_alt_count",filloutDF[m_chr][m_tumor_sample_name].get_value(filloutIndexTumor,"t_alt_count"))
                mafDF_copy.set_value(i_index,"n_depth",filloutDF[m_chr][m_normal_sample_name].get_value(filloutIndexNormal,"t_total_count"))
                mafDF_copy.set_value(i_index,"n_ref_count",filloutDF[m_chr][m_normal_sample_name].get_value(filloutIndexNormal,"t_ref_count"))
                mafDF_copy.set_value(i_index,"n_alt_count",filloutDF[m_chr][m_normal_sample_name].get_value(filloutIndexNormal,"t_alt_count"))
            else:
                logger.warn("replace_allele_counts: Skipping as could not find index in fillout as conversion of vcf to maf failed: %s, %d, %d, %s, %s", m_chr, m_start, m_end, m_ref, m_alt)
                continue
    q.put(mafDF_copy)
    return True 

def write_output(args,frames):
    if(args.outdir):
        outFile = os.path.join(args.outdir,args.outputMaf)
    else:
        outFile = os.path.join(os.getcwd(),args.outputMaf)
    pd.concat(frames).to_csv(outFile, sep='\t', index=False)
    
if __name__ == "__main__":
    logger = logging.getLogger('replace_allele_counts')
    logger.propagate=False
    logger.setLevel(logging.DEBUG)
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    ch=logging.StreamHandler()
    ch.setLevel(logging.DEBUG)
    ch.setFormatter(formatter)
    logger.addHandler(ch)

    start_time = time.time()  
    main()
    end_time = time.time()
    totaltime = end_time - start_time
    logger.info("replace_allele_counts: Elapsed time was %g seconds", totaltime)
    sys.exit(0)
