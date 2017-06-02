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

logging.basicConfig(
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        datefmt='%m/%d/%Y %I:%M:%S %p',
        level=logging.DEBUG)
logger = logging.getLogger('replace_allele_counts')
try:
    import coloredlogs
    coloredlogs.install(level='DEBUG')
except ImportError:
    logger.warning("replace_allele_counts: coloredlogs is not installed, please install it if you wish to see color in logs on standard out.")
    pass
try:
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
   if(args.verbose):
       logger.info("replace_allele_counts: Started the run for replacing allele counts.")
   (cvDF,mafDF) =read_maf(args)
   (filloutDF) =read_fillout(args)
   (cleanDF) = replace_allele_count(args, mafDF, filloutDF)
   write_output(args,cleanDF)
   if(args.verbose):
       logger.info("replace_allele_counts: Finished the run for replacing allele counts.")

def read_maf(args):
    dataDF = pd.read_table(args.inputMaf, comment="#", low_memory=False)
    complex_variant_dataDF = dataDF.loc[dataDF['TYPE'] == "Complex"]
    return(complex_variant_dataDF,dataDF)

def read_fillout(args):
    dataDF = pd.read_table(args.fillout, comment="#", low_memory=False)
    return(dataDF)

def get_index_for_fillout(filloutDF,m_chr,m_start,m_end,m_ref,m_alt, m_sample_name):
    index = None
    index = filloutDF[(filloutDF['Chromosome'] == m_chr) & 
                           (filloutDF['Tumor_Sample_Barcode'].str.contains(m_sample_name)) & 
                           (filloutDF['Start_Position'] == m_start) & 
                           (filloutDF['End_Position'] == m_end) & 
                           (filloutDF['Reference_Allele'] == m_ref) & 
                           (filloutDF['Tumor_Seq_Allele1'] == m_alt)].index.tolist()[0]
    return(index)
    

def replace_allele_count(args,mafDF,filloutDF):
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
            logging.warn("replace_allele_counts: Skipping as complex event: %s, %d, %d, %s, %s", m_chr, m_start, m_end, m_ref, m_alt)
            continue
        else:
            filloutIndexTumor = get_index_for_fillout(filloutDF, m_chr, m_start, m_end, m_ref, m_alt, m_tumor_sample_name)
            filloutIndexNormal = get_index_for_fillout(filloutDF, m_chr, m_start, m_end, m_ref, m_alt, m_normal_sample_name)
            if(filloutIndexTumor != None and filloutIndexNormal != None):
                #tumor_recordToReplace = filloutDF.iloc[[filloutIndexTumor]]
                #normal_recordToReplace = filloutDF.iloc[[filloutIndexNormal]]
                mafDF_copy.set_value(i_index,"t_depth",filloutDF.get_value(filloutIndexTumor,"t_total_count"))
                mafDF_copy.set_value(i_index,"t_ref_count",filloutDF.get_value(filloutIndexTumor,"t_ref_count"))
                mafDF_copy.set_value(i_index,"t_alt_count",filloutDF.get_value(filloutIndexTumor,"t_alt_count"))
                mafDF_copy.set_value(i_index,"n_depth",filloutDF.get_value(filloutIndexNormal,"t_total_count"))
                mafDF_copy.set_value(i_index,"n_ref_count",filloutDF.get_value(filloutIndexNormal,"t_ref_count"))
                mafDF_copy.set_value(i_index,"n_alt_count",filloutDF.get_value(filloutIndexNormal,"t_alt_count"))
            else:
                logging.warn("replace_allele_counts: Skipping as could not find index in fillout as conversion of vcf to maf failed: %s, %d, %d, %s, %s", m_chr, m_start, m_end, m_ref, m_alt)
                continue
    return(mafDF_copy)

def write_output(args,output_DF):
    if(args.outdir):
        outFile = os.path.join(args.outdir,args.outputMaf)
    else:
        outFile = os.path.join(os.getcwd(),args.outputMaf)
    output_DF.to_csv(outFile, sep='\t', index=False)
    
if __name__ == "__main__":
    start_time = time.time()  
    main()
    end_time = time.time()
    totaltime = end_time - start_time
    logging.info("replace_allele_counts: Elapsed time was %g seconds", totaltime)
    sys.exit(0)
