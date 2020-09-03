"""
Description: Script to annotate breakpoints with repeat regions and poor mappability regions
Author: Riccha Sethi
"""

import argparse
import pandas as pd
import numpy as np
import multiprocessing
import pysam
from Tree_Rep import *
import os

def pileup_pos(Line, samfile, STRING):
    ### calculates pileup in a window of 200bp having atleast 10 as base quality
    if STRING=="breakpoint1":
        chrom=Line['chrom1']
        pos=Line['pos1']
    else:
        chrom=Line['chrom2']
        pos=Line['pos2']
    a=[pileupcolumn.n for pileupcolumn in samfile.pileup(chrom, int(pos)-200, int(pos)+200, stepper='all', min_base_quality=10)]
    allreads=sum(a)
    return float(allreads)/400.0

def check_mappability(Line, mappabilityTrack, STRING):
    if STRING=="breakpoint1":
        chromo=Line['chrom1']
        pos=Line['pos1']
    else:
        chromo=Line['chrom2']
        pos=Line['pos2']
    NewDF= mappabilityTrack[(mappabilityTrack['chrom']==chromo) & (mappabilityTrack['pos1']<=pos) & (mappabilityTrack['pos2']>=pos)]
    if NewDF.shape[0]!=0:
        return "Unique"
    else:
        return "Not-unique"

def repeatAnnotation(bt, Line, STRING):
    if STRING=="breakpoint1":
        rep=bt.overlap(Line['chrom1'], Line['pos1'], 2, "+")
    else:
        rep=bt.overlap(Line['chrom2'], Line['pos2'], 2, "+")
    if rep:
        return rep[0][0]
    else:
        return "NA"

def getOverlap(a,b):
    return max(0,min(a[1],b[1])-max(a[0],b[0])+1)

def RepeatClass(chromo, reppos, RepMasker):
    tmpRep=RepMasker[(RepMasker['genoName']==chromo) & (RepMasker['repName']==reppos)]
    if tmpRep.shape[0]==0:
        return 'NA'
    else:
        return tmpRep['repClass'].values.tolist()[0]

def annotate(DF, RepMasker, MapFile, segmentalFile, blackList, BAMT, BAML, bt):
    BAM=pysam.AlignmentFile(BAMT, 'rb')
    BAMLinked=pysam.AlignmentFile(BAML, 'rb')
    DF['LocalCoverage_Pos1_SR'], DF['LocalCoverage_Pos1_LR'], DF['LocalCoverage_Pos2_SR'], DF['LocalCoverage_Pos2_LR'],DF['RepPos1'], DF['RepPos2'],DF['Mappable_Pos1'],DF['Mappable_Pos2'],DF['RepPos1_Class'], DF['RepPos2_Class']= 'NA','NA','NA','NA','NA','NA','NA','NA','NA','NA'
    NewDF= DF.copy(deep=True)
    NewDF['chrom2']=NewDF['chrom2:pos2'].apply(lambda x: x.split(":")[0])
    NewDF['pos2']=NewDF['chrom2:pos2'].apply(lambda x: int(x.split(":")[1]))
    for index,rows in NewDF.iterrows():
        tmp_a=pileup_pos(rows, BAM, "breakpoint1")
        tmp_b=pileup_pos(rows, BAMLinked, "breakpoint1")
        tmp_c=pileup_pos(rows, BAM, "breakpoint2")
        tmp_d=pileup_pos(rows, BAMLinked, "breakpoint2")
        tmp_e=repeatAnnotation(bt, rows, "breakpoint1")
        tmp_f=repeatAnnotation(bt, rows, "breakpoint2")
        tmp_g=check_mappability(rows, MapFile, "breakpoint1")
        tmp_h=check_mappability(rows, MapFile, "breakpoint2")
        NewDF.at[index, 'LocalCoverage_Pos1_SR']=tmp_a
        NewDF.at[index, 'LocalCoverage_Pos1_LR']=tmp_b
        NewDF.at[index, 'LocalCoverage_Pos2_SR']=tmp_c
        NewDF.at[index, 'LocalCoverage_Pos2_LR']=tmp_d
        NewDF.at[index, 'RepPos1']=tmp_e
        NewDF.at[index, 'RepPos2']=tmp_f
        NewDF.at[index, 'Mappable_Pos1']=tmp_g
        NewDF.at[index, 'Mappable_Pos2']=tmp_h
        NewDF.at[index,'RepPos1_Class']= RepeatClass(rows['chrom1'], tmp_e, RepMasker)
        NewDF.at[index,'RepPos2_Class']= RepeatClass(rows['chrom2'], tmp_f, RepMasker)
    BAM.close()
    BAMLinked.close()
    return NewDF

def split(dfm, chunk_size):
# divides dfm pandas dataframe into chunk_size
    rows=dfm.shape[0]
    CHUNK=[]
    for i in range(chunk_size):
        if i != chunk_size-1:
            CHUNK.append(dfm.iloc[int(rows/chunk_size)*i : int(rows/chunk_size)*(i+1), : ])
        else:
            CHUNK.append(dfm.iloc[int(rows/chunk_size)*i :, : ])
    return CHUNK

if __name__=='__main__':
    parser=argparse.ArgumentParser(description="Annotate breakpoints with repetitive regions and poor mappability regions and include local coverage around breakpoints")
    parser.add_argument('-repeatMasker','--repeats', help="enter repeatMasker file (Required)")
    parser.add_argument('-mappability','--mappability', help="enter mappability track file (Required)")
    parser.add_argument('-File','--File', help="enter .tsv file with structural variations from GEM quantification (Required)")
    parser.add_argument('-BAM','--BAMshort',help="enter BAM file with aligned short-reads (Required)")
    parser.add_argument('-BAMLinked','--BAMLinked', help="enter BAM file with aligned linked-reads (Default=None)", default=None)
    parser.add_argument('-Threads','--processes', help="Enter number of cores (default=1)", default=1, type=int)
    parser.add_argument('-outdir','--outdir', help="Enter output directory (Default=current directory)", default=".")
    args=parser.parse_args()
    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)

    RepMasker= pd.read_csv(args.repeats, sep='\t')
    RepMasker=RepMasker.sort_values(by=['genoName', 'genoStart', 'genoEnd'])
    RepFile = RepMasker[['genoName', 'genoStart', 'genoEnd', 'repName', 'swScore', 'strand']]
    MapFile= pd.read_csv(args.mappability, sep='\t', header=None, skiprows=1, names=["chrom","pos1","pos2","kmer","unique","strand"])

    CombinedFile= pd.read_csv(args.File, sep='\t')
    bt=bed_tree(RepFile, 2)

    nthreads= int(args.processes)
    p= multiprocessing.Pool(nthreads)
    chunks=split(CombinedFile, nthreads)
    workers=[p.apply_async(annotate, args=(i, RepMasker, MapFile, segmentalFile, blackList, args.BAMshort, args.BAMLinked, bt,)) for i in chunks]
    final_result=[worker.get() for worker in workers]
    p.close()
    p.join()
    FinalCombinedFile= pd.concat(results for results in final_result)
    FinalCombinedFile.to_csv(args.outdir+ '/Combined_SR_LR_500_annotated.tsv', sep='\t', index=False)
