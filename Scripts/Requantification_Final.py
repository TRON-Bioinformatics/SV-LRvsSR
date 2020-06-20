"""
Description: Script to perform requantification for SVs
Author: Riccha Sethi
"""

### Load packages

import argparse
import pandas as pd
import os
import itertools
import numpy as np
import time
import math
import multiprocessing
from template_final import *
pd.options.mode.chained_assignment = None
import subprocess
import glob
from Tree_Rep import *

def NumberTools(WholeLine):
    tools=0
    ListNull= WholeLine.isnull()
    if ListNull['DellyChrom1']==False:
        tools+=1
    if ListNull['LumpyChrom1']==False:
        tools+=1
    if ListNull['SvabaChrom1']==False:
        tools+=1
    return tools

def NumberToolsLinked(WholeLine):
    tools=0
    ListNull= WholeLine.isnull()
    if ListNull['Linked_chrom1']==False:
        tools+=1
    if ListNull['NAIBR_chrom1']==False:
        tools+=1
    if ListNull['GROCSV_chrom1']==False:
        tools+=1
    return tools

def Add_features(CombinedFile, refBit, areaRequant, tmpdir):
    ID=id_generator()
    totalTools, totalTools_LR=[], []
    FusedSeqT_fa=open(str(tmpdir)+"/TumorSV_"+str(ID)+".fa", 'w')
    requantificationBP=open(str(tmpdir)+"/requantBP_"+str(ID)+".txt",'w')
    for i in range(CombinedFile.shape[0]):
        totalTools.append(NumberTools(CombinedFile.iloc[i]))
        totalTools_LR.append(NumberToolsLinked(CombinedFile.iloc[i]))
        Requant, areaRequantT=find_requant(CombinedFile.iloc[i], refBit, areaRequant)
        requantificationBP.write(str(CombinedFile.iloc[i]['chrom1'])+":"+ str(CombinedFile.iloc[i]['pos1']) + "_"+ str(CombinedFile.iloc[i]['chrom2:pos2']) + "_"+ str(CombinedFile.iloc[i]['SVType'])+"_"+ str(CombinedFile.iloc[i]['Orientation'])+ '\t'+str(areaRequantT)+ '\n')
        FusedSeqT_fa.write(">"+str(CombinedFile.iloc[i]['chrom1'])+":"+ str(CombinedFile.iloc[i]['pos1']) + "_"+ str(CombinedFile.iloc[i]['chrom2:pos2']) + "_"+ str(CombinedFile.iloc[i]['SVType'])+"_"+ str(CombinedFile.iloc[i]['Orientation'])+'\n' + Requant+'\n')
    FusedSeqT_fa.close()
    requantificationBP.close()
    CombinedFile['NumberTools_SR']=pd.Series(totalTools, index=CombinedFile.index)
    CombinedFile['NumberTools_LR']=pd.Series(totalTools_LR, index=CombinedFile.index)
    return CombinedFile

def calculateRequant(outdir, CombinedFile, SR_R1, SR_R2, tmpdir, junc_cutoff, nthreads, LR_R1, LR_R2):
    ### reading bp position for requant fused synthetic template
    requantificationBP={}
    # combine all synthtic genomic templates for all SVs
    os.chdir(tmpdir)
    requantBPFiles=glob.glob('./requantBP*.txt')
    commandCombRequant=['cat']+requantBPFiles
    with open(str(outdir)+"/BP_All.txt",'w') as F:
        combineReq=subprocess.Popen(commandCombRequant, stdout=F)

    files=os.listdir(tmpdir)
    FaFiles=['cat']+[i for i in files if i.endswith(".fa")]
    with open(str(outdir)+"/FinalFasta.fa",'w') as COMB:
        combineFasta=subprocess.Popen(FaFiles, stdout=COMB)
    combineReq.wait()
    combineFasta.wait()
    F.close()
    COMB.close()
    with open(str(outdir)+"/BP_All.txt",'r') as f:
        for line in f:
            line_list=line.split('\t')
            requantificationBP[line_list[0]]=int(line_list[1])

    print "aligning cWGS reads to genomic templates (requantification)"
    align(str(outdir)+"/FinalFasta.fa",SR_R1, SR_R2, outdir, tmpdir, nthreads, "SR")
    COUNTS_SR=calculate_Junc_Span(str(outdir)+"/FinalFasta.fa", str(tmpdir)+"/filtered_SR.sorted.bam", junc_cutoff, outdir, requantificationBP, "SR")
    align(str(outdir)+"/FinalFasta.fa",LR_R1, LR_R2, outdir, tmpdir, nthreads, "LR")
    COUNTS_LR=calculate_Junc_Span(str(outdir)+"/FinalFasta.fa", str(tmpdir)+"/filtered_LR.sorted.bam", junc_cutoff, outdir, requantificationBP, "LR")
    print "Calculating junction and spanning reads"
    JUNC_READS_SR, SPAN_READS_SR=[],[]
    JUNC_READS_LR, SPAN_READS_LR=[],[]
    Ref_Ids=COUNTS_SR.keys()
    for i in range(CombinedFile.shape[0]):
        refid=str(CombinedFile.iloc[i]['chrom1'])+":"+ str(CombinedFile.iloc[i]['pos1']) + "_"+ str(CombinedFile.iloc[i]['chrom2:pos2']) + "_"+ str(CombinedFile.iloc[i]['SVType'])+"_"+ str(CombinedFile.iloc[i]['Orientation'])
        if refid in Ref_Ids:
            JUNC_READS_SR.append(COUNTS_SR[refid][0])
            SPAN_READS_SR.append(COUNTS_SR[refid][1])
            JUNC_READS_LR.append(COUNTS_LR[refid][0])
            SPAN_READS_LR.append(COUNTS_LR[refid][1])
        else:
            print "missing SV for requantification",refid
    CombinedFile['JR_SR']=pd.Series(JUNC_READS_SR, index=CombinedFile.index)
    CombinedFile['SP_SR']=pd.Series(SPAN_READS_SR, index=CombinedFile.index)
    CombinedFile['JR_LR']=pd.Series(JUNC_READS_LR, index=CombinedFile.index)
    CombinedFile['SP_LR']=pd.Series(SPAN_READS_LR, index=CombinedFile.index)
    return CombinedFile

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

def main():
    parser=argparse.ArgumentParser(description= "A python script to perform requantification for SVs using read-pairs from cWGS and 10XWGS technology")
    parser.add_argument('-inputFile','--inputFile', help="enter file containing combined calls")
    parser.add_argument('-out','--output', help="enter final csv file to which requantification counts would be added")
    parser.add_argument('-n','--processes', help="enter number of cores to run the process", default=1, type=int)
    parser.add_argument('-refBit','--twoBit', help="enter two bit file of reference genome")
    parser.add_argument('-area','--areaRequant',help="enter area around breakpoints for synthetic genomic template", default=500, type=int)
    parser.add_argument('-Read1_SR','--Read1_SR',help="enter read1.fq from cWGS")
    parser.add_argument('-Read2_SR','--Read2_SR',help="enter read2.fq from cWGS")
    parser.add_argument('-Read1_LR','--Read1_LR',help="enter read1.fq from 10XWGS")
    parser.add_argument('-Read2_LR','--Read2_LR',help="enter read2.fq from 10XWGS")
    parser.add_argument('-outdir','--outdir',help="enter output directory", default=".")
    parser.add_argument('-tmpdir','--tmpdir',help="enter path of temporary directory", default="")
    parser.add_argument('-cutoff','--junc_cutoff', type=int, help="cutoff for junction reads. a read has to overlap breakpoint with atleast <INT> bases", default=10)
    parser.add_argument('-lengths','--lenChrom',help="enter file contaning length of chromosomes (Required)")
    return parser.parse_args()

if __name__=='__main__':
    args=main()
    CombinedFile1=pd.read_csv(args.inputFile)
    CombinedFile=CombinedFile1.sample(frac=1).reset_index(drop=True)
    lengths=dict()
    with open(args.lenChrom) as f:
        for line in f:
            line_list=line.split('\t')
            lengths[line_list[0]]=int(line_list[-1])
    nthreads=int(args.processes)
    p= multiprocessing.Pool(nthreads)
    chunks=split(CombinedFile, nthreads)
    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)
    if not os.path.exists(args.tmpdir):
        os.makedirs(args.tmpdir)
    workers= [p.apply_async(Add_features, args=(i, args.twoBit, args.areaRequant, args.tmpdir,)) for i in chunks]
    final_result= [worker.get() for worker in workers]
    FinalCombinedFile_1=pd.concat(results for results in final_result)
    FinalCombinedFile_2=calculateRequant(args.outdir, FinalCombinedFile_1, args.Read1_SR, args.Read2_SR, args.tmpdir, args.junc_cutoff, nthreads, args.Read1_LR, args.Read2_LR)
    FinalCombinedFile_2.to_csv(str(args.outdir)+"/"+str(args.output), sep="\t", index=False)
