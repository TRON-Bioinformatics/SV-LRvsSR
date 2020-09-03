"""
Description: Script to count GEMs supporting a particular type and orientation of SV
Author: Riccha Sethi
"""
import pandas as pd
import argparse
import pysam
import numpy as np
import multiprocessing
import math
pd.options.mode.chained_assignment = None

def count_barcodes(bamF, region, chrom1, pos1, chrom2, pos2, SVType, orient):
    barcodes, molecules=set(),set()
    if region:
        for read in bamF.fetch(region[0], max(0,min(region[1],region[2])), max(region[1],region[2])):
            if read.has_tag("BX") and not read.is_duplicate and read.flag<255 and read.is_paired:
                if read.next_reference_name==region[3] and read.next_reference_start in range(min(region[4], region[5]), max(region[4], region[5])):
                    bar_id=read.get_tag("BX")
                    barcodes.add(bar_id)
                    if read.has_tag("MI"):
                        mol_id=read.get_tag("MI")
                        molecules.add(mol_id)
        for read in bamF.fetch(region[3], max(0,min(region[4],region[5])), max(region[4],region[5])):
            if read.has_tag("BX") and not read.is_duplicate and read.flag<255 and read.is_paired:
                if read.next_reference_name==region[0] and read.next_reference_start in range(min(region[1], region[2]), max(region[1], region[2])):
                    bar_id=read.get_tag("BX")
                    barcodes.add(bar_id)
                    if read.has_tag("MI"):
                        mol_id=read.get_tag("MI")
                        molecules.add(mol_id)
    return barcodes,molecules

def  split(dfm, chunk_size):
    rows=dfm.shape[0]
    CHUNK=[]
    for i in range(chunk_size):
        if i != chunk_size-1:
            CHUNK.append(dfm.iloc[int(rows/chunk_size)*i : int(rows/chunk_size)*(i+1), : ])
        else:
            CHUNK.append(dfm.iloc[int(rows/chunk_size)*i :, : ])
    return CHUNK

def find_GEM(df, bamF, window):
    bamfile= pysam.AlignmentFile(bamF,'rb')
    num_commonbar, num_commonmol=[],[]
    sv_select= df[['chrom1','pos1','chrom2:pos2', 'SVType','Orientation']].values.tolist()
    for (chrom1, pos1, chrom2_pos2, SVType, Orientation) in sv_select:
        chromA, posA, chromB, posB, svtype, orient= str(chrom1), int(pos1), str(chrom2_pos2.split(":")[0]), int(chrom2_pos2.split(":")[1]), SVType, Orientation
        if svtype=="Dels" or (svtype=="Trans" and orient=="3to5"):
            region=[chromA, posA-int(window), posA, chromB, posB, posB+int(window)]
        elif svtype=="Dups" or (svtype=="Trans" and orient=="5to3"):
            region=[chromB, posB-int(window), posB, chromA, posA, posA+int(window)]
        elif (svtype=="Invs" and orient=="3to3") or (svtype=="Trans" and orient=="3to3"):
            region=[chromA, posA-int(window), posA, chromB, posB-int(window), posB]
        elif (svtype=="Invs" and orient=="5to5") or (svtype=="Trans" and orient=="5to5"):
            region=[chromB, posB, posB+int(window), chromA, posA, posA+int(window)]
        else:
            region=None
        pos_all_bar, pos_all_mol= count_barcodes(bamfile, region, chromA, posA, chromB, posB, svtype, orient)
        common_bar=set(pos_all_bar)
        common_mol=set(pos_all_mol)
        num_commonbar.append(len(common_bar))
        num_commonmol.append(len(common_mol))
    bamfile.close()
    df['Count_common_bar']=pd.Series(num_commonbar, index=df.index)
    df['Count_common_mol']=pd.Series(num_commonmol, index=df.index)
    return df

if __name__=='__main__':
    parser=argparse.ArgumentParser(description= "A python script to calculate number of GEMs/barcodes(tag:BX) and molecules(tag:MI) support a particular type and orientation of structural variation")
    parser.add_argument('-inputFile','--inputFile', help="enter .tsv file containing combined calls")
    parser.add_argument('-bam','--bamFile', help="enter coordinate sorted bam file from Long Ranger (Lariat) (Required)")
    parser.add_argument('-n','--n', help="enter number of cores to be used (Default=1)", default=1, type=int)
    parser.add_argument('-out','--output', help="enter output .tsv file name (Required)")
    parser.add_argument('-w','--window', help="window size in <INT> bases to calculate supporting GEMs(Default=1000)", default=1000, type=int)
    args=parser.parse_args()
    All_SV= pd.read_csv(args.inputFile, sep='\t', header=(0))
    p= multiprocessing.Pool(int(args.n))
    chunks=split(All_SV, int(args.n))
    workers=[p.apply_async(find_GEM, args=(i,args.bamFile, args.window,)) for i in chunks]
    final_result=[worker.get() for worker in workers]
    CombinedGEM=pd.concat(results for results in final_result)
    mainCol=list(All_SV.columns.values)
    mainCol=mainCol+ ['Count_common_bar','Count_common_mol']
    CombinedGEM.columns=mainCol
    CombinedGEM.to_csv(args.output, sep='\t', index=False)
