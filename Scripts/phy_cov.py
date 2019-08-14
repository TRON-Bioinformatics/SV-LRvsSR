### calculates number of molecules at the position that covers the position

import argparse
import pysam
import multiprocessing
import math
import pandas as pd
pd.options.mode.chained_assignment = None

def Mol_depth(L, bamF, chromo):
    bamFile=pysam.AlignmentFile(bamF, 'rb')
    chromosomes, Positions, N_base, depth, Molecules=[],[],[],[],[]
    for i in range(L[0],L[1]+1):
        tmp_mol=[]
        for read in bamFile.fetch(chromo, i, i+1):
            if read.has_tag("MI") and not read.is_duplicate and read.flag<255:
                if read.get_tag("MI") not in tmp_mol:
                    tmp_mol.append(read.get_tag("MI"))
        chromosomes.append(chromo)
        Positions.append(i+1)
        N_base.append("N")
        depth.append(len(tmp_mol))
#        Molecules.append(tmp_mol)
    bamFile.close()
    return pd.concat([pd.Series(chromosomes),pd.Series(Positions),pd.Series(N_base),pd.Series(depth)], axis=1)

def divide_coord(L,n):
    pos1,pos2=L[0],L[1]
    index=(pos2-pos1)/(n)
    tmp=[]
    for i in range(pos1,pos2,index):
        tmp.append([i,i+index])
    tmp.pop()
    tmp[-1][-1]=pos2+1
    return tmp

if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('-bam','--bamfile',help="enter linked-reads generated BAM file")
    parser.add_argument('-chrom','--chrom', help="enter chromosome for molecule pileup")
    parser.add_argument('-out','--output', help="output tab-delimited file with pileup values")
    parser.add_argument('-n','--n',help="enter number of cores to be used")
    args=parser.parse_args()
    size_dict=dict()
    size_dict['chr1']=[9992,248946430]
    size_dict['chr2']=[10007,242183531]
    size_dict['chr3']=[9994,198235457]
    size_dict['chr4']=[9998,190204558]
    size_dict['chr5']=[9992,181478241]
    size_dict['chr6']=[60117,170745945]
    size_dict['chr7']=[9994,159335976]
    size_dict['chr8']=[60653,145076004]
    size_dict['chr9']=[9998,138334475]
    size_dict['chr10']=[9992,133787427]
    size_dict['chr11']=[60351,135076624]
    size_dict['chr12']=[9991,133265313]
    size_dict['chr13']=[18171258,114354334]
    size_dict['chr14']=[18223522,106883717]
    size_dict['chr15']=[19775307,101981192]
    size_dict['chr16']=[9998,90228338]
    size_dict['chr17']=[60185,83247444]
    size_dict['chr18']=[9997,80263286]
    size_dict['chr19']=[60004,58607620]
    size_dict['chr20']=[60000,64334154]
    size_dict['chr21']=[5010085,46699986]
    size_dict['chr22']=[10510020,50808472]
    size_dict['chrX']=[9996,156030894]
    Lists=divide_coord(size_dict[args.chrom], int(args.n))
    p= multiprocessing.Pool(int(args.n))
    workers=[p.apply_async(Mol_depth, args=(i, args.bamfile, args.chrom)) for i in Lists]
    final_result=[worker.get() for worker in workers]
    CombinedMol=pd.concat(results for results in final_result)
    CombinedMol.to_csv(args.output, sep='\t', index=False)
