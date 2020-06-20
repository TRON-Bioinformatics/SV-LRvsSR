
"""
Description: Script to combine SV calls from short-reads and linked-reads
Author: Riccha Sethi
"""

### Load packages
import argparse, __main__ as main
import pandas as pd
import os, sys
import numpy as np
import itertools
pd.options.mode.chained_assignment = None
CurVersion=1.0

def get_orient(orient):
    if orient=="+-":
        return "3to5"
    elif orient=="-+":
        return "5to3"
    elif orient=="++":
        return "3to3"
    elif orient=="--":
        return "5to5"

def get_SVType(orient, chrom1, chrom2):
    if chrom1==chrom2:
        if orient=="+-":
            return "Dels"
        elif orient=="-+":
            return "Dups"
        elif orient=="++":
            return "Invs"
        elif orient=="--":
            return "Invs"
    else:
        return "Trans"

def oppOrient(orient):
    if orient=="3to5":
        return "5to3"
    elif orient=="5to3":
        return "3to5"
    elif orient=="3to3":
        return "3to3"
    elif orient=="5to5":
        return "5to5"

def check_common(Line):
    if Line=="Only cWGS":
        return "Common"
    elif Line=="Only 10XWGS":
        return "Only 10XWGS"
    elif Line=="Common":
        return "Common"
    else:
        return "Only 10XWGS"

def getSize(x, tool):
    if tool=="NAIBR":
        if x['NAIBR_SVType'] in ['Dels', 'Dups','Invs']:
            return x['NAIBR_Pos2']-x['NAIBR_Pos1']+1
        else:
            return 0
    else:
        if x['GROCSV_SVType'] in ['Dels','Dups','Invs']:
            return x['GROCSV_pos2']-x['GROCSV_pos1']+1
        else:
            return 0

def check_category(mate1,mate2):
    chrB=str(mate1.split(":")[0].split("]")[-1].split("[")[-1])
    pos2=int(mate1.split(":")[1].split("[")[0].split("]")[0])
    chrA=str(mate2.split(":")[0].split("]")[-1].split("[")[-1])
    pos1=int(mate2.split(":")[1].split("[")[0].split("]")[0])
    if chrA==chrB and "[" in mate1 and "]" in mate2:
        if pos2>pos1:
            return "Dels", chrA, pos1, chrA, pos2, "3to5"
        else:
            return "Dels", chrA, pos2, chrA, pos1, "3to5"
    elif chrA==chrB and "]" in mate1 and "[" in mate2:
        if pos2>pos1:
            return "Dups",chrA, pos1, chrA, pos2, "5to3"
        else:
            return "Dups",chrA,pos2,chrA,pos1,"5to3"
    elif chrA==chrB and ("]" in mate1 and "]" in mate2):
        if pos2>pos1:
            return "Invs",chrA, pos1, chrA, pos2,"3to3"
        else:
            return "Invs",chrA, pos2, chrA, pos1,"3to3"
    elif chrA==chrB and ("[" in mate1 and "[" in mate2):
        if pos2>pos1:
            return "Invs", chrA, pos1, chrA, pos2,"5to5"
        else:
            return "Invs", chrA, pos2, chrA, pos1,"5to5"
    else:
        if "[" in mate1 and "]" in mate2:
            return "Trans",chrA, pos1, chrB, pos2,"3to5"
        elif "]" in mate1 and "[" in mate2:
            return "Trans",chrA, pos1, chrB, pos2,"5to3"
        elif "]" in mate1 and "]" in mate2:
            return "Trans",chrA, pos1, chrB, pos2,"3to3"
        elif "[" in mate1 and "[" in mate2:
            return "Trans",chrA, pos1, chrB, pos2,"5to5"

def read_GROCSV(File):
    grocsv={'GROCSV_chrom1':[], 'GROCSV_pos1':[], 'GROCSV_chrom2':[], 'GROCSV_pos2':[], 'GROCSV_SVType':[], 'GROCSV_Orient':[], 'GROCSV_Filter':[], 'GROCSV_SupportBarcode':[]}
    with open(File, 'r') as f:
        for line1,line2 in itertools.izip_longest(f, f, fillvalue=''):
            line1_list=line1.rstrip('\n').split('\t')
            line2_list=line2.rstrip('\n').split('\t')
            info1, info2= line1_list[7].split(";"), line2_list[7].split(";")
            if int(line1_list[2].split(":")[1])< int(line2_list[2].split(":")[1]):
                mate1,mate2=line1_list[4],line2_list[4]
            else:
                mate1,mate2=line2_list[4],line1_list[4]
            FILTER=line1_list[6]
            bs=int(line1_list[9].split(":")[1])
            cat=check_category(mate1, mate2)
            if cat:
                grocsv['GROCSV_chrom1'].append(cat[1])
                grocsv['GROCSV_pos1'].append(int(cat[2]))
                grocsv['GROCSV_chrom2'].append(cat[3])
                grocsv['GROCSV_pos2'].append(int(cat[4]))
                grocsv['GROCSV_SVType'].append(cat[0])
                grocsv['GROCSV_Orient'].append(cat[5])
                grocsv['GROCSV_Filter'].append(FILTER)
                grocsv['GROCSV_SupportBarcode'].append(bs)
            else:
                print "missing category\n", line1
    return pd.DataFrame(grocsv)

def add_calls(calls, CombFile, window, svtype, tool):
    all1_chrom1, all1_chrom2, all1_pos1, all1_pos2, all1_orient, all1_SV=CombFile.chrom1.values, CombFile.chrom2.values, CombFile.pos1.values, CombFile.pos2.values, CombFile.Orientation.values, CombFile.SVType.values
    if tool=="NAIBR":
        all2_chrom1, all2_chrom2, all2_pos1, all2_pos2, all2_orient, all2_SV=calls.NAIBR_chrom1.values, calls.NAIBR_chrom2.values, calls.NAIBR_Pos1.values, calls.NAIBR_Pos2.values, calls.NAIBR_Orient.values, calls.NAIBR_SVType.values
    else:
        all2_chrom1, all2_chrom2, all2_pos1, all2_pos2, all2_orient, all2_SV=calls.GROCSV_chrom1.values, calls.GROCSV_chrom2.values, calls.GROCSV_pos1.values, calls.GROCSV_pos2.values, calls.GROCSV_Orient.values, calls.GROCSV_SVType.values
    all2_pos1_l=all2_pos1-window
    all2_pos1_h=all2_pos1+window
    all2_pos2_l=all2_pos2-window
    all2_pos2_h=all2_pos2+window
    if svtype=="Trans":
        i, j =np.where(((all1_chrom1[:, None]==all2_chrom1) & (all1_chrom2[:,None]==all2_chrom2) & (all1_pos1[:,None] >= all2_pos1_l) & (all1_pos1[:,None] <= all2_pos1_h) & (all1_pos2[:,None] >= all2_pos2_l) & (all1_pos2[:,None] <= all2_pos2_h) & (all1_orient[:,None]==all2_orient) & (all1_SV[:,None]==all2_SV)) | ((all1_chrom1[:,None]==all2_chrom2) & (all1_chrom2[:,None]==all2_chrom1) & (all1_pos1[:,None] >= all2_pos2_l) & (all1_pos1[:,None] <= all2_pos2_h) & (all1_pos2[:,None] >= all2_pos1_l) & (all1_pos2[:,None] <= all2_pos1_h) & (all1_orient[:,None]==list(map(oppOrient, all2_orient))) & (all1_SV[:,None]==all2_SV)))
    else:
        i, j =np.where((all1_chrom1[:, None]==all2_chrom1) & (all1_chrom2[:,None]==all2_chrom2) & (all1_pos1[:,None] >= all2_pos1_l) & (all1_pos1[:,None] <= all2_pos1_h) & (all1_pos2[:,None] >= all2_pos2_l) & (all1_pos2[:,None] <= all2_pos2_h) & (all1_orient[:,None]==all2_orient) & (all1_SV[:,None]==all2_SV))
    tmpDF=calls.copy(deep=True)
    Commons=pd.DataFrame(np.column_stack([CombFile.values[i], tmpDF.values[j]]), columns= CombFile.columns.append(tmpDF.columns))
    Commons['Predicted_by']=Commons['Predicted_by'].apply(lambda x: check_common(x))
    OnlyCalls= calls[~np.in1d(np.arange(calls.shape[0]), np.unique(j))]
    if tool=="NAIBR":
        OnlyCalls['chrom1']=OnlyCalls['NAIBR_chrom1']
        OnlyCalls['pos1']=OnlyCalls['NAIBR_Pos1']
        OnlyCalls['chrom2:pos2']=OnlyCalls['NAIBR_chrom2']+":"+OnlyCalls['NAIBR_Pos2'].astype(str)
        OnlyCalls['SVType']=OnlyCalls['NAIBR_SVType']
        OnlyCalls['Orientation']=OnlyCalls['NAIBR_Orient']
    else:
        OnlyCalls['chrom1']=OnlyCalls['GROCSV_chrom1']
        OnlyCalls['pos1']=OnlyCalls['GROCSV_pos1']
        OnlyCalls['chrom2:pos2']=OnlyCalls['GROCSV_chrom2']+":"+OnlyCalls['GROCSV_pos2'].astype(str)
        OnlyCalls['SVType']=OnlyCalls['GROCSV_SVType']
        OnlyCalls['Orientation']=OnlyCalls['GROCSV_Orient']
    OnlyCalls['Size']=OnlyCalls.apply(lambda x: getSize(x, tool), axis=1)
    OnlyCalls['Predicted_by']="Only 10XWGS"
    AllDF=pd.concat([Commons, CombFile[~np.in1d(np.arange(CombFile.shape[0]), np.unique(i))], OnlyCalls], ignore_index=True, sort=False)
    return AllDF

if __name__=='__main__':
    parser= argparse.ArgumentParser(description = "A python script to add NAIBR and GROC-SV calls")
    parser.add_argument('-NAIBR','--NAIBR', help="enter .bed file containing SV calls from NAIBR (REQUIRED)")
    parser.add_argument('-w', '--window', help="window to match SV breakpoints from short-reads and longranger\n (DEFAULT=100bp)", default=int(100))
    parser.add_argument('-GROCSV','--GROCSV', help="enter .vcf file containing SV calls from GROC-SV (REQUIRED)")
    parser.add_argument('-File','--CombinedFile', help="enter .csv file containing calls from short-reads and longranger (REQUIRED)")
    args= parser.parse_args()

    ### reading combined file
    CombFile=pd.read_csv(args.CombinedFile, sep=',')
    CombFile['chrom2']=CombFile['chrom2:pos2'].apply(lambda x: x.split(":")[0])
    CombFile['pos2']=CombFile['chrom2:pos2'].apply(lambda x: int(x.split(":")[1]))

    ### adding NAIBR calls
    NAIBR_File= pd.read_csv(args.NAIBR, sep='\t')
    NAIBR_File['NAIBR_Orient']= NAIBR_File['Orientation'].apply(lambda x: get_orient(x))
    NAIBR_File= NAIBR_File[(NAIBR_File['Chr1'].isin(["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y"]) & (NAIBR_File['Chr2'].isin(["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y"])))]
    NAIBR_File['Chr1']=NAIBR_File['Chr1'].astype(str).apply(lambda x: "chr"+x)
    NAIBR_File['Chr2']=NAIBR_File['Chr2'].astype(str).apply(lambda x: "chr"+x)
    NAIBR_File['NAIBR_SVType']= NAIBR_File.apply(lambda x: get_SVType(x['Orientation'], x['Chr1'], x['Chr2']), axis=1)
    NAIBR_File= NAIBR_File.sort_values(by=['Chr1','Break1','Chr2', 'Break2'])
    NAIBR_File= NAIBR_File.rename(columns={"Chr1":"NAIBR_chrom1", "Break1":"NAIBR_Pos1", "Chr2": "NAIBR_chrom2", "Break2":"NAIBR_Pos2", "Split molecules":"NAIBR_Split_molecules", "Discordant reads":"NAIBR_Discordant_reads", "Orientation":"NAIBR_Orientation", "Haplotype":"NAIBR_Haplotype", "Score":"NAIBR_Score", "Pass filter":"NAIBR_Filter"})
    print "dels"
    AllDels= add_calls(NAIBR_File[(NAIBR_File['NAIBR_SVType']=="Dels")], CombFile[CombFile['SVType']=="Dels"], int(args.window), "Dels", "NAIBR")
    print "dups"
    AllDups= add_calls(NAIBR_File[(NAIBR_File['NAIBR_SVType']=="Dups")], CombFile[CombFile['SVType']=="Dups"], int(args.window), "Dups", "NAIBR")
    print "invs"
    AllInvs= add_calls(NAIBR_File[(NAIBR_File['NAIBR_SVType']=="Invs")], CombFile[CombFile['SVType']=="Invs"], int(args.window), "Invs", "NAIBR")
    print "trans"
    AllTrans= add_calls(NAIBR_File[(NAIBR_File['NAIBR_SVType']=="Trans")], CombFile[CombFile['SVType']=="Trans"], int(args.window), "Trans", "NAIBR")
    print "done"
    FinalFile1= pd.concat([AllDels, AllDups, AllInvs, AllTrans], ignore_index=True, sort=False)
    FinalFile1= FinalFile1.sort_values(by=['chrom1','pos1','chrom2:pos2'])

    ### adding GROC-SV calls
    grocCalls= read_GROCSV(args.GROCSV)
    print "dels"
    AllDels= add_calls(grocCalls[(grocCalls['GROCSV_SVType']=="Dels")], FinalFile1[FinalFile1['SVType']=="Dels"], int(args.window), "Dels", "GROCSV")
    print "dups"
    AllDups= add_calls(grocCalls[(grocCalls['GROCSV_SVType']=="Dups")], FinalFile1[FinalFile1['SVType']=="Dups"], int(args.window), "Dups", "GROCSV")
    print "invs"
    AllInvs= add_calls(grocCalls[(grocCalls['GROCSV_SVType']=="Invs")], FinalFile1[FinalFile1['SVType']=="Invs"], int(args.window), "Invs", "GROCSV")
    print "trans"
    AllTrans= add_calls(grocCalls[(grocCalls['GROCSV_SVType']=="Trans")], FinalFile1[FinalFile1['SVType']=="Trans"], int(args.window), "Trans", "GROCSV")
    FinalFile= pd.concat([AllDels, AllDups, AllInvs, AllTrans], ignore_index=True, sort=False)
    FinalFile= FinalFile.sort_values(by=['chrom1','pos1','chrom2:pos2'])

    FinalFile.to_csv('Combined_SR_LR_500_NAIBR_GROCSV.csv', sep=',', index=False)
