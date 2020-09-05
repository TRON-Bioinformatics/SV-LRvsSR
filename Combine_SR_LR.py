#!/usr/bin/env python

"""
Description: Script to combine SV calls from short-reads from Integrate.py and Long Ranger, NAIBR and GROC-SV from linked-reads
Author: Riccha Sethi
"""

### Load packages
import argparse, __main__ as main
import pandas as pd
import gzip
import os, sys
import numpy as np
import itertools
pd.options.mode.chained_assignment = None
CurVersion=1.0

### for linked-reads vcf files
def find_orient(mate1, mate2):
    chrom2=str(mate1[4].split(":")[0].split("]")[-1].split("[")[-1])
    pos2=int(mate1[4].split(":")[1].split("[")[0].split("]")[0])
    chrom1=str(mate2[4].split(":")[0].split("]")[-1].split("[")[-1])
    pos1=int(mate2[4].split(":")[1].split("[")[0].split("]")[0])
    if chrom1==chrom2 and "[" in mate1[4] and "]" in mate2[4]:
        return "Dels", chrom1, pos1, chrom1, pos2, "3to5", pos2-pos1+1
    elif chrom1==chrom2 and "]" in mate1[4] and "[" in mate2[4]:
        return "Dups",chrom1, pos1, chrom1, pos2, "5to3", pos2-pos1+1
    elif chrom1==chrom2 and ("]" in mate1[4] and "]" in mate2[4]):
        return "Invs",chrom1, pos1, chrom1, pos2,"3to3", pos2-pos1+1
    elif chrom1==chrom2 and ("[" in mate1[4] and "[" in mate2[4]):
        return "Invs", chrom1, pos1, chrom1, pos2,"5to5", pos2-pos1+1
    else:
        if "[" in mate1[4] and "]" in mate2[4]:
            return "Trans",chrom1, pos1, chrom2, pos2,"3to5", 0
        elif "]" in mate1[4] and "[" in mate2[4]:
            return "Trans",chrom1, pos1, chrom2, pos2,"5to3", 0
        elif "]" in mate1[4] and "]" in mate2[4]:
            return "Trans",chrom1, pos1, chrom2, pos2,"3to3",0
        elif "[" in mate1[4] and "[" in mate2[4]:
            return "Trans",chrom1, pos1, chrom2, pos2,"5to5",0

def add_bnd(bndList, linked, done):
    bndListN=sorted(bndList.keys())
    for mate1 in bndListN:
        if mate1 not in done:
            mate1Value=bndList[mate1]
            mate2=[i for i in mate1Value[7].split(";") if i.startswith("MATEID=")][0].split("=")[1]
            mate2Value=bndList[mate2]
            done.extend([x for x in [mate1,mate2]])
            svtype,chrom1,pos1,chrom2,pos2,orient,size= find_orient(mate1Value, mate2Value)
            if svtype not in linked.keys():
                linked[svtype]={}
            if orient not in linked[svtype].keys():
                linked[svtype][orient]={}
            if chrom1 not in linked[svtype][orient].keys():
                linked[svtype][orient][chrom1]={}
            linked[svtype][orient][chrom1][pos1]=[chrom2+":"+str(pos2), mate1Value[6], [i for i in mate1Value[7].split(";") if i.startswith("SOURCE=")][0].split("=")[1], size, int(mate1Value[5]), [i for i in mate1Value[7].split(";") if i.startswith("HAP_ALLELIC_FRAC")][0].split("=")[1], [i for i in mate1Value[7].split(";") if i.startswith("ALLELIC_FRAC")][0].split("=")[1],[i for i in mate1Value[7].split(";") if i.startswith("PAIRS=")][0].split("=")[1], [i for i in mate1Value[7].split(";") if i.startswith("SPLIT=")][0].split("=")[1] ]
    return linked,done

def read_linked(file1, linked):
    bndList,done={},[]
    with gzip.open(file1,'r') as f:
        for line in f:
            if line.startswith("#"):
                continue
            else:
                line_list=line.rstrip().split('\t')
                chrom1, pos1, info, infoAll= line_list[0], int(line_list[1]), line_list[6], line_list[7].split(";")
                svtype=[i for i in infoAll if i.startswith('SVTYPE=')][0].split("=")[1]
                source=[i for i in infoAll if i.startswith('SOURCE=')][0].split("=")[1]
                if svtype !='BND':
                    pos2=int([i for i in infoAll if i.startswith('END=')][0].split("=")[1])
                else:
                    bndList[line_list[2]]=line_list
                if svtype=="DEL" and source!='CNV':
                    if 'Dels' not in linked:
                        linked['Dels']={}
                        linked['Dels']['3to5']={}
                    if chrom1 not in linked['Dels']['3to5']:
                        linked['Dels']['3to5'][chrom1]={}
                    linked['Dels']['3to5'][chrom1][pos1]=[chrom1+":"+str(pos2), info, [i for i in infoAll if i.startswith("SOURCE=")][0].split("=")[1], pos2-pos1+1, line_list[5], [i for i in infoAll if i.startswith("HAP_ALLELIC_FRAC")][0].split("=")[1], [i for i in infoAll if i.startswith("ALLELIC_FRAC")][0].split("=")[1],[i for i in infoAll if i.startswith("PAIRS=")][0].split("=")[1], [i for i in infoAll if i.startswith("SPLIT=")][0].split("=")[1]]
                elif svtype=="DUP" and source!='CNV':
                    if 'Dups' not in linked:
                        linked['Dups']={}
                        linked['Dups']['5to3']={}
                    if chrom1 not in linked['Dups']['5to3']:
                        linked['Dups']['5to3'][chrom1]={}
                    linked['Dups']['5to3'][chrom1][pos1]=[chrom1+":"+str(pos2), info, [i for i in infoAll if i.startswith("SOURCE=")][0].split("=")[1], pos2-pos1+1, line_list[5], [i for i in infoAll if i.startswith("HAP_ALLELIC_FRAC")][0].split("=")[1], [i for i in infoAll if i.startswith("ALLELIC_FRAC")][0].split("=")[1],[i for i in infoAll if i.startswith("PAIRS=")][0].split("=")[1], [i for i in infoAll if i.startswith("SPLIT=")][0].split("=")[1]]
                elif svtype=="INV":
                    if 'Invs' not in linked:
                        linked['Invs']=dict()
                    if '5to5' not in linked['Invs']:
                        linked['Invs']['5to5']=dict()
                    if '3to3' not in linked['Invs']:
                        linked['Invs']['3to3']=dict()
                    if chrom1 not in linked['Invs']['5to5'].keys():
                        linked['Invs']['5to5'][chrom1]={}
                    if chrom1 not in linked['Invs']['3to3'].keys():
                        linked['Invs']['3to3'][chrom1]={}
                    linked['Invs']['5to5'][chrom1][pos1]=[chrom1+":"+str(pos2), info,[i for i in infoAll if i.startswith('SOURCE=')][0].split("=")[1], pos2-pos1+1, line_list[5], [i for i in infoAll if i.startswith("HAP_ALLELIC_FRAC")][0].split("=")[1], [i for i in infoAll if i.startswith("ALLELIC_FRAC")][0].split("=")[1],[i for i in infoAll if i.startswith("PAIRS=")][0].split("=")[1], [i for i in infoAll if i.startswith("SPLIT=")][0].split("=")[1]]
                    linked['Invs']['3to3'][chrom1][pos1]=[chrom1+":"+str(pos2), info,[i for i in infoAll if i.startswith('SOURCE=')][0].split("=")[1], pos2-pos1+1, line_list[5], [i for i in infoAll if i.startswith("HAP_ALLELIC_FRAC")][0].split("=")[1], [i for i in infoAll if i.startswith("ALLELIC_FRAC")][0].split("=")[1],[i for i in infoAll if i.startswith("PAIRS=")][0].split("=")[1], [i for i in infoAll if i.startswith("SPLIT=")][0].split("=")[1]]
                elif svtype=="UNK" and source!='CNV':
                    if 'UNK' not in linked:
                        linked['UNK']={}
                        linked['UNK']["."]={}
                    if chrom1 not in linked['UNK']['.']:
                        linked['UNK']['.'][chrom1]={}
                    linked['UNK']['.'][chrom1][pos1]=[chrom1+":"+str(pos2), info, [i for i in infoAll if i.startswith('SOURCE=')][0].split("=")[1], pos2-pos1+1, line_list[5], [i for i in infoAll if i.startswith("HAP_ALLELIC_FRAC")][0].split("=")[1], [i for i in infoAll if i.startswith("ALLELIC_FRAC")][0].split("=")[1],[i for i in infoAll if i.startswith("PAIRS=")][0].split("=")[1], [i for i in infoAll if i.startswith("SPLIT=")][0].split("=")[1]]
    finalDict, done=add_bnd(bndList, linked, done)
    return finalDict

def oppOrient(orient):
    if orient=='3to5':
        return '5to3'
    elif orient=='5to3':
        return '3to5'
    elif orient=='3to3':
        return '3to3'
    elif orient=='5to5':
        return '5to5'
    elif orient=='None' or orient==None:
        return 'None'
    else:
        return orient

def presentLinked(linked, chrom1, pos1, chrom2, pos2, svtype, orient, window, string):
    if string=="NotUNK":
        DICT1=linked[svtype][orient]
        if svtype!="Trans":
            if chrom1 in DICT1.keys():
                AllPOS1=sorted(DICT1[chrom1].keys())
                for position1 in AllPOS1:
                    if pos1 in range(position1-window, position1+window) and pos2 in range(int(DICT1[chrom1][position1][0].split(":")[1])-window, int(DICT1[chrom1][position1][0].split(":")[1])+window):
                        return [[chrom1, position1, DICT1[chrom1][position1][0], svtype, orient]+ DICT1[chrom1][position1][1:9], svtype, orient, chrom1, position1]
                    elif position1+window > pos1:
                        return None
                    else:
                        continue
            else:
                return None
        else:
            DICT2= linked[svtype][oppOrient(orient)]
            if chrom1 in DICT1.keys():
                AllPOS1=sorted(DICT1[chrom1].keys())
                for position1 in AllPOS1:
                    if pos1 in range(position1-window, position1+window) and pos2 in range(int(DICT1[chrom1][position1][0].split(":")[1])-window, int(DICT1[chrom1][position1][0].split(":")[1])+window) and DICT1[chrom1][position1][0].split(":")[0]==chrom2:
                        return [[chrom1, position1, DICT1[chrom1][position1][0], svtype, orient]+ DICT1[chrom1][position1][1:9], svtype, orient, chrom1, position1]
                    elif position1+window > pos1:
                        if chrom2 in DICT2.keys():
                            AllPOS2=sorted(DICT2[chrom2].keys())
                            for position2 in AllPOS2:
                                if pos2 in range(position2-window, position2+window) and pos1 in range(int(DICT2[chrom2][position2][0].split(":")[1])-window, int(DICT2[chrom2][position2][0].split(":")[1])+window) and DICT2[chrom2][position2][0].split(":")[0]==chrom1:
                                    return [[chrom1, DICT2[chrom2][position2][0].split(":")[1], str(chrom2)+":"+str(position2), svtype, orient]+ DICT2[chrom2][position2][1:9], svtype, oppOrient(orient), chrom2, position2]
                                elif position2+window > pos2:
                                    return None
                                else:
                                    continue
                    else:
                        if position1==AllPOS1[-1]:
                            if chrom2 in DICT2.keys():
                                AllPOS2=sorted(DICT2[chrom2].keys())
                                for position2 in AllPOS2:
                                    if pos2 in range(position2-window, position2+window) and pos1 in range(int(DICT2[chrom2][position2][0].split(":")[1])-window, int(DICT2[chrom2][position2][0].split(":")[1])+window) and DICT2[chrom2][position2][0].split(":")[0]==chrom1:
                                        return [[chrom1, DICT2[chrom2][position2][0].split(":")[1], str(chrom2)+":"+str(position2), svtype, orient]+ DICT2[chrom2][position2][1:9], svtype, oppOrient(orient), chrom2, position2]
                                    elif position2+window > pos2:
                                        return None
                                    else:
                                        continue
            else:
                # print "here"
                if chrom2 in DICT2.keys():
                    AllPOS2=sorted(DICT2[chrom2].keys())
                    for position2 in AllPOS2:
                        if pos2 in range(position2-window, position2+window) and pos1 in range(int(DICT2[chrom2][position2][0].split(":")[1])-window, int(DICT2[chrom2][position2][0].split(":")[1])+window) and DICT2[chrom2][position2][0].split(":")[0]==chrom1:
                            return [[chrom1, DICT2[chrom2][position2][0].split(":")[1], str(chrom2)+":"+str(position2), svtype, orient]+ DICT2[chrom2][position2][1:9], svtype, oppOrient(orient), chrom2, position2]
                        elif position2+window > pos2:
                            return None
                        else:
                            continue
    else:
        DICT1= linked['UNK']["."]
        if chrom1 in DICT1.keys():
            AllPOS1=sorted(DICT1[chrom1].keys())
            for position1 in AllPOS1:
                if pos1 in range(position1-window, position1+window) and pos2 in range(int(DICT1[chrom1][position1][0].split(":")[1])-window, int(DICT1[chrom1][position1][0].split(":")[1])+window):
                    return [[chrom1, position1, DICT1[chrom1][position1][0], "UNK", "."]+ DICT1[chrom1][position1][1:9], 'UNK','.',chrom1, position1]
                elif position1+window > pos1:
                    if chrom2!=chrom1:
                        if chrom2 in DICT1.keys():
                            AllPOS2=sorted(DICT1[chrom2].keys())
                            for position2 in AllPOS2:
                                if pos2 in range(position2-window, position2+window) and pos1 in range(int(DICT1[chrom2][position2][0].split(":")[1])-window, int(DICT1[chrom2][position2][0].split(":")[1])+window):
                                    return [[chrom2, position2, DICT1[chrom2][position2][0], "UNK", "."]+DICT1[chrom2][positions2][1:9], 'UNK', '.', chrom2, position2]
                                elif position2+window > pos2:
                                    return None
                                else:
                                    continue
                else:
                    continue
        else:
            return None

def compare_SR_LR(linked,short, w):
    list1=[]
    Category, Linked_chrom1, Linked_pos1, Linked_chrom2_pos2, Linked_SVType, Linked_Orient, Linked_Filter, Linked_Source, Linked_Size, Linked_Qual, Linked_HapAllelicFrac, Linked_AllelicFrac,  Linked_Pairs, Linked_Split =list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list()
    TofillList=[Linked_chrom1, Linked_pos1, Linked_chrom2_pos2, Linked_SVType, Linked_Orient, Linked_Filter, Linked_Source, Linked_Size, Linked_Qual, Linked_HapAllelicFrac, Linked_AllelicFrac,  Linked_Pairs, Linked_Split]
    short=short.reindex(range(0, short.shape[0])) # important to re-index so don't exclude it
    for index,row in short.iterrows():
        check1= presentLinked(linked, row["chrom1"], row['pos1'], row['chrom2:pos2'].split(":")[0], int(row['chrom2:pos2'].split(":")[1]), row['SVType'], row['Orientation'], w, "NotUNK")
        if check1==None:
            check2= presentLinked(linked, row["chrom1"], row['pos1'], row['chrom2:pos2'].split(":")[0], int(row['chrom2:pos2'].split(":")[1]), row['SVType'], row['Orientation'], w, "UNK")
            if check2==None:
                for k in TofillList:
                    k.append("NA")
                Category.append("Only cWGS")
            else:
                a=linked[check2[1]][check2[2]][check2[3]].pop(check2[4],None)
                for k in range(len(TofillList)):
                    TofillList[k].append(check2[0][k])
                Category.append("Common")
        else:
            b=linked[check1[1]][check1[2]][check1[3]].pop(check1[4],None)
            for k in range(len(TofillList)):
                TofillList[k].append(check1[0][k])
            Category.append("Common")
    short['Category']=Category
    short['Linked_chrom1']=Linked_chrom1
    short['Linked_pos1']=Linked_pos1
    short['Linked_chrom2_pos2']=Linked_chrom2_pos2
    short['Linked_SVType']=Linked_SVType
    short['Linked_Orient']=Linked_Orient
    short['Linked_Filter']=Linked_Filter
    short['Linked_Source']=Linked_Source
    short['Linked_Size']=Linked_Size
    short['Linked_Qual']=Linked_Qual
    short['Linked_HapAllelicFrac']=Linked_HapAllelicFrac
    short['Linked_AllelicFrac']=Linked_AllelicFrac
    short['Linked_Pairs']=Linked_Pairs
    short['Linked_Split']=Linked_Split
    # add remaining linked read seq SV to common dataframe
    for svtype in linked.keys():
        for orient in linked[svtype].keys():
            for chrom1 in linked[svtype][orient].keys():
                for pos1 in linked[svtype][orient][chrom1].keys():
                    tmp=[]
                    tmp=tmp+[chrom1, pos1, linked[svtype][orient][chrom1][pos1][0], svtype, linked[svtype][orient][chrom1][pos1][3], orient]
                    tmp.append(chrom1)
                    tmp.append(pos1)
                    tmp.append(linked[svtype][orient][chrom1][pos1][0])
                    tmp.append(svtype)
                    tmp.append(orient)
                    for i in range(5,13):
                        tmp.append(linked[svtype][orient][chrom1][pos1][i-4])
                    tmp.append("Only 10XWGS")
                    list1.append(tmp)
    header=['chrom1','pos1','chrom2:pos2','SVType', 'Size','Orientation','Linked_chrom1', 'Linked_pos1', 'Linked_chrom2_pos2', 'Linked_SVType', 'Linked_Orient', 'Linked_Filter', 'Linked_Source', 'Linked_Size', 'Linked_Qual', 'Linked_HapAllelicFrac', 'Linked_AllelicFrac',  'Linked_Pairs', 'Linked_Split', 'Category']
    df2=pd.DataFrame(list1 , columns=header)
    finalDF=pd.concat([short, df2], sort=False, ignore_index=True)
    return finalDF
    # finalDF.to_csv("Combined_SR_LR_"+str(w)+".csv",index=False)

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
    Commons['Category']=Commons['Category'].apply(lambda x: check_common(x))
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
    OnlyCalls['Category']="Only 10XWGS"
    AllDF=pd.concat([Commons, CombFile[~np.in1d(np.arange(CombFile.shape[0]), np.unique(i))], OnlyCalls], ignore_index=True, sort=False)
    return AllDF

### parse the arguments

def parse_arguments():
    parser= argparse.ArgumentParser(description = "A Python script to combine results from short-reads and linked-reads structural variation calls")
    parser.add_argument('-SR','--shortreads', help="enter .tsv file containing SV calls from short-reads (REQUIRED)")
    parser.add_argument('-LR1','--linkedDels', help="enter .vcf file containing small and mid sized deletion calls by longranger (REQUIRED)")
    parser.add_argument('-LR2','--linkedLarge', help="enter .vcf file containing large sized SV calls by longranger (REQUIRED)")
    parser.add_argument('-NAIBR','--NAIBR', help="enter .bed file containing SV calls from NAIBR (REQUIRED)")
    parser.add_argument('-GROCSV','--GROCSV', help="enter .vcf file containing SV calls from GROC-SV (REQUIRED)")
    parser.add_argument('-w','--window', help="window to match SV breakpoints from short-reads and linked-reads\n(DEFAULT=500bp)", default=int(500))
    parser.add_argument('-outdir','--outdir', help="enter path of output directory (DEFAULT=current directory)", type=str, default=".")
    return parser.parse_args()

if __name__=='__main__':
    args=parse_arguments()
    if (not args.shortreads or not args.linkedDels or not args.linkedLarge or not args.NAIBR or not args.GROCSV):
        print os.path.basename(main.__file__) + " missing file\n"
        sys.exit(1)
    else:
        linked, bndList={}, {}
        linked1= read_linked(args.linkedDels, linked)
        linked2= read_linked(args.linkedLarge, linked1)
        SV_short=pd.read_csv(args.shortreads, sep='\t')
        SV_short=SV_short.sort_values(['chrom1','pos1'])
        ### adding longranger calls
        CombFile=compare_SR_LR(linked2, SV_short, int(args.window))
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

        AllDels= add_calls(NAIBR_File[(NAIBR_File['NAIBR_SVType']=="Dels")], CombFile[CombFile['SVType']=="Dels"], int(args.window), "Dels", "NAIBR")
        AllDups= add_calls(NAIBR_File[(NAIBR_File['NAIBR_SVType']=="Dups")], CombFile[CombFile['SVType']=="Dups"], int(args.window), "Dups", "NAIBR")
        AllInvs= add_calls(NAIBR_File[(NAIBR_File['NAIBR_SVType']=="Invs")], CombFile[CombFile['SVType']=="Invs"], int(args.window), "Invs", "NAIBR")
        AllTrans= add_calls(NAIBR_File[(NAIBR_File['NAIBR_SVType']=="Trans")], CombFile[CombFile['SVType']=="Trans"], int(args.window), "Trans", "NAIBR")
        FinalFile1= pd.concat([AllDels, AllDups, AllInvs, AllTrans], ignore_index=True, sort=False)
        FinalFile1= FinalFile1.sort_values(by=['chrom1','pos1','chrom2:pos2'])

        ### add GROC-SV calls
        grocCalls= read_GROCSV(args.GROCSV)
        AllDels= add_calls(grocCalls[(grocCalls['GROCSV_SVType']=="Dels")], FinalFile1[FinalFile1['SVType']=="Dels"], int(args.window), "Dels", "GROCSV")
        AllDups= add_calls(grocCalls[(grocCalls['GROCSV_SVType']=="Dups")], FinalFile1[FinalFile1['SVType']=="Dups"], int(args.window), "Dups", "GROCSV")
        AllInvs= add_calls(grocCalls[(grocCalls['GROCSV_SVType']=="Invs")], FinalFile1[FinalFile1['SVType']=="Invs"], int(args.window), "Invs", "GROCSV")
        AllTrans= add_calls(grocCalls[(grocCalls['GROCSV_SVType']=="Trans")], FinalFile1[FinalFile1['SVType']=="Trans"], int(args.window), "Trans", "GROCSV")
        FinalFile= pd.concat([AllDels, AllDups, AllInvs, AllTrans], ignore_index=True, sort=False)
        FinalFile= FinalFile.sort_values(by=['chrom1','pos1','chrom2:pos2'])

        FinalFile.to_csv(str(args.outdir)+"/Combined_SR_LR.tsv", sep='\t', index=False)
