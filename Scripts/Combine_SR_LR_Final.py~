
"""
Description: Script to combine SV calls from short-reads and linked-reads from LongRanger
Author: Riccha Sethi
"""

### Load packages
import argparse, __main__ as main
import pandas as pd
import gzip
import os, sys
pd.options.mode.chained_assignment = None
CurVersion=1.0

### read linked-reads vcf files
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
    Predicted_by, Linked_chrom1, Linked_pos1, Linked_chrom2_pos2, Linked_SVType, Linked_Orient, Linked_Filter, Linked_Source, Linked_Size, Linked_Qual, Linked_HapAllelicFrac, Linked_AllelicFrac,  Linked_Pairs, Linked_Split =list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list()
    TofillList=[Linked_chrom1, Linked_pos1, Linked_chrom2_pos2, Linked_SVType, Linked_Orient, Linked_Filter, Linked_Source, Linked_Size, Linked_Qual, Linked_HapAllelicFrac, Linked_AllelicFrac,  Linked_Pairs, Linked_Split]
    short=short.reindex(range(0, short.shape[0])) # important to re-index so don't exclude it
    for index,row in short.iterrows():
        check1= presentLinked(linked, row["chrom1"], row['pos1'], row['chrom2:pos2'].split(":")[0], int(row['chrom2:pos2'].split(":")[1]), row['SVType'], row['Orientation'], w, "NotUNK")
        if check1==None:
            check2= presentLinked(linked, row["chrom1"], row['pos1'], row['chrom2:pos2'].split(":")[0], int(row['chrom2:pos2'].split(":")[1]), row['SVType'], row['Orientation'], w, "UNK")
            if check2==None:
                for k in TofillList:
                    k.append("NA")
                Predicted_by.append("Only cWGS")
            else:
                a=linked[check2[1]][check2[2]][check2[3]].pop(check2[4],None)
                for k in range(len(TofillList)):
                    TofillList[k].append(check2[0][k])
                Predicted_by.append("Common")
        else:
            b=linked[check1[1]][check1[2]][check1[3]].pop(check1[4],None)
            for k in range(len(TofillList)):
                TofillList[k].append(check1[0][k])
            Predicted_by.append("Common")
    short['Predicted_by']=Predicted_by
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
    header=['chrom1','pos1','chrom2:pos2','SVType', 'Size','Orientation','Linked_chrom1', 'Linked_pos1', 'Linked_chrom2_pos2', 'Linked_SVType', 'Linked_Orient', 'Linked_Filter', 'Linked_Source', 'Linked_Size', 'Linked_Qual', 'Linked_HapAllelicFrac', 'Linked_AllelicFrac',  'Linked_Pairs', 'Linked_Split', 'Predicted_by']
    df2=pd.DataFrame(list1 , columns=header)
    finalDF=pd.concat([short, df2], sort=False, ignore_index=True)
    finalDF.to_csv("Combined_SR_LR_"+str(w)+".csv",index=False)

### parse the arguments

def parse_arguments():
    parser= argparse.ArgumentParser(description = "A Python script to combine results from short-reads and longranger vcf SV calls")
    parser.add_argument('-SR','--shortreads',help="enter .csv file containing SV calls from short-reads seq (REQUIRED)")
    parser.add_argument('-LR1','--linkedDels', help="enter vcf file containing small and mid sized deletion calls by longranger (REQUIRED)")
    parser.add_argument('-LR2','--linkedLarge', help="enter vcf file containing large sized SV calls by longranger (REQUIRED)")
    parser.add_argument('-w','--window', help="window to match SV breakpoints from short-reads and longranger\n(DEFAULT=100bp)", default=int(100))
    parser.add_argument('--version', action='version', version='%(prog)s '+str(CurVersion))
    return parser.parse_args()

if __name__=='__main__':
    args=parse_arguments()
    if (not args.shortreads or not args.linkedDels or not args.linkedLarge):
        print os.path.basename(main.__file__) + " missing file\n"
        sys.exit(1)
    else:
        linked, bndList={}, {}
        linked1= read_linked(args.linkedDels, linked)
        linked2= read_linked(args.linkedLarge, linked1)
        SV_short=pd.read_csv(args.shortreads)
        SV_short=SV_short.sort_values(['chrom1','pos1'])
        comparedDF=compare_SR_LR(linked2, SV_short, int(args.window))
