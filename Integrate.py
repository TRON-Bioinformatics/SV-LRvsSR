#!/usr/bin/env python

"""
Author: Riccha Sethi
Description: Integrating SV calls from delly, lumpy and SvABA
"""

import argparse
import itertools
import pandas as pd

def extractSV_Delly(File,SVType, allchromo):
    DICT={}
    with open(File,'r') as f:
        for line in f:
            if line.startswith("#"):
                continue
            else:
                line_list=line.rstrip('\n').split()
                info=line_list[7].split(";")
                chrom1,pos1,chrom2,pos2=line_list[0],int(line_list[1]),[i for i in info if i.startswith('CHR2=')][0].split("=")[1],int([i for i in info if i.startswith('END=')][0].split("=")[1])
                CT=[i for i in info if i.startswith('CT=')][0].split("=")[1]
                size=int(pos2)-int(pos1)
                if SVType in ["Dels","Dups"] and size>50 and chrom1 in all_chromo and chrom2 in all_chromo:
                    if chrom1 not in DICT:
                        DICT[chrom1]={}
                    if CT not in DICT[chrom1]:
                        DICT[chrom1][CT]={}
                    DICT[chrom1][CT][pos1]=[pos1,pos2,chrom2,CT,line_list[6],"D",SVType]
                elif SVType in ["Trans","Invs"] and chrom1 in all_chromo and chrom2 in all_chromo:
                    if CT not in DICT:
                        DICT[CT]={}
                    if chrom1 not in DICT[CT]:
                        DICT[CT][chrom1]={}
                    if chrom2 not in DICT[CT][chrom1]:
                        DICT[CT][chrom1][chrom2]={}
                    DICT[CT][chrom1][chrom2][pos1]=[pos1,pos2,chrom2,CT,line_list[6],"D",SVType]
    return DICT

def read_delly(file_del,file_dup,file_invs,file_trans, all_chromo):
    dels=extractSV_Delly(file_del,"Dels", all_chromo)
    dups=extractSV_Delly(file_dup,"Dups", all_chromo)
    invs=extractSV_Delly(file_invs,"Invs", all_chromo)
    trans=extractSV_Delly(file_trans,"Trans", all_chromo)
    return [dels,dups,invs,trans]

def add_entry(chrom1,pos1,chrom2,pos2,svtype,CTF,dels,dups,invs,trans, all_chromo):
    size=int(pos2)-int(pos1)
    if svtype=="DEL" and size>50 and chrom1 in all_chromo and chrom2 in all_chromo:
        if chrom1 not in dels:
            dels[chrom1]={}
        if CTF not in dels[chrom1]:
            dels[chrom1][CTF]={}
        dels[chrom1][CTF][pos1]=[pos1,pos2,chrom2,CTF,"L","Dels"]
    elif svtype=="DUP" and size>50 and chrom1 in all_chromo and chrom2 in all_chromo:
        if chrom1 not in dups:
            dups[chrom1]={}
        if CTF not in dups[chrom1]:
            dups[chrom1][CTF]={}
        dups[chrom1][CTF][pos1]=[pos1,pos2,chrom2,CTF,"L","Dups"]
    elif svtype=="INV" and size>50 and chrom1 in all_chromo and chrom2 in all_chromo:
        if CTF not in invs:
            invs[CTF]={}
        if chrom1 not in invs[CTF]:
            invs[CTF][chrom1]={}
        if chrom2 not in invs[CTF][chrom1]:
            invs[CTF][chrom1][chrom2]={}
        invs[CTF][chrom1][chrom2][pos1]=[pos1,pos2,chrom2,CTF,"L","Invs"]
    elif svtype=="TRA" and chrom1 in all_chromo and chrom2 in all_chromo:
        if CTF not in trans:
            trans[CTF]={}
        if chrom1 not in trans[CTF]:
            trans[CTF][chrom1]={}
        if chrom2 not in trans[CTF][chrom1]:
            trans[CTF][chrom1][chrom2]={}
        trans[CTF][chrom1][chrom2][pos1]=[pos1,pos2,chrom2,CTF,"L","Trans"]
    return dels,dups,invs,trans

def read_lumpy(fileAll, all_chromo):
    dels,dups,invs,trans={},{},{},{}
    with open(fileAll,'r') as f:
        for line in f:
            if line.startswith("#"):
                continue
            else:
                line_list=line.rstrip('\n').split()
                chrom1,pos1, info=line_list[0],int(line_list[1]), line_list[7].split(";")
                svtype=[i for i in info if i.startswith("SVTYPE=")][0].split("=")[1]
                if svtype!="BND":
                    CT=[i for i in info if i.startswith("STRANDS=")][0].split("=")[1].split(":")[0]
                    chrom2=line_list[0]
                    pos2=int([i for i in info if i.startswith("END=")][0].split("=")[1])
                    if CT=="++":
                        CTF="3to3"
                    elif CT=="--":
                        CTF="5to5"
                    elif CT=="+-":
                        CTF="3to5"
                    elif CT=="-+":
                        CTF="5to3"
                    dels,dups,invs,trans= add_entry(chrom1,pos1,chrom2,pos2,svtype,CTF,dels,dups,invs,trans, all_chromo)
                else:
                    if int(line_list[2].split("_")[1])==1:
                        CT=[i for i in info if i.startswith("STRANDS=")][0].split("=")[1].split(":")
                        chrom2=str(line_list[4].split(":")[0].split("[")[-1].split("]")[-1])
                        pos2=int(line_list[4].split(":")[1].split("[")[0].split("]")[0])
                        if "++" in CT:
                            CTF="3to3"
                            if chrom1==chrom2:
                                svtype="INV"
                            else:
                                svtype="TRA"
                        elif "--" in CT:
                            CTF="5to5"
                            if chrom1==chrom2:
                                svtype="INV"
                            else:
                                svtype="TRA"
                        elif "+-" in CT:
                            CTF="3to5"
                            if chrom1==chrom2:
                                svtype="DEL"
                            else:
                                svtype="TRA"
                        elif "-+" in CT:
                            CTF="5to3"
                            if chrom1==chrom2:
                                svtype="DUP"
                            else:
                                svtype="TRA"
                        dels,dups,invs,trans= add_entry(chrom1,pos1,chrom2,pos2,svtype,CTF,dels,dups,invs,trans, all_chromo)
    return [dels,dups,invs,trans]

def check_category(mate1,mate2):
    chrB=str(mate1.split(":")[0].split("]")[-1].split("[")[-1])
    pos2=int(mate1.split(":")[1].split("[")[0].split("]")[0])
    chrA=str(mate2.split(":")[0].split("]")[-1].split("[")[-1])
    pos1=int(mate2.split(":")[1].split("[")[0].split("]")[0])
    if chrA==chrB and "[" in mate1 and "]" in mate2:
        return "category1","Dels", chrA, pos1, chrA, pos2, "3to5"
    elif chrA==chrB and "]" in mate1 and "[" in mate2:
        return "category2","Dups",chrA, pos1, chrA, pos2, "5to3"
    elif chrA==chrB and ("]" in mate1 and "]" in mate2):
        return "category3","Invs",chrA, pos1, chrA, pos2,"3to3"
    elif chrA==chrB and ("[" in mate1 and "[" in mate2):
        return "category3","Invs", chrA, pos1, chrA, pos2,"5to5"
    else:
        if "[" in mate1 and "]" in mate2:
            return "category4","Trans",chrA, pos1, chrB, pos2,"3to5"
        elif "]" in mate1 and "[" in mate2:
            return "category4","Trans",chrA, pos1, chrB, pos2,"5to3"
        elif "]" in mate1 and "]" in mate2:
            return "category4","Trans",chrA, pos1, chrB, pos2,"3to3"
        elif "[" in mate1 and "[" in mate2:
            return "category4","Trans",chrA, pos1, chrB, pos2,"5to5"

def read_svaba(File,mode,all_chromo):
    dels,dups,invs,trans={},{},{},{}
    with open(File,'r') as f:
        for line1,line2 in itertools.izip_longest(f, f, fillvalue=''):
            line1_list=line1.rstrip('\n').split('\t')
            line2_list=line2.rstrip('\n').split('\t')
            if mode=="germline":
                feature_list=line1_list[9].split(":")
            else:
                feature_list=line1_list[10].split(":")
            info1, info2=line1_list[7].split(";"), line2_list[7].split(";")
            if int(line1_list[2].split(":")[1])==1:
                mate1,mate2=line1_list[4],line2_list[4]
            else:
                mate1,mate2=line2_list[4],line1_list[4]
            cat=check_category(mate1,mate2)
            chrom1,pos1,chrom2,pos2,CT=cat[2],cat[3],cat[4],cat[5],cat[6]
            if pos2-pos1>50 and chrom1 in all_chromo and chrom2 in all_chromo and "Trans" not in cat:
                if "Dels" in cat:
                    if chrom1 not in dels:
                        dels[chrom1]={}
                    if CT not in dels[chrom1]:
                        dels[chrom1][CT]={}
                    dels[chrom1][CT][pos1]=[pos1,pos2,chrom2,CT,line1_list[6],"S","Dels"]
                elif "Dups" in cat:
                    if chrom1 not in dups:
                        dups[chrom1]={}
                    if CT not in dups[chrom1]:
                        dups[chrom1][CT]={}
                    dups[chrom1][CT][pos1]=[pos1,pos2,chrom2,CT,line1_list[6],"S","Dups"]
                elif "Invs" in cat:
                    if CT not in invs:
                        invs[CT]={}
                    if chrom1 not in invs[CT]:
                        invs[CT][chrom1]={}
                    if chrom2 not in invs[CT][chrom1]:
                        invs[CT][chrom1][chrom2]={}
                    invs[CT][chrom1][chrom2][pos1]=[pos1,pos2,chrom2,CT,line1_list[6],"S","Invs"]
            if "Trans" in cat and chrom1 in all_chromo and chrom2 in all_chromo:
                if CT not in trans:
                    trans[CT]={}
                if chrom1 not in trans[CT]:
                    trans[CT][chrom1]={}
                if chrom2 not in trans[CT][chrom1]:
                    trans[CT][chrom1][chrom2]={}
                trans[CT][chrom1][chrom2][pos1]=[pos1,pos2,chrom2,CT,line1_list[6],"S","Trans"]
    return [dels,dups,invs,trans]

def find_oppOrient(orient):
    if orient=="3to5":
        return "5to3"
    elif orient=="5to3":
        return "3to5"
    elif orient=="3to3":
        return "3to3"
    elif orient=="5to5":
        return "5to5"

def check_opp(chrom1,posA,chrom2,posB, orient, window, DictOpp, done):
    count=0
    if chrom2 in DictOpp:
        DictOrientOpp_chrom2=DictOpp[chrom2]
    else:
        DictOrientOpp_chrom2={}
    if chrom1 in DictOrientOpp_chrom2:
        DictOrientOpp=DictOrientOpp_chrom2[chrom1]
    else:
        DictOrientOpp={}
    if DictOrientOpp:
        PosBList=sorted(DictOrientOpp.keys())
        LEN=len(PosBList)
        while PosBList:
            PosBCompare=PosBList.pop(0)
            count+=1
            value=DictOrientOpp[PosBCompare]
            if count==1:
                tmpDelly, tmpLumpy, tmpSvaba=['NA']*6, ['NA']*5, ['NA']*6
            if PosBCompare < posB-window:
                continue
            elif PosBCompare in range(posB-window, posB+window) and chrom1==value[2] and value[1] in range(posA-window, posA+window):
                done.append([chrom2]+value)
                if 'D' in value and tmpDelly[0]=='NA':
                    tmpDelly=[chrom2]+value[0:-2]
                elif 'D' in value and tmpDelly[0]!='NA':
                    tmpDelly[0]=str(tmpDelly[0])+";"+str(chrom2)
                    tmpDelly[1]=str(tmpDelly[1])+";"+str(value[0])
                    tmpDelly[2]=str(tmpDelly[2])+";"+str(value[1])
                    tmpDelly[3]=str(tmpDelly[3])+";"+str(value[2])
                    tmpDelly[4]=str(tmpDelly[4])+";"+str(value[3])
                    if 'PASS' in tmpDelly[5] or 'PASS' in value:
                        tmpDelly[5]="PASS"
                if 'L' in value and tmpLumpy[0]=='NA':
                    tmpLumpy=[chrom2]+value[0:-2]
                elif 'L' in value and tmpLumpy[0]!='NA':
                    tmpLumpy[0]=str(tmpLumpy[0])+";"+str(chrom2)
                    tmpLumpy[1]=str(tmpLumpy[1])+";"+str(value[0])
                    tmpLumpy[2]=str(tmpLumpy[2])+";"+str(value[1])
                    tmpLumpy[3]=str(tmpLumpy[3])+";"+str(value[2])
                    tmpLumpy[4]=str(tmpLumpy[4])+";"+str(value[3])
                if 'S' in value and tmpSvaba[0]=='NA':
                    tmpSvaba=[chrom2]+value[0:-2]
                elif 'S' in value and tmpSvaba[0]!='NA':
                    tmpSvaba[0]=str(tmpSvaba[0])+";"+str(chrom2)
                    tmpSvaba[1]=str(tmpSvaba[1])+";"+str(value[0])
                    tmpSvaba[2]=str(tmpSvaba[2])+";"+str(value[1])
                    tmpSvaba[3]=str(tmpSvaba[3])+";"+str(value[2])
                    tmpSvaba[4]=str(tmpSvaba[4])+";"+str(value[3])
                    if 'PASS' in tmpSvaba[5] or 'PASS' in value:
                        tmpSvaba[5]="PASS"
                if LEN==count:
                    return tmpDelly,tmpLumpy,tmpSvaba,done
            elif PosBCompare > posB+window:
                return tmpDelly,tmpLumpy,tmpSvaba,done

def joinCalls_trans(List, chromo1, chromo2, orient1, DictOpp_delly, DictOpp_lumpy, DictOpp_svaba, window, done, svtype):
    MergedList=[]
    for i in range(len(List)):
        pos1,pos2,chrom1,chrom2, orient=List[i][0], List[i][1], chromo1,chromo2, orient1
        call=List[i]
        if [chrom1]+call not in done:
            done.append([chrom1]+call)
            if svtype=="Trans":
                tmpExtras=[chromo1, pos1, str(chrom2)+":"+str(pos2), call[-1], 0, call[3]]
            else:
                tmpExtras=[chromo1, pos1, str(chrom2)+":"+str(pos2), call[-1], int(call[1])-int(call[0])+1, call[3]]
            tmpDelly, tmpLumpy, tmpSvaba=['NA']*6, ['NA']*5, ['NA']*6
            if 'D' in call:
                tmpDelly=[chromo1]+call[0:-2]
            elif 'L' in call:
                tmpLumpy=[chromo1]+call[0:-2]
            elif 'S' in call:
                tmpSvaba=[chromo1]+call[0:-2]
            if i+1==len(List):
                MergedList= MergedList+[tmpExtras+tmpDelly+tmpLumpy+tmpSvaba] #MergedList.append(tmpExtras+tmpDelly+tmpLumpy+tmpSvaba)
            else:
                for j in range(i+1, len(List)):
                    callN=List[j]
                    pos1N,pos2N,chrom1N,chrom2N, orientN= callN[0], callN[1], chromo1, callN[2], callN[3]
                    if pos1 in range(pos1N-window, pos1N+window) and pos2 in range(pos2N-window, pos2N+window):
                        if [chrom1N]+callN not in done:
                            done.append([chrom1N]+callN)
                            if 'D' in callN and tmpDelly[0]=='NA':
                                tmpDelly=[chrom1N]+callN[0:-2]
                            elif 'D' in callN and tmpDelly[0]!='NA':
                                tmpDelly[0]=str(tmpDelly[0])+";"+str(chrom1N)
                                tmpDelly[1]=str(tmpDelly[1])+";"+str(callN[0])
                                tmpDelly[2]=str(tmpDelly[2])+";"+str(callN[1])
                                tmpDelly[3]=str(tmpDelly[3])+";"+str(callN[2])
                                tmpDelly[4]=str(tmpDelly[4])+";"+str(callN[3])
                                if 'PASS' in tmpDelly[5] or 'PASS' in callN:
                                    tmpDelly[5]="PASS"
                            if 'L' in callN and tmpLumpy[0]=='NA':
                                tmpLumpy=[chrom1N]+callN[0:-2]
                            elif 'L' in callN and tmpLumpy[0]!='NA':
                                tmpLumpy[0]=str(tmpLumpy[0])+";"+str(chrom1N)
                                tmpLumpy[1]=str(tmpLumpy[1])+";"+str(callN[0])
                                tmpLumpy[2]=str(tmpLumpy[2])+";"+str(callN[1])
                                tmpLumpy[3]=str(tmpLumpy[3])+";"+str(callN[2])
                                tmpLumpy[4]=str(tmpLumpy[4])+";"+str(callN[3])
                            if 'S' in callN and tmpSvaba[0]=='NA':
                                tmpSvaba=[chrom1N]+callN[0:-2]
                            elif 'S' in callN and tmpSvaba[0]!='NA':
                                tmpSvaba[0]=str(tmpSvaba[0])+";"+str(chrom1N)
                                tmpSvaba[1]=str(tmpSvaba[1])+";"+str(callN[0])
                                tmpSvaba[2]=str(tmpSvaba[2])+";"+str(callN[1])
                                tmpSvaba[3]=str(tmpSvaba[3])+";"+str(callN[2])
                                tmpSvaba[4]=str(tmpSvaba[4])+";"+str(callN[3])
                                if 'PASS' in tmpSvaba[5] or 'PASS' in callN:
                                    tmpSvaba[5]="PASS"
                            if svtype=="Trans":
                                OppD=check_opp(chrom1, pos1, chrom2, pos2, orient, int(window), DictOpp_delly, done)
                                OppL=check_opp(chrom1, pos1, chrom2, pos2, orient, int(window), DictOpp_lumpy, done)
                                OppS=check_opp(chrom1, pos1, chrom2, pos2, orient, int(window), DictOpp_svaba, done)
                                if OppD!=None:
                                    done=OppD[3]
                                if OppL!=None:
                                    done=OppL[3]
                                if OppS!=None:
                                    done=OppS[3]
                                if OppD!=None and OppD[0][1]!='NA' and tmpDelly[1] != 'NA':
                                    tmpDelly[0]=str(tmpDelly[0])+";"+str(OppD[0][0])
                                    tmpDelly[1]=str(tmpDelly[1])+";"+str(OppD[0][1])
                                    tmpDelly[2]=str(tmpDelly[2])+";"+str(OppD[0][2])
                                    tmpDelly[3]=str(tmpDelly[3])+";"+str(OppD[0][3])
                                    tmpDelly[4]=str(tmpDelly[4])+";"+str(OppD[0][4])
                                    if 'PASS' in tmpDelly[5] or 'PASS' in OppD[0][5]:
                                        tmpDelly[5]="PASS"
                                elif OppD!=None and OppD[0][1]!='NA' and tmpDelly[1]=='NA':
                                    tmpDelly=OppD[0]
                                if OppL!=None and OppL[1][1]!='NA' and tmpLumpy[1]!='NA':
                                    tmpLumpy[0]=str(tmpLumpy[0])+";"+str(OppL[1][0])
                                    tmpLumpy[1]=str(tmpLumpy[1])+";"+str(OppL[1][1])
                                    tmpLumpy[2]=str(tmpLumpy[2])+";"+str(OppL[1][2])
                                    tmpLumpy[3]=str(tmpLumpy[3])+";"+str(OppL[1][3])
                                    tmpLumpy[4]=str(tmpLumpy[4])+";"+str(OppL[1][4])
                                elif OppL!=None and OppL[1][1]!='NA' and tmpLumpy[1]=='NA':
                                    tmpLumpy=OppL[1]
                                if OppS!=None and OppS[2][0]!='NA' and tmpSvaba[1]!='NA':
                                    tmpSvaba[0]=str(tmpSvaba[0])+";"+str(OppS[2][0])
                                    tmpSvaba[1]=str(tmpSvaba[1])+";"+str(OppS[2][1])
                                    tmpSvaba[2]=str(tmpSvaba[2])+";"+str(OppS[2][2])
                                    tmpSvaba[3]=str(tmpSvaba[3])+";"+str(OppS[2][3])
                                    tmpSvaba[4]=str(tmpSvaba[4])+";"+str(OppS[2][4])
                                    if 'PASS' in tmpSvaba[5] or 'PASS' in OppS[2][5]:
                                        tmpSvaba[5]="PASS"
                                elif OppS!=None and OppS[2][1]!='NA' and tmpSvaba[1]=='NA':
                                    tmpSvaba=OppS[2]
                            if j+1==len(List):
                                MergedList= MergedList+[tmpExtras+tmpDelly+tmpLumpy+tmpSvaba]
                    elif pos1 in range(pos1N -window, pos1N+window) and pos2< pos2N-window:
                        continue
                    else:
                        if svtype=="Trans":
                            OppD=check_opp(chrom1, pos1, chrom2, pos2, orient, int(window), DictOpp_delly, done)
                            OppL=check_opp(chrom1, pos1, chrom2, pos2, orient, int(window), DictOpp_lumpy, done)
                            OppS=check_opp(chrom1, pos1, chrom2, pos2, orient, int(window), DictOpp_svaba, done)
                            if OppD!=None:
                                done=OppD[3]
                            if OppL!=None:
                                done=OppL[3]
                            if OppS!=None:
                                done=OppS[3]
                            if OppD!=None and OppD[0][1]!='NA' and tmpDelly[1] != 'NA':
                                tmpDelly[0]=str(tmpDelly[0])+";"+str(OppD[0][0])
                                tmpDelly[1]=str(tmpDelly[1])+";"+str(OppD[0][1])
                                tmpDelly[2]=str(tmpDelly[2])+";"+str(OppD[0][2])
                                tmpDelly[3]=str(tmpDelly[3])+";"+str(OppD[0][3])
                                tmpDelly[4]=str(tmpDelly[4])+";"+str(OppD[0][4])
                                if 'PASS' in tmpDelly[5] or 'PASS' in OppD[0][5]:
                                    tmpDelly[5]="PASS"
                            elif OppD!=None and OppD[0][1]!='NA' and tmpDelly[1]=='NA':
                                tmpDelly=OppD[0]
                            if OppL!=None and OppL[1][1]!='NA' and tmpLumpy[1]!='NA':
                                tmpLumpy[0]=str(tmpLumpy[0])+";"+str(OppL[1][0])
                                tmpLumpy[1]=str(tmpLumpy[1])+";"+str(OppL[1][1])
                                tmpLumpy[2]=str(tmpLumpy[2])+";"+str(OppL[1][2])
                                tmpLumpy[3]=str(tmpLumpy[3])+";"+str(OppL[1][3])
                                tmpLumpy[4]=str(tmpLumpy[4])+";"+str(OppL[1][4])
                            elif OppL!=None and OppL[1][1]!='NA' and tmpLumpy[1]=='NA':
                                tmpLumpy=OppL[1]
                            if OppS!=None and OppS[2][0]!='NA' and tmpSvaba[1]!='NA':
                                tmpSvaba[0]=str(tmpSvaba[0])+";"+str(OppS[2][0])
                                tmpSvaba[1]=str(tmpSvaba[1])+";"+str(OppS[2][1])
                                tmpSvaba[2]=str(tmpSvaba[2])+";"+str(OppS[2][2])
                                tmpSvaba[3]=str(tmpSvaba[3])+";"+str(OppS[2][3])
                                tmpSvaba[4]=str(tmpSvaba[4])+";"+str(OppS[2][4])
                                if 'PASS' in tmpSvaba[5] or 'PASS' in OppS[2][5]:
                                    tmpSvaba[5]="PASS"
                            elif OppS!=None and OppS[2][1]!='NA' and tmpSvaba[1]=='NA':
                                tmpSvaba=OppS[2]
                        MergedList= MergedList+[tmpExtras+tmpDelly+tmpLumpy+tmpSvaba]
                        break
    return MergedList, done

def joinCalls(List, chromo, window):
    MergedList, done=[],[]
    for i in range(len(List)):
        start,end,chrom2=List[i][0],List[i][1],List[i][2]
        call=List[i]
        if call not in done:
            done.append(call)
            size=int(call[1])-int(call[0])+1
            tmpExtras=[chromo, call[0], str(call[2])+":"+str(call[1]), call[-1], size, call[3]]
            tmpDelly, tmpLumpy, tmpSvaba=['NA']*6, ['NA']*5, ['NA']*5
            if 'D' in call:
                tmpDelly=[chromo]+call[0:-2]
            elif 'L' in call:
                tmpLumpy=[chromo]+call[0:-2]
            elif 'S' in call:
                tmpSvaba=[chromo]+call[0:-2]
            if (i+1)==len(List):
                MergedList= MergedList+[tmpExtras+tmpDelly+tmpLumpy+tmpSvaba]
            else:
                for j in range(i+1, len(List)):
                    callN=List[j]
                    startN,endN=callN[0], callN[1]
                    if start in range(startN-window, startN+window+1) and end in range(endN-window, endN+window+1):
                        if callN not in done:
                            done.append(callN)
                            if 'D' in callN and tmpDelly[0]=='NA':
                                tmpDelly=[chromo]+callN[0:-2]
                            elif 'D' in callN and tmpDelly[0]!='NA':
                                tmpDelly[1]=str(tmpDelly[1])+";"+str(callN[0])
                                tmpDelly[2]=str(tmpDelly[2])+";"+str(callN[1])
                                if 'PASS' in tmpDelly[5] or 'PASS' in callN:
                                    tmpDelly[5]="PASS"
                            if 'L' in callN and tmpLumpy[0]=='NA':
                                tmpLumpy=[chromo]+callN[0:-2]
                            elif 'L' in callN and tmpLumpy[0]!='NA':
                                tmpLumpy[1]=str(tmpLumpy[1])+";"+str(callN[0])
                                tmpLumpy[2]=str(tmpLumpy[2])+";"+str(callN[1])
                            if 'S' in callN and tmpSvaba[0]=='NA':
                                tmpSvaba=[chromo]+callN[0:-2]
                            elif 'S' in callN and tmpSvaba[0]!='NA':
                                tmpSvaba[1]=str(tmpSvaba[1])+";"+str(callN[0])
                                tmpSvaba[2]=str(tmpSvaba[2])+";"+str(callN[1])
                            if j+1==len(List):
                                MergedList= MergedList+[tmpExtras+tmpDelly+tmpLumpy+tmpSvaba]
                    elif start in range(startN-window, startN+window) and end< endN-window:
                        continue
                    else:
                        MergedList= MergedList+[tmpExtras+tmpDelly+tmpLumpy+tmpSvaba]
                        break
    return MergedList

def not_empty(DICT, chromo):
    if chromo in DICT:
        return DICT[chromo]
    else:
        return {}

def combine_all(svtype, delly, lumpy, svaba, window):
    MergedList=[]
    if svtype in ["Dels","Dups"]:
        for chromo in list(set(delly.keys()+lumpy.keys()+svaba.keys())):
            Ddict,Ldict,Sdict=not_empty(delly,chromo), not_empty(lumpy,chromo), not_empty(svaba,chromo)
            for orient in list(set(Ddict.keys() +Ldict.keys() + Sdict.keys())):
                DStartList, LStartList, SStartList= not_empty(Ddict,orient), not_empty(Ldict, orient), not_empty(Sdict, orient)
                combinedSorted=sorted(DStartList.values()+LStartList.values()+SStartList.values(), key=  lambda x: (int(x[0]), int(x[1])))
                Combined=joinCalls(combinedSorted, chromo, int(window))
                MergedList=MergedList+list(Combined)
    elif svtype=="Invs":
        done=[]
        for orient in ['3to3','5to5']:
            Ddict,Ldict,Sdict= not_empty(delly,orient), not_empty(lumpy,orient), not_empty(svaba,orient)
            for chrom1 in list(set(Ddict.keys() +Ldict.keys() +Sdict.keys())):
                DChrom2List, LChrom2List, SChrom2List= not_empty(Ddict, chrom1), not_empty(Ldict, chrom1), not_empty(Sdict, chrom1)
                for chrom2 in list(set(DChrom2List.keys() +LChrom2List.keys() +SChrom2List.keys())):
                    DStartList,LStartList,SStartList= not_empty(DChrom2List,chrom2), not_empty(LChrom2List,chrom2), not_empty(SChrom2List,chrom2)
                    combinedSorted=sorted(DStartList.values()+ LStartList.values()+ SStartList.values(), key=lambda x:(int(x[0]),int(x[1])))
                    Combined, done=joinCalls_trans(combinedSorted, chrom1, chrom2, orient,  not_empty(delly,find_oppOrient(orient)), not_empty(lumpy,find_oppOrient(orient)), not_empty(svaba,find_oppOrient(orient)), int(window), done, svtype)
                    MergedList=MergedList+list(Combined)
    elif svtype=="Trans":
        done=[]
        for orient in ['3to5','5to3','3to3','5to5']:
            Ddict,Ldict,Sdict= not_empty(delly,orient), not_empty(lumpy,orient), not_empty(svaba,orient)
            for chrom1 in list(set(Ddict.keys() +Ldict.keys() +Sdict.keys())):
                DChrom2List, LChrom2List, SChrom2List= not_empty(Ddict, chrom1), not_empty(Ldict, chrom1), not_empty(Sdict, chrom1)
                for chrom2 in list(set(DChrom2List.keys() +LChrom2List.keys() + SChrom2List.keys())):
                    DStartList,LStartList,SStartList= not_empty(DChrom2List,chrom2), not_empty(LChrom2List,chrom2), not_empty(SChrom2List,chrom2)
                    combinedSorted=sorted(DStartList.values()+ LStartList.values()+ SStartList.values(), key=lambda x:(int(x[0]),int(x[1])))
                    Combined, done=joinCalls_trans(combinedSorted, chrom1, chrom2, orient,  not_empty(delly,find_oppOrient(orient)), not_empty(lumpy,find_oppOrient(orient)), not_empty(svaba,find_oppOrient(orient)), int(window), done, svtype)
                    MergedList=MergedList+list(Combined)
    return MergedList

def merging_All(dellyL, lumpyL, SvabaL, window, lengthChrom):
    MergedDels=combine_all("Dels",dellyL[0], lumpyL[0],SvabaL[0],window)
    MergedDups=combine_all("Dups",dellyL[1], lumpyL[1],SvabaL[1],window)
    MergedInvs=combine_all("Invs",dellyL[2], lumpyL[2],SvabaL[2],window)
    MergedTrans=combine_all("Trans",dellyL[3], lumpyL[3],SvabaL[3],window)
    MergedListAll=MergedDels+ MergedDups+ MergedInvs+ MergedTrans
    return MergedListAll

def read_len(file1):
    lengthsChrom=dict()
    with open(file1,'r') as f:
        for line in f:
            line_list=line.split('\t')
            lengthsChrom[line_list[0]]=int(line_list[1])
    return lengthsChrom

if __name__=='__main__':
    parser=argparse.ArgumentParser(description = "Python script to integrate SV calls from Delly, Lumpy and SvABA")
    parser.add_argument('-chromLengths','--lengthChrom',help="enter tab separated file with first column as chromosome name and second column as length of chromosomes (Required)")
    parser.add_argument('-delly','--delly',help="enter calls from delly in order:deletion,duplication,inversion,translocation (Required). \n Example: DEL.vcf,DUP.vcf,INV.vcf,TRA.vcf")
    parser.add_argument('-lumpy','--lumpy',help="enter calls from lumpy (Required)")
    parser.add_argument('-svaba','--svaba',help="enter calls from SvABA (Required)")
    parser.add_argument('-frag','--frag',help="enter size of window in which breakpoint would be integrated (default=500)", type=int, default=500)
    parser.add_argument('-outdir','--outdir',help="enter path of output directory (default=current directory)", default=".", type=str)
    parser.add_argument('-mode','--mode',help="germline or somatic (default=germline)", default="germline", type=str)
    parser.add_argument('-finalFile','--finalFile', help="Output file name (Required)")
    args=parser.parse_args()
    ### reading files ###
    lenChrom=read_len(args.lengthChrom)
    all_chromo=["chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY"]
    dellyR=read_delly(args.delly.split(",")[0],args.delly.split(",")[1],args.delly.split(",")[2],args.delly.split(",")[3], all_chromo)
    lumpyR=read_lumpy(args.lumpy, all_chromo)
    svabaR=read_svaba(args.svaba,args.mode, all_chromo)
    CombinedSV=merging_All(dellyR,lumpyR,svabaR,args.frag, lenChrom)
    extras=['chrom1','pos1','chrom2:pos2','SVType','Size', 'Orientation']
    columnsD=['DellyChrom1','DellyStart','DellyEnd','DellyChromo2','DellyOrientation','DellyFilter']
    columnsL=['LumpyChrom1','LumpyStart','LumpyEnd','LumpyChromo2','LumpyOrientation']
    columnsS=['SvabaChrom1','SvabaStart','SvabaEnd','SvabaChromo2','SvabaOrientation', 'SvabaFilter']
    columnsAll=extras+ columnsD+ columnsL + columnsS
    CombinedDF=pd.DataFrame(CombinedSV,columns= columnsAll)
    CombinedDF.to_csv(args.finalFile, sep=",", index=False)
