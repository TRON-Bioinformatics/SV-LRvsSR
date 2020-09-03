### script for requantification
import os
import random
import string
import pysam
import subprocess

def id_generator(size=10, chars=string.ascii_uppercase + string.digits):
    return ''.join(random.choice(chars) for _ in range(size))

def find_seq(chromo, pos1, pos2, strand, reference):
    result=subprocess.Popen(['python','./Scripts/twobit.py',str(reference),str(chromo),str(max(pos1,1)),str(pos2),str(strand)], stdout=subprocess.PIPE)
    out=result.communicate()[0].rstrip('\n')
    result.stdout.close()
    return out

def find_requant(WholeLine, refBit, areaRequant):
    chromA, pos1, chromB, pos2=WholeLine['chrom1'], WholeLine['pos1'], WholeLine['chrom2:pos2'].split(":")[0], WholeLine['chrom2:pos2'].split(":")[1]
    SVType, orient, size=WholeLine['SVType'], WholeLine['Orientation'], int(WholeLine['Size'])
    if SVType=="Dels":
        if size<= 300:
            areaRequant=150
        elif 300<size<= 1000:
            areaRequant=min(500,size)
        area1=find_seq(chromA, int(pos1)-int(areaRequant), int(pos1), "+", refBit)
        area2=find_seq(chromA, int(pos2), int(pos2)+int(areaRequant), "+", refBit)
    elif SVType=="Dups":
        area1=find_seq(chromA, int(pos2)-int(areaRequant), int(pos2), "+", refBit)
        area2=find_seq(chromA, int(pos1), int(pos1)+int(areaRequant), "+", refBit)
    elif SVType=="Invs":
        if orient=="3to3":
            area1=find_seq(chromA, int(pos1)-int(areaRequant), int(pos1), "+", refBit)
            area2=find_seq(chromA, int(pos2)-int(areaRequant), int(pos2), "-", refBit)
        elif orient=="5to5":
            area1=find_seq(chromA, int(pos1), int(pos1)+int(areaRequant), "-", refBit)
            area2=find_seq(chromA, int(pos2), int(pos2)+int(areaRequant), "+", refBit)
    elif SVType=="Trans":
        if orient=="3to5":
            area1=find_seq(chromA, int(pos1)-int(areaRequant), int(pos1), "+", refBit)
            area2=find_seq(chromB, int(pos2), int(pos2)+int(areaRequant), "+", refBit)
        elif orient=="5to3":
            area1=find_seq(chromB, int(pos2)-int(areaRequant), int(pos2),"+", refBit)
            area2=find_seq(chromA, int(pos1), int(pos1)+int(areaRequant),"+", refBit)
        elif orient=='3to3':
            area1=find_seq(chromA, int(pos1)-int(areaRequant), int(pos1), "+", refBit)
            area2=find_seq(chromB, int(pos2)-int(areaRequant), int(pos2), "-", refBit)
        elif orient=='5to5':
            area1=find_seq(chromA, int(pos1), int(pos1)+int(areaRequant), "-", refBit)
            area2=find_seq(chromB, int(pos2), int(pos2)+int(areaRequant), "+", refBit)
    else:
        area1='N'*499
        area2='N'*499
    combined=str(area1)+str(area2)
    return combined, areaRequant

def align(FaFile, Read1, Read2, outdir,tmpdir, nthreads , typeof):
    # create index of fasta file written
    os.chdir(tmpdir)
    p_index=subprocess.Popen(['bwa','index',FaFile])
    p_index.wait()
    with open(str(tmpdir)+"/aln_"+typeof+"_read1.sai", 'w') as alignmentR1:
        p_align1=subprocess.Popen(['bwa','aln','-t',str(int(nthreads)/2),str(FaFile),str(Read1)], stdout=alignmentR1)
    with open(str(tmpdir)+"/aln_"+typeof+"_read2.sai", 'w') as alignmentR2:
        p_align2=subprocess.Popen(['bwa','aln','-t',str(int(nthreads)/2),str(FaFile),str(Read2)], stdout=alignmentR2)
    p_align1.communicate(), p_align1.wait(), p_align2.communicate(), p_align2.wait()
    p_main=subprocess.Popen(['bwa','sampe',str(FaFile),str(tmpdir)+"/aln_"+typeof+"_read1.sai",str(tmpdir)+"/aln_"+typeof+"_read2.sai", str(Read1), str(Read2)], stdout=subprocess.PIPE)

    with open(str(tmpdir)+"/filtered_"+typeof+".bam",'w') as outstream:
        p_filter=subprocess.Popen(['samtools','view','-bS','-F','4'], stdin=p_main.stdout, stdout=outstream)
    p_filter.communicate(), p_filter.wait()
    print tmpdir, str(tmpdir)+"/filtered_"+typeof+".bam"
    p_sort=subprocess.Popen(['sambamba','sort','-t',str(nthreads),'--tmpdir=',str(tmpdir), str(tmpdir)+"/filtered_"+typeof+".bam"])
    p_sort.wait()

def findPercentMapping(cigar_read):
    totalLen=0.0
    matches=0.0
    for i in cigar_read:
        tmp=list(i)
        if tmp[0]==0:
            matches+=float(tmp[1])
            totalLen+=float(tmp[1])
        else:
            totalLen+=float(tmp[1])
    return float((matches/totalLen)*100)

def calculate_Junc_Span(FaFile, bamFile, junc_cutoff, outdir, requantbp, typeof):
    samfile=pysam.AlignmentFile(bamFile,'rb')
    ref_ids=samfile.references
    counts={}
    for ref_id in ref_ids:
        pos1=int(ref_id.split("_")[0].split(":")[1])
        pos2=int(ref_id.split("_")[1].split(":")[1])
        svtype=ref_id.split("_")[2]
        reads_on_ref_l = {}
        reads_on_ref_r = {}
        for read in samfile.fetch(ref_id):
            if read.flag > 255:
                continue
            query_name= read.query_name
            start= read.reference_start
            end= read.reference_end
            cigar=read.cigartuples
            if read.is_read1:
                reads_on_ref_l[query_name]= (start,end,cigar)
            else:
                reads_on_ref_r[query_name]= (start,end,cigar)
        bp_pos=int(requantbp[ref_id])
        junc_reads=0
        span_pairs=0
        for query_name in list(set(list(reads_on_ref_l.keys())+list(reads_on_ref_r.keys()))):
            if query_name in reads_on_ref_l:
                (start_r1, end_r1, cigar_r1)= reads_on_ref_l[query_name]
            else:
                (start_r1, end_r1, cigar_r1)=(None,None, None)
            if query_name in reads_on_ref_r:
                (start_r2, end_r2, cigar_r2)= reads_on_ref_r[query_name]
            else:
                (start_r2, end_r2, cigar_r2)= (None,None,None)
            if (end_r1!=None and end_r2 != None) and ((start_r1 < bp_pos and end_r1 < bp_pos and start_r2 > bp_pos and end_r2 > bp_pos) or (start_r1 > bp_pos and end_r1 > bp_pos and start_r2 < bp_pos and end_r2 < bp_pos)):
                Percent_r1=findPercentMapping(cigar_r1)
                Percent_r2=findPercentMapping(cigar_r2)
                if Percent_r1 >= 70 and Percent_r2 >= 70:
                    span_pairs +=1
            if (end_r1 != None) and ((start_r1 + junc_cutoff) <= bp_pos <= (end_r1 - junc_cutoff)):
                Percent_r1=findPercentMapping(cigar_r1)
                if Percent_r1 >= 70:
                    junc_reads +=1
            if (end_r2 != None) and ((start_r2 + junc_cutoff) <= bp_pos <= (end_r2 - junc_cutoff)):
                Percent_r2=findPercentMapping(cigar_r2)
                if Percent_r2 >= 70:
                    junc_reads +=1
        counts[ref_id]= [junc_reads, span_pairs]
    countFile=open(str(outdir)+"/Counts_"+typeof+".txt",'w')
    countFile.write("{}\t{}\t{}\n".format("ref_id", "junc_reads", "span_pairs"))
    for ref_id in ref_ids:
        if ref_id in counts:
            [junc_reads, span_pairs] = counts[ref_id]
        else:
            [junc_reads, span_pairs] = [0,0]
        countFile.write("{}\t{}\t{}\n".format(ref_id, junc_reads, span_pairs))
    countFile.close()
    return counts
