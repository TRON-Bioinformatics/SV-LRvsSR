import operator
import math
import sys

# -> 34 min for 76M reads

# runtime m log n
# with  m = number of reads reads
# and   n = number of exons


class binary_node(object):
  """ Class representing a node in a binary tree. A node has a numeric value determining the contents on its left or right child. If a node is a leaf,, the value may be of any type """
  left = None
  right = None
  value = None
  is_root = True
  is_leaf = False
  def __init__(self,left,right,value,root,leaf):
    self.left = left
    self.right = right
    self.value = value
    self.is_root = root
    self.is_leaf = leaf
  def search(self,value):
    """ Searches this node for value and returns the corresponding child or the own value if this node is a leaf """
    if self.is_leaf:
      return self.value
    if value >= self.value:
      return self.right.search(value)
    else:
      return self.left.search(value)
  def dfs(self,l):
    if self.is_leaf:
      l.append(self)
      return
    else:
      self.left.dfs(l)
      self.right.dfs(l)

class bed_tree(object):
  """ This class represents a binary tree build out of a bed file. The input bed file is expected to have the full 11 columns. The read_length argument defines the minimum distance between the interval sets in the leafs of the search tree. For internal use, it is multiplied by three."""
  def __init__(self,bedfile,read_length,force_small_bed_line=False,mls=50,coordinate_modifier=0):
    t,e,nt,self.ne = read_reference_transcripts(bedfile,force_small_bed_line,coordinate_modifier)
    self.treed = {}
    self.n_leafs = 0
    self.read_length = read_length
    self.t = t
    self.e = e
    for chr in e.keys():
       self.treed[chr],n = construct_tree(e[chr],mls,read_length*3)
       self.n_leafs += n
  def _calc_overlaps(self,l,r_start,r_stop):
    l_overlap = []
    for element in l:
      e_start,e_stop,e_info = element
      if (r_start >= e_start and r_start <= e_stop) or (r_stop >= e_start and r_stop <= e_stop):
        x = sorted([r_start,r_stop,e_start,e_stop])
        d = x[2] - x[1] + 1
        transcript,exon_nr,strand = e_info
        # might check also strand here
        l_overlap.append((transcript,exon_nr,strand,d,x))
    return l_overlap
  def overlap(self,chr,r_start,r_length,r_strand):
    """ Searches the bed_tree for all intervals (=exons) which overlap with the input interval (=read, defined by chromosome, 0-based start and length; the strand argument is not used currently). The method returns a list of 4-tuples. Each tuple represents an exon with the values transcript-id, 0-based exon number in transcript (numbered in 5'-3' direction of mRNA), strand of transcript and number of overlapping bases between exon and read. """
    r_stop = r_start + r_length - 1
    l = self.treed[chr].search(r_stop)
    return self._calc_overlaps(l,r_start,r_stop)
  def test(self):
    s = 0
    s2 = 0
    for k in self.treed.keys():
      l = []
      self.treed[k].dfs(l)
      ll = []
      for i in l:
        starts,stops,data = zip(*i.value)
        ll.append((min(starts),max(stops)))
        s += len(i.value)
      s2 += len(l)
      for i in range(1,len(ll)):
        if ll[i][0] - ll[i-1][1] < (self.read_length*3):
          print "[TEST] error",ll[i],ll[i-1]
    print >> sys.stderr, "[TEST] counted number of exons:", s
    print >> sys.stderr, "[TEST] counted number of leaves:", s2

#def main_repeat(bedfile, chromo, pos):
#  return bt.overlap(chromo,int(pos), read_length, "+")


def construct_tree(interval_list, max_leaf_size, min_difference):
  """ Entry function for recursive algorithm which creates a binary tree. interval_list is a list of triples with the elemets interval start (0-based), interval end (0-based) and some data of any type. The interval represented by the start and stop is expected to be a closed interval. THe interval_list is expected to be sorted by the interval starts. The algorithm guarantees that a maximum of intervals is at one leaf of the tree (as given by argument max_leaf_size) UNLESS any interval in any other node would overlap with an interval of the current leaf OR the distance between neigbouring leafs would be smaller than min_difference. These side conditions gurantee that a search will return all overlapping elements (among some non-overlapping candiates). This function returns the root of the tree, an object of class binary_node. """
  return _construct_tree(interval_list,max_leaf_size,min_difference,0)

def _construct_tree(interval_list, max_leaf_size, min_difference, dp):
  """ Recursive function for binary tree generation. """
  is_root = (dp == 0)
  if len(interval_list) <= max_leaf_size:
    return binary_node(None,None,interval_list,is_root,True),1
  starts,stops,data = zip(*interval_list)
  split = int(math.floor(len(interval_list)/2))
  d = 1; sign = -1

  while True:
    if stops[split-1] >= (starts[split] - min_difference):
      split = split + (d * sign)
      d += 1; sign = sign * (-1)
    else:
      break
    if (split < 1) or (split == len(interval_list)):
      #print >> sys.stderr, "no disjunct set possible, adding leaf with size {0}".format(len(interval_list))
      return binary_node(None,None,interval_list,is_root,True),1
  #if (stops[split-1] - starts[0] < min_difference) or (stops[-1] - starts[split] < min_difference):
  if (min(starts[split:]) - max(stops[:split]) < min_difference):
    #print >> sys.stderr, "split produces small distance between leafs, adding leaf with size {0}".format(len(interval_list))
    #print >> sys.stderr,min(starts[split:]), max(stops[:split]), min(starts[split:]) - max(stops[:split])
    #print >> sys.stderr, interval_list[:10]
    #sys.exit()
    return binary_node(None,None,interval_list,is_root,True),1

  l,n_left = _construct_tree(interval_list[:split],max_leaf_size,min_difference,dp + 1)
  r,n_right = _construct_tree(interval_list[split:],max_leaf_size,min_difference,dp + 1)
  m = binary_node(l,r,starts[split],is_root,False)
  return m, n_left + n_right



##################################################### BED file handling
def __convert_bed_line_to_exons(words):
  """Returns a list of exons (chr,start,stop,transcript name,exon number,strand). Input is a list representing an eleven field bed line. The start and stop coordinates are returned as zero base closed interval."""
  exon_list = []
  (chrom,start,end,name) = words[0:4]
  start = int(start)
  strand = words[5]
  exons = int(words[9])
  exon_starts = words[11].strip(",").split(",")
  exon_sizes = words[10].strip(",").split(",")
  exon_startsI = [start + int(x) for x in exon_starts]
  exon_sizesI = [int(x) for x in exon_sizes]
  r = range(exons)
  if strand == "-":
    r = range(exons-1,-1,-1)
  for ne in r:
    abs_e_end = exon_startsI[ne] + exon_sizesI[ne] - 1
    exon_number_correct_order = (exons - 1 - ne) if strand == "-" else ne
    exon_words = [chrom,exon_startsI[ne],abs_e_end,name,exon_number_correct_order,strand]
    #exon_words = map(str,exon_words)
    exon_list.append(exon_words)
  return exon_list

def __convert_small_bed_line_to_exons(words):
  """Returns a list with a single element, which represents a gene/transcript in a seven field bed file. The start and stop coordinates are returned as zero base closed interval."""
  (chrom,start,end,name) = words[0:4]
  strand = words[5]
  w = [chrom,int(start),int(end),name,-1,strand]
  return [w]

def read_reference_transcripts(bedfile,force_small_bed_line=False,coordinate_modifier=0):
  exons = {}
  transcripts = {}
  numt = 0

  # read bed file
  with open(bedfile,"r") as f:
    for line in f:
      words = line[:-1].split("\t")
      # store bed entry
      l = transcripts.get(words[3],[])
      l.append(words)
      transcripts[words[3]] = l
      # convert to exons
      if (len(words) == 12) and (not force_small_bed_line):
        es = __convert_bed_line_to_exons(words)
      else:
        es = __convert_small_bed_line_to_exons(words)
      for e in es:
        l = exons.get(e[0],[])
        l.append((e[1],e[2],(e[3],e[4],e[5])))
        exons[words[0]] = l
      numt += 1


  nume = 0
  for i in exons.keys():
      if coordinate_modifier != 0:
        for j in range(len(exons[i])):
          w = list(exons[i][j])
          w[0] -= coordinate_modifier
          w[1] += coordinate_modifier
          if w[0] < 0:
            w[0] == 0
          exons[i][j] = tuple(w)
      exons[i].sort(key=operator.itemgetter(1))
      exons[i].sort(key=operator.itemgetter(0))
      nume += len(exons[i])

  return transcripts,exons,numt,nume
