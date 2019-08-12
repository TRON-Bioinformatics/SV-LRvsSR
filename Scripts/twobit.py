import struct
import sys
import math
from exceptions import IndexError

class dna_2bit:
	num_entries = 0
	header = []
	offsets = {}
	starts = {}
	file = ""
	def __init__(self, filename):
    		self.file = filename
		self.header = []
		f = open(filename, "rb")
		s = f.read(16)
		header = struct.unpack("<LLLL", s)
		self.header.append("%X" % header[0])
		self.header.append(header[1:4])
		self.num_entries = header[2]
		self.offsets = {}
		for i in range(self.num_entries):
			s = f.read(1)
			name_len = int(struct.unpack("<B", s)[0])
			s = f.read(name_len)
			name_seq = struct.unpack("<"+str(name_len)+"s", s)
			s = f.read(4)
			offset = struct.unpack("<L", s)
			self.offsets[name_seq[0]] = offset[0]
		f.close()
		for i in self.offsets.keys():
			self.starts[i] = self.get_seq_start(i)
	
	def get_seq_start(self,name):
		f = open(self.file, "rb")
		f.seek(self.offsets[name], 0)
		s = f.read(4)
		dnaSize = struct.unpack("<L", s)[0]
		s = f.read(4)
		nBlockCount = struct.unpack("<L", s)[0]
		s = f.read(4*nBlockCount)
		nBlockStarts = struct.unpack("<"+str(nBlockCount)+"L", s)
		s = f.read(4*nBlockCount)
		nBlockSizes = struct.unpack("<"+str(nBlockCount)+"L", s)
		#f.seek(8*nBlockCount,1)
		s = f.read(4)
		maskBlockCount = struct.unpack("<L", s)[0]
		s = f.read(4*maskBlockCount)
		maskBlockStarts = struct.unpack("<"+str(maskBlockCount)+"L", s)
		s = f.read(4*maskBlockCount)
		maskBlockSizes = struct.unpack("<"+str(maskBlockCount)+"L", s)
		s = f.read(4)
		reserved = struct.unpack("<L", s)[0]
		#f.seek(8*maskBlockCount,1)
		bits_needed = dnaSize * 2
		bits_used = ((bits_needed/32)+1) * 32
		bit_padding = bits_used - bits_needed
		extra_nucleotides = bit_padding / 2
		packedDna_offset = f.tell()
		f.close()
		return (packedDna_offset,extra_nucleotides,dnaSize)

	def toBase2(self,number):
    		"""toBase2(number): given an int/long, converts it to
    		a string containing the number in base 2."""
		_nibbles = {"0":"0000", "1":"0001", "2":"0010", "3":"0011",
        	"4":"0100", "5":"0101", "6":"0110", "7":"0111",
            	"8":"1000", "9":"1001", "A":"1010", "B":"1011",
            	"C":"1100", "D":"1101", "E":"1110", "F":"1111",
            	"-":"-"}
		result = [_nibbles[nibble] for nibble in "%X"%number]
    		return "".join(result)

	def unsigned_byte_to_nucleotides(self,B,strand):
		r = ""
		n = ""
		if strand == "+":
			n = "TCAG"
		else:
			n = "AGTC"
		for i in [64,16,4,1]:
			r = r + (n[B/i])
			B = B % i
		return r
	
	def get_seq_range_bytes(self,name,start,stop):
		pos,padding,size = self.starts[name]
		if (stop >= size):
                        stop = size - 1
			#raise IndexError("Stop index too large: " + str(stop))
		if (start > stop):
			raise IndexError("Start index larger or equal to stop index: " + str(start) + ", " + str(stop))
		if (start < 0):
			raise IndexError("Indices must be >= zero: " + str(start))	
		start_byte = start/4
		end_byte = stop/4
		byte_len = end_byte - start_byte + 1
		f = open(self.file, "rb")
		f.seek(pos+start_byte, 0)
		s = f.read(byte_len)
		f.close()
		bytes = list(struct.unpack(">"+str(byte_len)+"B", s))
		return bytes

	def get_slice_from_byte_string(self,byte_str,start,stop,factor=1):
		start_corr = factor * start%4 
		end_corr = factor * ((stop%4)-3) 
		if end_corr == 0:
			byte_str = byte_str[start_corr:]
		else:
			byte_str = byte_str[start_corr:end_corr]
		return byte_str

	def get_seq_range_bits(self,name,start,stop):
		bytes = self.get_seq_range_bytes(name,start,stop)
		bits = []
		for i in bytes:
			bits.append(self.toBase2(i))
		bits = "".join(bits)
		bits = self.get_slice_from_byte_string(bits,start,stop,factor=2)
		return bits
        # gets the closed nucleotide interval, e.g. len(get_seq_range_nucleotides(chr1,10,20)) = 11
	def get_seq_range_nucleotides(self,name,start,stop,strand="+"):
		bytes = self.get_seq_range_bytes(name,start,stop)
		nucl = []
		for i in bytes:
			nucl.append(self.unsigned_byte_to_nucleotides(i,strand))
		nucl = "".join(nucl)
		nucl = self.get_slice_from_byte_string(nucl,start,stop)
		if strand != "+":
			nucl = list(nucl)
			nucl.reverse()
			nucl = "".join(nucl)
		return nucl
		
if __name__ == '__main__':
	if len(sys.argv) > 1:
		x = dna_2bit(sys.argv[1])
		print x.get_seq_range_nucleotides(sys.argv[2],int(sys.argv[3]),int(sys.argv[4]),sys.argv[5])
	else:
		x = dna_2bit("/var/www/gbdb/hg18/nib/hg18.2bit")
		print x.starts
	        print x.get_seq_range_bits("chr10",52000,52031)
		print x.get_seq_range_nucleotides("chr10",52000,52031)
		print x.get_seq_range_nucleotides("chr10",52000,52031,"-")
		print x.get_seq_range_nucleotides("chr10",52000,52011,"-")
		try:
			print x.get_seq_range_nucleotides("chr10",1,1000000000)
		except IndexError as e:
	     		print "Error: " + str(e)

		"""chrX    150831651       150844298       uc004fez.1      0       +       150842792       150843746       0       3       19,72,1571,     0,10929,11076,"""
		
		"""chrX    151086289       151370487       uc010ntk.1      0       -       151087355       151283698       0       10      1402,212,153,144,83,221,68,122,166,214,       0,22568,30471,40839,57601,88616,117506,178419,197269,283984,"""
		print "@exon_+_uc004fez.1"
		y=x.get_seq_range_nucleotides("chrX", 150842757,150842817)
		print y
		print "+\n" + "h"*len(y)

		print "@exon_-_uc010ntk.1"
		y=x.get_seq_range_nucleotides("chrX", 151086319,151086379,"-")
		print y
		print "+\n" + "h"*len(y)
		
		print "@junction_+"
		y=x.get_seq_range_nucleotides("chrX",  150832037,  150832096) + x.get_seq_range_nucleotides("chrX", 150842580, 150842621)
		print y
		print "+\n" + "h"*len(y)
	
		print "@junction_-"
		y=x.get_seq_range_nucleotides("chrX",151108857,151108887,"-")+x.get_seq_range_nucleotides("chrX",151087660,151087690,"-")  
		print y
		print "+\n" + "h"*len(y)
	
		print "@intron_+"
		y=x.get_seq_range_nucleotides("chrX", 150832297,150832357)
		print y
		print "+\n" + "h"*len(y)
	
		print "@intron_-"
		y=x.get_seq_range_nucleotides("chrX", 151087789,151087849,"-")
		print y
		print "+\n" + "h"*len(y)
	
		print "@extra_+"
		y=x.get_seq_range_nucleotides("chrX", 150831651-500,150831651-440)
		print y
		print "+\n" + "h"*len(y)
	
		print "@extra_-"
		y=x.get_seq_range_nucleotides("chrX", 151086289-1000,151086289-940,"-")
		print y
		print "+\n" + "h"*len(y)
	
		print "@nil"
		y="GATAAAGCAGGAATTACTACTGCTTGTTTACGAATTAAATCGAAGTGGACTGCTGGCGGAAAATGAGAAA"
		print y
		print "+\n" + "h"*len(y)

