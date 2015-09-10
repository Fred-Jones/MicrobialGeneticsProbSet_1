from Bio import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqUtils import GC

i = Seq.Seq("TTCCCGCTCAGCCTGTCATGGCTGAGCGATTTTTTTTTTACGT"
            ,IUPAC.unambiguous_dna)
i.complement = i.complement()
i.complement.perc_GC = GC(i.complement)
i.__name__ = 'i'

ii = Seq.Seq("TGATCGGATTCGTAGAACGTAGTTAGTGTTTTTTTTTTTAGCT"
            ,IUPAC.unambiguous_dna)
ii.complement = ii.complement()
ii.complement.perc_GC = GC(ii.complement)
ii.__name__ = 'ii'

seqs = list([i,ii])
GCList = []
print "~*--------------------------*~"
for item in seqs:
    item.complement.__name__ = item.__name__
    GCList.append(item.complement)
    print "Percent GC of the complement of sequence {0}: {1}\n ".format(item.__name__, GC(item.complement))
maxGC = {"max":max(map(GC, GCList)), "__name__": None}
for ss in seqs:
    if GC(ss.complement) == maxGC["max"]:
        maxGC["__name__"] = ss.__name__
nnx = ["__name__", "max"]
print "Strand with max GC is {0} at {1}%".format(maxGC[nnx[0]], maxGC[nnx[1]])

print "~*--------------------------*~"
sig_H = "ATACCTTCTCACCGTTGAATACGCGCAGTTACGTTCCGGATTTACGGTCCGCATTTAGGCCTTTTTTATGGAAGAGTGGCAACTTATGCGCGTCAATGCAAGGCCTAAATGCCAGGCGTAAATCCGGAAAAA"
print "sig_H template", sig_H
TCTC = sig_H.find("TCTCACC")
print "Most likely Sig_H -35 region is at +", TCTC, " relative to 5' end of strand"

AUG = ["AUG", "GUA", str(Seq.Seq("AUG").complement()), str(Seq.Seq("AUG").reverse_complement())]
print "Any AUG codons? \n{0} \n{1}".format(AUG, map(sig_H.find, AUG))
print "AT pair aligned at {0} and TA pair aligned at {1}".format(sig_H[31:len(sig_H)+1].find("AT"),
                                                                 sig_H[31:len(sig_H)+1].find("TA"))

import re
book_seq_neg10 = "T?TA"
print "TNTA alignments: ", re.search(sig_H[31:len(sig_H)+1], book_seq_neg10)
print "~*--------------------------*~"
##Worksheet says 'this is a 'real' sequence'
##So I blast it
# from Bio.Blast import NCBIWWW, NCBIXML
# res_handle = NCBIWWW.qblast("blastn", "nt", sig_H)
# result = res_handle.read()
# save_f = open('result.xml', 'wr')
# xml_result, http_result = save_f.write(result), result
# blast_record = NCBIXML.BlastParser().parse(http_result)
#
# print http_result
# print blast_record
# #print blast_record.alignments[0]._def
# #print blast_record.alignments[1]._def

#Find longest ORF
STOPS = ["UAA", "UAG", "UGA"]
org_seq_0 = "AGGCCUACUCGAUCUAAGGGCUUCCCAGUCAGGUCAUAAAGCUAGGUAGGGC"
org_seq_1 = org_seq_0[1:]
org_seq_2 = org_seq_1[1:]
orf_0 = map(org_seq_0.find, STOPS)
orf_1 = map(org_seq_1.find, STOPS)
orf_2 = map(org_seq_2.find, STOPS)
# print orf_0, orf_1, orf_2
# print len(org_seq_0)

##Refactored code from above for my sanity
seq_x = [org_seq_0, org_seq_1, org_seq_2]
stops_in_seq_x = []
for seqx in seq_x:
    stops_in_seq_x.append(map(seqx.find, STOPS))
print "Original sequence", org_seq_0
print "Stops in open reading frames: ", stops_in_seq_x
print "Max: ", max(stops_in_seq_x)
print "Longest ORF: ", max(max(stops_in_seq_x))
