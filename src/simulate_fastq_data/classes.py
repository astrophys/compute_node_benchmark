import sys
import re
import datetime
import random
import copy
import operator
import numpy
#********************************************************************************
#********************************************************************************
#*******************************  CLASSES  **************************************
#********************************************************************************
#********************************************************************************
class GTF_ENTRY:
    """
    Container to hold each GTF entry. 

    See Example : 
    /reference/homo_sapiens/GRCh38/ensembl/Annotation/Genes/gtf/Homo_sapiens.GRCh38.83.gtf

    See Specification:
    http://useast.ensembl.org/info/website/upload/gff.html
    http://www.gencodegenes.org/gencodeformat.html
    """
    def __init__(self, Chromosome = None, Source = None, EntryType= None, Start = None,
                 Stop = None, Score = None, Strand = None, Frame = None, Attribute = None):
        """
        ARGS:
            Chromosome  : Chromosome, can only be 1,2,...,X,Y
            Source      : Database or Project name for entry
            EntryType   : exon, transcript, CDS, gene, etc
            Start       : start position on chromosome
            Stop        : stop position on chromosome
            Score       : UNKNOWN
            Strand      : '+' (forward) or '-' (reverse)
            Frame       : 0 indexed position of first base of codon
            Attribute   : semicolon sep list of tag-value pairs
            
        RETURN:
            NONE : Initializes GTF_ENTRY

        DESCRIPTION:
    
        DEBUG:
            Tested by reading in gtf file and printing out list of GTF_ENTRYs.
            compared using full Homo_sapiens.GRCh38.83.gtf and the output 
            was _identical_.

        FUTURE: 
        """
        self.chrm    = None        # str, Chromosome
        self.src     = None        # str, Source of data
        self.etype   = None        # str, exon, transcript, CDS, gene
        self.start   = None        # int, Start position on chrm
        self.stop    = None        # int, End position on chrm
        self.score   = None        # str, expects empty field, ie : '.'
        self.strand  = None        # str, '+' (forward) or '-' (reverse)
        self.frame   = None        # str, 0 indexed position of first base of codon
        self.attribute = None      # str, semicolon sep list of tag-value pairs
        # below vars are parsed from attribute list
        self.geneID  = None        # str, Gene ID ... starts with ENSG..
        self.geneName= None        # str, common gene name
        self.transID = None        # str, transcript ID ... starts with ENST...
        self.transName= None       # str, transcript Name.. = gene common name + num
        self.exonID  = None        # str, exon ID ... starts with ENSE
        self.exonNum = None        # int, exon number
        self.biotype = None        # str, gene_biotype

        if(Chromosome is not None):
            self.chrm = Chromosome
        else:
            sys.stderr.write("ERROR! Chromosome not specified!\n")
            sys.exit(1)
        if(Source is not None):
            self.src = Source
        else:
            sys.stderr.write("ERROR! Source not specified!\n")
            sys.exit(1)
        if(EntryType is not None):
            self.etype = EntryType
        else:
            sys.stderr.write("ERROR! EntryType not specified!\n")
            sys.exit(1)
        if(Start is not None):
            self.start = int(Start)
        else:
            sys.stderr.write("ERROR! Start not specified!\n")
            sys.exit(1)
        if(Stop is not None):
            self.stop = int(Stop)
        else:
            sys.stderr.write("ERROR! Stop not specified!\n")
            sys.exit(1)
        if(Score is not None):
            self.score = Score
            # Check for empty field
            if(self.score != '.'):
                sys.stderr.write("ERROR! Score = %s. Expected empty field\n"%(self.score))
                sys.exit(1)
        else:
            sys.stderr.write("ERROR! Score not specified!\n")
            sys.exit(1)
        if(Strand is not None):
            self.strand = Strand
            if(self.strand != '+' and self.strand != '-'):
                sys.stderr.write("ERROR! Strand = %s, wrong format\n"%(self.strand))
                sys.exit(1)
        else:
            sys.stderr.write("ERROR! Strand not specified!\n")
            sys.exit(1)
        if(Frame is not None):
            self.frame = Frame
        else:
            sys.stderr.write("ERROR! Frame not specified!\n")
            sys.exit(1)
        if(Attribute is not None):
            self.attribute = Attribute.split("\n")[0]

            # Now parse attribute list for relevant fields. Create dictionary for easy lookup
            attributeDict = {}
            attributeSplit = self.attribute.split(";")[:-1]     # lop off last due to ; at end
            for attrEntry in attributeSplit:
                attrEntry = attrEntry.strip()                   # Eliminate white space on ends
                key   = attrEntry.split(" ")[0]
                value = attrEntry.split(" ")[1]
                attributeDict[key] = value
            # Every key is not necessarily in the attributeDict, handle accordingly
            try:  self.geneID   = attributeDict["gene_id"]
            except KeyError:  pass
            try:  self.geneName = attributeDict["gene_name"]
            except KeyError:  pass
            try:  self.transID  = attributeDict["transcript_id"]
            except KeyError:  pass
            try:  self.transName= attributeDict["transcript_name"]
            except KeyError:  pass
            try:  self.exonID   = attributeDict["exon_id"]
            except KeyError:  pass
            try:  self.exonNum  = int((attributeDict["exon_number"])[1:-1])
            except KeyError:  pass
            try:  self.biotype  = attributeDict["gene_biotype"]
            except KeyError:  pass
        else:
            sys.stderr.write("ERROR! Attribute not specified!\n")
            sys.exit(1)


    def print_entry(self):
        print("%s\t%s\t%s\t%i\t%i\t%s\t%s\t%s\t%s"%(self.chrm, self.src, self.etype,
              self.start, self.stop, self.score, self.strand, self.frame,
              self.attribute))


class FASTQ_READ:
    def __init__(self, Transcript = None, Insert = None, ReadLength = None, MetaData = None, exonList = None, R1_R2 = None):
        """
        ARGS:
            Transcript = a TRANSCRIPT instance
            Insert = an INSERT instance
            ReadLength = length of desired read.
        RETURN:
            NONE : Initializes GENE

        DESCRIPTION:
            Currently, synthetic reads are forced to be completely contained within
            the transcript. i.e, I am not permitting them to read into the
            adapters, we'll save that for future work

        DEBUG:
            I've only done some spot debugging. I _really_ _really_ need to 
            invest thoroughly debugging this

            To ensure that I am placing reads where I think 
            they should be and getting the number of reads.

            Using ENST00000350501, I printed 10000 reads and found that it 
            NEVER gave me read positions such that the read would step out of
            the Transcript.seq bounds. 

            See debug section for create_fastq_file(). CONCLUSION is that
            the start positions w/r/t genomic coordinates is correct. Also 
            this means that the exonIDs written to the header are also correct.

        FUTURE: 
        """
        self.seq  = None
        self.qual = None
        self.readlen = None
        self.metadata = None
        self.read = R1_R2           # whether the fastq read is read 1 or read 2
        exonsSpanned = []           # List of exons spanned by read
        startWrtGenome = 0          # coord relative to genome / chromosome start
        stopWrtGenome = 0           # coord relative to genome / chromosome start
        #exonCoordsWrtTrans = {}    # dict of lists : relative exon coords w/r/t transcript start
        curTransPos  = 0            # current transcipt coordinate position.

        #insertList = []             # List of inserts from this transcript

        # type check 
        if(not isinstance(Insert, INSERT)):
            sys.stderr.write("ERROR! Insert is not of class type INSERT\n")
            sys.exit(1) 

        if(ReadLength is not None):
            # ReadLength = int(ReadLength)
            self.readlen = ReadLength
        else:
            sys.stderr.write("ERROR! ReadLength not specified!\n")
            sys.exit(1)

        self.get_qual()
        # Currently, reads are forced to be completely contained within
        # the transcript. i.e, I am not permitting them to read into the
        # adapters, we'll save that for future work

        # ALSO : recall that Transcript.seq is always in the direction that
        # transcripts are transcribed.

        if (self.read == "R1"):
            self.seq = Insert.seq[0:ReadLength]
        else:
            seqLen = len (Insert.seq)
            self.seq = reverse_complement(Insert.seq[seqLen-ReadLength:seqLen])


        # Add useful information to reads
        # shouldn't be so complicated. will neet to revise Transcript.exonList to simplify this part
        if(MetaData is not None):
            exonsSpanned = []
            curTransPos  = 0                                 # Find exon start pos w/r/t trans start

            # Get read start position w/r/t genome.
            # Also, get exons spanned by read.

            #sys.stderr.write("Transcript:\n%s\n\n"%Transcript.transID)
            #sys.stderr.write("Insert start %i\t stop:%i\n"%(Insert.t_start,Insert.t_stop))


            for exonIdx in Transcript.exonIdxList:
                # below transcript coordinates are 0 indexed, and the coordinates are _inclusive_
                exon = exonList[exonIdx]                # get actual EXON instance
                exonLen = len(exon.seq)
                # exonCoordsWrtTrans[exonIdx] = [curTransPos, curTransPos + exonLen - 1]

                #sys.stderr.write("exonIdx:\t%i\tstart: %i\tstop: %i\texonLen: %i\t"%(exonIdx,exon.start,exon.stop,exonLen))
    

                # The start of the insert is to the right of the current exon pointer
                #   Note the pointer always point to the start position of the exon
                if(Insert.t_start >= curTransPos):
                    # The start of the insert is to the right of the current exon, 
                    #   so break the loop and  check next exon
                    if (Insert.t_start > curTransPos + exonLen - 1):
                        curTransPos += exonLen
                        #sys.stderr.write("\tcurTransPos: %i\n"%curTransPos)
                        continue

                    # The start of the insert falls into this exon 
                    #   then we define the position of the insert start w/r/t genome
                    #   exon.start and exon.stop is 1 indexed
                    if (Transcript.strand == '+'):
                        startWrtGenome = exon.start + (Insert.t_start - curTransPos)
                    else:
                        startWrtGenome = exon.stop - (Insert.t_start - curTransPos)
                    #sys.stderr.write("startWrtGenome defined!\t")

                # The start of the insert is to the left of the current exon
                #   Now we check the stop position of the insert
                #   If the stop of the insert is in this exon, then we define the position of t
                #      the insert stop w/r/t the genome and break the loop
                if (Insert.t_stop <= curTransPos + exonLen - 1):
                    if(Transcript.strand == '+'):
                        stopWrtGenome = exon.start + (Insert.t_stop - curTransPos)
                    else:
                        stopWrtGenome = exon.stop - (Insert.t_start - curTransPos)
                    #sys.stderr.write("stopWrtGenome defined!\t")
                    
                    exonsSpanned.append(exon.exonID)
                    #sys.stderr.write("\texonsSpanned: %i\tcurTransPos: %i\n"%(len(exonsSpanned),curTransPos))
                    # break the loop
                    break
                else: # The stop of the insert is still on the horizon, increment the variables and 
                      # continue
                    exonsSpanned.append(exon.exonID)
                    curTransPos += exonLen

                #sys.stderr.write("\texonsSpanned: %i\tcurTransPos: %i\n"%(len(exonsSpanned),curTransPos))

            if(len(exonsSpanned) == 0):
                sys.stderr.write("ERROR! Read does _not_ span any exons!\n")
                sys.exit(1)

            exonsSpanned = list(set(exonsSpanned))
            self.metadata = "%s:trans:%s:start:%i:exons"%(MetaData, Transcript.transID, startWrtGenome)
            for exonName in exonsSpanned:
                self.metadata = "%s:%s"%(self.metadata, exonName)
        else:
            sys.stderr.write("ERROR! MetaData not specified!\n")
            sys.exit(1)




    def get_qual(self):
        """
        ARGS:
            NONE 
        RETURN:
            quality scores based on Illumina, Phread+33 between chars '!' (33)
            and 'J' (74), inclusive

        DESCRIPTION:
            Initializes quality information

        DEBUG:
            Should assign only G's. It does.

        FUTURE: 
            Perhaps add some randomness to this later or distribution to this later.
        """
        qual = ""
        for i in range(self.readlen):
            qual += 'G'
        self.qual = qual




class CHROMOSOME:
    def __init__(self, Chromosome = None, Sequence = None):
        """
        ARGS:
            Chromosome = chromosome string
            Sequence   = chromosome sequence

        RETURN:
            NONE : Initializes CHROMOSOME

        DESCRIPTION:
            Gets chromosome sequence so that we can use it to get the exon
            sequences.

        DEBUG:
            For a shortened whole genome fasta file, it correctly reads it in.
            See read_genome() for debugging code and futher details.

        FUTURE: 
        """
        self.chrm= None         # str, Chromosome
        self.seq     = []           # str, sequence

        if(Chromosome is not None):
            self.chrm = Chromosome
        else:
            sys.stderr.write("ERROR! Chromosome is not specified!\n")
            sys.exit(1)
        if(Sequence is not None):
            self.seq = Sequence
        else:
            sys.stderr.write("ERROR! Sequence is not specified!\n")
            sys.exit(1)




class GENE:
    def __init__(self, GtfEntry):
        """
        ARGS:
            GtfEntry = a single GTF_ENTRY class element

        RETURN:
            NONE : Initializes GENE

        DESCRIPTION:

        DEBUG:

        FUTURE: 
        """
        self.chrm   = None        # str, Chromosome
        self.start  = None        # int, Start position on chrm
        self.stop   = None        # int, End position on chrm
        self.strand = None        # str, '+' (forward) or '-' (reverse)
        self.geneID = None        # str, nominal geneID, but may belong to multiple genes
        self.geneName= None        # str, common gene name
        self.transList    = []    # transcript names that are part of this gene
        self.transIdxList = []    # list, use to quickly map to AllTransList
        transIdx    = 0
        
        # type check
        if(not isinstance(GtfEntry, GTF_ENTRY)):
            sys.stderr.write("ERROR! GtfEntry is not of class type GTF_ENTRY\n")
            sys.exit(1)

        if(GtfEntry.chrm is not None):
            self.chrm = GtfEntry.chrm
        else:
            sys.stderr.write("ERROR! GtfEntry.chrm is None\n")
            sys.exit(1)
        if(GtfEntry.start is not None):
            self.start = int(GtfEntry.start)
        else:
            sys.stderr.write("ERROR! GtfEntry.start is None\n")
            sys.exit(1)
        if(GtfEntry.stop is not None):
            self.stop = int(GtfEntry.stop)
        else:
            sys.stderr.write("ERROR! GtfEntry.stop is None\n")
            sys.exit(1)
        if(GtfEntry.strand is not None):
            self.strand = GtfEntry.strand
        else:
            sys.stderr.write("ERROR! GtfEntry.strand is None\n")
            sys.exit(1)
        if(GtfEntry.geneID is not None):
            self.geneID = GtfEntry.geneID
        else:
            sys.stderr.write("ERROR! GtfEntry.geneID is None\n")
            sys.exit(1)
        if(GtfEntry.geneName is not None):
            self.geneName = GtfEntry.geneName
        else:
            sys.stderr.write("ERROR! GtfEntry.geneName is None\n")
            sys.exit(1)

        # Get all the exons and their indices associated with this transcript
        # This is _HIGHLY_ inefficient. Scales as O(N**2)
        #for trans in AllTransList:
        #    if(trans.geneID == self.geneID):
        #        self.transList.append(trans.transID)
        #        self.transIdxList.append(transIdx)
        #    transIdx = transIdx + 1


class TRANSCRIPT:
    def __init__(self, GtfEntry):
        """
        ARGS:
            GtfEntry = a single GTF_ENTRY class element

        RETURN:
            NONE : Initializes TRANSCRIPT

        DESCRIPTION:

        FUTURE: 
        """
        self.seq    = None        # str, sequence
        self.chrm   = None        # str, Chromosome
        self.start  = None        # int, Start position on chrm
        self.stop   = None        # int, End position on chrm
        self.strand = None        # str, '+' (forward) or '-' (reverse)
        self.geneID = None        # str, nominal geneID, but may belong to multiple genes
        self.transID = None       # str, nominal transcript ID, may belong to mult. trans.
        self.transNum= None       # str, only number portion of transID for sorting by transcript num
        self.exonList= []         # exon names that are part of this transcript
        self.exonIdxList = []     # list, use to quickly map to AllExonList
        exonIdx     = 0

        self.copy = 1               # copy number of this particular transcript in the transcriptome
                                    # It is set to 1 at the moment.

        # type check
        if(not isinstance(GtfEntry, GTF_ENTRY)):
            sys.stderr.write("ERROR! GtfEntry is not of class type GTF_ENTRY\n")
            sys.exit(1)

        if(GtfEntry.chrm is not None):
            self.chrm = GtfEntry.chrm
        else:
            sys.stderr.write("ERROR! GtfEntry.chrm is None\n")
            sys.exit(1)
        if(GtfEntry.start is not None):
            self.start = int(GtfEntry.start)
        else:
            sys.stderr.write("ERROR! GtfEntry.start is None\n")
            sys.exit(1)
        if(GtfEntry.stop is not None):
            self.stop = int(GtfEntry.stop)
        else:
            sys.stderr.write("ERROR! GtfEntry.stop is None\n")
            sys.exit(1)
        if(GtfEntry.strand is not None):
            self.strand = GtfEntry.strand
        else:
            sys.stderr.write("ERROR! GtfEntry.strand is None\n")
            sys.exit(1)
        if(GtfEntry.geneID is not None):
            self.geneID = GtfEntry.geneID
        else:
            sys.stderr.write("ERROR! GtfEntry.geneID is None\n")
            sys.exit(1)
        if(GtfEntry.transID is not None):
            self.transID = GtfEntry.transID
            self.transNum = self.transID[4:]
        else     :
            sys.stderr.write("ERROR! GtfEntry.transcriptID is None\n")
            sys.exit(1)


class INSERT:
    def __init__(self, Transcript = None, start = 0, stop = 0):
        """
        ARGS:
            Transcript = a TRANSCRIPT instance
            start = the start position on chrm
            stop = the stop position on chrm
            sequence = the sequence of the insert

        RETURN:
            NONE : Initializes INSERT

        DESCRIPTION:

        FUTURE: 
        """
        self.seq    = None        # str, sequence
        self.chrm   = None        # str, Chromosome
        self.t_start  = start        # int, Start position on the transcript (starts from 0)
        self.t_stop   = stop        # int, End position on the transcript
        self.start = None       # int, start position w/r/t the chromosome
        self.stop = None        # int, start position w/r/t the chromosome
        self.strand = None        # str, '+' (forward) or '-' (reverse)
        self.geneID = None        # str, nominal geneID, but may belong to multiple genes
        self.transID = None       # str, nominal transcript ID, may belong to mult. trans.
        self.transNum= None       # str, only number portion of transID for sorting by transcript num
       
        #self.insertIdx = idx       # int, '0' indexed insert position w/r/t other inserts


        # type check
        if(not isinstance(Transcript, TRANSCRIPT)):
            sys.stderr.write("ERROR! Transcript is not of class type TRANSCRIPT 2\n")
            sys.exit(1)

        if(Transcript.chrm is not None):
            self.chrm = Transcript.chrm
        else:
            sys.stderr.write("ERROR! Transcript.chrm is None\n")
            sys.exit(1)
        if(Transcript.strand is not None):
            self.strand = Transcript.strand
        else:
            sys.stderr.write("ERROR! Transcript.strand is None\n")
            sys.exit(1)
        if(Transcript.geneID is not None):
            self.geneID = Transcript.geneID
        else:
            sys.stderr.write("ERROR! Transcript.geneID is None\n")
            sys.exit(1)
        if(Transcript.transID is not None):
            self.transID = Transcript.transID
            self.transNum = self.transID[4:]
        else:
            sys.stderr.write("ERROR! Transcript.transcriptID is None\n")
            sys.exit(1)

        # get the sequence of the insert
        self.seq = Transcript.seq[start:stop]

        # get the start and stop position relative to the chromosome


class EXON:
    def __init__(self, GtfEntry):
        """
        ARGS:
            GtfEntry = a single GTF_ENTRY class element

        RETURN:
            NONE : Initializes EXON 

        DESCRIPTION:

        DEBUG:

        FUTURE: 
        """
        self.seq    = None        # str, sequence
        self.chrm   = None        # str, Chromosome
        self.start  = None        # int, Start position on chrm
        self.stop   = None        # int, End position on chrm
        self.strand = None        # str, '+' (forward) or '-' (reverse)
        self.exonNum= None        # int, '1' indexed exon position w/r/t other exons
        self.exonID = None        # str, exon ID ... starts with ENSE
        self.geneID = None        # str, nominal geneID, but may belong to multiple genes
        self.transID = None  # str, nominal transcript ID, may belong to mult. trans.

        # type check
        if(type(GtfEntry) is GTF_ENTRY):
            sys.stderr.write("ERROR! GtfEntry is not of class type GTF_ENTRY\n")
            sys.exit(1)
            
        if(GtfEntry.chrm is not None):
            self.chrm = GtfEntry.chrm
        else:
            sys.stderr.write("ERROR! GtfEntry.chrm is None\n")
            sys.exit(1)
        if(GtfEntry.start is not None):
            self.start = int(GtfEntry.start)
        else:
            sys.stderr.write("ERROR! GtfEntry.start is None\n")
            sys.exit(1)
        if(GtfEntry.stop is not None):
            self.stop = int(GtfEntry.stop)
        else:
            sys.stderr.write("ERROR! GtfEntry.stop is None\n")
            sys.exit(1)
        if(GtfEntry.strand is not None):
            self.strand = GtfEntry.strand
        else:
            sys.stderr.write("ERROR! GtfEntry.strand is None\n")
            sys.exit(1)
        if(GtfEntry.exonNum is not None):
            self.exonNum = GtfEntry.exonNum
        else:
            sys.stderr.write("ERROR! GtfEntry.exonNum is None\n")
            sys.exit(1)
        if(GtfEntry.exonID is not None):
            self.exonID = GtfEntry.exonID
        else:
            sys.stderr.write("ERROR! GtfEntry.exonID is None\n")
            sys.exit(1)
        if(GtfEntry.geneID is not None):
            self.geneID = GtfEntry.geneID
        else:
            sys.stderr.write("ERROR! GtfEntry.geneID is None\n")
            sys.exit(1)
        if(GtfEntry.transID is not None):
            self.transID = GtfEntry.transID
        else:
            sys.stderr.write("ERROR! GtfEntry.transcriptID is None\n")
            sys.exit(1)
            


def reverse_complement(seq):
    """
    ARGS:
        seq : sequence with _only_ A, T, C or G (case sensitive)

    RETURN:
        rcSeq : reverse complement of sequenced passed to it.

    DESCRIPTION:

    DEBUG: 
        Compared several sequences.  Is working.

    FUTURE: 
    """ 
    rcSeq = ""            # Reverse Complement sequence
    # Complement
    for char in seq:
        if(char == 'A' ):
            rcSeq += 'T' 
            continue
        if(char == 'T' ):
            rcSeq += 'A' 
            continue
        if(char == 'G' ):
            rcSeq += 'C' 
            continue
        if(char == 'C' ):
            rcSeq += 'G' 
            continue
        if(char == 'N' ):
            rcSeq += 'N' 
            continue

        if(char not in "ATCGN"):
            sys.stderr.write("ERROR! char %s is not a valid sequencing character!\n"%(char))
            sys.exit(1)
    # Revese
    rcSeq = rcSeq[::-1]
    return rcSeq

            



