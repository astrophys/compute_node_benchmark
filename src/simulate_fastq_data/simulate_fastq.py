#!/gpfs0/export/opt/anaconda-2.3.0/bin/python
# 
# Author  : Ali Snedden
# Date    : 5/27/16
# License : MIT
# 
# PURPOSE :
#   Read in a GTF file and create a set of simulated fastq data that 
#   can be aligned by tophat to GRCH38. The purpose is to be able to 
#   control directly the relative ratio of different transcripts 
#   (ie isoforms) expressed. This will give us a way to test 
#   MISO, rMATS and other such programs.
#   
# VERSIONS:
#
#
#
#
# TESTING:
#   1. GTF is read in correctly and parsed by GTF_ENTRY.__init__
#      -> Tested by diffing output w/ original gtf file
#   2. Scales like crap ...O(N**2). This needs sped up b/c it 
#      cannot complete the whole genome and gtf in 12+ hours, which is 
#      unacceptable. Get's bogged down in finding exons for transcripts.
#      --> in link_exons_trans_and_genes() I fixed this issue.
#          After reading in exons, transcripts and genes, it does one final
#          pass through the list of GTF_ENTRYs and _assuming_ that 
#          all the transcripts are directly _after_ the gene but _before_
#          the next gene. Same goes for exons between transcripts.
#
# FUTURE: 
#   1. If I was a clever object oriented programmer, I'd appropriately use
#      inheritence in the parts of the classes, GENE, EXON and TRANSCRIPT
#
#   2. Improve algorithm so that I don't need to loop through Gtf 3 times
#      to read in all the exons, transcripts and genes.
#
#   3. Handle CDS, stop_codon, etc. and other features.  Currently I only
#      handle genes, transcripts and exons.
#
#   4. Add strandedness options. fr-firststrand, fr-secondstrand, fr-unstraned
#      See https://dbrg77.wordpress.com/2015/03/20/library-type-option-in-the-tuxedo-suite/
#      for discussion.
#      Use IGV for seeing if I got it correct.
#
#   5. Paired end reads.
#
#   6. Fix bug, for instance with 100 nt reads starting at 930330 for trans
#      ENST00000616125, it spans 3 exons, but I only report 2
#
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
            


            


#********************************************************************************
#********************************************************************************
#******************************  FUNCTIONS **************************************
#********************************************************************************
#********************************************************************************
def print_help():
    """This prints the Help"""
    print("\nUSAGE: ./simulate_fastq.py pathToGtf pathToSequenc pathToConfig numReads pathToFastq\n"
            "pathToGtf      = Path to the _ensembl_ GTF file\n"
            "pathToSequence = Path to the _ensembl_ whole genome file, fasta format\n"
            "pathToConfig   = Configuration file. Comment lines begin with '#'.\n"
            "                 contains list of transcripts and their relative abundance \n"
            "                 (an integer) and a line with ReadLength\n"
            "numReads       = the number of reads generated\n"
            "pathToFastq    = the output file name prefix")



def read_gtf(PathToGtf = None):
    """
    ARGS:
        PathToGtf = path to gene transfer file

    RETURN:
        A list of GTF_ENTRYs

    DESCRIPTION:
        Reads in gtf file.

    DEBUG:
        Can reproduce input gtf file. Works as expected.

    FUTURE: 
    """
    if(PathToGtf[-3:] != "gtf"):
        sys.stderr.write("ERROR! You did not pass a file with the .gtf extention\n")
        sys.exit(1)

    gtfFile = open(PathToGtf, 'r')
    gtfList = []
    gtfEntry = None
    timeBegin = datetime.datetime.now()

    for line in gtfFile:
        if(line[0] == "#"):
            continue
        line = line.split("\t")
        # Format check
        if(len(line) != 9):
            sys.stderr.write("ERROR! There should be 9 tab separated columns in a GTF file\n"
                             "You only have %i\n"%(len(line)))
            sys.exit(1)
        gtfEntry = GTF_ENTRY(Chromosome = line[0], Source = line[1], EntryType = line[2],
                             Start = line[3], Stop = line[4], Score = line[5], Strand = line[6],
                             Frame = line[7], Attribute = line[8])
        gtfList.append(gtfEntry)

    gtfFile.close()
    timeEnd = datetime.datetime.now()
    # Debug Code
    #for gtfEntry in gtfList:
    #    gtfEntry.print_entry()
    print("read_gtf()      run time = %s"%(timeEnd - timeBegin))
    return gtfList


def get_exon_list(gtfList):
    """
    ARGS:
        PathToGtf = path to gene transfer file

    RETURN:
        A list of GTF_ENTRYs

    DESCRIPTION:
        Note : this list is not unique. I.e. there are exons that will be duplicated.
        I make no effort to correct / eliminate duplicates b/c it would complicate
        mapping trans.exonIdxList for all the TRANSCRIPTs

    DEBUG:

    FUTURE: 
    """
    exonList = []
    timeBegin = datetime.datetime.now()
    for gtfEntry in gtfList:
        if(gtfEntry.etype == "exon"):
            exon = EXON(gtfEntry)
            exonList.append(exon)
    timeEnd = datetime.datetime.now()
    print("get_exon_list() run time = %s"%(timeEnd - timeBegin))
    return exonList


def get_transcript_list(gtfList, allExonList):
    """
    ARGS:
        GtfList = list of GTF_ENTRYs
        AllExonList = list of EXON entries

    RETURN:
        A list of TRANSCRIPTS

    DESCRIPTION:
        Scales like crap. O(N**2), which is too slow for the whole genome + gtf

    DEBUG: 
        Appears to correctly read in transcripts and associated exons.
        I checked by eye on first 87 lines of grch38.83.gtf

    FUTURE: 
    """
    transList = []
    timeBegin = datetime.datetime.now()
    for gtfEntry in gtfList:
        if(gtfEntry.etype == "transcript"):
            trans = TRANSCRIPT(gtfEntry)
            transList.append(trans)
    timeEnd = datetime.datetime.now()
    print("get_transcript_list() run time = %s"%(timeEnd - timeBegin))
    return transList



def get_gene_list(gtfList, allTransList):
    """
    ARGS:
        GtfList = list of GTF_ENTRYs
        AllTransList = list of TRANSCRIPT entries

    RETURN:
        A list of GENES

    DESCRIPTION:
        Scales like crap. O(N**2), which is too slow for the whole genome + gtf

    DEBUG: 
        Appears to correctly read in genes and associated transcripts.
        I checked by eye on first 87 lines of grch38.83.gtf

    FUTURE: 
    """
    geneList = []
    timeBegin = datetime.datetime.now()
    for gtfEntry in gtfList:
        if(gtfEntry.etype == "gene"):
            gene = GENE(gtfEntry)
            geneList.append(gene)
    timeEnd = datetime.datetime.now()
    print("get_gene_list() run time = %s"%(timeEnd - timeBegin))
    return geneList





def read_genome(seqPath):
    """
    ARGS:

    RETURN:
        A list of Chromosomes

    DESCRIPTION:
        Reads in a _whole_ genome fasta file and returns a list of Chromosomes. It
        follows the ensemble file format. Read in all sequences as upper case.

    DEBUG: 
        For a shortened whole genome fasta file, it correctly reads it in.

    FUTURE: 
    """
    seqFile = open(seqPath, "r")    # 1st fasta file 
    chrmList = []
    chrm = None
    seq = ""
    timeBegin = datetime.datetime.now()

    # Read file
    for line in seqFile:
        if(line[0] == "#"):                  # Skip ali's comments in shortened file
            continue
        if(line[0] == '>' or line[0] == '\n'):
            # Create chromosome, reset seq
            if(chrm is not None):
                #sys.stdout.write("chromosome : %s\n"%(chrm))
                chrmList.append(CHROMOSOME(chrm, seq))
                seq = ""
            #chrm = line[1:]            # Skip '>' for chrm name
            line = re.split(' ', line)      # only want chrm name
            chrm = (line[0])[1:]        # Skip '>' for chrm name
            continue
        else:
            line = re.split('\n', line)     # line is now a list of strings
            seq += line[0].upper()
    chrmList.append(CHROMOSOME(chrm, seq))  # pick up final chromosome.
    #sys.stdout.write("chromosome : %s\n"%(chrm))

    # Short debugging section, check that genome.fa duplicates original
    #for chrm in chrmList:
    #    sys.stdout.write("\n%s\n"%chrm.chrm)
    #    for i in range(0,len(chrm.seq)):
    #        sys.stdout.write("%c"%(chrm.seq[i]))
    #        if((i+1) % 60 == 0 and i != 0):
    #            sys.stdout.write("\n")

    seqFile.close()
    timeEnd = datetime.datetime.now()
    print("read_genome()   run time = %s"%(timeEnd - timeBegin))
    return chrmList




def get_exon_seq(exonList, chrmList):
    """
    ARGS:

    RETURN:
        None

    DESCRIPTION:
        Gets sequences for exons from chromosome list

    DEBUG: 
        Spot checked 3 exons, are all ok. More testing needed, however it is challenging
        to get a list of all the exons (incl. seqs) in a single file

    FUTURE: 
    """
    timeBegin = datetime.datetime.now()
    for exon in exonList:
        for chrm in chrmList:
            chrmLen = len(chrm.seq)

            if(chrm.chrm == exon.chrm):
                start = exon.start - 1        # -1 b/c python is 0 indexed, gtf file isnot
                end   = exon.stop
                if(start >= chrmLen or end >= chrmLen):
                    sys.stderr.write("ERROR!! start (%i) or stop (%i) Position > chromosome length (%i)\n"%
                                     (start, end, chrmLen))
                    sys.exit(1)

                if(exon.strand == '+'):
                    exon.seq = chrm.seq[start:end]
                elif(exon.strand == '-'):
                    exon.seq = reverse_complement(chrm.seq[start:end])
                else:
                    sys.stderr.write("ERROR! strand char = %s is invalid", exon.strand)
                    sys.exit(1)
    timeEnd = datetime.datetime.now()
    print("get_exon_seq()  run time = %s"%(timeEnd - timeBegin))
    


def get_trans_seq(transList, exonList):
    """
    ARGS:

    RETURN:
        None

    DESCRIPTION:
        Gets sequences for transcripts from chromosome list

    DEBUG: 
        Tested on 2 transcripts, more testing required. Getting a transcript file
        with the transcripts and sequences is challenging though

    FUTURE: 
    """
    timeBegin = datetime.datetime.now()
    for trans in transList:
        exonNum = 0         # use to check that indexs are loaded in order
        prevExonNum = 0
        for exon in exonList:
            if(exon.transID == trans.transID):
                exonNum = int(exon.exonNum)
                if(exonNum - prevExonNum != 1):
                    sys.stderr.write("ERROR! exon numbers for %s are loaded out of order!\n"%
                                    (trans.transID))
                    sys.exit(1)
                if(trans.seq is None):
                    trans.seq = exon.seq
                else:
                    trans.seq += exon.seq
                prevExonNum = exonNum
    timeEnd = datetime.datetime.now()
    print("get_trans_seq() run time = %s"%(timeEnd - timeBegin))


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



def get_list_of_unique_gtf_features(gtfList):
    """
    ARGS:
        gtfList : list of all GTF_ENTRYs

    RETURN:
        uniqueFeatureList : list of unique features

    DESCRIPTION:
        Finds all the unique features (ie column 3) of the GTF file

    DEBUG: 

    FUTURE: 
    """ 
    featureList = []
    uniqueFeatureList = []
    for gtfEntry in gtfList:
        featureList.append(gtfEntry.etype)
    for feature in featureList:
        if feature not in uniqueFeatureList:
            uniqueFeatureList.append(feature)
    return uniqueFeatureList



def link_exons_trans_and_genes(gtfList, exonList, transList, geneList):
    """
    ARGS:
        gtfList  : list of all GTF_ENTRYs
        exonList : list of all EXONS 
        transList: list of all TRANSCRIPTS
        geneList : list of all GENES

    RETURN:

    DESCRIPTION:
        Loops through gtfList and captures the indices of exons in exonList and passes it to 
        the transcripts in transList.  Also captures indices of transcripts in transList and
        passes it to genes in geneList.
        
        Does this in one pass through gtfList and scales roughly O(N). Should be faster
        than previous versions.

    DEBUG: 
        I validated by using print_transcripts_with_seqs() and comparing against the biomart
        download for chromosome 1. My data file was _identical_ to biomart's. 

        For how this was done, see the debug comment in print_transcripts_with_seqs() 
    FUTURE: 
    """ 
    gIdx = 0            # Gene index, for geneList
    tIdx = 0            # Transcript index, for transList
    eIdx = 0            # Exon index, for exonList
    gtfListLen = len(gtfList)
    timeBegin = datetime.datetime.now()
    
    # Ugly, and non-pythonic b/c i cant find cleaner way of accessing the next gtfEntry in the list
    for i in range(len(gtfList)):
        if(gtfList[i].etype == "gene"):
            # Check that genes in geneList are same order as gtfList 
            if(gtfList[i].geneID != geneList[gIdx].geneID):
                sys.stderr.write("ERROR! gtfList[%i].geneID = %s and geneList[%i].geneID = %s"%(
                                 i, gtfEntry.geneID, gIdx, geneList[gIdx].geneID))
                sys.exit(1)
            j = i + 1

            # Get constituent transcripts between gene entries
            while(gtfList[j].etype != "gene"):
                if(gtfList[j].etype == "transcript"):

                    # Check that transcripts in transList are same order as gtfList 
                    # Checking transcripts after gene in gtf _actually_ are members of the gene
                    # Add trans info to appropriate geneList[]
                    if(gtfList[j].transID == transList[tIdx].transID and
                      gtfList[i].geneID == transList[tIdx].geneID    and
                      gtfList[i].geneID == geneList[gIdx].geneID):
                        geneList[gIdx].transList.append(transList[tIdx].transID)
                        geneList[gIdx].transIdxList.append(tIdx)
                        k = j + 1

                        # Get constituent exons between transcript entries
                        while(gtfList[k].etype != "transcript"):
                            if(gtfList[k].etype == "exon"):
                                # Check that exons in exonList are same order as gtfList 
                                # Checking exons after transcript in gtf are members of the transcript 
                                # Add exon info to appropriate transList[]
                                if(gtfList[k].transID == exonList[eIdx].transID and
                                  gtfList[i].geneID == exonList[eIdx].geneID   and
                                  gtfList[i].geneID == geneList[gIdx].geneID):
                                    transList[tIdx].exonList.append(exonList[eIdx].exonID)
                                    transList[tIdx].exonIdxList.append(eIdx)
                                    eIdx += 1
                                else:
                                    sys.stderr.write("ERROR! gtfList[%i].transID = %s and exonList[%i]."
                                                  "transID = %s\n\tgtfList[%i].geneID = %s and "
                                                  "transList[%i].geneID = "
                                                  "%s\n\tand geneList[%i].geneID = %s\n"%
                                                  (k, gtfList[k].transID, eIdx, exonList[eIdx].transID,
                                                  k, gtfList[k].geneID, tIdx, transList[tIdx].geneID,
                                                  gIdx, geneList[gIdx].geneID))
                                    sys.exit(1)
                            k += 1
                            if(k == gtfListLen):
                                break
                        tIdx += 1   
                    else:
                        sys.stderr.write("ERROR! gtfList[%i].transID = %s and transList[%i].transID = "
                                         "%s\n\tgtfList[%i].geneID = %s and transList[%i].geneID = "
                                         "%s\n\tand geneList[%i].geneID = %s\n"%
                                         (j, gtfList[j].transID, tIdx, transList[tIdx].transID,
                                         j, gtfList[j].geneID, tIdx, transList[tIdx].geneID,
                                         gIdx, geneList[gIdx].geneID))
                        sys.exit(1)
                j += 1
                if(j == gtfListLen):
                    break
            gIdx += 1
            
    # Now get transcript sequences.
    for trans in transList:
        trans.seq = ""
        for eIdx in trans.exonIdxList:
            trans.seq += exonList[eIdx].seq

    timeEnd = datetime.datetime.now()
    print("link_exons_trans_and_genes() run time = %s"%(timeEnd - timeBegin))

    

def create_gene_and_trans_lookup_dict(geneList, transList):
    """
    ARGS:
        geneList  : list of GENEs
        transList : list of TRANSCRIPTS

    RETURN:
        geneDict  : keys = geneID, value geneList indices
        transDict : keys = transID, value transList indices

    DESCRIPTION:
        Dictionaries which are associative arrays with the
        geneID and transID as the keys and the associated geneList and transList
        indices as the values.

    DEBUG: 
        Spot checked the resulting dictionaries, appears to be correct

    FUTURE: 
    """ 
    geneDict = {}
    transDict = {}
    timeBegin = datetime.datetime.now()
    
    for gIdx in range(len(geneList)):
        geneDict[(geneList[gIdx].geneID)[1:-1]] = gIdx      # [1:-1] to strip leading and trailing "
    for tIdx in range(len(transList)):
        transDict[(transList[tIdx].transID)[1:-1]] = tIdx

    timeEnd = datetime.datetime.now()
    print("create_gene_and_trans_lookup_dict() run time = %s"%(timeEnd - timeBegin))
    return geneDict, transDict



def print_gtf_statistics(exonList, transList, geneList):
    """
    ARGS:
        exonList : list of EXONs
        transList : list of TRANSCRIPTs
        geneList  : list of GENEs

    RETURN:

    DESCRIPTION:
        Prints interesting statistics
        

    DEBUG: 
        I spot checked the results of transPairsDiffFile.txt and it appears to 
        return sets of transcripts that actually differ by 1 exon.  I checked 
        by searching in ensembl.

        Using a short.gtf, I checked the using Excel (identical):
            meanExonLen
            sigmaExonLen
            meanTransLen
            sigmaTransLen
            maxTransLen
            minTransLen
            minExonLen
            maxExonLen
            exonWMaxLen
            exonWMinLen
            transWMaxLen
            transWMinLen
        
    FUTURE: 
    """ 
    # exons
    meanExonLen = 0
    sigmaExonLen= 0
    maxExonLen  = 0
    minExonLen  = 10000
    exonWMaxLen = 0
    exonWMinLen = 10000
    meanExonNum = 0
    sigmaExonLen = 0
    # transcripts
    meanTransLen= 0
    maxTransLen = 0
    minTransLen = 10000
    transWMaxLen = ""
    transWMinLen = ""
    sigmaTransLen = 0
    # gene
    meanTransNum = 0
    timeBegin = datetime.datetime.now()

    # genes 
    for gene in geneList:
        meanTransNum += len(gene.transIdxList)
    
    # transcripts
    for trans in transList:
        transLen = len(trans.seq)
        meanTransLen += transLen
        # Try max and min later...
        if(transLen < minTransLen):
            minTransLen = transLen
            transWMinLen = trans.transID
        if(transLen > maxTransLen):
            maxTransLen = transLen
            transWMaxLen = trans.transID
        meanExonNum += len(trans.exonIdxList)

    # Exons
    for exon in exonList:
        exonLen     = len(exon.seq)
        meanExonLen += exonLen
        #Find max and min
        if(exonLen < minExonLen):
            minExonLen = exonLen
            exonWMinLen = exon.exonID
        if(exonLen > maxExonLen):
            maxExonLen = exonLen
            exonWMaxLen = exon.exonID

    # Averages
    meanTransLen  = meanTransLen / float(len(transList))
    meanTransNum  = meanTransNum / float(len(geneList))
    meanExonNum   = meanExonNum  / float(len(transList))
    meanExonLen   = meanExonLen  / float(len(exonList))
    
    # Standard deviations 
    # Transcripts
    for trans in transList:
        transLen = len(trans.seq)
        sigmaTransLen += (transLen - meanTransLen)**2
    sigmaTransLen = (sigmaTransLen / float(len(transList) - 1))**0.5
    # Exons
    for exon in exonList:
        sigmaExonLen += (len(exon.seq) - meanExonLen)**2
    sigmaExonLen = (sigmaExonLen / float(len(exonList) - 1))**0.5
    

    print("\nMean exon length : %f"%(meanExonLen))
    print("Std. Dev. exon length : %f"%(sigmaExonLen))
    print("Exon w/ Min length : %s, len : %i"%(exonWMinLen, minExonLen))
    print("Exon w/ Max length : %s, len : %i\n"%(exonWMaxLen, maxExonLen))
    print("Mean trans length : %f"%(meanTransLen))
    print("Std. Dev. trans Length : %s"%(sigmaTransLen))
    print("Trans w/ Min length : %s, len : %i"%(transWMinLen, minTransLen))
    print("Trans w/ Max length : %s, len : %i"%(transWMaxLen, maxTransLen))
    print("Mean num of exons per trans : %f\n"%(meanExonNum))

    timeEnd = datetime.datetime.now()
    print("print_gtf_statistics() run time = %s"%(timeEnd - timeBegin))
                        
                    

def find_trans_that_differ_by_1_exon(geneList, transList):
    """
    ARGS:
        transList : list of TRANSCRIPTs
        geneList  : list of GENEs

    RETURN:

    DESCRIPTION:
        Prints interesting statistics
        

    DEBUG: 
        I spot checked the results of transPairsDiffFile.txt and it appears to 
        return sets of transcripts that actually differ by 1 exon.  I checked 
        by searching in ensembl.
    FUTURE: 
    """ 
    transPairsDiffFile = open("transPairsDiffFile.txt", "w+")
    transPairsDiffFile.write("# This file contains a list of transcripts that differ only by 1 exon\n")
    transPairsDiffFile.write("# exon transcript1 transcript2 exonDiff\n")

    # Find transcripts that vary by 1 exon. Write to transPairsDiffFile.txt
    # This for loop should probably be its own separate program, but since
    # I already have all the functions in this program
    # This loop could be commented out in the final version of the code.
    #
    # Find transcripts that differ by 1 exon
    for gene in geneList:
        for tidx1 in range(len(gene.transIdxList)):              # transcriptIdx1
            trans1 = transList[gene.transIdxList[tidx1]]
            tidx2 = tidx1 + 1
            while tidx2 < len(gene.transIdxList):
                trans2 = transList[gene.transIdxList[tidx2]]
                
                # Follow Differ by one exon only
                if((len(trans1.exonList) > len(trans2.exonList) + 1)   or 
                   (len(trans1.exonList) < len(trans2.exonList) - 1)):
                    tidx2 += 1              # without this, it will loop infinitely
                    continue    
                else:
                    # below two operations return the elements in argument 1 that are not in
                    # argument 2. Hence why we do it twice, exonDiff and exonDiff2 won't be the
                    # same when there is a difference of 1 exon.
                    exonDiff = list(set(trans1.exonList) - set(trans2.exonList))
                    exonDiff2 = list(set(trans2.exonList) - set(trans1.exonList))
                    #exonDiff = [exn for exn in trans1.exonList if exn not in trans2.exonList]
                    if(len(exonDiff) == 1 and len(exonDiff2) == 0):
                        transPairsDiffFile.write("%s %s %s %s\n"%(gene.geneID, trans1.transID,
                                                 trans2.transID, exonDiff))
                    if(len(exonDiff) == 0 and len(exonDiff2) == 1):
                        transPairsDiffFile.write("%s %s %s %s\n"%(gene.geneID, trans1.transID,
                                                 trans2.transID, exonDiff2))
                tidx2 += 1
    transPairsDiffFile.close()





def read_config(pathToConfig):
    """
    ARGS:
        pathToConfig : str, path to configuration file

    RETURN:
        readLength      = readlength desired
        desiredTransList= List of transcripts to use
        abundanceList   = list of relative abundances of transcripts

    DESCRIPTION:
        Config file format :
            1. Comment lines begin with '#'
            2. first non-header line begins with 'ReadLength'
            3. All subsequent lines must be transcripts with relative abundance
               The relative abundance can be any integer. Scaling is done 
               automatically. 
               E.g.
                    ENST00000488147 10
                    ENST00000473358 5
    DEBUG: 
        For small config files it reads in all the fields correctly.

    FUTURE: 
    """ 
    desiredTransList = []
    abundanceList = []              # integers used to get relative abundance of transcripts
    readLength = 0
    numOfReads = 0

    configFile = open(pathToConfig, 'r')
    for line in configFile:
        if(line[0] == "#"):
            continue
        line = (line.split("\n"))[0]        # remove trailing \n

        # Check for tabs, only spaces permitted
        if(re.search('\t',line)):
            sys.stderr.write("ERROR! Tabs not permitted in config file!\n")
            sys.exit(1)
        line = line.split(" ")

        # ReadLength
        if(line[0] == "ReadLength"):
            if(readLength == 0):
                readLength = int(line[1])
                continue
            else:
                sys.stderr.write("ERROR! multiple instances of ReadLength in config file\n")
                sys.exit(1)

        # NumberOfReads
        if(line[0] == "NumberOfReads"):
            if(numOfReads == 0):
                numOfReads = int(line[1])
                continue
            else:
                sys.stderr.write("ERROR! multiple instances of ReadLength in config file\n")
                sys.exit(1)
        
        # Transcripts
        if(re.search('ENST', line[0])):
            desiredTransList.append(line[0])
            abundanceList.append(int(line[1]))
        else:
            sys.stderr.write("ERROR! Incorrect transcript entry : %s\n"
                             " All entries should begin with 'ENST'\n"%(line))
            sys.exit(1)
            
    if(readLength == 0 or numOfReads == 0):
        sys.stderr.write("ERROR! ReadLength or NumberOfReads not specified in config.txt\n")
        sys.exit(1)

    print("Config File Parameters : \nReadLength : %i\nNumberOfReads : %i"%(readLength, numOfReads))
    i = 0
    for trans in desiredTransList:
        print("%s %i"%(trans,abundanceList[i]))
        i += 1
    print("\n")

    return readLength, desiredTransList, abundanceList, numOfReads



def print_transcripts_with_seqs(transList):
    """
    ARGS:
        transList

    RETURN:
        None

    DESCRIPTION:
        This takes in a transList, makes a deep copy b/c I sort it.  I do not want 
        to mess up the actual order of transList, since the indices stored by 
        GENEs in transIdxList[] depends on the maintaining the order in transList.
    
        This function sorts all the transcripts numerically.
        This function is primarily used to debug link_exons_trans_and_genes()

    DEBUG: 
        1. Created chr1.fa and chr1.gtf, passed as args to this program
        2. Ran this function to print out _sorted_ transcript fasta file, called
           trans_list_and_seq.fa
        3. Went to http://useast.ensembl.org/biomart/martview, selected 
           chr1, sequences->cdna and downloaded the file. Downloaded as 
           biomart_export_cdna.fa 
        4. Created a program, validate/sort_biomart.py to read in biomart_export_cdna.fa 
           and generated sorted_biomart_export_cdna.fa
        5. ran 'diff trans_list_and_seq.fa validate/sorted_biomart_export_cdna.fa', it
           was identical.
           
           CONCLUSION:
               I correctly handle mapping the sequences to exons, and exons to 
               transcripts. If these weren't done correctly, (e.g. not handling 
               reverse complements), there is no way that this would work correctly.
        
               The following functions _MUST_ be correctly working (at least w/r/t
               transcripts):
                        GTF_ENTRY.__init__() 
                        CHROMOSOME.__init__()
                        TRANSCRIPT.__init__()
                        EXON.__init__()
                        read_gtf()  
                        get_exon_list()
                        get_transcript_list()
                        read_genome()
                        get_exon_seq()
                        get_trans_seq()
                        reverse_complement()
                        ---->  link_exons_trans_and_genes() <----
    FUTURE: 
    """ 
    transCopyList = copy.deepcopy(transList)
    outFile = open("trans_list_and_seq.fa", "w+")
    transCopyList.sort(key=operator.attrgetter('transNum'))
    
    for trans in transCopyList:
        outFile.write(">%s|%s\n"%(trans.geneID[1:-1], trans.transID[1:-1]))
        i = 0 

        for i in range(len(trans.seq)):
            outFile.write("%s"%(trans.seq[i]))
            if((i + 1) % 100 == 0):
                outFile.write("\n")

        outFile.write("\n")

    outFile.close()
    

def create_insert(Transcript, readLength, mu, sigma):
    """
    ARGS:
        Transcript = a TRANSCRIPT instance
        mu = the mean of fragment length distribution
        sigma = the standard deviation of fragment length distribution

    RETURN:
        AN INSERT of length n, where n fall in a distribution of rnorm(mu,sigma)

    DESCRIPTION:

    DEBUG: 

    FUTURE: 
    """
    start = 0
    stop = 0
    timeBegin = datetime.datetime.now()
    transLength = len(Transcript.seq)

    # type check
    if(not isinstance(Transcript, TRANSCRIPT)):
        sys.stderr.write("ERROR! Transcript is not of class type TRANSCRIPT 1\n")
        sys.exit(1)

    insertLength = 0

    while ( insertLength < readLength ):
        start = random.randint (0, transLength -1 )    
        stop = start + int (numpy.random.normal(mu, sigma))
        if (stop > transLength - 1):
            insert = INSERT(Transcript, start, transLength-1)
        else:  
            insert = INSERT(Transcript, start, stop)
        insertLength = len(insert.seq)        

    timeEnd = datetime.datetime.now()
    # print("get_insert_list() run time = %s"%(timeEnd - timeBegin))
    return insert


def create_fastq_file(pathToFastq, desiredTransList, abundanceList, numOfReads, readLength,
                      transDict, transList, exonList):
    """
    ARGS:
        pathToFastq      : Path to output fastq file
        desiredTransList : transcripts read from config file
        abundanceList    : list of integers that sum to a value used to normalize 
                           the number of reads. 
                                E.g. trans1 has 5 and trans2 has 10, 
                                     the ratio of trans1 : trans2 = 1 : 2
        numOfReads       : Rough number of desired reads, the ratios from abundanceList
                           is maintained at the expense of numOfReads. 
                                E.g from the above example if numOfReads = 10, 
                                    the actual number of reads would be 
                                    3 for trans1, 6 for trans2
        readLength       : length of reads
        transDict        : Dictionary used map transID quickly to the correct transcript in 
                           transList
        transList        : List of TRANSCRIPTs. Contains sequence to pass to instance of 
                           FASTQ_READ()
        exonList         : List of EXONs. Passed to FASTQ_READ() to get metadata for each
                           fastq read. E.g. the start position and exons a read spans.
    RETURN:
        None. File written

    DESCRIPTION:
        Writes fastq file.

    DEBUG: 
        Blasted against ensembl database, spot checked a couple of transcripts.
        Need further debugging. 

        Took synthetic_sample.fastq and operated in vim on transcript ENST00000473358:
            Exons are ordered as : ENSE00001947070 ENSE00001922571 ENSE00001827679
            Copied synthetic_sample.fastq to poop.fq

            ****** in vim ******
            %s/^@Read.\{1,1000\}:start:\([0-9]\{1,100\}\):exons:\(.\{1,100\}\)\n\(^[A-Z]\{50\}\)\n^+\n
                                \(^.\{50\}\)/\1\t\3\t\2/gc
            %s/:/\t/g   # remove colon sep exon names
            %s/"//g     # remove " around exon names

            ****** in shell, want exon reads start at (see order above) ******
            ****** Avoid grepping enties with start positions on the exon prior ******
            grep ENSE00001947070 poop.fq &> ENSE00001947070.txt
            grep ENSE00001922571 poop.fq | grep -v ENSE00001947070 &> ENSE00001922571.txt
            grep ENSE00001827679 poop.fq | grep -v ENSE00001922571 &> ENSE00001827679.txt
            awk '{print $1 "\t" $2}' ENSE00001947070.txt &> ENSE00001947070_1and2.txt
            awk '{print $1 "\t" $2}' ENSE00001922571.txt &> ENSE00001922571_1and2.txt
            awk '{print $1 "\t" $2}' ENSE00001827679.txt &> ENSE00001827679_1and2.txt
            awk '{print $2}' ENSE00001947070.txt | xargs -I{} grep -aob {} ENST00000473358.txt | 
                    awk 'BEGIN{FS=":"}{start = $1 + 29554; print start "\t" $2}' &> awk_out.txt
            awk '{print $2}' ENSE00001922571.txt | xargs -I{} grep -aob {} ENST00000473358.txt | 
                    awk 'BEGIN{FS=":"}{start = $1 + 30564 - 486; print start "\t" $2}' &> awk_out.txt
            awk '{print $2}' ENSE00001827679.txt | xargs -I{} grep -aob {} ENST00000473358.txt | 
                    awk 'BEGIN{FS=":"}{start = $1 + 30976 - 486 - 104; print start "\t" $2}' &> awk_out.txt
            Used diff to compare all the awk_out.txt to ENSE*_1and2.txt files.
                    CONCLUSION : they are identical. Therefor I get the correct start position from the
                                 correct sequences.
                    THEREFOR : I believe that create_fastq_file and FASTQ_READ() are working as expected.

    FUTURE: 
        Include more error checking for goofy parameters, e.g. not enough reads for
        the ratios, etc.
    """

    pathToFastqR1 = pathToFastq + "-R1.fq"
    pathToFastqR2 = pathToFastq + "-R2.fq"
    fastqFileR1 = open(pathToFastqR1, "w+")
    fastqFileR2 = open(pathToFastqR2, "w+")
    fastqListR1 = []
    fastqListR2 = []
    abundanceSum = 0
    transIdx = 0
    readIdx = 0
    
    for abundance in abundanceList:
        abundanceSum += abundance
    #abundanceNormalization = abundanceNormalization / len(abundanceList)    # integer division
    if(abundanceSum < 1):
        sys.stderr.write("ERROR! abundanceSum = %i\n"
                         "Please enter abundance values > 1\n"%(abundanceNormalization))
        sys.exit(1)
    
    for transName in desiredTransList:
        try:
            trans = transList[transDict[transName]]
        except KeyError:
            sys.stderr.write("ERROR! %s is not a transcript annotated in your gtf file\n"%(transName))
            sys.exit(1)
        for i in range( int(float(abundanceList[transIdx]) / float(abundanceSum) * numOfReads)):
            insert = create_insert (trans, readLength, 150, 15 )
            fastqEntry = FASTQ_READ(trans, insert, readLength, "@Read_num:%i"%(readIdx), exonList, "R1")
            fastqListR1.append(fastqEntry)

            fastqEntry = FASTQ_READ(trans, insert, readLength, "@Read_num:%i"%(readIdx), exonList, "R2")
            fastqListR2.append(fastqEntry)
            readIdx += 1
        transIdx += 1

    for fastqEntry in fastqListR1:
        fastqFileR1.write("%s\n%s\n+\n%s\n"%(fastqEntry.metadata, fastqEntry.seq, fastqEntry.qual))
        
    for fastqEntry in fastqListR2:
        fastqFileR2.write("%s\n%s\n+\n%s\n"%(fastqEntry.metadata, fastqEntry.seq, fastqEntry.qual))

    fastqFileR1.close()
    fastqFileR2.close()



#********************************************************************************
#********************************************************************************
#********************************  MAIN  ****************************************
#********************************************************************************
#********************************************************************************
def main():
    timeBegin = datetime.datetime.now()
    if(len(sys.argv) != 6):
        if(len(sys.argv) > 1 and (sys.argv[1] == "--help" or sys.argv[1] == "-h")):
            print_help()
            sys.exit(0)
        else:
            print_help()
            sys.exit(1)
    #pathToGtf = "/reference/homo_sapiens/GRCh38/ensembl/Annotation/Genes/gtf/Homo_sapiens.GRCh38.83.gtf"
    #pathToSeq = "/reference/homo_sapiens/GRCh38/ensembl/Sequence/WholeGenomeFasta/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
    pathToGtf = sys.argv[1]
    pathToSeq = sys.argv[2]
    pathToConfig = sys.argv[3]
    pathToFastq = sys.argv[5]
    #pathToFastq = "synth_samp.fastq"

    gtfList   = read_gtf(pathToGtf)
    exonList  = get_exon_list(gtfList)
    transList = get_transcript_list(gtfList, exonList)
    geneList  = get_gene_list(gtfList, transList)
    chrmList  = read_genome(pathToSeq)
    uniqueFeatureList = get_list_of_unique_gtf_features(gtfList)
    get_exon_seq(exonList, chrmList)
    link_exons_trans_and_genes(gtfList, exonList, transList, geneList)
    # print_transcripts_with_seqs(transList)      # Debug link_exons_trans_and_genes()

    geneDict, transDict = create_gene_and_trans_lookup_dict(geneList, transList)
    print_gtf_statistics(exonList, transList, geneList)
    # find_trans_that_differ_by_1_exon(geneList, transList) # Uncomment this if you want a complete list
    readLength, desiredTransList, abundanceList, numOfReads = read_config(pathToConfig)
    #random.seed(42)

    numOfReads = int(sys.argv[4])    

    for i in range(0,100):
        pathToFastq_2 = pathToFastq + "-" + str(i) + "-1to1"
        random.seed()
        create_fastq_file(pathToFastq_2, desiredTransList, abundanceList, numOfReads, readLength,
                      transDict, transList, exonList)
       

    print("Unique features in Gtf : ")
    for feature in uniqueFeatureList:
        print("\t%s"%(feature))
    timeEnd = datetime.datetime.now()
    print("Running time = %s"%(timeEnd - timeBegin))
    sys.exit(0)


if __name__ == "__main__":
    main()

