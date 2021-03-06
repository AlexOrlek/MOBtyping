import sys
sys.path.append('./mobtypingmodules')
from mobmods import extractfromfasta, writeFastaSeqs, runsubprocess
from Bio.Seq import Seq #use biopython to translate sequence in all 6 frames
from Bio.Alphabet import generic_dna


#if there is not a fasta file representing untyped plasmids to be added to the plasmid database, convert from genbank file

from Bio import SeqIO
try:
    fileObj=open("./untyped_plasmids/%s.fa" %sys.argv[1]) #sys.argv[1] is the name of the untyped plasmid(s), used to label fasta file
    fileObj.close()
except:
    try:
        SeqIO.convert("./untyped_plasmids/%s.gb" %sys.argv[1], "genbank", "./untyped_plasmids/%s.fa" %sys.argv[1], "fasta")
    except:
        print "there is no fasta file or genbank file in the untyped_plasmids directory; check file extensions are .fa or .gb"


#extract plasmid dna sequences in order to translate to protein

fastafilepath="./untyped_plasmids/%s.fa" %sys.argv[1]
accessionseq=extractfromfasta(fastafilepath)
accessions=accessionseq[0]
sequences=accessionseq[1]


#translate in 6 frames

table = 11
translation_dict={}

for j in range(len(sequences)):
    print j
    translation_dict[accessions[j]]=[] 
    dna=str(sequences[j])
    dna=dna.strip()
    dna=Seq(dna,generic_dna)

    #check orfs in all 6 frames #source code: http://biopython-documentation-chinese-translate.googlecode.com/git-history/82ea42d22864986a277cd0647e2039c3e2c8880e/_build/html/ch16.html
    for strand, nuc in [(+1, dna), (-1, dna.reverse_complement())]:  #strand, nuc are +1 dna , then -1 rev.compl
        for frame in range(3):
            translation_dict[accessions[j]].append(str(nuc[frame:].translate(table))+'|'+'strand'+str(strand)+'|'+'frame'+str(frame))



#write to fasta

comments=[]
sequences=[]

for i in range(len(accessions)):
    accession=accessions[i]
    accessionid=accession.split(' ')[0]
    for j in range(len(translation_dict[accession])):
        data=translation_dict[accession][j].split('|')
        protein=data[0]
        strand=data[1]
        frame=data[2]
        proteinlength=len(protein)
        comments.append(str(accessionid)+';'+str(strand)+';'+str(frame)+';'+str(proteinlength))
        sequences.append(protein)


writeFastaSeqs(comments, sequences, "./plasmiddatabase/plasmidproteins_%s_untyped.fa" %sys.argv[1])


#concatenate plasmidproteins_original, plasmidproteins_untyped

runsubprocess(['cat ./plasmiddatabase/translatedproteinseq.fa | tr "|" ";" > ./plasmiddatabase/translatedproteinseq_v2.fa'],shell=True)

file1="./plasmiddatabase/translatedproteinseq_v2.fa" #curated plasmid database (proteins translated from nucleotide sequences); the file is on Figshare (https://figshare.com/s/18de8bdcbba47dbaba41)
file2="./plasmiddatabase/plasmidproteins_%s_untyped.fa" %sys.argv[1]
writefilepath="./plasmiddatabase/%s_plasmidproteins.fa" %sys.argv[1] 
args=['cat', file1, file2]
runsubprocess(args, writefile=writefilepath)
