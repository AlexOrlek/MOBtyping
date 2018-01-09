import sys
sys.path.append('./mobtypingmodules')
from mobmods import unlist


#make accession : (MOB type(s), accession;strand;frame(s)) dictionary from psiblast output

psiblast_dict={}

blast_file = open("./output_final/%s_BLASTplasmidproteins_%s_%s_best_hits.tsv" %(sys.argv[1],sys.argv[2], sys.argv[3])) 

for i, line in enumerate(blast_file):
    data=line.split('\t')
    strandframe=';'.join(data[1].split(';')[1:]) #strand;frame;plasmidlength
    accession=data[1].split(';')[0]
    mobtype=data[0].split('|')[0]
    sstart=data[8]
    send=data[9]
    if accession in psiblast_dict:
        psiblast_dict[accession].append((mobtype,strandframe,sstart,send))
    else:
        psiblast_dict[accession]=[]  
        psiblast_dict[accession].append((mobtype,strandframe,sstart,send))
        
blast_file.close()
print strandframe

#get list of untyped plasmid accessions

untypedaccessions=[]
fileObj=open('./plasmiddatabase/plasmidproteins_%s_untyped.fa'%sys.argv[1])
for line in fileObj:
    if line.startswith('>'):
        line=line.strip()
        data=line.split(';')
        accession=data[0].lstrip('>')
        if accession in untypedaccessions:
            pass
        else:
            untypedaccessions.append(accession)
    else:
        continue
fileObj.close()


#get list of original plasmid accessions


originalaccessions=[]

fileObj=open('./plasmiddatabase/%s_plasmidproteins.fa'%sys.argv[1])
for line in fileObj:
    if line.startswith('>'):
        data=line.split(';')
        accession=data[0].lstrip('>')
        if accession in untypedaccessions:
            continue
        elif accession in originalaccessions:
            pass
        else:
            originalaccessions.append(accession)
    else:
        continue

fileObj.close()


accessions=originalaccessions+untypedaccessions

#######################################


#TYPING
import operator

fileObj=open("./output_final/%s_FINAL_OUTPUT_mobtyped_plasmids_originaldatabase_%s_%s.tsv" %(sys.argv[1],sys.argv[2],sys.argv[3]), "w")
fileObj2=open("./output_final/%s_FINAL_OUTPUT_mobtyped_plasmids_%s_%s.tsv" %(sys.argv[1],sys.argv[2],sys.argv[3]), "w")

fileheader=["Plasmid accession id", "MOB type(s)","Strand;Frame;Plasmid_length_(nucleotide_bp)","Start position in subject","End position in subject"]

psiblast_keys=list(psiblast_dict.keys())

for indx, accession in enumerate(accessions):
    if indx==0:
        fileObj.write('%s\n' %'\t'.join(fileheader))
        fileObj2.write('%s\n' %'\t'.join(fileheader))
        
    if accession in psiblast_keys:
        mytuplelist=psiblast_dict[accession] #list of tuples of mobtype(s),accession;strand;frame(s)
        mobtypes=unlist(sorted(map(operator.itemgetter(0), mytuplelist)))
        strandframe=unlist(sorted(map(operator.itemgetter(1), mytuplelist)))
        sstart=unlist(sorted(map(operator.itemgetter(2), mytuplelist)))
        send=unlist(sorted(map(operator.itemgetter(3), mytuplelist)))
    else:
        mobtypes='-'
        strandframe='-'
        sstart='-'
        send='-'
    if accession in originalaccessions:
        fileObj.write('%s\t%s\t%s\t%s\t%s\n' % (accession, mobtypes, strandframe,sstart,send))
    elif accession in untypedaccessions:
        fileObj2.write('%s\t%s\t%s\t%s\t%s\n' % (accession, mobtypes,strandframe,sstart,send))
    else:
        print 'error: accession not recognised'
        sys.exit()
print "finished mobtyping accessions"
