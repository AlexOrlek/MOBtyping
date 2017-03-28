import sys
sys.path.append('./mobtypingmodules')
from mobmods import unlist


#make accession : MOB type(s) dictionary from psiblast output

psiblast_dict={}

blast_file = open("./output_final/%s_BLASTplasmidproteins_%s_%s_best_hits.tsv" %(sys.argv[1],sys.argv[2], sys.argv[3])) 

for i, line in enumerate(blast_file):
    data=line.split('\t')
    key=data[1].split('|')
    key=str(key[0])
    mobtype=data[0].split('|')
    mobtype=str(mobtype[0])
    if key in psiblast_dict:
        psiblast_dict[key].append(mobtype)
    else:
        psiblast_dict[key]=[]  
        psiblast_dict[key].append(mobtype)
        
blast_file.close()


#get list of all plasmid accessions

accessions=[]

fileObj=open('./plasmiddatabase/%s_plasmidproteins.fa'%sys.argv[1])
for line in fileObj:
    if line.startswith('>'):
        data=line.split('|')
        accession=data[0].lstrip('>')
        if accession in accessions:
            pass
        else:
            accessions.append(accession)
    else:
        continue




#######################################


#TYPING

fileObj=open("./output_final/FINAL_OUTPUT_mobtyped_plasmids_%s_%s.tsv" %(sys.argv[2],sys.argv[3]), "w")

fileheader=["Plasmid accession id", "MOB type"]

psiblast_keys=list(psiblast_dict.keys())

for indx, accession in enumerate(accessions):
    if indx==0:
        fileObj.write('%s\n' %'\t'.join(fileheader))
    if accession in psiblast_keys:
        mobtypes=psiblast_dict[accession] #list of mobtype(s)
        mobtypes=sorted(mobtypes) #sort in alphabetical order
        mobtypes=str(unlist(mobtypes))
    else:
        mobtypes="-"
    fileObj.write('%s\t%s\n' % (accession, mobtypes))

print "finished mobtyping accessions"
