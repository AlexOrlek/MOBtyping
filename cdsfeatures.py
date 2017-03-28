import sys, pickle
sys.path.append('./mobtypingmodules')
from mobmods import unlist
from Bio import Entrez, SeqIO


#this script extracts details on CDS features of untyped plasmid(s) from a Genbank file, or from a tsv file; then appends these details to pickle.CDSfeatures_original dictionary

CDSfeatures_dict={}
accessions=[]

#check file availability
try:
    fileObj=open("./untyped_plasmids/%s_CDSfeatures.tsv" %sys.argv[1]) #sys.argv[1] is the name of the untyped plasmid(s)
    fileObj.close()
    availablefile="tsv"
except:
    try:
        fileObj=open("./untyped_plasmids/%s.gb" %sys.argv[1])
        fileObj.close()
        availablefile="gb"
    except:
        availablefile="none"


#extract CDS feature information and create pickle.CDSfeatures
if availablefile=="tsv":  #if there is a tsv file, extract CDS feature information from there
    fileObj=open("./untyped_plasmids/%s_CDSfeatures.tsv" %sys.argv[1]) #sys.argv[1] is the name of the untyped plasmid(s)
    for line in fileObj:
        line=line.strip()
        data=line.split('\t')
        accession=data[0]
        if accession in accessions:
            pass
        else:
            accessions.append(accession)
        start=data[1]
        start=start.strip('<')
        start=start.strip('>') 
        end=data[2]
        end=end.strip('<')
        end=end.strip('>')    
        annotation=data[3]
        plasmidlength=data[4]

        difference=(int(start)-int(end))*int(-1)  #exclude misannotated start-ends (if CDS start-end length is close to plasmid length)
        if int(difference)>(int(plasmidlength)-int(10)):
            continue
        else:
            CDSrange=(start,end, annotation)
            #append to dictionary
            if accession not in CDSfeatures_dict:
                CDSfeatures_dict[accession]=[CDSrange]
            else:
                CDSfeatures_dict[accession].append(CDSrange)

    fileObj.close()

    picklefileObj=open("./plasmiddatabase/pickle.CDSfeatures_original", "rb")
    CDSfeatures_original=pickle.load(picklefileObj)
    num_accessions_original=len(CDSfeatures_original)
    picklefileObj.close()

    CDSfeatures_original.update(CDSfeatures_dict) #add CDSfeatures info from untyped plasmid(s) to CDSfeatures_origial dictionary
        
    picklefileObj=open("./plasmiddatabase/pickle.CDSfeatures_%s" %sys.argv[1], "wb")
    pickle.dump(CDSfeatures_original, picklefileObj)
    picklefileObj.close()
    print "extracted CDS feature information from tsv file"
        
else:
    
    if availablefile=="gb":  #if there is a genbank file, extract CDS feature information from there
        fileObj=open("./untyped_plasmids/%s.gb" %sys.argv[1])
        fileObj2=open("./untyped_plasmids/%s_CDSfeatures.tsv" %sys.argv[1],"w")
        for record in SeqIO.parse(fileObj, "genbank"):
            accession=record.id
            accessions.append(accession)
            plasmidlength=len(record.seq)
            for feature in record.features:
                start=str(feature.location.start)
                start=start.strip('<')
                start=start.strip('>') 
                end=str(feature.location.end)
                end=end.strip('<')
                end=end.strip('>')
                
                try:
                    annotation=unlist(feature.qualifiers["product"])
                except KeyError:
                    annotation="-"

                difference=(int(start)-int(end))*int(-1)  #exclude misannotated start-ends (if CDS start-end length is close to plasmid length)
                if int(difference)>(int(plasmidlength)-int(10)):
                    continue
                else:
                    CDSrange=(start,end, annotation)
                    #append to dictionary
                    if accession not in CDSfeatures_dict:
                        CDSfeatures_dict[accession]=[CDSrange]
                    else:
                        CDSfeatures_dict[accession].append(CDSrange)
                    #write to tsv file
                    fileObj2.write('%s\t%s\t%s\t%s\t%s\n' %(accession, start, end, annotation, plasmidlength))
                    
        fileObj.close()
        fileObj2.close()
            
        picklefileObj=open("./plasmiddatabase/pickle.CDSfeatures_original", "rb")
        CDSfeatures_original=pickle.load(picklefileObj)
        num_accessions_original=len(CDSfeatures_original)
        picklefileObj.close()

        CDSfeatures_original.update(CDSfeatures_dict) #add CDSfeatures info from untyped plasmid(s) to CDSfeatures_origial dictionary
        
        picklefileObj=open("./plasmiddatabase/pickle.CDSfeatures_%s" %sys.argv[1], "wb")
        pickle.dump(CDSfeatures_original, picklefileObj)
        picklefileObj.close()
        print "extracted CDS feature information from gb file"
                        
                                                            
    elif availablefile=="none":  #if CDS feature details are not provided as gb or tsv, just use pickle.CDSfeatures_original (CDS feature information from untyped plasmid(s) will remain missing)
        print "CDS features of untyped plasmid(s) not provided; check file extensions are correct (.gb or _CDSfeatures.tsv); pickle.CDSfeatures_original will be used"

        picklefileObj=open("./plasmiddatabase/pickle.CDSfeatures_original", "rb")
        CDSfeatures_original=pickle.load(picklefileObj)
        num_accessions_original=len(CDSfeatures_original)
        picklefileObj.close()
        
        picklefileObj=open("./plasmiddatabase/pickle.CDSfeatures_%s" %sys.argv[1], "wb")
        pickle.dump(CDSfeatures_original, picklefileObj)
        picklefileObj.close()
    else:
        print "code error"


print num_accessions_original, "number of accessions included in pickle.CDSfeatures_original"
print len(CDSfeatures_original), "number of accessions included in pickle.CDSfeatures"
print accessions, "untyped plasmid accession ids"
if availablefile!="none":
    for accession in accessions:
        print CDSfeatures_original[accession], "CDSfeatures dictionary contents for accession: %s" %accession
