import sys, pickle
sys.path.append('./mobtypingmodules')
from mobmods import runsubprocess

#sys.argv[1] dataprefix; sys.argv[2]/[3] 0 0 (no id/coverage thresholds); sys.argv[4] evalue; sys.argv[5] maxiter

#extract sequence ids, alignment lengths and gapopen from blast table
fileObj=open('./output_intermediate/%s_BLASTplasmidproteins_%s_%s_iterationadded.tsv' %(sys.argv[1],sys.argv[4],sys.argv[5]))
qseqids=[] #query sequence id (i.e. mob protein id)
alignmentlengths=[]
gapopens=[]
pids=[]
for i, line in enumerate(fileObj):
    data=line.split('\t')
    qseqids.append(data[0]) #used to extract mob protein lengths from id table so they can be appended to blast table  
    alignmentlengths.append(data[3]) #used in coverage calculation
    gapopens.append(data[5]) #used in coverage calculation (subtract gapopens from alignment length otherwise they inflate coverage calculation)
    pids.append(data[2])
fileObj.close()
print len(qseqids), "number of mob query alignments (prior to selection of best hits)"


#extract mob protein lengths and names from idtable
fileObj2=open('./mobqueries/idtable_mobproteins.tsv')
moblengths=[]
mobnames=[]
for line in fileObj2:
    data=line.split('\t')
    moblengths.append(data[0])
    mobnames.append(data[2].strip())
fileObj2.close()


#extract strand/frame and calculate nucleotide qstart/qend
fileObj=open('./output_intermediate/%s_BLASTplasmidproteins_%s_%s_iterationadded.tsv' %(sys.argv[1],sys.argv[4],sys.argv[5]))

nuclstarts=[]
nuclends=[]
for i, line in enumerate(fileObj):
    data=line.split('\t')
    identifier=data[1].split(';')
    strand=str(identifier[1])  #strand1 or strand-1
    subjectstart=int(data[8])
    subjectend=int(data[9])
    if strand=='strand1':
        nuclstart=int((subjectstart-1)*3)+1
        nuclend=int(subjectend*3)
        nuclstarts.append(nuclstart)
        nuclends.append(nuclend)
    elif strand=='strand-1':
        plasmidlength=int(identifier[3])
        nuclstart=int(plasmidlength*3)-(int(subjectend*3)-1)
        nuclend=int(plasmidlength*3)-(int(int(subjectstart-1)*3))
        nuclstarts.append(nuclstart)
        nuclends.append(nuclend)
    else:
        print 'error: nuclstart/end not found', data, strand
    

fileObj.close()
print len(nuclstarts), 'nuclstarts' ###testing
print len(nuclends), 'nuclends'


#filter mob lengths to match qseqids
filteredmoblengths=[]
for i in range(len(qseqids)):
    filteredmoblengths.append(moblengths[mobnames.index(qseqids[i])])


#add mob length and coverage proportion columns + hitscore
data=open('./output_intermediate/%s_BLASTplasmidproteins_%s_%s_iterationadded.tsv' %(sys.argv[1],sys.argv[4],sys.argv[5])).readlines()
fileObj=open('./output_intermediate/%s_BLASTplasmidproteins2_%s_%s.tsv' %(sys.argv[1],sys.argv[4],sys.argv[5]), "w")
for indx, line in enumerate(data):
    moblength=filteredmoblengths[indx]
    pid=float(pids[indx])/float(100)
    coverage=(float(alignmentlengths[indx])-float(gapopens[indx]))/float(filteredmoblengths[indx])
    hitscore=coverage*pid
    fileObj.write("%s\t%s\t%s\t%s\n" %(line.strip(), moblength, coverage, hitscore))
fileObj.close()
    

#add nucl start and end columns
fileObj=open('./output_intermediate/%s_BLASTplasmidproteins2_%s_%s.tsv' %(sys.argv[1],sys.argv[4],sys.argv[5]))
fileObj2=open('./output_intermediate/%s_BLASTplasmidproteins2.1_%s_%s.tsv' %(sys.argv[1],sys.argv[4],sys.argv[5]), "w")
for i, line in enumerate(fileObj):
    line=line.strip()
    data=line.split('\t')
    data.append(str(nuclstarts[i]))
    data.append(str(nuclends[i]))
    fileObj2.write('%s\n' % '\t'.join(data))
fileObj.close()
fileObj2.close()



#filter blast table by PID and coverage (not necessary to set filtering thresholds i.e. sys.argv[2]/[3] == 0)
#also add accession column (to allow sorting of hits for each plasmid accession)

fileObj=open('./output_intermediate/%s_BLASTplasmidproteins3_%s_%s.tsv' %(sys.argv[1],sys.argv[4],sys.argv[5]), 'w')
for line in open('./output_intermediate/%s_BLASTplasmidproteins2.1_%s_%s.tsv' %(sys.argv[1],sys.argv[4],sys.argv[5])):
    line=line.strip()
    data=line.split('\t')
    accession=data[1].split(';')[0]
    data.append(accession)
    if float(data[2])>int(sys.argv[2]) and float(data[14])>float(sys.argv[3]):  #for no filtering use pid 0 and coverage 0 
        fileObj.write('%s\n' %'\t'.join(data))
fileObj.close()




######################################write new file with best hits


#order hits by subject plasmid accession, then by iteration, and then by hit score
filteredfile='./output_intermediate/%s_BLASTplasmidproteins3_%s_%s.tsv' %(sys.argv[1],sys.argv[4],sys.argv[5])
writefile='./output_intermediate/%s_BLASTplasmidproteins3_%s_%s_sorted_initial.tsv' %(sys.argv[1],sys.argv[4],sys.argv[5])
args=['sort', filteredfile, '-k19,19V', '-k13,13n', '-k16,16nr', '-t\t']  #sort by col19 (accession) alphanumerically (V); then sort col13 (iteration) numerically; then sort col16 (score) reverse numerically (i.e. favouring higher scores)
runsubprocess(args, writefile=writefile)


#re-write file, without accession info
fileObj=open('./output_intermediate/%s_BLASTplasmidproteins3_%s_%s_sorted_initial.tsv' %(sys.argv[1],sys.argv[4],sys.argv[5]))
fileObj2=open('./output_intermediate/%s_BLASTplasmidproteins3_%s_%s_sorted.tsv' %(sys.argv[1],sys.argv[4],sys.argv[5]),'w')
for line in fileObj:
    line=line.strip()
    data=line.split('\t')
    data=data[:-1]
    fileObj2.write('%s\n' %'\t'.join(data))
fileObj.close()
fileObj2.close()


#make dictionary based on sorted file, with subjects as keys and sorted hits as nested list of values                                            
sortedfile='./output_intermediate/%s_BLASTplasmidproteins3_%s_%s_sorted.tsv' %(sys.argv[1],sys.argv[4],sys.argv[5])

uniquesseqids=[]
subjectdict={}
fileObj=open(sortedfile)
for indx, line in enumerate(fileObj):
    line=line.strip()
    data=line.split('\t')
    subjectname=data[1].split(';')[0]
    if subjectname in subjectdict:
        subjectdict[subjectname].append(data)
    else:
        subjectdict[subjectname]=[]
        subjectdict[subjectname].append(data)
        uniquesseqids.append(subjectname)

print "finished making subject:hit dict"


#rangeslist: for each subject (plasmid) accession, get list of hit ranges associated with that accession

rangeslist=[]
for j in range(len(uniquesseqids)):
    rangeslistnest=[]
    for i in range(len(subjectdict[uniquesseqids[j]])):
        querylocus=[]                              
        querylocus.append(int((subjectdict[uniquesseqids[j]][i][16])))                               
        querylocus.append(int((subjectdict[uniquesseqids[j]][i][17])))                                
        querymax=max(querylocus)
        querymin=min(querylocus)
        rangeslistnest.append(range(int(querymin),int(querymax)+int(1)))                  
    rangeslist.append(rangeslistnest)


#find indices of CDSs intersected by hits  e.g. [[14,14,14,14],[...]] would mean that for the first accession, associated nucl start-end hit ranges all intersect at the same CDS (index 14)
#also recording cdsrange, and cds annotation as well as indices

picklefileObj= open('./plasmiddatabase/pickle.CDSfeatures_%s' %sys.argv[1])
CDSrange_dict=pickle.load(picklefileObj)
picklefileObj.close()


cdsindices=[] #listed list of indices representing the CDS that is intersected by a hit
cdsranges=[]  #listed list of nucl start-ends of the CDS
cdsannotations=[]
for j in range(len(uniquesseqids)):
    cdsindicesnestfinal=[]
    cdsrangesnestfinal=[]
    cdsannotationsnestfinal=[]
    for queryhitrange in rangeslist[j]: #queryhitrange is a list of consecutive numbers repersenting a range between nuclstart/nuclend
        queryhitrange=list(queryhitrange)
        cdsindicesnest=[]
        cdsrangesnest=[]
        cdsannotationsnest=[]
        cdsintersections=[] #if a hit intersects with multiple CDSs, choose the one it intersects most with
        try:
            for indx, cdsrange in enumerate(CDSrange_dict[uniquesseqids[j]]): #iterates through tuples of cds ranges; which cds does a hit range intersect with?
                cdsstart=cdsrange[0]
                cdsend=cdsrange[1]
                cdsannotation=cdsrange[2]
                cdsrangelist=list(range(int(float(cdsstart)), int(float(cdsend))))
                cdsintersection=len(set(queryhitrange).intersection(set(cdsrangelist)))
                if cdsintersection>0: #if there is any overlap with a cds, add the cds indx to cdsindicesnest
                    cdsindicesnest.append(indx)
                    cdsrangesnest.append(str(cdsrange[0])+'-'+str(cdsrange[1]))
                    cdsannotationsnest.append(cdsannotation)
                    cdsintersections.append(cdsintersection)
                    if len(cdsindicesnest)>1: #each hit should only intersect with one CDS
                        cdsindicesnest=[cdsindicesnest[cdsintersections.index(max(cdsintersections))]]
                        cdsrangesnest=[cdsrangesnest[cdsintersections.index(max(cdsintersections))]]
                        cdsannotationsnest=[cdsannotationsnest[cdsintersections.index(max(cdsintersections))]]
                        cdsintersections=[cdsintersections[cdsintersections.index(max(cdsintersections))]]
                    else:
                        pass
                else:
                    continue
        except:
            print "warning: accession not found in cdsrange dict - if features of untyped plasmids were added to pickle.CDSfeatures, all accessions should be included"
        if len(cdsindicesnest)==1:
            cdsindicesnestfinal.append(cdsindicesnest[0]) #ideally each hit range from a given plasmid accession should intersect with one cds range; if not, cdsindicesnest not appended and length discrepancy between number of intersects and number of hit ranges will produce warning (below)
            cdsrangesnestfinal.append(cdsrangesnest[0])
            cdsannotationsnestfinal.append(cdsannotationsnest[0])
        else: #i.e len(cdsindicesnest)==0
            pass
    if len(cdsindicesnestfinal) != len(rangeslist[j]):
        print "warning: num alignments on plasmid != num intersected CDSs on plasmid; best hit will be selected by overlap method rather than intersection at the same locus. Accession; num alignments:", uniquesseqids[j], len(rangeslist[j]) #this may occur if a cds region is not annotated e.g. NC_008382.1, of 13 hits to the plasmid, there are 4 cases where there is no intersection between the hit on the plasmid nucleotide and a cds on the plasmid
        cdsindices.append(None) #label as None - overlap method will be used for best hit selection
        cdsranges.append(None)
        cdsannotations.append(None)
    else:
        cdsindices.append(cdsindicesnestfinal)
        cdsranges.append(cdsrangesnestfinal)
        cdsannotations.append(cdsannotationsnestfinal)



#re-write file with intersect information (cdsindices/cdsranges/annotations or None - None means overlap method)

fileObj=open('./output_final/%s_BLASTplasmidproteins_%s_%s_all_hits.tsv' %(sys.argv[1],sys.argv[4],sys.argv[5]), 'w')

fileheader=["MOB query","Subject plasmid","Percent identity","Alignment length","Number of mismatches","Number of gap opens","Start position in query","End position in query","Start position in subject","End position in subject","E-value","Bit score","Iteration of hit retrieval","MOB query length","Coverage score","Overall score","Start position on plasmid nucleotide sequence","End position on plasmid nucleotide sequence","CDS index (arbitrary)", "Nucleotide position of plasmid CDS intersected by alignment","CDS annotation"]

cdsindices2=list(cdsindices) #cdsindices2 is duplicate but with None as string for each hit where cds intersect method can't be used
for indxa, accession in enumerate(uniquesseqids):
    if indxa==0:
        fileObj.write('%s\n' %'\t'.join(fileheader))
    if cdsindices2[indxa]==None: #change so there is a "None" for each hit
        cdsindices2[indxa]=["None"]*len(subjectdict[accession])
        cdsranges[indxa]=["None"]*len(subjectdict[accession])
        cdsannotations[indxa]=["None"]*len(subjectdict[accession])
    for indxb, data in enumerate(subjectdict[accession]):
        if len(subjectdict[accession])==len(cdsindices2[indxa]) and len(subjectdict[accession])==len(cdsranges[indxa]) and len(subjectdict[accession])==len(cdsannotations[indxa]):
            pass
        else:
            print "error - number of intersects associated with accession doesn't match number of hits", accession, subjectdict[accession], cdsindices2[indxa]
            sys.exit()
        fileObj.write('%s\t%s\t%s\t%s\n'% ('\t'.join(data), cdsindices2[indxa][indxb], cdsranges[indxa][indxb], cdsannotations[indxa][indxb]))        
    
fileObj.close()



#####################use a mixture of intersect approach and overlap approach (in few cases where intersects didn't match cds ranges) to filter



#first create dictionary of accessions: hits as before (now with intersect information); then write to file with filtering

fileObj=open('./output_final/%s_BLASTplasmidproteins_%s_%s_all_hits.tsv' %(sys.argv[1],sys.argv[4],sys.argv[5])) 
subjectdict2={}
for indx, line in enumerate(fileObj):
    if line.startswith("MOB query"):
        continue
    line=line.strip()
    data=line.split('\t')
    subjectname=data[1].split(';')[0]
    if subjectname in subjectdict2:
        subjectdict2[subjectname].append(data)
    else:
        subjectdict2[subjectname]=[]
        subjectdict2[subjectname].append(data)

print "finished making subjectdict2"

fileObj.close()



fileObj=open('./output_final/%s_BLASTplasmidproteins_%s_%s_best_hits.tsv' %(sys.argv[1],sys.argv[4],sys.argv[5]), 'w')

for indxa, accession in enumerate(uniquesseqids):
    if indxa==0:
        fileObj.write('%s\n' %'\t'.join(fileheader))
    if cdsindices[indxa]==None: #OVERLAP APPROACH
        print indxa, "indxa, blastfilter, overlap method"
        hitrangemoblengths=[] #list of tuples of included hits and the associated mob lengths                                                
        for indxb, data in enumerate(subjectdict2[accession]): #running through all hits associated with a given query           
            if indxb==0:  #for the highest scoring hit, include, irrespecitve of overlaps                                                      
                ranges=[int(data[16]),int(data[17])]
                moblength=int(data[13])
                minrange=ranges[ranges.index(min(ranges))]
                maxrange=ranges[ranges.index(max(ranges))]
                hitrange=range(minrange,maxrange+int(1))   #this is a list of consecutive numbers spanning the hit range                         
                hitrangemoblengths.append((hitrange, moblength))
                fileObj.write('%s\n' % '\t'.join(data))
            if indxb>0:
                ranges=[int(data[16]),int(data[17])]
                moblength=int(data[13])
                minrange=ranges[ranges.index(min(ranges))]
                maxrange=ranges[ranges.index(max(ranges))]
                hitrange=range(minrange,maxrange+int(1))
                #if hit range intersects, dont include; otherwise, include and add hitrange to ranges.                                            
                for indxc, includedhit in enumerate(hitrangemoblengths):  #checking non-best hit against all inlcuded hits; there will be at least one hit (the best hit) in hit ranges
                    includedhitrange=includedhit[0]
                    includedmoblength=includedhit[1]
                    moblengths=[moblength, includedmoblength]
                    minmoblength=moblengths[moblengths.index(min(moblengths))]
                    intersectlength=len(set(hitrange).intersection(set(includedhitrange)))
                    pairwiseoverlap=float(float(intersectlength)/float(minmoblength))            
                    if pairwiseoverlap<float(0.001): #0.001 is arbitrary threshold
                        add=True
                        continue
                    else:
                        add=False
                        break
                if add==True:
                    hitrangemoblengths.append((hitrange, moblength))
                    fileObj.write('%s\n'% '\t'.join(data))
                else:
                    pass #if add is false i.e. there was an overlap with a better hit, don't write to file, don't append to hitrangemoblengths


    else: #INTERSECT APPROACH (NEW)
        indices=[]
        for indx, cdsindex in enumerate(cdsindices[indxa]):
            if indx==0:
                indices.append(cdsindex)
                data=subjectdict2[accession][indx]
                fileObj.write('%s\n' % '\t'.join(data))
            else:
                if cdsindex in indices:
                    continue
                else:
                    indices.append(cdsindex)
                    data=subjectdict2[accession][indx]
                    fileObj.write('%s\n' % '\t'.join(data))



#remove all intermediate files

args=['find ./output_intermediate/ -maxdepth 1 ! -name ".*" -type f -delete']
runsubprocess(args, shell=True)
