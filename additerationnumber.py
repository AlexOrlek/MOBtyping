import sys
sys.path.append('./mobtypingmodules')
from mobmods import runsubprocess,getcol


dataprefix=sys.argv[1]
evalue=sys.argv[2]
maxiter=int(sys.argv[3])
mobtypes=getcol('./mobtyping_evalues.tsv', 0) 
evalues=getcol('./mobtyping_evalues.tsv', evalue)

convergencedict={}  #mob : maxmobiter
with open('./output_intermediate/convergence_%s_%s_%s.tsv'%(dataprefix,evalue,maxiter)) as f:
    for line in f:
        data=line.strip().split('\t')
        mob=data[0]
        iteration=data[1]
        converged=data[2]
        convergencedict[mob]=iteration


hitcountdict={} #mob : [iter 1 count, iter 2 count ...]
for (mob,e) in zip(mobtypes,evalues):
    maxmobiter=int(convergencedict[mob])
    iterations=list(range(1,(maxmobiter+1)))
    for iteration in iterations:
        hitcount=runsubprocess(['cat ./output_intermediate/%s_BLASTplasmidproteins_%s_%s_%s.tsv | sed "/^\s*$/d" | sed "/^Search has CONVERGED!*$/d" | wc -l'%(dataprefix, e,iteration,mob)],shell=True) #remove blank lines/search has converged message
        hitcount=hitcount.strip()
        if mob in hitcountdict:
            hitcountdict[mob].append(hitcount)
        else:
            hitcountdict[mob]=[]
            hitcountdict[mob].append(hitcount)

#print hitcountdict


#write a file: concatenate final iteration files for each mob type, annotating with iteration

f2=open('./output_intermediate/%s_BLASTplasmidproteins_%s_%s_iterationadded.tsv'%(dataprefix, evalue, maxiter),'w')
for (mob,e) in zip(mobtypes,evalues):
    maxmobiter=int(convergencedict[mob])
    with open('./output_intermediate/%s_BLASTplasmidproteins_%s_%s_%s.tsv'%(dataprefix,e,maxmobiter,mob)) as f:
        mobindx=int(0)
        hitcounts=hitcountdict[mob]
        #print hitcounts,'hitcounts'
        for indx, line in enumerate(f):
            line=line.strip()
            if not line.strip():
                continue
            if line.startswith('Search has CONVERGED!'):
                continue
            data=line.split('\t')
            #print indx, mobindx, 'indx,mobindx'
            if indx < int(hitcounts[mobindx]):
                data.append(str(mobindx+1))
                f2.write('%s\n'%'\t'.join(data))
            else:
                mobindx=mobindx+1
                data.append(str(mobindx+1))
                f2.write('%s\n'%'\t'.join(data))


    
