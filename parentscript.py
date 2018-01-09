import sys, subprocess, time
sys.path.append('./mobtypingmodules')
from mobmods import runsubprocess, getcol

#e.g. python parentscript untypedplasmids 14 10    #i.e. run to 14 iterations, using 10 threads

dataprefix=str(sys.argv[1])

runsubprocess(['python', './plasmidtranslate.py', dataprefix])
runsubprocess(['python', './cdsfeatures.py', dataprefix])
runsubprocess(['python', './makeblastdb.py', dataprefix])


evalue=1 #this is the column index in mobtype_evalues.tsv from which evalues for evalues are extracted for each MOB query - default is column 1, but additional columns could be added to test other sets of evalue thresholds

mobtypes=getcol('./mobtyping_evalues.tsv', 0) #col0 is mobtypes
evalues=getcol('./mobtyping_evalues.tsv', evalue)

maxiter=int(sys.argv[2])
threads=int(sys.argv[3])
iterations=list(range(1,(maxiter+1)))

f=open('./output_intermediate/convergence_%s_%s_%s.tsv'%(dataprefix,evalue,maxiter),'w')
convergedproteins=[]
for iteration in iterations:
    for (mob, e) in zip(mobtypes, evalues): #running through mob protein/evalue threshold pairs for a given iteration
        if mob in convergedproteins:
            continue
        runsubprocess(['python', './psiblast.py', dataprefix, str(e), str(iteration), str(mob), str(threads)])
        checkconvergence=runsubprocess(['cat ./output_intermediate/%s_BLASTplasmidproteins_%s_%s_%s.tsv | tail'%(dataprefix,e,iteration,mob)],shell=True)
        if 'Search has CONVERGED!' in checkconvergence:
            print '%s search has converged at iteration %s'%(mob, iteration)
            f.write('%s\t%s\t%s\n'%(mob,iteration,'converged'))
            convergedproteins.append(mob)
        elif int(iteration)==int(maxiter):
            f.write('%s\t%s\t%s\n'%(mob,iteration,'notconverged'))
f.close()

runsubprocess(['python', './additerationnumber.py', dataprefix, str(evalue), str(maxiter)])
runsubprocess(['python', './blastfilter.py', dataprefix, '0', '0', str(evalue), str(maxiter)]) #(0 0 means no pid/coverage cutoffs by default)    
runsubprocess(['python', './mobtyping.py', dataprefix, str(evalue), str(maxiter)])
