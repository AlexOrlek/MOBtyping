import sys, subprocess, time
sys.path.append('./mobtypingmodules')
from mobmods import runsubprocess, getcol

dataprefix=str(sys.argv[1])

runsubprocess(['python', './plasmidtranslate.py', dataprefix])
runsubprocess(['python', './cdsfeatures.py', dataprefix])
runsubprocess(['python', './makeblastdb.py', dataprefix])


evalue=1 #this is the column index in mobtype_evalues.tsv from which evalues for evalues are extracted for each MOB query - default is column 1, but additional columns could be added to test other sets of evalue thresholds

mobtypes=getcol('./mobtyping_evalues.tsv', 0) #col0 is mobtypes
evalues=getcol('./mobtyping_evalues.tsv', evalue)

maxiter=int(sys.argv[2])
iterations=list(range(1,(maxiter+1)))

for iteration in iterations:
    ps=[]
    for (mob, e) in zip(mobtypes, evalues): #running through mob protein/evalue threshold pairs for a given iteration
        p = subprocess.Popen(['python', './psiblast.py', dataprefix, str(e), str(iteration), str(mob), '&'])
        ps.append(p)
    while True:  #once a for loop for a given iteration has completed, call indefinite while loop; break (if subprocesses are complete) or wait
        ps_status=[p.poll() for p in ps]
        if all([x is not None for x in ps_status]):
            print "completed iteration:", iteration
            break
        else:
            time.sleep(5)
            print "waiting for subprocesses to finish; iteration:", iteration


runsubprocess(['python', './additerationnumber.py', dataprefix, str(evalue), str(maxiter)])
runsubprocess(['python', './blastfilter.py', dataprefix, '0', '0', str(evalue), str(maxiter)]) #(0 0 means no pid/coverage cutoffs by default)    
runsubprocess(['python', './mobtyping.py', dataprefix, str(evalue), str(maxiter)])
