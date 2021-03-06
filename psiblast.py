import sys
sys.path.append('./mobtypingmodules')
from mobmods import runpsiblast

#query mob proteins are PSI-BLASTed against the database of 6-frame translated plasmids

evalue=sys.argv[2]
iteration=sys.argv[3]
mobprotein=sys.argv[4]
threads=sys.argv[5]

query="./mobqueries/mobdownload_%s.fa" %mobprotein
database="./plasmiddatabase/%s_plasmidproteinsdb" %sys.argv[1]
writefilepath='./output_intermediate/%s_BLASTplasmidproteins_%s_%s_%s.tsv' %(sys.argv[1], evalue, iteration, mobprotein)

runpsiblast(query, database, writefilepath, iteration, evalue, num_threads=str(threads))
