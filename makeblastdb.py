import sys
sys.path.append('./mobtypingmodules')
from mobmods import makeBLASTdb


fastafile="./plasmiddatabase/%s_plasmidproteins.fa" %sys.argv[1]  
databasename="./plasmiddatabase/%s_plasmidproteinsdb" %sys.argv[1]

makeBLASTdb(fastafile, databasename, 'prot', parse_seqids=True)


