"""This module contains functions for general bioinformatic tasks such as handling fasta files, running BLAST searches, and using the built-in subprocess module flexibly and safely"""


#HANDLING FASTA FILES

def extractfromfasta(fastafilepath):
    """extract FASTA headers ('comments') and sequences from FASTA file"""
    fileObj=open(fastafilepath)                                                                                                        

    comments=[]
    sequences=[]
    seq=[]

    for idx, line in enumerate(fileObj):
        data=line.strip()
        if data.count('>')>0 and idx==0:
            data=data.lstrip('>')
            comments.append(data)                                                                                                                      
        
        elif data.count('>')>0 and idx>0:
            sequences.append(''.join(seq))
            seq=[]
            data=data.lstrip('>')
            comments.append(data)
        
        else:                                        
            seq.append(data)
        
    sequences.append(''.join(seq))
    fileObj.close()
    print "finished extracting from fasta"
    return (comments, sequences)
    

def writeFastaSeqs(comments, sequences, fastafilepath, width=60):
    """from correspinding lists of comments and sequences, write FASTA file"""
    fileObj=open(fastafilepath, 'w')
    for i, seq in enumerate(sequences):
        seqlength=len(seq)
        numLines=1+(seqlength-1)//width
        seqLines=[seq[width*x:width*(x+1)] for x in range(numLines)]
        seq='\n'.join(seqLines)
        fileObj.write('>%s\n%s\n' % ((comments[i].strip()), (seq.strip())))

    fileObj.close()
    print "finished writing fasta"



#BLAST
    
def makeBLASTdb(fastafile, databasename, dbtype, parse_seqids=False): #dbtype can be 'nucl' or 'prot'
    """takes fastafile filepath, databasename filepath and dbtype args"""
    import subprocess
    if parse_seqids==False:
        cmdArgs=['makeblastdb', '-dbtype', dbtype, '-in', fastafile, '-out', databasename]
    else:
        cmdArgs=['makeblastdb', '-dbtype', dbtype, '-in', fastafile, '-out', databasename, '-parse_seqids']
    subprocess.call(cmdArgs)

    
def runpsiblast(query, database, writefilepath, iteration, evalue):
    import subprocess
    cmdArgs=['psiblast', '-query', query, '-db', database, '-num_iterations', '%s' %iteration, '-num_threads', '1', '-out', writefilepath, '-outfmt', '6', '-evalue', '%s' %evalue]
    subprocess.call(cmdArgs)


    
#OTHER FUNCTIONS

def unlist(listed):
    """takes a list and converts to comma-separated string"""
    unlisted=(",".join(a for a in listed))
    return unlisted


def getcol(filepath, col, strip=True):
    """gets a column from a file and assigns to a variable"""
    collist=[]
    fileObj=open(filepath)
    for line in fileObj:
        if strip==True:
            line=line.strip()
        else:
            pass
        data=line.split('\t')
        collist.append(data[int(col)])
    return collist
    fileObj.close()

    
def runsubprocess(args, stderrpath=None, stdoutpath=None, writefile=None, shell=False):
    import subprocess,sys,thread
    """takes a subprocess argument list and runs Popen/communicate(); by default, both output and error are printed to screen; stderrpath and stdoutpath for saving output can be optionally set; a redirect can be optionally set (writefile argument); subthread error handling using keyboard interrupt; function can be used fruitfully since stdout is returned"""
    if shell==True:
        processname=args[0]
        processname=processname.split() 
    else:
        processname=args
    processname=("-".join(a for a in args))

    if stderrpath==None:
        pass
    else:
        if stderrpath.endswith('stderr.txt'): #ensure file ends with non-duplicated 'stderr.txt'
            stderrpath=str(stderrpath[:-10]).strip()
        stderrstrip=stderrpath.split('/')[-1]
        if stderrstrip=='': #there was nothing to strip after / i.e. was just /stderr.txt or stderr.txt
            pass
        else:
            stderrpath=stderrpath[:-(len(stderrstrip))]

        stderrpath=stderrpath+processname+'_'+stderrstrip+'stderr.txt'

    if stdoutpath==None:
        pass
    else:
        if stdoutpath.endswith('stdout.txt'): 
            stdoutpath=str(stdoutpath[:-10]).strip()
        stdoutstrip=stdoutpath.split('/')[-1]
        if stdoutstrip=='': 
            pass
        else:
            stdoutpath=stdoutpath[:-(len(stdoutstrip))]

        stdoutpath=stdoutpath+processname+'_'+stdoutstrip+'stdout.txt'

    
    print processname, "processname"
    
    try:

        if writefile==None:
            if shell==False:
                p=subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            if shell==True:
                p=subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
            stdout, stderr= p.communicate()
            print stdout, "stdout"
            print stderr, "stderr"
            
            if stdoutpath==None:
                pass
            else:
                with open(stdoutpath,'w') as stdoutfile:
                    stdoutfile.write(stdout)
                
            if stderrpath==None:
                pass
            else:
                with open(stderrpath,'w') as stderrfile:
                    stderrfile.write(stderr)

        else:
            with open(writefile,'w') as stdout:
                if shell==False:
                    p=subprocess.Popen(args,stdout=stdout, stderr=subprocess.PIPE)
                if shell==True:
                    p=subprocess.Popen(args,stdout=stdout, stderr=subprocess.PIPE, shell=True)
                stdout, stderr= p.communicate()
                print stdout, "stdout"
                print stderr, "stderr"
                
                if stderrpath==None:
                    pass
                else:
                    with open(stderrpath,'w') as stderrfile:
                        stderrfile.write(stderr)
        
                
        if p.returncode==0:
            print processname, "code has run successfully"
        else:
            print "source code fail"
        
    except:
        print "runsubprocess function code error"
        sys.exit()
        
    if p.returncode!=0:
        if 'driverscript' in str(sys.argv[0]): #if error is triggered in the main thread (driverscript.py) invoke sys.exit(); otherwise (when error triggered in subthread) invoke keyboard interrrupt
            sys.exit()
        else:
            thread.interrupt_main() #Raise a KeyboardInterrupt exception in the main thread. A subthread can use this function to interrupt the main thread.
    else:
        return stdout
