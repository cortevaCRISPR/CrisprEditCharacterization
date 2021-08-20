###Script by Ananta Acharya, Corteva Agriscience
###Process bam to characterise edits



import sys, os
import pysam
import re
from collections import Counter 
import itertools
from collections import namedtuple, defaultdict
AlignPos = namedtuple('AlignPos', ('qpos', 'qbase', 'rpos', 'rbase'))

def get_pairs(aln):

    qseq = aln.query_sequence
    pairs = (AlignPos(qpos=qp,
                      qbase=qseq[qp] if qp is not None else '-',
                      rpos=rp,
                      rbase=rb if rp is not None else '-'
                      )
             for qp, rp, rb in aln.get_aligned_pairs(with_seq=True)
             )
    return pairs 

def processbyInterval(bam,interval,key,overlaplen=10, NMpar=5):
    samfile = pysam.AlignmentFile(bam, "rb")
    #outsam=pysam.AlignmentFile(bam+".filtered.bam", "wb", template=samfile)
    piledict=dict()

    mystring = ''
    mylist = []
    start=int(interval[1])
    end=int(interval[2])+1
    poslist=[]
    readlist=[]
    reflist=[]
    readdict=defaultdict(list)
    posdict=defaultdict(list)
    header=samfile.header
    log=open(bam+".log", "a") 
    for read in samfile.fetch(reference=interval[0], start=start, end=end):
      
        name=read.query_name
      
        if read.is_proper_pair and not read.is_unmapped:   
            #print(read)
            tags=dict(read.tags)
            AS=tags['AS']
            
            if read.has_tag('XS'):
                XS=tags["XS"]
                # print(AS,XS,"AS", "XS")
                # print(read)

                if AS<XS*0.99:
                    log.writelines("%s discarded for %s as AS less than 0.99 XS\n" %(name,interval))
                    #print("skipping this")
                    #print(name,interval)
                    continue
                
            if AS>1:
                 
                if read.query_alignment_end>read.query_alignment_start:
                    mstrand="+"
                else:
                    mstrand="-"
                namestrand=name+mstrand
                cigar=read.cigarstring
                
                ap=read.get_aligned_pairs(with_seq=True)
                ol=read.get_overlap(start=start, end=end)
                rp1=read.get_reference_positions()
                rp2=read.get_blocks()
               
                MD=tags['MD']
                
                NM=tags['NM']
               

                if (AS>overlaplen-1 and NM<NMpar) or "D" in cigar:   
                  
                    if ol>overlaplen-1 or "D" in cigar:
                                               
                        pairs=get_pairs(read)
                        pos=start
                        newread=""
                        ref=""
                        qb=""
                        rb=""
                        rpprev=start
                        for pair in itertools.dropwhile(lambda x: x.rpos is None or x.qpos is None, pairs):

                            if pair.rpos!=None:
                                rpprev=pair.rpos
                                
                            else:
                                rpprev=rpprev
                            ##print(pair.rpos,rpprev,start,end)
                            try:
                            
                                #break
                                ##print(pair)
                                if pair.rpos!=None and (pair.rpos>=start and pair.rpos<=end+1):
                                    pos=pair.rpos
                                    qb=pair.qbase
                                    rb=pair.rbase
                                    posdict[namestrand,pos]=(qb.upper(),rb.upper())
                                    newread+=pair.qbase
                                    ref+=pair.rbase
                                else:
                                    if rpprev>=start and rpprev<=end+1:
                                        qb=qb+pair.qbase
                                        rb=rb
                                        posdict[namestrand,pos]=(qb.upper(),rb.upper())
                                        newread+=pair.qbase
                                        ref+=pair.rbase
                                    
                                    #with pysam.AlignmentFile(outbamgood, "ab", header=header) as outf:
                                    
                            except Exception as e:
                              
                                pass
                          
                        ind=[x for x, v in enumerate(ref) if v == '-']
                        ##print(ind)
                        l=len(ref)-1
                        ##print(l)
                        if 0 in ind or l in ind:
                            ##print("true")
                            ref="".join([char for idx, char in enumerate(ref) if idx not in ind])
                            newread="".join([char for idx, char in enumerate(newread) if idx not in ind])
                            ##print(ref,read)
                            ###input()

                        readlist.append(newread.upper())
                        reflist.append(ref.upper())
                        #poslist.append(posdict)
                        #print(interval)
                        readdict[(key,"\t".join(interval))].append((newread.upper(),ref.upper()))
                        
   
                    else:
                   
                        log.writelines("%s discarded for %s as Qual \n" %(name,interval))

                        
                else:
  
                    log.writelines("%s discarded for %s as Qual \n" %(name,interval))

        else:
         
            log.writelines("%s discarded for %s as Qual \n" %(name,interval))


    samfile.close()

    log.close()
    return posdict,readlist,reflist, readdict
    
def processbam(bam, bed, offbase, sample,outdir):
    writeline=""
    window=5
    with open(bed,"r") as bedin:
        baseposdict=defaultdict(list)
        seqdict=defaultdict(list)
        refdict=defaultdict(list)
        for line in bedin:
            if not (line.startswith("#") or line.startswith("type")) and len(line.strip())>0:
                linelist=line.strip().split()
                strand=linelist[4]

                start=linelist[1]
                end=linelist[2]
                chrm=linelist[0]

                group=linelist[5]
                name=linelist[3]
                interval=getDSB([chrm,start,end],group,strand,window)

                
                key=name
                # print(interval)
                seqwin="_".join([interval[0],str(interval[1]),str(interval[2])])
                #print(interval2)
                #input()
                edit=processbyInterval(bam, interval,key)

                posdict=edit[0]
                windowedits=window2table(edit[3])
                with open(outdir+"/"+sample+"_"+offbase+".hap.txt","a") as outf:   
                    outf.writelines(windowedits)
               


def window2table(readdict):
    writeline=""
    for interval,reads in readdict.items():
        # print(interval,reads)
        # print("reads2table")
        # input()
        observed=[pair[0] for pair in reads]
        ref=list(set([pair[1] for pair in reads]))[0]
        #print(ref)
        counted=Counter(observed)
        for hap, count in counted.items():
            writeline+="\t".join([interval[0],interval[1],ref,hap,str(count),"\n"])
            #print(interval,hap,str(count),ref[0])
            #input()
    return writeline

def getDSB(interval, guide,strand,window):
    chr,start,end=interval[0],int(interval[1]),int(interval[2])
    if "Cas12f" in guide:
        if strand=="+":
            DSB=start+27
        elif strand=="-":
            DSB=start-4
        else:
            print("Not valid strand")
            DSB=0
    elif "Cas9" in guide:
        if strand=="+":
            DSB=start+18
        elif strand=="-":
            DSB=start+5
        else:
            print("Not valid strand")
            DSB=0
    else:
        print("Not valid Guide")
        DSB=0
    return(chr,str(DSB-5),str(DSB+5))


def main():
    bam=sys.argv[1]
    bed=sys.argv[2]

    outdir=sys.argv[3]

    sample=os.path.basename(bam)
    offbase=os.path.basename(bed)
   
    logb=bam+".log"
    hapb=bam+"_"+offbase+".hap.txt"
    if os.path.exists(hapb):

      os.remove(hapb)
    if os.path.exists(logb):

      os.remove(logb)



    #outdir=os.path.dirname(os.path.abspath(bam))
    print(bam, bed,offbase, sample,outdir)
    samfile = pysam.AlignmentFile(bam, "rb")
   
    
    processbam(bam, bed, offbase, sample, outdir)


main()