#!/data/SW/anaconda3/envs/xli_new/bin/python

import os,sys,operator, math, argparse
import gzip, sqlite3, glob, matplotlib
import numpy as np
import glob, itertools
from Bio import PDB
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Align.Applications import ClustalwCommandline
from Bio.SeqRecord import SeqRecord
from Bio import AlignIO
from Bio.Align import AlignInfo
from Bio.PDB import PDBIO
from Bio.PDB.Polypeptide import *
from difflib import SequenceMatcher
import timeit
from pylab import *
from matplotlib.ticker import*
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid import make_axes_locatable
from mpl_toolkits import *
from matplotlib import cm
from pylab import *
import pylab

#############################################################################################################################
#### XLinterpreter: a tool for fast structural annotation of restraints from structural proteomics (e.g. primarily MS/XL)####
########################################### Raimondi F. and Russell RB ######################################################

##############Requires modifications of the following modules:
#Bio/PDB/MMCIFParser.py
#Bio/SeqIO/PdbIO.py


AA=['K','R','H','E','D','N','Q','Y','S','T','C','A','G','P','F','I','L','V','M','W']


def similar(a, b):
    return SequenceMatcher(None, a, b).ratio()

class SqueezedNorm(matplotlib.colors.Normalize):
    def __init__(self, vmin=None, vmax=None, mid=0, s1=2, s2=2, clip=False):
        self.vmin = vmin # minimum value
        self.mid  = mid  # middle value
        self.vmax = vmax # maximum value
        self.s1=s1; self.s2=s2
        f = lambda x, zero,vmax,s: np.abs((x-zero)/(vmax-zero))**(1./s)*0.5
        self.g = lambda x, zero,vmin,vmax, s1,s2: f(x,zero,vmax,s1)*(x>=zero) - \
                                             f(x,zero,vmin,s2)*(x<zero)+0.5
        matplotlib.colors.Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        r = self.g(value, self.mid,self.vmin,self.vmax, self.s1,self.s2)
        return np.ma.masked_array(r)


def GetPdbRes(struct, ch):
  resmaps={}
  for m in struct:
    for c in m:
      if c.id == ch:
        ii=1
        for r in c:
          resmaps[ii]=r.id[1]
          ii+=1
  
  return resmaps
  
def GetID(aa):
  labels=""
  resid=""
  for c in aa:
    if c.isdigit() == False:
      labels+=c
    elif c.isdigit() == True:
      resid+=c

  if int(resid)>0:
    return int(resid)



def FindPos(aa, start, seq):
  idx=0
  for s in seq:
    if int(start) == int(aa):
      return idx
    #if s != "-" and s != ".":
    if s != "-" and s != ".":# and s != "X":
      start+=1
    idx+=1

def CheckXRes(pos, refseq):
  aac=0
  Xc=0
  outstr=""
  for aa in refseq.split("\n")[1]:
    if aa == "X":
        Xc += 1
    outstr+=aa
    if int(aac) == int(pos):
      #print (refseq)
      #print (pos,aac, Xc,outstr)
      return pos-Xc
    aac += 1

def GetPos(idx,start,seq,seqchain):
  pos=0
  for s in seq:
    if int(pos) == int(idx):
      if s != "-" and s != "." and s != "X":
###I need to check the preceding number of X position occuring before that positions in the full length sequence    
        return CheckXRes(start, seqchain),s
        
      else:
        return 0,""
        
    #if s != "-" and s != ".":
    if s != "-" and s != ".":# and s != "X":
      start+=1
    pos+=1


def MapVariants(tstart,tseq,rstart,rseq,start_count,simcutoff,seqchain):
  offset=5
  #simcutoff=0
  varidx=FindPos(start_count,tstart,tseq)
  if start_count == 0:
    fampos,res=GetPos(varidx-1,rstart,rseq,seqchain)
  else:
    fampos,res=GetPos(varidx,rstart,rseq,seqchain)
  #Checking on a sequence window equal to offset if the query and reference sequences are identical or not
  #Some tolerance is applied through the similar function to allow small differences in the sequence
  ###Remember to remove the similarity cutoff when performing remote homology searches with either PSI BLAST or HHLIBITS!!!
  #return fampos, res  
  if similar(tseq[varidx:varidx+offset], rseq[varidx:varidx+offset]) >= simcutoff and tseq[varidx:varidx+offset] != "" and rseq[varidx:varidx+offset] != "":# or tseq[varidx:varidx-offset] == rseq[varidx:varidx-offset]:
    #print tseq[varidx:varidx+offset], rseq[varidx:varidx+offset] 
    return fampos, res
  elif  similar(tseq[varidx:varidx-offset], rseq[varidx:varidx-offset]) >= simcutoff  and tseq[varidx:varidx-offset] != "" and rseq[varidx:varidx-offset] != "":
    #print tseq[varidx:varidx-offset], rseq[varidx:varidx-offset] 
    return fampos, res
  else:
    #print tseq[varidx:varidx+offset], rseq[varidx:varidx+offset]
    return 0,""



def calc_dist(resA, resB, atomtype):
  ###Check the distance between specified atom types from two residues
  for atm1 in resA:
    if (atm1.fullname).strip() == atomtype:
      #print "OK", resA, atm1.fullname, atm1.coord
      atm_sele1=atm1
      
  for atm2 in resB:
    if (atm2.fullname).strip() == atomtype:
      #print "OK",resB, atm2.fullname, atm2.coord
      atm_sele2=atm2
  diff  = atm_sele1.coord - atm_sele2.coord
  return np.sqrt(np.sum(diff * diff))


def ExtractFasta(db, idslist, fastaout):
  alias1={}
  alias2={}
  alias3={}
  fastaoutfile=open(fastaout,"w")

  conn = sqlite3.connect(db)

  c = conn.cursor()
  ###We need to retrieve all the pdb2uniprot matches for the input XL
  sql_query=""
  #sql_query = 'select * from swissprot where ' + ' or '.join(("id_ac = '" + str(n)+"'"+" or id_id = '" + str(n)+"'"+" or id_gn = '" + str(n)+"'" for n in idslist)) 
  for n in idslist:
    ###This is a workaround I implemented to process as much data as possible from XlinkAnalyzer
    ###There are a lot of identifiers screwed up, i.e. they are lower case with no correspondence in the Uniprot
    ###while instead Uniprot cotains the corresponding uppercase identifier.
    #sql_query = 'select * from uniprot where id_ac = '+"'"+ str(n)+"'"+" or id_id = '" + str(n)+"'"+" or id_gn = '" + str(n)+"'"
    sql_query = 'select * from uniprot where id_ac = '+"'"+ str(n)+"'"
    c.execute(sql_query)
    data=c.fetchall()
    if len(data) == 0:
      N=n.upper()
      #sql_query = 'select * from uniprot where id_ac = '+"'"+ str(N)+"'"+" or id_id = '" + str(N)+"'"+" or id_gn = '" + str(N)+"'"
      sql_query = 'select * from uniprot where id_ac = '+"'"+ str(N)+"'"
      #c.execute(sql_query):
      #data=c.fetchall()
    print (sql_query)
    
    for line in c.execute(sql_query):
      uniac=line[0]
      unid=line[1]
      gn=line[2]
      fasta=line[3]
      fastaoutfile.write(fasta)
      alias1[uniac]=unid
      alias2[uniac]=gn
      alias3[n]=uniac
      continue

  conn.close()
  
  return alias1, alias2, alias3


def GetBlastOut(infile, evalcutoff,nriter,overlap_flag):
  blastout=open(infile,"r")
  pdb_matches={}
  flag1=0
  flag2=0
  flag3=0
  start=0
  nrround=0
  protein_matches=0
  buffer_list={}
  for line in blastout:
    if line.find("Results from round")!= -1:

      start=1
      nrround=int((line.split()[3]).strip("\n"))
      continue
      
    if start == 1 and len(line.split()) > 0 and line.split()[0] == "Query=":
      flag1=1
      uniac=line.split()[1]
      unid=uniac.split("|")[2]
      uniac=uniac.split("|")[1]
      matches=0

    if start == 1 and line[0] == ">":
      #chain=line.split("_")[1]
      #pdbid=line.split("_")[0]
      ###Using rsplit to handle cases where we have "_" in the name of the pdb file
      chain=line.rsplit("_",1)[1]
      pdbid=line.rsplit("_",1)[0]
      chain=chain.split()[0]
      pdbid=(pdbid.lstrip(">")).lstrip()
      flag2=1
      flag3=0
      cont=0
      query_start=0
      query_seq=""
      query_stop=0
      sbjc_start=0
      sbjc_seq=""
      sbjc_stop=0
      f1=0
      f2=0
      matches+=1
      continue

    if start == 1 and flag1 == 1 and flag2 == 1:
      if len(line.split()) > 0 and line.split()[0] == "Score" and line.split()[1] == "=":
        Evalue=(line.split()[7]).strip(",")
        if Evalue[0] == "e":
          Evalue="1"+Evalue
        Evalue=float(Evalue)          
      if len(line.split()) > 0 and line.split()[0] == "Identities" and line.split()[1] == "=":
        identity=(line.split()[3]).strip(",")
        identity=identity.strip("(")
        identity=identity.strip(")")
        identity=identity.strip("%")
        if Evalue < evalcutoff:
          #print uniac, unid, line
          flag3=1
          #flag2=0
          #continue
        #continue
      #continue

    if start == 1 and  len(line.split()) > 0 and line.split()[0] == "Query" and flag3 == 1:
      #print line
      query_start=line.split()[1]
      query_seq=line.split()[2]
      query_stop=line.split()[3]
      f1=1

    if start == 1 and len(line.split()) > 0 and line.split()[0] == "Sbjct" and flag3 == 1:
      #print line
      sbjc_start=line.split()[1]
      sbjc_seq=line.split()[2]
      sbjc_stop=line.split()[3]
      f2=1

    if start == 1 and flag3 == 1 and f1 == 1 and f2 == 1:
      ###Using targ_range and ref_range to dynamically monitor whether a given sequence alignment is overlapping with previously aligned regions
      ###If so, discard them to avoid confusing sequence-structure assignments which more often happen with repeat motives
      ###Currently the ali_overlap flag also controls the usage of single or multi-chain matches of a protein query to a pdb.
      ###This check shouldn't happen here, so FIX IT!!!!
      f1=0
      f2=0
      #print "Save stuff!"
      if pdbid not in buffer_list:
        buffer_list[pdbid]={}
        buffer_list[pdbid][uniac]=[]
        buffer_list[pdbid][uniac].append([Evalue,int(query_start), int(query_stop),query_seq,int(sbjc_start), int(sbjc_stop), sbjc_seq, chain, identity])
        targ_range=[]
        for ii in range(int(query_start), int(query_stop)):
          targ_range.append(ii)
        ref_range=[]
        for ii in range(int(sbjc_start), int(sbjc_stop)):
          ref_range.append(ii)
      elif pdbid in buffer_list and uniac not in buffer_list[pdbid]:
        buffer_list[pdbid][uniac]=[]
        buffer_list[pdbid][uniac].append([Evalue,int(query_start), int(query_stop),query_seq,int(sbjc_start), int(sbjc_stop), sbjc_seq, chain, identity])
        targ_range=[]
        for ii in range(int(query_start), int(query_stop)):
          targ_range.append(ii)
        ref_range=[]
        for ii in range(int(sbjc_start), int(sbjc_stop)):
          ref_range.append(ii)
      else:
        ###I probably have to change here to make overlap check chain specific
        if overlap_flag == 0:
          if int(query_start) not in targ_range and int(query_stop) not in targ_range and int(sbjc_start) not in ref_range and int(sbjc_stop) not in ref_range:
            buffer_list[pdbid][uniac].append([Evalue,int(query_start), int(query_stop),query_seq,int(sbjc_start), int(sbjc_stop), sbjc_seq, chain, identity])
            for ii in range(int(query_start), int(query_stop)):
              targ_range.append(ii)
            for ii in range(int(sbjc_start), int(sbjc_stop)):
              ref_range.append(ii)
        elif overlap_flag == 1:
          buffer_list[pdbid][uniac].append([Evalue,int(query_start), int(query_stop),query_seq,int(sbjc_start), int(sbjc_stop), sbjc_seq, chain, identity])
    if (start == 1 and line.find("Search has CONVERGED!") != -1) or (start == 1 and line.find("Effective search space used:") != -1 and nrround == nriter):
      for pid in buffer_list.keys():
        if pid not in pdb_matches:
          pdb_matches[pid]={}
        for uac in buffer_list[pid].keys():
          #print pdbid, uniac
          pdb_matches[pid][uac]=buffer_list[pid][uac]
      start=0
      #flag2=0
      buffer_list={}
      protein_matches+=1

  return pdb_matches

def aaLabel(resname):
  if resname in standard_aa_names:
    aaout=three_to_one(resname)
  else:
    aaout=resname

  return aaout

def GetResIdx(chain, residue):
  cc=0
  for rr in chain:
    if rr.id[1] == residue.id[1]:
      return cc  
    cc+=1

def GetHeader(headline):
  ###I also need to determine here the Field Separator type 
  count=0
  F1=0
  F2=0
  F3=0
  F4=0
  FS=""
  headline=headline.strip("\n")
  headline=headline.strip("\r")
  if headline.find("\t") != -1:
    FS="\t"
  elif headline.find(",") != -1:
    FS=","
  elif headline.find(";") != -1:
    FS=";"
  for field in headline.split(FS):
    if field.find("Protein1") != -1:
      F1=count
    if field.find("Protein2") != -1:
      F2=count
    if field.find("AbsPos1") != -1:
      F3=count
    if field.find("AbsPos2") != -1:
      F4=count
    count+=1
  return F1, F2, F3, F4, FS

def ReadRestr(restraint_list, LOGFILE):
  lc=0
  rest_list=[]
  prot_list=[]
  for l in open(restraint_list, "r"):
    if lc==0:
      f1,f2,f3,f4,fs=GetHeader(l)
      #LOGFILE.write("%s %s %s %s %s\n" % (f1,f2,f3,f4,fs))
    if len(l.split(fs)) >= 3 and lc>0 and l.split(fs)[f1] != "-" and l.split(fs)[f2] != "-": #let's not consider at this stage monolinks
      if l.split(fs)[f1].find("|") != -1:
        prot1=(l.split(fs)[f1]).split("|")[1]
      else:
        prot1=(l.split(fs)[f1])
      if l.split(fs)[f2].find("|") != -1:      
        prot2=(l.split(fs)[f2]).split("|")[1]
      else:
        prot2=(l.split(fs)[f2])
      pos1=l.split(fs)[f3]
      pos2=l.split(fs)[f4]
      if prot1 not in prot_list:
        prot_list.append(prot1)
      if prot2 not in prot_list:
        prot_list.append(prot2)
      rest_list.append([prot1, prot2, pos1, pos2])
    if l[0] != "#":
      lc+=1
  return rest_list, prot_list


def XLScore(distances, threshold, total, ac2gn, ac2xl):
  interfaces={}
  max_sc=len(distances.keys())
  tot_max_sc=total
  d_list=[]
  ###Getting the maximum distance for normalization of the penalty term
  for d in distances.keys():
    for val in distances[d]:
      d_list.append(val)
  max_d=max(d_list)
  sc=0
  count=0
  tot_xl={}
  ####Creating subset of unique combintations of chains and proteins
  for d in distances.keys():
    dis=distances[d]
    p1=d.split("|")[0]
    #c1=p1.split("/")[0]
    p1=p1.split("/")[0]
    p2=d.split("|")[1]
    #c2=p2.split("/")[0]
    p2=p2.split("/")[0]
    gn1=""
    gn2=""
    if p1 in ac2gn:
      gn1=ac2gn[p1]
    if p2 in ac2gn:
      gn2=ac2gn[p2]
    intf=gn1+"|"+gn2
    invintf=gn2+"|"+gn1
    if min(dis) <= threshold:
      part_sc=1
      count+=1
      if intf not in interfaces and invintf not in interfaces:
        interfaces[intf]=1
      elif intf in interfaces:
        interfaces[intf]+=1
      elif invintf in interfaces:
        interfaces[invintf]+=1
      
    else:
      dev=min(dis)-threshold
      #print dev, max_d, threshold
      part_sc=1-(float(dev)/float((max_d-threshold)))
    sc+=part_sc
    
    if p1 in ac2xl:
      if gn1 not in tot_xl:
        tot_xl[gn1]=ac2xl[p1]
    if p2 in ac2xl:
      if gn2 not in tot_xl:
        tot_xl[gn2]=ac2xl[p2]
    
  strout=""
  for intf in interfaces.keys():
    gn1=(intf.split("|")[0])
    gn2=(intf.split("|")[1])
    totxl=[]
    totxl=len(set(tot_xl[gn1]+tot_xl[gn2]))
    fraction=float(interfaces[intf])/float(totxl)
    strout+=intf+":("+str(interfaces[intf])+":"+str(totxl)+"=%3.2f)," % (fraction)

  return count, float(sc)/float(max_sc), float(sc)/float(tot_max_sc), interfaces, strout


def GenPym(xlinfo,outname,inpdb):
  for pdb in xlinfo.keys():
    print (outname, pdb)
    pymfile=open("%s_%s.pml" % (outname,pdb),"w")
    if inpdb != None:
      #os.popen("cp /net/home.isilon/ds-russell/pdb/%s/pdb%s.ent.gz ." % (pdb[1:3], pdb))
      outpdb=inpdb
      pymfile.write("load %s\n" % (inpdb))    
    else:
      os.popen("cp /net/home.isilon/ds-russell/pdb/%s/pdb%s.ent.gz ." % (pdb[1:3], pdb))
      #pymfile.write("load pdb%s.ent.gz\n" % (pdb))
      outpdb=pdb

    pymfile.write("load pdb%s.ent.gz\n" % (outpdb))

    outstring="""hide everything, all
    #set cartoon_smooth_loops=0
    #set cartoon_fancy_helices=1
    #set cartoon_cylindrical_helices = 1
    set transparency_mode=1
    set ray_shadows=0
    set ray_transparency_shadows=0
    #set ribbon_trace,1
    #set cartoon_transparency=0.5
    color white, all
    """
    
    pymfile.write(outstring)
    #pymfile.write("\nshow ribbon\n")
    pymfile.write("\nshow cartoon, %s\n" % (outpdb[:-4]))
    
    cc=0
    for xl in xlinfo[pdb]:
      c1=xl[0]
      c2=xl[1]
      r1=xl[2]
      r2=xl[3]
      idx1=xl[4]
      idx2=xl[5]
      satisfaction=xl[6]
      if inpdb != None:
        sele1="/%s//%s/%s/CA" % (outpdb[:-4],c1,idx1)
        sele2="/%s//%s/%s/CA" % (outpdb[:-4],c2,idx2)        
      else:
        sele1="/outpdb%s//%s/%s/CA" % (pdb,c1,idx1)
        sele2="/outpdb%s//%s/%s/CA" % (pdb,c2,idx2)
      
      if satisfaction == "satisfied":
        pymfile.write('cmd.distance("dist%s", "(%s)", "(%s)")\n' % (cc, sele1, sele2))
        pymfile.write('color green,dist%s\n' % (cc))
        pymfile.write('set label_outline_color, green, dist*\n')
        pymfile.write('set label_color, green, dist*\n')
        
      elif satisfaction == "unsatisfied":
        pymfile.write('cmd.distance("undist%s", "(%s)", "(%s)")\n' % (cc, sele1, sele2))
        pymfile.write('color red,undist%s\n' % (cc))
        pymfile.write('set label_outline_color, red, undist*\n')
        pymfile.write('set label_color, red, undist*\n')
      pymfile.write('show spheres, %s or %s\n' % (sele1,sele2))
      pymfile.write('label %s, "%s"\n' % (sele1, r1))
      pymfile.write('label %s, "%s"\n' % (sele2, r2))
      cc+=1

    pymfile.write('util.cbc outpdb%s\n' % (outpdb))
 
    pymfile.write("select Armadillo, %s and chain A and resi 62-668\n" % (outpdb[:-4]))
    pymfile.write("select Ankyrin, %s and chain A and resi 669-863\n" % (outpdb[:-4]))
    pymfile.write("select LRR, %s and chain A and resi 943-1314\n" % (outpdb[:-4]))
    pymfile.write("select Roc, %s and chain A and resi 1315-1519\n" % (outpdb[:-4]))
    pymfile.write("select COR, %s and chain A and resi 1520-1848\n" % (outpdb[:-4]))
    pymfile.write("select Kinase, %s and chain A and resi 1879-2138\n" % (outpdb[:-4]))
    pymfile.write("select WD40, %s and chain A and resi 2139-2527\n" % (outpdb[:-4]))
   
    outstring="""set dash_gap, 0
    color salmon, kinase
    color gold, COR
    color radium, Roc
    color argon, LRR
    color firebrick, WD40
    color nitrogen, Ankyrin
    color purple, Armadillo
    set dash_width, 1
    set float_labels, on
    set ortho=on
    set depth_cue=off
    set antialias=1

    cmd.bg_color('white')
    set ray_shadows=0
    set ray_transparency_shadows=1

set_view (\
     0.828978658,   -0.337360770,    0.446071982,\
     0.558496535,    0.541509449,   -0.628369451,\
    -0.029566927,    0.770034611,    0.637319565,\
    -0.000509888,   -0.000513971, -346.543640137,\
   221.254165649,  160.886825562,  197.022186279,\
   276.140258789,  417.041168213,   20.000000000 )
    """
    pymfile.write(outstring)    

  return


def GenBubblePlot(xlinfo,xllist,outname,pdb_sorted, Cutoff, Tots, Sat, ParScore,TotScore, uniacref):
  SystemSize=2527
  #fig=plt.figure(1,figsize=(20,20))
  #fig.subplots_adjust(left=0.05,bottom=0.10,right=0.85,top=0.90)
  #ax=fig.add_subplot(1,1,1)
  fig = plt.figure(figsize=(20,20))
  ax = fig.add_axes([0.2, 0.2,0.6,0.6])

  ax2=ax.twiny()
  ax3=ax2.twinx()
  ax4=ax.twinx()
  ax5 = fig.add_axes([0.2, 0.82,0.483,0.025])
  ax6 = fig.add_axes([0.16, 0.2,0.025,0.6])
  #ax5 = fig.add_subplot(1,2,2)
  #ax6 = fig.add_subplot()
  #ax5=fig.subplot(1,0.5, 1)

  pdb_sorted=sorted(pdb_sorted.items(), key=operator.itemgetter(1), reverse=True)
  xl_sorted=sorted(xllist, key=lambda a: int((a.split("|")[0]).split("/")[1]))

###Filling the distance matrix
  x_pos1_all=[]
  y_pos1_all=[]
  x_pos1=[]
  y_pos1=[]
  DistVal=[]
  for pdb,score in pdb_sorted:
    for xl in xl_sorted:
      site1=int((xl.split("|")[0]).split("/")[1])
      site2=int((xl.split("|")[1]).split("/")[1])
      invxl=xl.split("|")[1]+"|"+xl.split("|")[0]
      uac1=(xl.split("|")[0]).split("/")[0]
      uac2=(xl.split("|")[1]).split("/")[0]
      if uac1 == uniacref and uac2 == uniacref:
        if xl in xlinfo[pdb]:
          D=min(xlinfo[pdb][xl])
        elif invxl in xlinfo[pdb]:
          D=min(xlinfo[pdb][invxl])

        if site1 < site2:
          x_pos1_all.append(site1)
          y_pos1_all.append(site2)
          DistVal.append(D)
          if site1 not in x_pos1:
            x_pos1.append(site1)
          if site2 not in y_pos1:
            y_pos1.append(site2)
        if site1 > site2:
          x_pos1_all.append(site2)
          y_pos1_all.append(site1)
          DistVal.append(D)
          if site2 not in x_pos1:
            x_pos1.append(site2)
          if site1 not in y_pos1:
            y_pos1.append(site1)

  ###Working on final matrix 
  xDim=SystemSize
  yDim=SystemSize
  x=arange(xDim)
  y=arange(yDim)

  scFac=3
  norm=SqueezedNorm(vmin=0, vmax=np.nanmax(DistVal), mid=Cutoff, s1=1.7, s2=4)
  sc=ax.scatter(x_pos1_all,y_pos1_all,s=100,alpha=0.5, c=DistVal, edgecolors=None, cmap='RdYlGn_r',norm=norm, linewidths=None)
  
  ax4.yaxis.set_visible(False)
  ax.set_xlim((-1,xDim+1))
  ax.set_ylim((-1,yDim+1))
  ax2.set_xlim((-1,xDim+1))
  ax3.set_ylim((-1,yDim+1))
  ax5.set_xlim((-1,xDim+1))
  ax6.set_ylim((-1,yDim+1))
  
  ###Drawing crossing lines
  ###This function will be used to highlight XL mapping on structural models.
  
  for i in range(len(x_pos1_all)):
    #print x_pos1_all[i], y_pos1_all[i]
    l1=Line2D([0,x_pos1_all[i]],[y_pos1_all[i], y_pos1_all[i]],linestyle=':',linewidth=0.05,c='lightgray')
    l2=Line2D([x_pos1_all[i],x_pos1_all[i]],[y_pos1_all[i],SystemSize],linestyle=':',linewidth=0.05,c='lightgray')
    ax.add_line(l1)
    ax.add_line(l2)
  
  ax.set_xticks(x_pos1)
  ax.set_yticks(y_pos1)
  ax.yaxis.set_ticks_position('left')
  ax.xaxis.set_ticks_position('top')
  ax.set_xticklabels(x_pos1,rotation='vertical', multialignment='center',size=5, linespacing=2.0)
  ax.set_yticklabels(y_pos1,rotation='horizontal', multialignment='center',size=5, linespacing=2.0)
  ax2.set_xticklabels([])
  ax3.set_yticklabels([])

  fig.colorbar(sc, ax=ax)
  ax.text(10, 10, "Nr. of satisfied XL=%s (out of %s; %3.2f; %3.2f)" % (Sat, Tots, ParScore,TotScore), fontsize=15)

####This section is tailored to LRRK2 domain architectured  
  domains_block=[[62,655],[688,863],[943,1325],[1330,1519],[1529,1856],[1884,2135],[2141,2525]]
  domains=["ARM", "ANK", "LRR", "ROC", "COR", "KIN", "WD40"]
  colors_block=["darkviolet", "blue", "lightskyblue", "green", "gold", "coral","maroon"]


  
  
  hc=1
  for ii in range(len(domains_block)):
    b=domains_block[ii]
    c=colors_block[ii]
    start=b[0]
    end=(b[1]-b[0])+2
    middle=float(start)+float(end/2)-60
    #middle=float(start)+float(end/2)
    patch_hor=matplotlib.patches.FancyBboxPatch((start,0), end, 30, boxstyle="round,pad=0.1", color=c)
    patch_ver=matplotlib.patches.FancyBboxPatch((0,start), 30, end, boxstyle="round,pad=0.1", color=c)
    ax5.text(middle, 0.4,"%s" %(domains[ii]),size=20)
    ax5.add_patch(patch_hor)
    ax6.text(0.4,middle, "%s" %(domains[ii]),size=20,rotation=90)
    ax6.add_patch(patch_ver)
    hc+=1

  ax5.set_yticks([])
  ax6.set_xticks([])
  ax5.set_yticklabels([])
  ax6.set_xticklabels([])

  ax5.tick_params(axis = "both", which = "both", bottom = False, top = False, left=False, right=False, labelbottom = False, labeltop = False, labelleft=False, labelright=False)
  ax6.tick_params(axis = "both", which = "both", bottom = False, top = False, left=False, right=False, labelbottom = False, labeltop = False, labelleft=False, labelright=False)

  savefig('%s_mat.png' % (outname))
  savefig('%s_mat.eps' % (outname))  

  return

def GenHeatmap(xlinfo,xllist,outname,pdb_sorted, Cutoff):
  #import matplotlib
  
  pdb_sorted=sorted(pdb_sorted.items(), key=operator.itemgetter(1), reverse=True)
  xl_sorted=sorted(xllist, key=lambda a: int((a.split("|")[0]).split("/")[1]))
  print (xl_sorted)
###Filling the distance matrix
  mat=[]
  for pdb,score in pdb_sorted:
    row=[]
    for xl in xl_sorted:
      invxl=xl.split("|")[1]+"|"+xl.split("|")[0]
      if xl in xlinfo[pdb]:
        row.append(min(xlinfo[pdb][xl]))
      elif  invxl in xlinfo[pdb]:
        row.append(min(xlinfo[pdb][invxl]))
      else:
        row.append(float('NaN'))
    mat.append(row)
  
  mat=np.array(mat, dtype=float)
  print (mat)
  print (mat.shape)
  print (np.nanmax(mat))
###Drawing the matrix
  xlab_sorted=[]
  for xl in xl_sorted:
    site1=(xl.split("|")[0]).split("/")[1]
    site2=(xl.split("|")[1]).split("/")[1]
    xlab_sorted.append(site1+":"+site2)
  fig, ax = plt.subplots()
  norm=SqueezedNorm(vmin=0, vmax=np.nanmax(mat), mid=Cutoff, s1=1.7, s2=4)
  img=ax.imshow(mat, cmap='RdYlGn_r', norm=norm, interpolation='none', aspect=2)
  ax.set_xticks(np.arange(len(xl_sorted)))
  ax.set_yticks(np.arange(len(xlinfo.keys())))
  ax.set_yticklabels([pdb for pdb, score in pdb_sorted])
  ax.set_xticklabels(xlab_sorted, rotation=80)
  ax.tick_params(axis='x', which='major', labelsize=8)
  ax.tick_params(axis='y', which='major', labelsize=10)
  ax.xaxis.grid(which='major', color="lightgray", linestyle="--", linewidth=0.25)
  ax.yaxis.grid(which='major', color="lightgray", linestyle="--", linewidth=0.25)
  fig.colorbar(img, ax=ax)
  savefig('%s_distance_matrix.eps' % (outname))
  

  return



def GetNonRedPDB(pdb_site_list):
  ###1)Check the pdb with the highest xl coverage
  ###2)If the coverage between 2 pdb is the same , then take the one with best Evalue from PSI-BLAST
  
  residue_list={}
  eval_list={}
  for pdb in pdb_site_list.keys():
    for match in pdb_site_list[pdb]:
      p1=match[2]
      p2=match[3]
      r1=match[4]
      r2=match[5]
      E1=match[10]
      E2=match[12]
      pair=str(p1)+"|"+str(p2)
      invpair=str(p2)+"|"+str(p1)
      
      setinfo=[[p1,r1,E1],[p2,r2,E2]]
      for si in setinfo:
        p=si[0]
        r=si[1]
        E=si[2]
        #print p
        if pdb not in residue_list:
          residue_list[pdb]={}
          residue_list[pdb][p]=[r]
          eval_list[pdb]={}
          eval_list[pdb][p]=[E]
        elif pdb in residue_list and p not in residue_list[pdb]:
          residue_list[pdb][p]=[r]
          eval_list[pdb][p]=[E]
        else:
          if r not in residue_list[pdb][p]:
            residue_list[pdb][p].append(r)
          if E not in eval_list[pdb][p]:
            eval_list[pdb][p].append(E)
  
  pdb_matches={}
  for pdb in residue_list.keys():
    totres=0
    for p in residue_list[pdb].keys():
      totres+=len(residue_list[pdb][p])   
    #sorted_pdb=sorted(residue_list[p].items(), key=lambda a: len(a[1]), reverse=True)
    print ("INITIAL LIST",pdb, len(residue_list[pdb].keys()), totres, residue_list[pdb])
    pdb_matches[pdb]=totres
  sorted_pdb_matches=sorted(pdb_matches.items(), key=operator.itemgetter(1), reverse=True)

  ###Sorting here according to Evalue, so to consider first matches with better E-value
  ###Remember that in some cases, there might be template matches with same nr of XL satisfied but different E-value. In that case, the one with best E value has to be considered. 
  
  final_pdb=[]
  final_uniac={}  
  for pdb, cc in sorted_pdb_matches:
    for uniac in residue_list[pdb].keys():
      if uniac not in final_uniac:
        final_uniac[uniac]=[]
        for res in residue_list[pdb][uniac]:
          final_uniac[uniac].append(res)
        if pdb not in final_pdb:
          final_pdb.append(pdb)
      elif uniac in final_uniac:
        for res in residue_list[pdb][uniac]:
          if res not in final_uniac[uniac]:
            final_uniac[uniac].append(res)
            if pdb not in final_pdb:
              final_pdb.append(pdb)

  for pdb in final_pdb:
    print ("FINAL LIST", pdb, residue_list[pdb])
  
  return final_pdb


def TopologyFinder(SITES):
  outconflist={}
  for PDB in SITES.keys():
    outconflist[PDB]=[]
    conf_dict={}
    for inte in SITES[PDB]:
      c1=inte[0]
      c2=inte[1]
      p1=inte[2]
      p2=inte[3]
      id1=c1+"_"+p1
      id2=c2+"_"+p2
      flag=0    
      ###Creating a configuration dictionary for each chain_protein pair
      if id1 not in conf_dict and id2 not in conf_dict:
        conf_dict[id1]=[id2]
      elif id1 in conf_dict and id2 not in conf_dict[id1]:
        conf_dict[id1].append(id2)
      elif id2 in conf_dict and id1 not in conf_dict[id2]:
        conf_dict[id2].append(id1)
    ###Generating a list of compatible configurations
    #conf_list=[]
    #for id1 in conf_dict.keys():
    #  conf_list.append([id1])
    #  flag=0
    #  for id2 in conf_dict[id1]:
    #    #flag=0
    #    for idx in conf_list[-1]:
    #      cx=idx.split("_")[0]
    #      px=idx.split("_")[1]
    #      c2=id2.split("_")[0]
    #      p2=id2.split("_")[1]
    #      if (px == p2 and cx != c2) or (px != p2 and cx == c2):
    #        flag=1
    #    if flag == 0:
    #      conf_list[-1].append(id2)
    for id1 in conf_dict.keys():
      toremove=[]
      for id2 in conf_dict[id1]:
        c1=id1.split("_")[0]
        p1=id1.split("_")[1]
        c2=id2.split("_")[0]
        p2=id2.split("_")[1]
        if (p1 == p2 and c1 != c2) or (p1 != p2 and c1 == c2) and id2 not in toremove:
          toremove.append(id2)
      toappend=[x for x in conf_dict[id1] if x not in toremove]
      invtoappend=toappend.reverse()
      if toappend not in outconflist[PDB] and invtoappend not in outconflist[PDB]:
        outconflist[PDB].append(toappend)
    #if PDB == "6RIR" or PDB == "5SZJ":
    #  print (PDB, conf_dict)
    #  print (PDB, outconflist[PDB])

  return outconflist

def GetFasta(infile, outseqs):
  seq=""
  flag=0
  tracker=[]
  if infile[:-3] == ".gz":
    infile_stream=gzip.open(infile,"r")
  else:
    infile_stream=open(infile, "r")
  for l in infile_stream:
    if l[0] == ">":
      flag=0
      seq=""
      #gn=l.split()[0]
      seqid=l.split()[0].lstrip(">")
      outseqs[seqid]=l
    if l[0] != ">":
      outseqs[seqid]+=l.strip("\n")



#########################################################

def main():
  parser = argparse.ArgumentParser(description='Tool to map experimental restraints (usually MS/XL) on available 3D structures from the PDB | Basic usage: XL_interprets.py -restraint Input_XL_list.txt')
  parser.add_argument( '-input', action="store", dest='input', required=False, help='List of restraints', default="/net/home.isilon/ag-russell/bq_fraimondi/PROJECTs/XL_interpreter/XlinkAnalyzer/TFIIIC/xlinks/TFIIIC_allxlinks_without_tDNA.csv")
  parser.add_argument( '-inputpdb', action="store", dest='inputpdb', required=False, help='Input PDB: you must also provide a blast DB file with the fasta sequences from the PDB file', default=None)
  parser.add_argument( '-cutoff', action="store", dest='cutoff', required=False, help='Distance cutoff', type=float, default=35 )
  parser.add_argument( '-db', action="store", dest='db', required=False, help='SQL database with PDB2Uniprot blast alignments (default=/net/home.isilon/ds-russell/blastdb/pdbseq)', default="/net/home.isilon/ds-russell/blastdb/pdbseq")
  parser.add_argument( '-output', action="store", dest='output', required=False, help='Output filename', default="test")
  parser.add_argument( '-atom2monitor', action="store", dest='atom2monitor', required=False, help='Atom type whose distance will be monitored in order to model the restraints', default="CA")
  parser.add_argument( '-pymol', action="store_true", dest='pymol', required=False, help='Generates 3D representations of the mapped XL' )
  parser.add_argument( '-mode', action="store", dest='mode', required=False, help='Execution mode fast|thorough \n   -fast=only a non-redundant list of best pdb-query sequence matches is considered\n   -thorough=all pdb-query sequence matches are considered', default="fast" )
  parser.add_argument( '-blast_eval', action="store", dest='Eval', required=False, help='E-value cutoff for Blast results', type=float, default=2 )
  parser.add_argument( '-blast_v', action="store", dest='blast_v', required=False, help='Number of database sequences to show one-line descriptions', type=float, default=500)
  parser.add_argument( '-blast_b', action="store", dest='blast_b', required=False, help='Number of database sequence to show alignments for', type=float, default=250)
  parser.add_argument( '-psiblast_iterations', action="store", dest='psiblast_nriter', required=False, help='Number of PsiBlast iterations', type=int, default=2 )
  parser.add_argument( '-ali_overlap', action="store_true", dest='ali_overlap', required=False, help='Allows overlap between aligned sequence regions' )
  parser.add_argument( '-heatmap', action="store_true", dest='heatmap', required=False, help='Generates heatmap with XL-related distances' )
  parser.add_argument( '-bubbleplot', action="store", dest='bubbleplot', required=False, help='Generates bubble-plot with circles colored with XL-related distances', default="Q5S007" )
  parser.add_argument( '-assembly', action="store", dest='assembly', required=False, help='Type of PDB assembly:Asymmetric Unit|Bioassembly [AU|BioA]', default="BioA" ) 
  
  
  arguments = parser.parse_args()
  restr = arguments.input
  cutoff = float(arguments.cutoff)
  Eval = float(arguments.Eval)
  output = arguments.output
  outfile=open(output+"_XL_matches.txt", "w+t")
  atom2monitor=arguments.atom2monitor
  pymol=arguments.pymol
  heatmap=arguments.heatmap
  bubbleplot=arguments.bubbleplot
  blast_v=arguments.blast_v
  blast_b=arguments.blast_b
  psiblast_nriter=int(arguments.psiblast_nriter)
  mode=arguments.mode
  ali_overlap_flag=0
  if arguments.ali_overlap:
    ali_overlap_flag=1
  
  logfile=open(output+"_log.txt", "w+t")
  sysargv_str=""
  for sa in sys.argv:
    sysargv_str += sa + " "
  logfile.write("Command line: "+str(sysargv_str)+"\n")
  
  db=arguments.db
  inputpdb=arguments.inputpdb
  assembly=arguments.assembly 
  start_time=timeit.default_timer()
 
  print ("Got parameters") 
#### 1) Reading the experimental restraint input
  restr_list, id_list=ReadRestr(restr,logfile)
  print ("Read input XL") 
  
#### 2) Retrieving corresponding fasta sequences 
  fastaname=(output+"_inputseq.fasta")
  #pac2pid, pac2gn,origid2pac = ExtractFasta("/data/DB/uniprot/sprot_varsplic_seq_ids.db", id_list, fastaname)
  pac2pid, pac2gn,origid2pac = ExtractFasta("/data/DB/uniprot/uniprot_seq_ids_new.db", id_list, fastaname)
  print ("Retrieving Fasta sequences")
#### 3) PSIBLAST sequences on the PDB
  blastout_file="%s_blast.out" % (output)
#### Remember to uncomment this once the pevolution issue will be fixed and so the pdbseq
  #os.system("/net/home.isilon/ag-russell/install/CentOS-7.2.1511-x86_64/bin/blastall -p blastp -e %f -F F -m 0 -d /net/home.isilon/ds-russell/blastdb/pdbseq -i %s > %s" % (Eval, fastaname, blastout_file))
  os.system("/usr/local/bin/psiblast -evalue %f -db %s -query %s -num_iterations %d > %s" % (Eval, db, fastaname, psiblast_nriter, blastout_file))
###Reading PDB fasta sequences and storing them in a dictionary
  CASEQ={}
  if db[-3:] == ".fa":
    fastas=db
  else:
    fastas=db+".fa"
  GetFasta(fastas, CASEQ)
  print ("Blasting input sequences to PDB")
#### 4) Getting the Blast search output for sequences in the query against PDB reference sequences
  pdb2uniprot=GetBlastOut(blastout_file,Eval,psiblast_nriter,ali_overlap_flag)
  elapsed1 = timeit.default_timer() - start_time
#### 5) Determining how many XL sites are found within each PDB
  pdb_w_sites={}
  matched_xl=0
  #logfile.write(str(restr_list)+"\n")
  uniac2xl={}
  linked_chains={}
  unique_rstr=[]
  for xl in restr_list:
    if xl[0] in origid2pac:
      p1=origid2pac[xl[0]]
    else:
      p1=xl[0]
    if xl[1] in origid2pac:
      p2=origid2pac[xl[1]]
    else:
      p2=xl[1]
    r1=int(xl[2])
    r2=int(xl[3])
    xl=p1+"/"+xl[2]+"|"+p2+"/"+xl[3]
    invxl=p2+"/"+xl[3]+"|"+p1+"/"+xl[2]
    if xl not in unique_rstr and invxl not in unique_rstr:
      unique_rstr.append(xl)
    flag_match=0
    if p1 not in uniac2xl:
      uniac2xl[p1]=[xl]
    elif xl not in uniac2xl[p1] and invxl not in uniac2xl[p1]:
      uniac2xl[p1].append(xl)
    if p2 not in uniac2xl:
      uniac2xl[p2]=[xl]
    elif xl not in uniac2xl[p2] and invxl not in uniac2xl[p2]:
      uniac2xl[p2].append(xl)
    
    #print p1, p2, r1, r2
    for pdb in pdb2uniprot.keys():
      #print pdb
      if p1 in pdb2uniprot[pdb] and p2 in pdb2uniprot[pdb]:
        #print pdb, p1, p2
        ###Checking if the XL sites are within the blast match ranges
        flag1=0
        flag2=0
        m1=[]
        m2=[]
        ###First XL site
        for m in pdb2uniprot[pdb][p1]:
          if int(r1) in range(m[1], m[2]):
            #print p1, r1
            flag1=1
            m1.append(m)
            
        ###Second XL site
        for m in pdb2uniprot[pdb][p2]:
          if int(r2) in range(m[1], m[2]):
            #print p2, r2
            flag2=1
            m2.append(m)
            
        ### We don't know about the stechiometry of the XL, so we need to check all the possible chain combinations here
        if flag1 == 1 and flag2 == 1:
          for match1 in m1:
            for match2 in m2:
              E1=match1[0]
              I1=match1[8]
              E2=match2[0]
              I2=match2[8]
              chain1=match1[7]
              chain2=match2[7]
              seq_chain1=CASEQ[pdb+"_"+chain1]
              seq_chain2=CASEQ[pdb+"_"+chain2]
              idx1, aa1=MapVariants(match1[1],match1[3],match1[4],match1[6],int(r1),0.0,seq_chain1)
              idx2, aa2=MapVariants(match2[1],match2[3],match2[4],match2[6],int(r2),0.0,seq_chain2)
              print (r1,idx1,aa1,r2,idx2,aa2)
              flag_match=1
              #print "MATCHED_SITE:", pdb, match1[7], match2[7], p1, p2, r1, r2, idx1, idx2, aa1, aa2  
              if pdb not in pdb_w_sites:
                pdb_w_sites[pdb]=[[match1[7], match2[7], p1, p2, r1, r2, idx1, idx2, aa1, aa2, E1, I1, E2, I2]]
              else:
                pdb_w_sites[pdb].append([match1[7], match2[7], p1, p2, r1, r2, idx1, idx2, aa1, aa2, E1, I1, E2, I2])
      
  ####TopologyFinder: defining a list of plausible topologies, i.e. avoiding multi protein matches on the same chain
  pdb_sites_topology=TopologyFinder(pdb_w_sites)

###Make here pdb_w_sites non redundant if the mode flag is set to fast:
###With fast, we consider the pdb with highest number of xl mapped.
###Pdbs with xl already mapped in other pdb, are just skipped

  if mode == "fast":
    pdb_w_sites_final=GetNonRedPDB(pdb_w_sites)
  elif mode == "thorough":
    pdb_w_sites_final=pdb_w_sites.keys()

#### 6) Checking whether the XL sites found are compatible with XL theoretical distances (35A for DSS, 31A for DSG for example)
#### 7) Sorting pdb template according to restraint-satisfaction score

  outfile.write("MATCHED_SITE:\tPDB ID\tChain Site1\tChain Site2\tIdentity S1\tIdentity S2\tEval S1\tEval S2\tGene Name 1\t Gene Name 2\tUniprot Accession S1\tUniprot Accession S2\tS1 AA position (sequence)\tS2 AA position (sequence)\tS1 AA resseq PDB\tS2 AA resseq PDB\tS1 AA seqpos (structure sequence)\tS2 AA seqpos (structure sequence)\tS1 AA (structure)\tS2 AA (structure)\t%s-%s distance (Angstrom)\n" % (atom2monitor, atom2monitor))
  pdb_stat={}
  pdb4pym={}
  pdb_score={}
  #print (pdb_w_sites_final)
#### 6) Checking whether the XL sites found are compatible with XL theoretical distances (35A for DSS, 31A for DSG for example)
  
  XL_satisfied={}
  for pdb in pdb_w_sites_final:
###Calculating the distances on biological assembly only to speed up calculations
### Maybe I can speed up the process by first looping over PDB residues and internally checking
    if inputpdb != None:
      pdbpath=inputpdb
      Pdb=pdb
      pdb=pdb.lower()
    else:
      Pdb=pdb.upper()
      pdb=pdb.lower()
      pdbpath="/data/DB/pdb-mmCIF/%s/%s.cif.gz" % (pdb[1:3],pdb)

###Consider this only when we have local copies for biounit - need to check if these are provided for mmCIF format
    #if assembly == "BioA":
    #  pdbpath="/net/home.isilon/ds-russell/pdb-biounit/%s/%s.pdb*.gz" % (pdb[1:3],pdb)
    #else:
    #  pdbpath="/net/home.isilon/ds-russell/pdb/%s/pdb%s*gz" % (pdb[1:3],pdb)
    #print pdb, Pdb, pdbpath
    for pdbdir in glob.glob(pdbpath):
      #print pdbdir, pdb, Pdb
      #print pdb2uniprot[Pdb]
      if pdb not in pdb_stat:
        pdb_stat[pdb]={}
      if pdbdir[-3:] == ".gz":
        pdbpath=gzip.open(pdbdir, "r")
      else:
        #print pdbdir
        pdbpath=open(pdbdir, "r")
      if pdbdir[-3:] == "cif":
        parser = PDB.MMCIFParser(QUIET=True)
      elif pdbdir[-3:] == "pdb":
        parser = PDB.PDBParser(QUIET=True)
      struct = parser.get_structure('XXX',pdbpath)

      if os.path.isfile(pdbdir):
        print (pdb_w_sites[Pdb])
        for match in pdb_w_sites[Pdb]:
          #print (pdb_w_sites[Pdb])
          c1=match[0]
          c2=match[1]
          #print Pdb, c1, c2
          idx1=match[6]
          idx2=match[7]
          prot1=match[2]
          prot2=match[3]
          aa1=match[4]
          aa2=match[5]
          E1=match[10]
          I1=match[11]
          E2=match[12]
          I2=match[13]
          proto_id1=c1+"_"+prot1
          proto_id2=c2+"_"+prot2
          gn1=""
          gn2=""
          if prot1 in pac2gn:
            gn1=pac2gn[prot1]
          if prot2 in pac2gn:
            gn2=pac2gn[prot2]
          xl=prot1+"/"+str(aa1)+"|"+prot2+"/"+str(aa2)
          invxl=prot2+"/"+str(aa2)+"|"+prot1+"/"+str(aa1)
          chain1_check=0
          chain2_check=0
          topology_check=0
          topology_nr=0
          ####Topology check here to avoid conflicting matching of uniprot sequences on different pdb chains within the same configuration
          ###Eventually, generate MATCHED_SITE and SUMMARY output for each of these configurations
          for ii in range(len(pdb_sites_topology[Pdb])):
            topology=pdb_sites_topology[Pdb][ii]
            if proto_id1 in topology and proto_id2 in topology:
              topology_check=1
              topology_nr=ii+1
          model = struct[0]
          for chain in model:
            if chain.get_id() == c1:
              chain1_check = 1
            if chain.get_id() == c2:
              chain2_check = 1
          if chain1_check == 1 and chain2_check == 1 and topology_check == 1:# and final_check == 1:
            for ii, res1 in enumerate(model[c1]): 
              #print (ii, res1.get_full_id())
              if ii+1 == idx1:  ###Checking that pdb residue sequential number (i.e. from PDB Calpha sequence) matches
                for jj, res2 in enumerate(model[c2]):
                  if jj+1 == idx2:
                    distance = calc_dist(res1, res2, atom2monitor)
                    outfile.write("MATCHED_SITE:\t"+pdb+"\t"+c1+"\t"+c2+"\t"+str(I1)+"\t"+str(I2)+"\t"+str(E1)+"\t"+str(E2)+"\t"+gn1+"\t"+gn2+"\t"+str(match[2])+"\t"+str(match[3])+"\t"+str(match[4])+"\t"+str(match[5])+"\t"+str(res1.get_full_id()[3][1])+"\t"+str(res2.get_full_id()[3][1])+"\t"+str(idx1)+"\t"+str(idx2)+"\t"+str(match[8])+"\t"+str(match[9])+"\t"+str(distance)+"\n")
                    if pymol and distance < cutoff:
                      if pdb not in pdb4pym:
                        pdb4pym[pdb]=[[c1,c2,aa1,aa2,res1.get_full_id()[3][1],res2.get_full_id()[3][1], "satisfied",distance]]
                      else:
                        pdb4pym[pdb].append([c1,c2,aa1,aa2,res1.get_full_id()[3][1],res2.get_full_id()[3][1], "satisfied",distance])
                    else:
                      if pdb not in pdb4pym:
                        pdb4pym[pdb]=[[c1,c2,aa1,aa2,res1.get_full_id()[3][1],res2.get_full_id()[3][1], "unsatisfied",distance]]
                      else:
                        pdb4pym[pdb].append([c1,c2,aa1,aa2,res1.get_full_id()[3][1],res2.get_full_id()[3][1], "unsatisfied",distance])
                      
                    #pdb_stat[pdb].append(distance)
                    if xl not in XL_satisfied and invxl not in XL_satisfied:
                      XL_satisfied[xl]=[distance]
                    elif xl in XL_satisfied:
                      XL_satisfied[xl].append(distance)
                    elif xl != invxl and invxl in XL_satisfied:
                      XL_satisfied[invxl].append(distance)
                    if xl not in pdb_stat[pdb] and invxl not in pdb_stat[pdb]:
                      pdb_stat[pdb][xl]=[distance]
                    elif xl in pdb_stat[pdb]:
                      pdb_stat[pdb][xl].append(distance)
                    elif xl != invxl and invxl in pdb_stat[pdb]:
                      pdb_stat[pdb][invxl].append(distance)

          else:
            continue
          continue
  #print (pdb_stat)
  pdb_stat_sorted=sorted(pdb_stat.items(), key=lambda a: len(a[1]), reverse=True)
  outfile.write("SUMMARY:\tPDB ID\tNr. of total XL \tNr. of total XL mapped to the PDB\tNr. of XL mapped to the specific PDB\tNr. of XL distances satisfied\tPartial Score\tTotal Score\tRepresented Interactions(:Nr of Satisfied XL)\tXL sites mapped to PDB chains\n")
#  for pdb, dist in pdb_stat_sorted:
#    if len(dist) > 0:
#      satisfied, score=XLScore(dist, cutoff)
#      outfile.write("SUMMARY:\t"+pdb+"\t"+str(len(restr_list))+"\t"+str(len(XL_satisfied.keys()))+"\t"+str(len(dist))+"\t"+str(satisfied)+"\t"+str(score)+"\n")
  outali=open(output+".ali", "w+t")

  for pdb in pdb_stat.keys():
    dist=pdb_stat[pdb]
    uac_list=[]
    if len(dist) > 0:
      satisfied, partial_score, total_score, intfs, outstr=XLScore(dist, cutoff, len(unique_rstr), pac2gn, uniac2xl)
      outfile.write("SUMMARY:\t"+pdb+"\t"+str(len(unique_rstr))+"\t"+str(len(XL_satisfied.keys()))+"\t"+str(len(dist.keys()))+"\t"+str(satisfied)+"\t"+str(partial_score)+"\t"+str(total_score)+"\t"+outstr+"\t"+str(dist.keys())+"\n")
      pdb_score[pdb]=total_score
###Outputting here the multiple sequence alignments between template and query sequences

#[Evalue,int(query_start), int(query_stop),query_seq,int(sbj_start), int(sbj_stop), sbj_seq, chain, identity]
      #Pdb=pdb.upper()
      #Pdb=pdb
      for link in dist.keys():
        uac1=(link.split("|")[0]).split("/")[0]
        uac2=(link.split("|")[1]).split("/")[0]
        if uac1 not in uac_list:
          uac_list.append(uac1)
        if uac2 not in uac_list:
          uac_list.append(uac2)
      for uac in uac_list:
        print (pdb2uniprot[Pdb])
        evalue=pdb2uniprot[Pdb][uac][0][0]
        uac_start=pdb2uniprot[Pdb][uac][0][1]
        pdb_start=pdb2uniprot[Pdb][uac][0][4]
        head1=">P1;%s_%s\n" % (uac,evalue)
        head2=">P1;%s_%s\n" % (pdb,evalue)
        outstr1=""
        outstr2=""
        for lines in pdb2uniprot[Pdb][uac]:
          if lines[0] == evalue:
            uac_stop=lines[2]
            pdb_stop=lines[5]
            outstr1+=lines[3]
            outstr2+=lines[6]
            chain=lines[7]
          else:
            outstr1+="*\n"
            outstr2+="*\n\n"
            head1+="sequence:%s:%s:%s:%s:%s::::\n" % (uac, uac_start, chain, uac_stop, chain)
            head2+="structure:%s:%s:%s:%s:%s::::\n" % (pdb, pdb_start, chain, pdb_stop, chain)
            outali.write(head1+outstr1+head2+outstr2)
            evalue=lines[0]
            head1=">P1;%s_%s\n" % (uac,lines[0])
            head2=">P1;%s_%s\n" % (pdb,lines[0])
            outstr1=""
            outstr2=""
            uac_start=lines[1]
            pdb_start=lines[4]
            uac_stop=lines[2]
            pdb_stop=lines[5]
            outstr1+=lines[3]
            outstr2+=lines[6]
            chain=lines[7]
            
            
        outstr1+="*\n"
        outstr2+="*\n\n"
        head1+="sequence:%s:%s:%s:%s:%s::::\n" % (uac, uac_start, chain, uac_stop, chain)
        head2+="structure:%s:%s:%s:%s:%s::::\n" % (pdb, pdb_start, chain, pdb_stop, chain)
        outali.write(head1+outstr1+head2+outstr2)
      #print '3VKG', 'Q15276', pdb2uniprot['3VKG']['Q15276']


###Finding here the best combination of template structures: i.e. covering the highest fraction of sequence space and restraints
###In other words, I calculate the Jaccard score of XL sites in pairs of chains and sort the combinations with lowest violations 
###TO BE IMPLEMENTED IN THE FUTURE
### For the time being, an effective strategy is to define the individual best matching template structures (even if referred only to portions of the sequence/individual domains) and then using them as input for IMP
###IMP will then sort out the various pieces (even for the same protein)
  
  elapsed2 = timeit.default_timer() - start_time

  logfile.write("Total XL number: %d\n" % (len(restr_list))) 
  logfile.write("Nr. XL sites matched to the PDB: %d\n" % (len(XL_satisfied.keys()))) 
  xl_sat=0
  for k in XL_satisfied.keys():
    if min(XL_satisfied[k]) < cutoff:
      xl_sat+=1
  logfile.write("Nr. XL sites matched to the PDB and within distance threshold: %d\n" % ((xl_sat))) 
  logfile.write("Time to get fasta sequences and blast them on the PDB = %s\n" % (elapsed1))
  logfile.write("Time for the whole analysis = %s\n" % (elapsed2))

  #print pdb4pym
  if pymol:
    GenPym(pdb4pym, output,inputpdb)
  if heatmap:
    GenHeatmap(pdb_stat, XL_satisfied.keys(), output,pdb_score, cutoff)
  if bubbleplot:
    GenBubblePlot(pdb_stat, XL_satisfied.keys(), output,pdb_score, cutoff, len(dist.keys()), satisfied, partial_score, total_score, bubbleplot)

  sys.exit(0)
 
if __name__ == '__main__' : main()


