#!/usr/bin/python

import os
import sys
import math
import subprocess
import logging
import argparse
global blast_path
blast_path='/home/tang/Downloads/ncbi-blast-2.13.0+/bin/'

def get_blast_sim(fasta_file, threshold, outfile, num_cores=24):
    fasta_file=str(fasta_file)
    threshold=float(threshold)
    outfile=str(outfile)

    #logging.info("Process starts. protinfo file: {}".format(protinfo_file))
    logging.info("Similarity threshold: {} (minimum similarity).".format(threshold))
    blastdbname=fasta_file+"_blastdb"
    blast_result_file="blastresult.dat"
    blastdbcommand=blast_path+"makeblastdb -in {} -dbtype prot -out {}".format(fasta_file,blastdbname)
    output=subprocess.check_output(['bash','-c',blastdbcommand])
    blastpcommand=blast_path+"blastp -query {} -db {} -evalue 1e-5 -outfmt 6 -num_threads {} > {}".format(fasta_file,
         blastdbname, 
         num_cores, 
         blast_result_file)
    output=subprocess.check_output(['bash','-c',blastpcommand])
    calc_sim(blast_result_file, threshold, outfile)

    delete_file(blastdbname+".phr")
    delete_file(blastdbname+".pin")
    delete_file(blastdbname+".psq")
    delete_file(blast_result_file)
    logging.info("Process Done. BLAST similarity file= {}".format(outfile))
    logging.info("Similarity threshold: {} (minimum similarity).".format(threshold))
    
def delete_file(filename):
    delcommand="rm {}".format(filename)
    output=subprocess.check_output(['bash','-c',delcommand])

def get_protinfo(filename):
    idx2acc={}
    acc2idx={}
    with open(filename) as inf:
        next(inf)
        for line in inf:
            line=line.strip().split("\t")
            idx=str(line[0]).strip()
            acc=str(line[1]).strip()
            idx2acc[idx]=acc
            acc2idx[acc]=idx
    return idx2acc, acc2idx

def calc_sim(datfile, threshold, outfile):
    #calculate tanimoto similarity
    #query index start from 1+query_idx_pad (for splitted file)
    #score lower than threshold discarded
    threshold=float(threshold)
    outfile=str(outfile)
    fout=open(outfile,"a+")
    id2selfscore={}
    try:
        idx2acc, acc2idx=get_protinfo(protinfo_file)    
    except BaseException as e:
        print("Error occurred:\n"+str(e))
        print("Please provide valid protein info file.")
    with open(datfile,"r") as inf:
        for line in inf:
            line=line.strip().split('\t')
            qid=line[0].strip().split("|")[1]
            tid=line[1].strip().split("|")[1]
            score=float(line[-1])
            #qaccession=str(qids[1])
            #taccession=str(tids[1])
            #qidx=acc2idx[qaccession]
            #tidx=acc2idx[taccession]
            if qid==tid:
                id2selfscore[qid]=score
    with open(datfile,"r") as inf:
        for line in inf:
            line=line.strip().split('\t')
            qid=line[0].strip().split("|")[1]
            tid=line[1].strip().split("|")[1]
            score=float(line[-1])
            #qaccession=str(qids[1])
            #taccession=str(tids[1])
            #qidx=acc2idx[qaccession]
            #tidx=acc2idx[taccession]
            try:
              simscore=float(score / id2selfscore[qid])
              if simscore>=threshold:
                  fout.write("{}\t{}\t{}\n".format(qid,tid,simscore))
            except:
              print("{} may be too short for similarity search.".format(qid))
    fout.close()
if __name__ == '__main__':
    parser = argparse.ArgumentParser("Run BLAST similarity calculation")
    #parser.add_argument('--protinfo', type=str, default='', help="protein information tsv file")
    parser.add_argument('--input', type=str, default='',help='Input FASTA file')
    parser.add_argument('--threshold', type=float, default=0.5,help='Minimum similarity to be included')
    parser.add_argument('--num_threads', type=int, default=24,help='n_cpu per job')
    parser.add_argument('--output', type=str, default='output.csv',help='Input FASTA file')
    parser.add_argument('--log', type=str, default='info',help='logging level')
    opt = parser.parse_args()
    FORMAT = '%(asctime)-15s %(message)s'
    logging.basicConfig(format=FORMAT, level=getattr(logging, opt.log.upper()))
    logging.info(opt)
    get_blast_sim(opt.input,opt.threshold,opt.output,num_cores=opt.num_threads)
