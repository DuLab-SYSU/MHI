#!/usr/bin/python

import pandas as pd
import itertools
import time
import multiprocessing as mp
import sys, argparse

from indigo import *
indigo = Indigo()

def get_parser():
    description = 'Tanimoto 2D chemical similarity for drug pairs'
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-s', '--drug_smile_file', type=str, required=True, help='Drug SMILE')
    parser.add_argument('-o', '--out_file', type=str, required=True, help='Output file')
    return parser
    
def calculate_taniSim(d1, d2, drug_smile):
    s1 = drug_smile[d1]
    s2 = drug_smile[d2]
    m1 = indigo.loadMolecule(s1)
    m2 = indigo.loadMolecule(s2)
    
    # Calculate similarity between "similarity" fingerprints
    fp1 = m1.fingerprint('sim')
    fp2 = m2.fingerprint('sim')
    
    return indigo.similarity(fp1, fp2, 'tanimoto')
    
def drug_similarity(d1, d2, drug_smile, out_file):
    sim = calculate_taniSim(d1, d2, drug_smile)
    if sim >= 0.5:
        open(out_file, 'a').write("%s\t%s\t%f\n" % (d1, d2, sim))
        return
        
def get_drug_smile(drug_smile_file):
    smile_file = open(drug_smile_file)
    drug_smile = {}
    for row in itertools.islice(smile_file, 1, None):
        cols = row.strip().split('\t')
        drug_smile[cols[1]] = cols[3]
    return drug_smile
    
def run(args):
    drug_smile = get_drug_smile(args.drug_smile_file)

    start = time.time()
    combs = list(itertools.product([i for i in drug_smile.keys()], [i for i in drug_smile.keys()], [drug_smile], [args.out_file]))
    with mp.Pool(mp.cpu_count()) as pool:
        res = pool.starmap(drug_similarity, combs)
        pool.close()
        pool.join()
    end = time.time()
    print("The running time of this process is: {}h".format((end-start)/3600))
    

if __name__ == "__main__":
    run(get_parser().parse_args(sys.argv[1:]))
