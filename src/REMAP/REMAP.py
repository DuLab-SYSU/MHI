import numpy as np
import pandas as pd
import itertools
from scipy import sparse
from numba import jit

import sys
import argparse
import logging

def parse_args():
    parser=argparse.ArgumentParser(description="REMAP")
    parser.add_argument('--path',nargs='?', default='/datadir/',help='Input data path.')
    parser.add_argument('--R',nargs='?', default='chem_prot.csv')
    parser.add_argument('--chemsim',nargs='?', default='chem_chem.csv')
    parser.add_argument('--protsim',nargs='?', default='prot_prot.csv')
    parser.add_argument('--low_rank', type=int, default=100,
                        help='Low rank parameter.')
    parser.add_argument('--max_iter', type=int, default=100,
                        help='Maximum iteration.')
    parser.add_argument('--weight', type=float, default=0.1,
                        help='Global weight parameter.')
    parser.add_argument('--imp', type=float, default=0.1,
                        help='Global Imputation for unobserved associations.')
    parser.add_argument('--reg', type=float, default=0.1,
                        help='Regularization parameter.')
    parser.add_argument('--weight_chem', type=float, default=0.75,
                        help='Importance weight for chem-chem.')
    parser.add_argument('--weight_prot', type=float, default=0.25,
                        help='Importance weight for prot-prot.')
    parser.add_argument('--seed', type=int, default=1987,
                        help='Random seed.')
    parser.add_argument('--out_file', nargs='?', default='REMAP_pred.csv',
                        help='REMAP prediction')
    return parser.parse_args()

def REMAP(R,chem_chem,prot_prot,args):
    U,V=updateUV(R,chem_chem,prot_prot,args)
    return (U,V)

def updateUV(inputMatrix, chem_chem, prot_prot, args):
    lowrank=args.low_rank
    maxite=args.max_iter
    weight=args.weight
    imp=args.imp
    reg=args.reg
    pchem=args.weight_chem
    pprot=args.weight_prot
    m,n=inputMatrix.shape
    
    U0=np.asmatrix(np.random.rand(m,lowrank),dtype=np.float64)
    V0=np.asmatrix(np.random.rand(n,lowrank),dtype=np.float64)
    
    if (pchem==0) or (chem_chem is None):
        #gene sim not used
        updateU_=updateU_nosim
        Lu_plus=None
        Lu_minus=None
    else:
        Dm=sparse.diags(np.asarray(np.sum(chem_chem,1)).squeeze(),0,shape=(m,m))
        Lu=Dm-chem_chem
        Lu = sparse.csr_matrix(Lu)
        Lu_plus=(0.5*(np.abs(Lu)+Lu)).todense()
        Lu_minus=(0.5*(np.abs(Lu)-Lu)).todense()
        updateU_=updateU
        
    if (pprot==0) or (prot_prot is None):
        #cell sim not used
        updateV_=updateU_nosim
        Lv_plus=None
        Lv_minus=None
    else:
        Dn=sparse.diags(np.asarray(np.sum(prot_prot,1)).squeeze(),0,shape=(n,n))
        Lv=Dn-prot_prot
        Lv = sparse.csr_matrix(Lv)
        Lv_plus=(0.5*(np.abs(Lv)+Lv)).todense()
        Lv_minus=(0.5*(np.abs(Lv)-Lv)).todense()
        updateV_=updateU

    logging.debug("Preparing indicator matrix...")
    Nzs=np.zeros((m,n))
    row,col=np.nonzero(inputMatrix)
    for i in range(len(row)):
        Nzs[row[i],col[i]]=1
    logging.debug("updateUV started...")
    if inputMatrix.dtype != np.float64:
        inputMatrix=inputMatrix.astype(np.float64)
    if inputMatrix.__class__ in [sparse.coo_matrix,
                       sparse.csr_matrix,
                       sparse.csc_matrix,
                       sparse.bsr_matrix
                      ]:
        inputMatrix=inputMatrix.todense()
    for ite in range(0,maxite):
        UVT=get_UVT(Nzs,U0,V0)
        U0=fill_na(updateU_(inputMatrix,UVT,weight,imp,Lu_plus,Lu_minus,U0,V0,reg,pchem),val=0)
        V0=fill_na(updateV_(inputMatrix.T,UVT.T,weight,imp,Lv_plus,Lv_minus,V0,U0,reg,pprot),val=0)
        
    return U0, V0

def get_UVT(Nzs, U, V):
    return jit_mult(Nzs, jit_dot(U,V.T) )

def updateU(R,UVT,weight,imp,Lu_plus,Lu_minus,U0,V0,reg,importance):
    m,n=R.shape
    if imp>0:
        ua=(1-weight*imp)*jit_dot(R,V0)+jit_dot( (weight*imp*np.ones((m,n))), V0)+importance*jit_dot(Lu_minus,U0)
    else:
        ua=jit_dot(R,V0)+importance*jit_dot(Lu_minus,U0)
    ub=(1-weight)*jit_dot(UVT,V0)+ (weight* (jit_dot(U0, jit_dot(V0.T,V0) ))) + importance*jit_dot(Lu_plus,U0) + reg*U0
    u=jit_sqrt(jit_divide(ua,ub))
    U1=jit_mult(U0,u)
    return U1

def updateU_nosim(R,UVT,weight,imp,Lu_plus,Lu_minus,U0,V0,reg,importance):
    #high-efficiency version of updateU
    #in scRNA processing, cell-cell similarity may often be just identity matrix
    # and imputation score may not be used
    m,n=R.shape
    if imp>0:
        ua=(1-weight*imp)*jit_dot(R,V0)+jit_dot( (weight*imp*np.ones((m,n))), V0)
    else:
        ua=jit_dot(R,V0)
    ub=(1-weight)*jit_dot(UVT,V0)+ (weight* (jit_dot(U0, jit_dot(V0.T,V0) ))) + reg*U0
    u=jit_sqrt(jit_divide(ua,ub))
    U1=jit_mult(U0,u)
    return U1

def fill_na(A,val=0):
    A[np.isinf(A)]=val
    A[np.isnan(A)]=val
    return A

@jit(nopython=True)
def jit_dot(A,B):
    return np.dot(A,B)

@jit(nopython=True, parallel=True)
def jit_mult(A,B):
    return np.multiply(A,B)

@jit(nopython=True, parallel=True)
def jit_divide(A,B):
    return np.divide(A,B)

@jit(nopython=True, parallel=True)
def jit_sqrt(A):
    return np.sqrt(A)

@jit(nopython=True, parallel=True)
def jit_rownorm(A,ord=2):
    norm=[]
    for i in range(A.shape[0]):
        norm.append(np.linalg.norm(A[i,:],ord=ord))
    return np.array(norm,dtype=np.float32)

if __name__== '__main__':
    args=parse_args()
    lowrank=args.low_rank
    maxite=args.max_iter
    weight=args.weight
    imp=args.imp
    reg=args.reg

    protprot = pd.read_csv(args.path+args.protsim, header=None, delimiter='\t')
    protprot = protprot.sort_values(by = [0,1,2], ascending=[True, True, False])
    protprot = protprot.drop_duplicates(subset=[0, 1], keep='first', ignore_index=False)
    protprot = protprot.pivot(index = 0, columns = 1, values = 2)
    protprot = protprot.fillna(0)
    prots = protprot.index
    
    chemchem = pd.read_csv(args.path+args.chemsim, header=None, delimiter='\t')
    chemchem = chemchem.pivot(index = 0, columns = 1, values = 2)
    chemchem = chemchem.fillna(0)
    chems = chemchem.index
    
    combs = list(itertools.product([i for i in chemchem.index], [i for i in protprot.index]))
    chemprot = pd.DataFrame({'drug': [i[0] for i in combs], 
                             'gene': [i[1] for i in combs]})

    chemprot2 = pd.read_csv(args.path+args.R, delimiter='\t')
    chemprot2['score'] = [1] * chemprot2.shape[0]

    chemprot = pd.merge(chemprot, chemprot2, how='left',
                     on=['drug', 'gene'])
    chemprot = chemprot.pivot(index = 'drug', columns = 'gene', values = 'score')
    chemprot = chemprot.fillna(0)
    
    protprot = protprot.values
    chemchem = chemchem.values
    chemprot = chemprot.values
    
    
    U,V=REMAP(chemprot,chemchem,protprot,args) #run remap
    P=jit_dot(U,V.T) #prediction matrix
    P_df = pd.DataFrame(P)
    P_df.index = chems
    P_df.columns = prots
    P_df['drug'] = P_df.index
    P_df = P_df.melt(id_vars=['drug'], value_name='score', var_name='protein')
    
    #np.savetxt(args.out_file, P, delimiter = '\t') #save prediction matrix
    P_df.to_csv(path_or_buf=args.out_file, index = False, sep = '\t')
