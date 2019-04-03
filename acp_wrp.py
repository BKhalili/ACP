#*************************
#   Created on 30.04.2018
#   author: Bita Khalili
#*************************

import src.acp.acp_fn as acp_fn
import statistics
import pandas as pd
from optparse import OptionParser
import numpy as np
import sys
import os

def pseudospectrum(M,C,P,n):    #M is the initial matrix with feature vectore and C is the complete correlation matrix
    pseudospec_cr=[]
    for i in range(len(P)):
        cr_v=0.5*(C[P[i][0]]+C[P[i][1]])
        pseudospec_cr.append(cr_v)
    return(pseudospec_cr)

def compute_crscale(corr,n):
    iu1 = np.triu_indices(corr.shape[0], 1)
    z_scored=np.arctanh(corr[iu1])*np.sqrt(n-3)
    crscale=statistics.stdev(z_scored)
    #print("{0:.3f}".format(crscale))
    return("{0:.3f}".format(crscale))    

def main(input_data,outdir,numofps,dist,remNeigbPairsFlag):

    print("Flag for removing feature-feature pairs in a same neighboring distance=",remNeigbPairsFlag,"\nOff-diagonal distance=",dist)

    data = pd.read_csv(input_data,header=None)
    M = data.loc[1:,1:]
    M = M.values
    num_samples=M.shape[0]
    ppm = data.loc[0,:]
    ppm = ppm.drop(labels=0)
    ppm = ppm.values
    print('Data matrix loaded: '+\
              '{:d}'.format(M.shape[0])+'x'+'{:d}'.format(M.shape[1]))

    (C,P) = acp_fn.acp(M,ppm,dist,remNeigbPairsFlag) #Complete correlation matrix and the sorted indices of leading correlated features

    pseudospec_cr = pseudospectrum(M,C,P,num_samples)   #a list of pseusospecs each correspond to the respective pairs in P

    print("\n----------writing output files------------")
    #-------------Creating two files of pseoudospectrums from averaged correlation vectors and Fisher transformed of leading features
    output_cr=pd.DataFrame(pseudospec_cr)
    output_cr=output_cr.transpose()
    output_cr.insert(0,'shift',ppm)
    
    headers_cr=[]
    for i in range(output_cr.shape[1]):
        if(i==0):
            headers_cr.append('shift')
        else:
            headers_cr.append('cr/f'+'{:05}'.format(ppm[P[i-1,0]])+'f'+'{:05}'.format(ppm[P[i-1,1]]))

    fileName_cr=outdir+'ps.acp.'+items.split('.csv')[0]+'.pseudospectrum.tsv' 
    output_cr=output_cr.iloc[:,0:(numofps+1)]
    output_cr.to_csv(fileName_cr,index=False,header=headers_cr[0:(numofps+1)],sep='\t')

    #-------------Creating a parameters in and descreption files
    DescFile=open(outdir+'description.tsv','w')
    del headers_cr[0]
    tags=[x.replace('cr/','') for x in headers_cr]
    for i in range(len(tags)):
        DescFile.write(tags[i]+'\t')
        for j in range(P.shape[1]):
            if ((ppm[P[i,j]]*1e4%100)<50):
                DescFile.write('f'+str(round(ppm[P[i,j]],2))+'p')
            elif ((ppm[P[i,j]]*1e4%100)>50):
                DescFile.write('f'+str(round(ppm[P[i,j]],2))+'m')
            if j==0:
                DescFile.write('-')
        DescFile.write('\n')
    DescFile.close()

    #---------------producing a table with f-f labels and their leading corresponding correlation
    CorrV=np.ones((1,len(P)))
    for i in range(len(P)):
        CorrV[0,i]=C[P[i][0],P[i][1]]

    output=pd.DataFrame(CorrV)
    headers_corr=[x.replace('cr/','Corr/') for x in headers_cr]
    output=output.iloc[:,0:numofps]
    output.to_csv(outdir.rsplit('/',1)[0]+'/TableOfSortedHighlyCorrFeats.tsv',index=False,header=headers_corr[0:numofps],sep='\t')

    crscale=compute_crscale(C,num_samples)
    return (num_samples,crscale)

if __name__ == '__main__':

    print("----------Asigning parameters------------")
    usage = "usage: %prog [optional] arg"
    parser = OptionParser(usage)
    parser.add_option('-o','--rem',dest='remNeigbPairsFlag',default=True)
    parser.add_option('-d', '--dist',help='off-diagonal distance',dest='dist',type='float',default=0.1)
    parser.add_option('-n', '--numofps',help='number of ps',dest='numofps',type='int',default=179)

    (options, args) = parser.parse_args()

    source_dir=os.getcwd()
    infolder='/data/'

    if not os.path.exists(source_dir+infolder):
        print('Couldn not find the input folder!')
        sys.exit()
    for items in os.listdir(source_dir+infolder):
        if '.csv' in items:
            outdir='./data.out/ps.acp.'+items.split('.csv')[0]+'/'
            if not os.path.exists(outdir):
                os.makedirs(outdir)
            print(outdir)
            (numOfSamples,crScale)=main(source_dir+infolder+items,outdir,options.numofps,options.dist,options.remNeigbPairsFlag) 



