#*************************
#   Created on 30.04.2018
#   author: Bita Khalili
#*************************

import numpy as np
import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt
from copy import deepcopy
import sys


def acp(M,ppm,dist,remNeigbPairsFlag):
    print("\n----------Computing correlation matrix------------")
    M=correlation(M)
    print('Correlation matrix: '+\
          '{:d}'.format(M.shape[0])+'x'+'{:d}'.format(M.shape[1]))
        
    print("\n----------Computing off-diagonal distance------------")
    D=distanceMatrix(M.shape[0],ppm)
    # output=pd.DataFrame(D)
    # output.to_csv('distMatrix.csv',index=False,header=None)
    print('Distance matrix: '+\
        '{:d}'.format(D.shape[0])+'x'+'{:d}'.format(D.shape[1]))
          
    print("\n----------Removing feature-feature pairs in a same neighboring distance------------")
    # output=pd.DataFrame(M)
    # output.to_csv('CorrMatrix.csv',index=False,header=None)

    P=maxCorrPairs(deepcopy(M),D,dist,ppm)
    P=remove_Pairs(P,M,ppm,dist,remNeigbPairsFlag)
        
    return(M,P)

def correlation(M):
    C=np.corrcoef(M,rowvar=False)
    return(C)

def distanceMatrix(n,ppm):
    D = np.matmul(np.diag(ppm),np.ones((n,n)))   #a matrix with repeatition of ppms as columns
    D = abs(D-np.transpose(D))
    return (D)

def maxCorrPairs(CorrMat,DistMat,dist,ppm):
    P = np.zeros((CorrMat.shape[0] ,2),dtype=int)                  #list of pairs indices with max corr
    flag_feat =np.logical_not(DistMat<dist)
    CorrMat[np.invert(flag_feat)]=0
    max_array=np.max(CorrMat,1)
    P[:,1]=np.argmax(CorrMat,axis=0)
    P[:,0]=range(CorrMat.shape[0])
    P=P[np.argsort(max_array)][:]
    P=P[::-1][:]
    return (P)
    
def remove_Pairs(P,M,ppm,dist,remNeigbPairsFlag):    #removes repeated pairs or pairs that are whithin the diagonal distance
    removing_similar=np.ones(P.shape[0])
    for i in range(P.shape[0]):
        for j in range(i+1,P.shape[0]):
            if (P[i][0]==P[j][1] and P[i][1]==P[j][0]):       #removes repeated pairs
                removing_similar[j]=0

            #checking if any two pairs are within the dist threshold if yes the one with min correlation should be also eliminated
            elif remNeigbPairsFlag and (ppm[P[i][0]]<ppm[P[j][0]]+dist and ppm[P[i][0]]>ppm[P[j][0]]-dist and 
                ppm[P[i][1]]<ppm[P[j][1]]+dist and ppm[P[i][1]]>ppm[P[j][1]]-dist) or (ppm[P[i][0]]<ppm[P[j][1]]+dist and 
                ppm[P[i][0]]>ppm[P[j][1]]-dist and ppm[P[i][1]]<ppm[P[j][0]]+dist and ppm[P[i][1]]>ppm[P[j][0]]-dist):
                
                if M[P[i,0]][P[i,1]]>M[P[j,0]][P[j,1]]:
                    removing_similar[j]=0
                else:
                    removing_similar[i]=0
    removing_similar=removing_similar>0
    P=P[removing_similar]

    return (P)


