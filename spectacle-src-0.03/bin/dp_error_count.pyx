import sys
import os
import numpy
from numpy import inf

def editDistDP(str1, str2, m, n):
    # Create a table to store results of subproblems
    dp = [[0 for x in range(m+1)] for x in range(n+1)]
    pointer = [[(0,0) for x in range(m+1)] for x in range(n+1)]
    # Fill d[][] in bottom up manner
    for i in range(n+1):
        for j in range(m+1):
            if ( i==0 and j==0):
                dp[i][j] = 0
            # If first string is empty,
            # insert all characters of second string
            elif i == 0:
                dp[i][j] = -j   # Min. operations = j
                pointer[i][j] = (0,j-1)
            # If second string is empty,
            # remove all characters of second string
            elif j == 0:
                dp[i][j] = 0   # Min. operations = i
                pointer[i][j] = (i-1,0)
            # If last characters are same, ignore last char
            # and recur for remaining string
            elif str1[j-1] == str2[i-1]:
                dp[i][j] = dp[i-1][j-1]+1
                pointer[i][j] = (i-1,j-1)
            # If last characters are different, consider all
            # possibilities and find minimum
            else:
                matMax = max(dp[i][j-1], dp[i-1][j], dp[i-1][j-1]) ;
                if matMax == dp[i][j-1]:
                   dp[i][j] = dp[i][j-1]-1
                   pointer[i][j] = (i, j-1)
                elif matMax == dp[i-1][j]:
                   dp[i][j] = dp[i-1][j]-1
                   pointer[i][j] = (i-1, j)
                else:
                   dp[i][j] = dp[i-1][j-1]-1 
                   pointer[i][j] = (i-1,j-1)
    
    return dp, pointer

def predecessorSetDP(corr, ref, orig, dp, m, n, alignRO, alignOR) :
    predecessorSet = [[[] for x in range(m+1)] for x in range(n+1)] 
    error = [[[] for x in range(m+1)] for x in range(n+1)]
    
    #  print "pred set init"
    #  print predecessorSet
    for i in range(n+1):
        for j in range(m+1):
            if (i==0):
                predecessorSet[i][j].append(('-','-'))
            elif (j==0):
                predecessorSet[i][j].append(('-','-'))
            else:
             if (corr[j-1] == ref[i-1]) and (dp[i][j] == dp[i-1][j-1]+1):
             #   print ref[i -1], corr[j-1], dp[i][j],dp[i-1][j-1]
                predecessorSet[i][j].append((i-1,j-1))
             elif (corr[j-1] != ref[i-1]) and (dp[i][j] == dp[i-1][j-1]-1):
              #  print ref[i -1], corr[j-1], dp[i][j],dp[i-1][j-1]
                predecessorSet[i][j].append((i-1,j-1))
             if (dp[i][j] == dp[i-1][j] - 1):
                predecessorSet[i][j].append((i-1,j)) 
             if (dp[i][j] == dp[i][j-1] - 1):
                predecessorSet[i][j].append((i,j-1))
    for i in range(n+1):
        for j in range(m+1):
            if (i==0):
               predecessorSet[i][j][0] += ('-',)
            elif (j==0):
               predecessorSet[i][j][0] += ('-',)
            else:
                #print i,j  
                for k in range(len(predecessorSet[i][j])):
                    if (predecessorSet[i][j][k][0]==i-1 and predecessorSet[i][j][k][1]==j-1):
                       #check if it is a new error : YYN
                       if  (corr[j-1] != ref[i-1]) :
                         if (alignRO[j-1] == '-') or (ref[alignRO[j-1]] == orig[j-1]) :
                            predecessorSet[i][j][k] += (1,)
                         else:
                            predecessorSet[i][j][k] += (0,)
                       else:
                         predecessorSet[i][j][k] += (0,)

                    elif (predecessorSet[i][j][k][0] == i-1 and predecessorSet[i][j][k][1] == j):
                       if(alignOR[i-1] != '-'):
                         predecessorSet[i][j][k] += (1,)
                       else: 
                         predecessorSet[i][j][k] += (0,)
                    elif (predecessorSet[i][j][k][0] == i and predecessorSet[i][j][k][1] == j-1):
                       if(alignRO[j-1] != '-'):
                         predecessorSet[i][j][k] += (1,)
                       else:
                         predecessorSet[i][j][k] += (0,)
    return predecessorSet

def minNewErrors(predecessorSet, m, n):
    newErrors = [[[] for x in range(m+1)] for x in range(n+1)] 
    for i in range(n+1):
        for j in range(m+1):
            if (i==0):
              newErrors[i][j] = 0 ;
            elif (j==0):
              newErrors[i][j] = 0 ;
            else:
              minError = 2000 ;
              for k in range(len(predecessorSet[i][j])):
                # print i, j, k 
                currError = predecessorSet[i][j][k][2] + newErrors[predecessorSet[i][j][k][0]][predecessorSet[i][j][k][1]]
                if (currError <= minError):
                    minError = currError;
              newErrors[i][j] = minError;
    return newErrors


def predecessorSetErrors(predecessorSet, newErrors, m, n, corr, orig, ref, alignRO, alignOR):
     
    predecessorSetErr = [[[] for x in range(m+1)] for x in range(n+1)] 
    maxError = [[[] for x in range(m+1)] for x in range(n+1)]
    
    #  print "pred set init"
    #  print predecessorSet
    for i in range(n+1):
        for j in range(m+1):
            if (i==0):
                predecessorSetErr[i][j].append(('-','-'))
            elif (j==0):
                predecessorSetErr[i][j].append(('-','-'))
            else:
             for k in predecessorSet[i][j]:
                 if (k[0] == i-1 and k[1] == j-1) and (newErrors[i][j] == newErrors[i-1][j-1] + k[2]):
                    predecessorSetErr[i][j].append((i-1,j-1))
                 if (k[0] == i and k[1] == j-1) and (newErrors[i][j] == newErrors[i][j-1] + k[2]):
                    predecessorSetErr[i][j].append((i,j-1)) 
                 if (k[0] == i-1 and k[1] == j) and (newErrors[i][j] == newErrors[i-1][j] + k[2]):
                    predecessorSetErr[i][j].append((i-1,j))
    
    for i in range(n+1):
        for j in range(m+1):
            if (i==0):
               predecessorSetErr[i][j][0] += ('-','-','-')
            elif (j==0):
               predecessorSetErr[i][j][0] += ('-','-','-')
            else: 
                for k in range(len(predecessorSetErr[i][j])):
                    
                    if (predecessorSetErr[i][j][k][0]==i-1 and predecessorSetErr[i][j][k][1]==j-1):
                       #check if is a corrected error : NYY
                       if ( (ref[i-1] == corr[j-1]) and (alignRO[j-1] != '-') and (ref[alignRO[j-1]] != orig[j-1]) ):
                         predecessorSetErr[i][j][k] += (1,)
                       
                       else:
                         predecessorSetErr[i][j][k] += (0,)
                         
                      #check if it is a new error at an incorrect location (modified base) : NYN
                       if ((corr[j-1] != ref[i-1]) and (alignRO[j-1]!='-') and(ref[alignRO[j-1]] != orig[j-1]) and (corr[j-1] != orig[j-1])):
                         predecessorSetErr[i][j][k] += (1,)
                       else:
                         predecessorSetErr[i][j][k] += (0,)

                        #check if it is an unmodified base : NNN
                       if ((alignRO[j-1]!='-') and (ref[alignRO[j-1]] != orig[j-1]) and (corr[j-1] == orig[j-1])):
                         predecessorSetErr[i][j][k] += (1,)
                       else:
                         predecessorSetErr[i][j][k] += (0,)

                    else:
                       predecessorSetErr[i][j][k] += (0,0,0)


    return predecessorSetErr
def correctedErrors(predecessorSetErr, m, n):
    corrErrors = [[[] for x in range(m+1)] for x in range(n+1)] 
    nynErrors = [[[] for x in range(m+1)] for x in range(n+1)] 
    nnnErrors = [[[] for x in range(m+1)] for x in range(n+1)] 

    for i in range(n+1):
        for j in range(m+1):
            if (i==0):
              corrErrors[i][j] = 0 ;
              nynErrors[i][j] = 0 ;
              nnnErrors[i][j] = 0;
            elif (j==0):
              corrErrors[i][j] = 0 ;
              nynErrors[i][j] = 0 ;
              nnnErrors[i][j] = 0;
            else:
              maxError = -1 ;
              for k in range(len(predecessorSetErr[i][j])):
                
                currError = predecessorSetErr[i][j][k][2] + corrErrors[predecessorSetErr[i][j][k][0]][predecessorSetErr[i][j][k][1]]
                if (currError >= maxError):
                    maxError = currError;
                    nynError_tmp = predecessorSetErr[i][j][k][3] + nynErrors[predecessorSetErr[i][j][k][0]][predecessorSetErr[i][j][k][1]]
                    nnnError_tmp = predecessorSetErr[i][j][k][4] + nnnErrors[predecessorSetErr[i][j][k][0]][predecessorSetErr[i][j][k][1]]

              corrErrors[i][j] = maxError;
              nynErrors[i][j] = nynError_tmp
              nnnErrors[i][j] = nnnError_tmp

    return corrErrors, nynErrors, nnnErrors
   

def backtraceDP(dp, pointer, m, n ) :
    backtrace = [];
    ind_j = m;
    maxVal = -inf
    for i in range(0,n+1):
        if (dp[i][m] >= maxVal):
           maxVal = dp[i][m]
           ind_i = i
    
    i = ind_i;
    j = ind_j;
    backtrace.append((i,j)) 
    while ( j !=0):
            id_i = pointer[i][j][0]
            id_j = pointer[i][j][1]
            backtrace.append((id_i,id_j));
            i = id_i
            j = id_j
    return backtrace[::-1]

def getAlignment(backtrace,n,m):
    alignRO = [ '-' for x in range(n)];
    alignOR = [ '-' for x in range(m)];
    for x in range(len(backtrace)-1):
        i = backtrace[x][0]
        #print i 
        j = backtrace[x][1]
        #print j 
        i_plus = backtrace[x+1][0]
        #print i_plus
        j_plus = backtrace[x+1][1]
        #print j_plus
        if ( i == i_plus - 1) and (j == j_plus - 1):
           alignRO[j] = i
           alignOR[i] = j
        elif ( i == i_plus - 1) and (j == j_plus):
           alignOR[i] = '-'
        else:
           alignRO[j] = '-'


#    for i in range(len(backVtrace)-1):

#        if ( (backtrace[i][0] == backtrace[i+1][0]-1) and (backtrace[i][1] == backtrace[i+1][1]-1) ): 
#             alignmentRO[i]= backtrace[i][1]
#             alignmentOR[j]= bactrace[i][0]
#        elif ( (backtrace[i][1]==backtrace[i+1][1]) and (backtrace[i][0] == backtrace[i+1][0]-1)):
#             alignment[i] = 'D'

#    alignment = [x for x in alignment if x != "D"]
    return alignRO, alignOR

##main
##process inputs and call functions

def main_func(char* ref, char*  orig, char* corrected):
    #ref = sys.argv[1]
    #orig = sys.argv[2]
    #corrected = sys.argv[3]

    #print "Reference sequence"
    #print ref
    #print len(ref)
    #print "Original Read"
    #print orig
    #print len(orig)
    #print "Corrected Read"
    #print corrected

    dpOrig, pointerOrig = editDistDP( orig, ref,  len(orig), len(ref))
    backtraceOrig = backtraceDP(dpOrig, pointerOrig, len(orig), len(ref))
    alignRO,alignOR = getAlignment(backtraceOrig,len(orig),len(ref))

    #dpOrig2, pointerOrig2 = editDistDP( ref, orig,  len(ref), len(orig))
    #backtraceOrig2 = backtraceDP(pointerOrig2, len(ref), len(orig))
    #alignOR = getAlignment(backtraceOrig2)


    dpCorr, pointerCorr = editDistDP( corrected, ref, len(corrected), len(ref))
    backtraceCorr = backtraceDP(dpCorr, pointerCorr, len(corrected), len(ref))
    #alignRC,alignCR = getAlignment(backtraceCorr,len(corrected),len(ref))

    #dpCorr2, pointerCorr2 = editDistDP( ref, corrected, len(ref), len(corrected))
    #backtraceCorr2 = backtraceDP(pointerCorr2, len(ref), len(corrected))
    #alignCR = getAlignment(backtraceCorr2)

  #  print str("Edit distance matrix") 
   # print dpOrig
    #print dpCorr

   # print str("Pointers") 
   # print backtraceOrig
   # print pointerOrig
    #print backtraceOrig

    #print str("alignRO[i] = j means O[i] is aligned to R[j]") 

    # print alignRO  
    #print str("alignOR[i] = j means R[i] is aligned to O[j]")
    #print alignOR

    #print alignRC
    #print alignCR
    #print str("alignRC[i] = j means C[i] is aligned to R[j]") 

    #print alignRC  
    #print str("alignCR[i] = j means R[i] is aligned to C[j]")
    #print alignCR


    predecessorSetCorr = predecessorSetDP(corrected, ref, orig, dpCorr, len(corrected), len(ref), alignRO, alignOR)
 #   print "predSetCorr"
 #   print predecessorSetCorr
    newErrors = minNewErrors(predecessorSetCorr, len(corrected), len(ref))
 #   print "newErrors"
 #   print newErrors
    # nynErrors = minNYNErrors(predecessorSetCorr, len(corrected), len(ref))
    # print "nynErrors"
    # print nynErrors
    predecessorSetErr = predecessorSetErrors(predecessorSetCorr, newErrors, len(corrected), len(ref), corrected, orig, ref, alignRO, alignOR);
    #print "new predecessorSet"
    #print predecessorSetErr 
    corrErrors, nynError, nnnError = correctedErrors(predecessorSetErr, len(corrected), len(ref))
    #print "corr errors set"
    #print corrErrors
    #print "nyn errors set"
    #print nynError
    #print "nnn errors set"
    #print nnnError
    numInsertions = 0 ;
    numDeletions = 0 ;
    i_start = backtraceCorr[-1][0]
    j_start = backtraceCorr[-1][1]
    for i in range(len(alignRO)):
        if(alignRO[i] == '-'):
           numInsertions = numInsertions + 1
    for i in range(len(alignOR)):
        if(alignOR[i] == '-'):
           numDeletions = numDeletions  + 1

    #print "i_start"
    #print i_start
    #print "j_start"
    #print j_start
    yny = len(corrected) - corrErrors[i_start][j_start]  - newErrors[i_start][j_start] - nynError[i_start][j_start] - nnnError[i_start][j_start]
        
    #print "New Errors" 
    #print newErrors[len(corrected)][len(ref)]
    #print "Corrected Errors"
    #print corrErrors[len(corrected)][len(ref)]
    return corrErrors[i_start][j_start], newErrors[i_start][j_start], nynError[i_start][j_start], nnnError[i_start][j_start], yny

#nyy, yyn, nyn,nnn, yny= main_func()

#print "yyn"
#print yyn
#print "yny"
#print yny
#print "nyy"
#print nyy
#print "nyn"
#print nyn
#print "nnn"
#print nnn



