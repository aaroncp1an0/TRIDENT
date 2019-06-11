#!/usr/bin/python3

import numpy as np
#import str as str
import sys

import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib import colors

np.set_printoptions(suppress=True)


########################################
########################################
#start of module functions


#FUNCTIONS FOR TAKING TWO OVERLAPPING READS AND COMBINING THEIR Q-scores in the region of overlapp
#input: FASTQ-bowtie paired reads
#output: SINGLE FASTQ READ in bowtie output format

intab = b"ATCG"
outtab = b"TAGC"
tabCOMP = bytes.maketrans(intab, outtab) #makes sequence complement
    
intab = b"ATCGN"
outtab = b"01234"
tabHOT = bytes.maketrans(intab, outtab) #translate a sequence to one-hot encoding

#FUNCTION FOR COMPUTING SEQUENCE/Q SCORE MATRIX for two FASTQ files
def SeqScoreMAKER(filename1, filename2, filename3, qscoremesh=[20,25,30,33,35,37,40], cutoff=5, seqlength=2300, verbose=False):
    #FUNCTION FOR COMPUTING SEQUENCE/SCORE MATRIX for two FASTQ files
    #calls 1. open files
    #            file1=fq1.bowtie aligned; file2=fq2.bowtie aligned; file3=fq1 original non-aligned
    #      2. plan cycling and comparison of paired/unpaired
    #          -single/paired function call
    #          -sequence/score matrix increment
    #      3. next cycle
    
    print("warning: remember to edit code if there is substantial paired-end overlap. like for a 300x300 library!")
    ##################################
    #INITIATE VARIABLES
    
    #Qmap: 0-20 -> 0; 20:25(1) 27-40(2) 31-40(3) 33-40(4) 35-40(5) 37-40(6) 39,40(7)
    #Qscoremesh should be entered into here
    Qtab=bytes.maketrans(b"!\"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJ", b"000000000000000000001111122222334456788999")
                         # 0 0000000000000000000111112222233445678899
    #initializes 2300x4 matrix, where SSM[0] references position 1, with [A,T,C,G] indexes
    SSM=np.zeros((seqlength,5))
    
    #files are new!
    endfile1, endfile2 = False, False
    
    count=0
    
    ##################################
    #BEGIN MAIN FUNCTION EXECUTION
    
    #open filehooks
    f1, f2, f3 = readout(filename1, filename2, filename3)
    
    #print(f1, f2)
    FQstring1 = returnparsed(f1.readline(), Qtab)
    
    #read first line of each FASTQ bowtie paired file
    try: FQstring1 = returnparsed(f1.readline(), Qtab)
    except: print('empty files')
    try: FQstring2 = returnparsed(f2.readline(), Qtab)
    except: print('empty files')


    ###########################ALTERANTIVE FUNCTION (NON DUAL USAGE, NON OVERLAPPING)
    #this will only enter strings which are overlapping
    if True:
        print('USING NON-OVERLAP PROCESSING, largest possible dataset post-processing')
        Nover, Nexcept = 0, 0
        #
        #initialize a counter for making a histogram of how many mutations there are per string
        mutHIST=np.zeros(40)
        #
        for line in f1:
            ###updating SSM
            
            #print(FQstring1)
            
            try: SSM = updateSSmatrix(FQstring1[1][1],FQstring1[1][2],SSM,FQstring1[1][0],cutoff=cutoff)
            except: print('!! 1 instance of SSM update fail exception !!')
                
            #increment f1 to the next line in the file
            try: FQstring1 = returnparsed(f1.readline(), Qtab)
            except: endfile1=True
        
        for line in f2:
            ###updating SSM
            
            try: SSM = updateSSmatrix(FQstring2[1][1],FQstring2[1][2],SSM,FQstring2[1][0],cutoff=cutoff)
            except: print('!! 1 instance of SSM update fail exception !!')
                
            #increment f1 to the next line in the file
            try: FQstring2 = returnparsed(f2.readline(), Qtab)
            except: endfile2=True  
                
        if endfile1 and endfile2: 
            print('exiting')
            f1.close()
            f2.close()
            f3.close()
            return SSM
    #
    #in case we might want to record # mutations per read....
    #print( np.round( 100*mutHIST[:7]/np.sum(mutHIST), 2 ) )    
    #print("{0:.3f}".format( np.round( 100*mutHIST[:7]/np.sum(mutHIST), 2 ) ) )
    return SSM
      
    #GENERAL PROCESSOR - CAN BE PAIRED ENDS OR NOT
    ######################3
    for line in f3: #read a single line from f3, will compare f1 and f2 indexe to see if they match?
        #this function will now search for paired reads, or act on unpaired reads which lack quality pairingix
               
        #print('using NON-DualStrict processing')
        #process out the read number information using split
        FQstring3 = line.split()[0].split(':')[5:]
        
        #IF READ INDEX MATCHES MAIN FILE (F3) in either f1:reference or f2:reference, proceed to computations
        if FQstring1[0] == FQstring3 or FQstring2[0] == FQstring3:
            
            if verbose: print('updating')
            if verbose: print(FQstring1)
            
                    
            #does f1=f2?
            #yes? - execute overlap function and pairing and read next lines of f1,2,3
            if True:
                #execute single function, and update f1
                
                ###UPDATE COUNTER HERE
                try: SSM = updateSSmatrix(FQstring1[1][1],FQstring1[1][2],SSM,FQstring1[1][0],cutoff=cutoff)
                except: print('!! 1 instance of SSM update fail exception !!')
                
                #increment f1 to the next line in the file
                try: FQstring1 = returnparsed(f1.readline(), Qtab)
                except: endfile1=True
                
                ###execute single function, and update f2 and UPDATE COUNTER HERE
                try: SSM = updateSSmatrix(FQstring2[1][1],FQstring2[1][2],SSM,FQstring2[1][0],cutoff=cutoff)
                except: print('!! 1 instance of SSM update fail exception !!')
                
                #increment f2 to the next line in the file
                try: FQstring2 = returnparsed(f2.readline(), Qtab)
                except: endfile2=True                

        
            elif FQstring1[0] == FQstring3:
                #execute single function, and update f1
                
                 ###UPDATE COUNTER HERE
                try: SSM = updateSSmatrix(FQstring1[1][1],FQstring1[1][2],SSM,FQstring1[1][0],cutoff=cutoff)
                except: print('!! 1 instance of SSM update fail exception !!')
                
                #increment f1 to the next line in the file
                try: FQstring1 = returnparsed(f1.readline(), Qtab)
                except: endfile1=True
                
            elif FQstring2[0] == FQstring3:
                #execute single function, and update f2
                
                ###UPDATE COUNTER HERE
                try: SSM = updateSSmatrix(FQstring2[1][1],FQstring2[1][2],SSM,FQstring2[1][0],cutoff=cutoff)
                except: print('!! 1 instance of SSM update fail exception !!')
                
                #increment f2 to the next line in the file
                try: FQstring2 = returnparsed(f2.readline(), Qtab)
                except: endfile2=True                
                
            #just increment f3 at the begining of the next for loop
            else: pass
            
            
        #just increment f3 at the beginning of the next for loop
        else: pass
        
        if endfile1 and endfile2: 
            print('exiting')
            f1.close()
            f2.close()
            f3.close()
            return SSM
        
    return SSM
        
        
def SeqLength(filename1, filename2, filename3, qscoremesh=[20,25,30,33,35,37,40], cutoff=5, seqlength=2300, verbose=False):
    #FUNCTION FOR COMPUTING SEQUENCE/SCORE MATRIX for two FASTQ files
    #calls 1. open files
    #            file1=fq1.bowtie; file2=fq2.bowtie; file3=fq1 original
    #      2. plan cycling and comparison of paired/unpaired
    #          -single/paired function call
    #          -sequence/score matrix increment
    #      3. next cycle
    
    print("warning: remember to edit code if there is substantial paired-end overlap. like for a 300x300 library!")
    ##################################
    #INITIATE VARIABLES
    
    #Qmap: 0-20 -> 0; 20:25(1) 27-40(2) 31-40(3) 33-40(4) 35-40(5) 37-40(6) 39,40(7)
    #Qscoremesh should be entered into here
    Qtab=bytes.maketrans(b"!\"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHI", b"00000000000000000000111112222334455667788")
    
    #initializes 2300x4 matrix, where SSM[0] references position 1, with [A,T,C,G] indexes
    SSM=np.zeros((seqlength,5))
    
    Lrecord=np.zeros(5000)
    
    #files are new!
    endfile1, endfile2 = False, False
    
    count=0
    
    ##################################
    #BEGIN MAIN FUNCTION EXECUTION
    
    #open filehooks
    f1, f2, f3 = readout(filename1, filename2, filename3)
    
    #read first line of each FASTQ bowtie paired file
    try: FQstring1 = returnparsed(f1.readline(), Qtab)
    except: print('empty files')
    try: FQstring2 = returnparsed(f2.readline(), Qtab)
    except: print('empty files')

    for line in f3: #read a single line from f3, will compare f1 and f2 indexe to see if they match?
        #this function will now search for paired reads, or act on unpaired reads which lack quality pairingix
               
        #process out the read number information using split
        FQstring3 = line.split()[0].split(':')[5:]
        
        #IF READ INDEX MATCHES MAIN FILE (F3) in either f1:reference or f2:reference, proceed to computations
        if FQstring1[0] == FQstring3 or FQstring2[0] == FQstring3:
            
            if verbose: print('updating')
            if verbose: print(FQstring1)
            
            if FQstring1[0] == FQstring2[0]:
                #execute overlap function and pairing
                #update read f1 and f2 (f3 is taken care of by the main counter)
                
                #####NEED TO RETURN foverlap, in the same format as FQstring!
                
                ###UPDATE COUNTER HERE
                #SSM = updateSSmatrix(foverlap[1,1],foverlap[1,2],SSM,foverlap[1,0],cutoff)
                
                #all fixed
                
                x=int(np.abs(int(FQstring1[1][0])-int(FQstring2[1][0])))
                
                try: Lrecord[x]+=1 
                except: print('!! 1 instance of SSM dual update fail exception !!')
                    
                #lenthrec[]+=1
                
                #increment f1 and f2 to the next line in the file
                #exception method for boolean testing of file end
                try: FQstring1 = returnparsed(f1.readline(), Qtab)
                except: endfile1=True
                try: FQstring2 = returnparsed(f2.readline(), Qtab)
                except: endfile2=True
                    
            #does f1=f2?
            #yes? - execute overlap function and pairing and read next lines of f1,2,3
        
            elif FQstring1[0] == FQstring3:
                #execute single function, and update f1
                
                 ###UPDATE COUNTER HERE
                
                #increment f1 to the next line in the file
                try: FQstring1 = returnparsed(f1.readline(), Qtab)
                except: endfile1=True
                
            elif FQstring2[0] == FQstring3:
                #execute single function, and update f2
                
                ###UPDATE COUNTER HERE
                
                #increment f2 to the next line in the file
                try: FQstring2 = returnparsed(f2.readline(), Qtab)
                except: endfile2=True                
                
            #just increment f3 at the begining of the next for loop
            else: pass
            
            
        #just increment f3 at the beginning of the next for loop
        else: pass
        
        if endfile1 and endfile2: 
            print('exiting')
            f1.close()
            f2.close()
            f3.close()
            return Lrecord
        
    return Lrecord


#plan out end-logic and exception states
        
        
#PART OF THE PARSING FUNCTION      
def returnparsed(fastq, Qtab):
    y=fastq.split()
    #Bowtie already flips reads to fwd position
    #if y[2] == '-':
    #    y[5]=y[5].translate(tab)[::-1]
    #    y[6]=y[6][::-1]
    
    y[6]=y[6].translate(Qtab)
        
    #output: [['23140', '15797'],['1994','CATCAACACAGCAGATA','558888777555570']]
    return [y[0].split(':')[-2:], y[4:7]]     

            
#PART OF THE PARSING FUNCTION
#FUNCTION FOR OUTPUTING LINES FROM TWO FILES:
def readout(filename1, filename2, filename3):
    #input, filenames
    #output, fileholders
    
    #open files with 'holders'
    #options to extract filenames from shell command: filename1=str(sys.argv[1])
    try: filehold1=open(filename1,'r')
    except: print(filename1+str(' failed to open'))
    try: filehold2=open(filename2,'r')
    except: print(filename2+str(' failed to open'))
    try: filehold3=open(filename3,'r')
    except: print(filename3+str(' failed to open'))
    
    return filehold1, filehold2, filehold3

#PART OF THE PARSING FUNCTION
def updateSSmatrix(read,scores,SSM,pos=0,cutoff=3, tabhot=None, verbose=False): #FUNCTION takes read/scores and increments SSM-matrix at bases which meet the cutoff criteria
    #output, returns SSM matrix
    
    if verbose: print(read)
    if verbose: print(scores)
    
    readHOT = np.array(list(read.translate(tabHOT)), dtype=int)
    scoresH = np.array(list(scores), dtype=int)
    
    ####
    ####STRICT NEIGHBORING SCORES... means we require the +1 and -1 positions to have high Qscores
    if False:
        scoresH+=np.roll(scoresH, 1)
        scoresH+=np.roll(scoresH, -1)
        scoresM = np.zeros(len(scoresH)) #scoresM is the 0,1 multiplier
        scoresM[scoresH>cutoff*3]=1 #set scoresM=1 where scoresH (one-hot) meets the cutoff
    
    else:
        
        scoresM = np.zeros(len(scoresH)) #scoresM is the 0,1 multiplier
        scoresM[scoresH>cutoff]=1 #set scoresM=1 where scoresH (one-hot) meets the cutoff
    
    #pulls out index by position, then only those indexes referenced in A,T,C,G.
    #Pulled out indexes are then incremented by +1 according to scores, so only bases passing filter are accounted for
    SSM[np.arange(int(pos),int(pos)+len(read),1),readHOT]+=scoresM
    
    #returns incremented SSM matrix
    return SSM

#PART OF THE PARSING FUNCTION
def updateSSmatrixDUAL(read1,read2,scores1,scores2,SSM,pos1=0,pos2=0,cutoff=3, tabhot=None, verbose=False): #FUNCTION takes read/scores and increments SSM-matrix at bases which meet the cutoff criteria
    #output, returns SSM matrix
    
    #print('in SSDUAL update')
    
    if verbose: print(read1)
    if verbose: print(scores)
        
    pos1, pos2 = int(pos1), int(pos2)
        
    if abs(pos1-pos2)>=max(len(read1), len(read2)):
        SSM = updateSSmatrix(read1,scores1,SSM,pos1,cutoff,tabhot,verbose)
        SSM = updateSSmatrix(read2,scores2,SSM,pos2,cutoff,tabhot,verbose)
        return SSM
    
    
    #print(readHOT1, readHOT2)
    #print(scoresH1, scoresH2)
    
        
    elif pos1 < pos2 and abs(pos1-pos2)<len(read1):
        
        readHOT1 = np.array(list(read1.translate(tabHOT)), dtype=int)
        readHOT2 = np.array(list(read2.translate(tabHOT)), dtype=int)
    
        scoresH1 = np.array(list(scores1), dtype=int)
        scoresH2 = np.array(list(scores2), dtype=int)
        
        readNEW = np.zeros(len(read2)+pos2-pos1, dtype=int) #new read with terminal at pos1 going out to pos2+read2 length
        scoresN = np.zeros(len(read2)+pos2-pos1, dtype=int)
        
        #assign 1) read1 to the NEW sequences
        #assign 2) read2 to the NEW sequences: not the overlap region will have 2x the value as it normally should
        readNEW[:len(readHOT1)], scoresN[:len(readHOT1)] = readHOT1, scoresH1
        #print(readNEW)
        if verbose: print(readNEW, readHOT2, pos2, pos1, len(readHOT1), len(readNEW), len(readHOT2))
        
        readNEW[len(readHOT1):] += readHOT2[len(readHOT1)-pos2+pos1:]
        #print(readNEW)
        scoresN[pos2-pos1:] += scoresH2
        
        #if there is NOT agreement in the overlap between read1 and read2, set score to 0
        #print(len(scoresN[pos2:]))
        #print(readNEW)
        #print(scoresN)
        scoresN[pos2-pos1:][readNEW[pos2-pos1:]!=readHOT2] = 0
        
        #for the regions where there are agreement, set the score to 60 - we basically know this is the base
        #formally, we could add the two scores to get the actual number
        #scoresN[pos2-pos1:][readNEW[pos2-pos1:]==readHOT2] = 60
        
        #print(scoresN, readNEW)
        
        #0,1 multiplier scoresM; assign 1 depending on cutoff, pull out and add one to all proper positions
        scoresM = np.zeros(len(scoresN), dtype=int)
        scoresM[scoresN>cutoff]=1
        SSM[np.arange(int(pos1),int(pos1)+len(readNEW),1),readNEW]+=scoresM
        
        
    elif pos1 > pos2 and abs(pos1-pos2)<len(read2):
        
        SSM = updateSSmatrixDUAL(read2,read1,scores2,scores1,SSM,pos2,pos1,cutoff, tabhot=None, verbose=False)
        
    return SSM

############
############FILE HANDLING TOOLS

def load_files(filelist,pathA='./ip_M',pathB='',only4=True):
    #pathA/filelist#/pathB
    #ip_Ms4_S1
    
    M = []
    
    for name in filelist:
        M.append(np.load(pathA+name+pathB))
        
    return M

def moving_average(a, n=3) :
    #a[0,a[0,:]>.95]=0
    ret = np.cumsum(a, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return ret[n - 1:] / n

def M_norm(M):
    
    if len(np.shape(M)) > 1:
        norm=np.sum(M,axis=1,dtype=float)
        norm[norm==0]=1.0
        Mnorm=(1.0/norm)*M.T
        
    else:
        norm=np.sum(M)
        Mnorm=(1.0/norm)*M
    
    return Mnorm

def plot_i_average(Mnormarray,i=0,N=5,log=True,limitsy=[.000001,.01],limitsx=[750,1150],cutoff=.3):
    
    for array in Mnormarray:
        y=array[i,array[i,:]<cutoff]
        m=moving_average(y,n=N)
        #print(len(m.T))
        plt.plot(range(len(m)),m.T)
        
    plt.ylim(limitsy)
    plt.xlim(limitsx)
    
    if log: plt.yscale('log')
    else: pass
    
def plot_counts(Marray):
    #for plotting the total read count arrays, sanity checks
    
    try:
        for array in Marray:
            plt.plot(range(len(array)), np.sum(array, axis=1))
    except:
        plt.plot(range(len(Marray)), np.sum(Marray, axis=1))

        
def plot_NtoN(Mnormarray,baseStart='A',baseEnd='T',N=20,log=True,limitsy=[.000001,.01],limitsx=[750,1150],cutoff=.3, lims=True, ASScolor=True, colors=['gold'], linewidth=1, remove_out=True):
    
    lookup={'A':0,'T':1,'C':2,'G':3}
    start=lookup[baseStart]
    end=lookup[baseEnd]
    
    i=0
    
    for array in Mnormarray:
        #print(Mnormarray)
        xA = array[start,:]>.8
        x, y = np.arange(0,len(xA))[xA], array[end,xA]
        
        #print(len(m.T))
        
        ###OPTIONAL ANALYSIS FOR OUTLIER REMOVAL
        ###calculating mean and stddev, use these to remove points > 3stddev from mean (calculated w/o zero points)
        ###we then repeat this processing on the new data points in order to correct for re-calculation of the mean/stddev
        ###this processes only changes the moving average curve by removing extreme outliers; mostly correcting at the T-rich promoter
        
        if remove_out: #remove_outl
            
            ##PARAMETERS of 
            #number of standard deviations*N + mean to allow outliers in
            #reAveraging substitution window (what value to swap with)
            NUMstd, AVE = 3, 30
            
            #FIRST PASS OF REPLACEMENTS
            #calculate outlier index positions, outLys1
            mean = np.average(y[20:-20])
            stdD = np.std(y[y!=0][20:-20])
            outLys1 = y>NUMstd*stdD+mean
            
            #replace outLys1
            y[0:1-AVE][outLys1[0:1-AVE]]=moving_average(y, n=AVE)[outLys1[0:1-AVE]]
            
            #SECOND PASS OF REPLACEMENTS, accounting for new mean; stdDev
            #calculate outlier index positions, outLys2
            mean = np.average(y[20:-20])
            stdD = np.std(y[y!=0][20:-20])
            outLys2 = y>NUMstd*stdD+mean
            
            #replace outLys2
            y[0:1-AVE][outLys2[0:1-AVE]]=moving_average(y, n=AVE)[outLys2[0:1-AVE]]
            
            #re-replace outLys1, accounting for new averaging values
            y[0:1-AVE][outLys1[0:1-AVE]]=moving_average(y, n=AVE)[outLys1[0:1-AVE]]
            
            #re-replace outLys2, accounting for new averaging values
            y[0:1-AVE][outLys2[0:1-AVE]]=moving_average(y, n=AVE)[outLys2[0:1-AVE]]
         
        m = moving_average(y, n=N)
        
        #m = RunningMedian(y, M=N)
        
        #print(np.shape(x))
        if ASScolor: plt.plot(x[:len(m)],m, color=colors[i], linewidth=linewidth)
        else: plt.plot(x[:len(m)],m, linewidth=linewidth)
        
        #update color index
        i+=1
       
    if lims:
        plt.ylim(limitsy)
        
    plt.xlim(limitsx)
    
    if log: plt.yscale('log')
    else: pass
    
    return plt
    
    
def plot_ATCG_freq(Mnorm_single, base='A', log=True, area=.3, limitsy=[.0000001,.1],limitsx=[750,1150]):
    #input single base matrix array
    #output single plot w/ ATCG output mutations for a single starting base
    lookup = {'A':0,'T':1,'C':2,'G':3}
    i=lookup[base]    
    #np.max(Mnorm_single[:,:],axis=1)
    xA = Mnorm_single[i,:]>.93 #pull out all i values 
    plt.scatter(range(len(Mnorm_single.T)), Mnorm_single.T[:,0]*xA, s=area, color='red')
    plt.scatter(range(len(Mnorm_single.T)), Mnorm_single.T[:,1]*xA, s=area, color='green')
    plt.scatter(range(len(Mnorm_single.T)), Mnorm_single.T[:,2]*xA, s=area, color='blue')
    plt.scatter(range(len(Mnorm_single.T)), Mnorm_single.T[:,3]*xA, s=area, color='purple')
    
    plt.ylim(limitsy)
    plt.xlim(limitsx)
    if log: plt.yscale('log')

        
        
def runA_heatmap(Mctrlnorm, Msamplenorm, maskCT=False, maskTC=False, 
                 samplename='temp', filename='temp', sumX1X2=[900,1100], nozeros=True, 
                 SAVEFIG=False, normalize=True, directory='2017_analysis/'):
    
    #EXAMPLE RUN:
    #runA_heatmap(M1norm, M3norm, maskCT=False, maskTC=False,
    #        samplename='EXO1', filename='EXO1minusM2', sumX1X2=[900,1100], nozeros=True,
    #       SAVEFIG=True, directory='20170703_s4s10_analysis/')
    
    
    ##############################
    #calculate mutation counts
    try: freq_ctrl   = mutation_types(Mctrlnorm,x1=sumX1X2[0],x2=sumX1X2[1])
    except: print('fail ctrl mutation types')
    
    try: freq_sample = mutation_types(Msamplenorm,x1=sumX1X2[0],x2=sumX1X2[1])
    except: print('fail ctrl mtuation types')
        
    #print(freq_sample)
    
    if normalize: freq_difference = freq_sample - freq_ctrl
    else: freq_difference = freq_sample
    
    ##############################
    #zero out any negative values; these are not logical values to be realized bc/ of the noise threshold
    if nozeros: freq_difference[freq_difference<0] = 0.0
        
    if maskCT: freq_difference[9] =0
    if maskTC: freq_difference[12]=0
        
    ##############################
    #calculate total mutation rate
    
    #METHOD WHICH IGNORES BASE COMPOSITION OF SEQuENCE (each gets 1/4 of representation)
    
    #####METHOD WHICH REFLECTS BASE COMPOSITION OF SEQUENCE... not as general
    #xN1, xN2 = Mctrlnorm>.50, Msamplenorm>.50
    #totalR1 = np.max(Mctrlnorm.T,axis=1)[sumX1X2[0]:sumX1X2[1]] 
    #totalR2 = np.max(Msamplenorm.T,axis=1)[sumX1X2[0]:sumX1X2[1]]
    
    #if maskCT:
    #    totalR1 = np.max(Mctrlnorm.T,axis=1)[sumX1X2[0]:sumX1X2[1]]
    #if maskTC:
    #
    #totalR1, totalR2 = 1-np.average(totalR1), 1-np.average(totalR2)
    #
    #if normalize: totalrate = totalR2 - totalR1
    #else: totalrate = totalR2
    
    totalrate = calculate_rate(freq_difference, Msamplenorm, x1=sumX1X2[0], x2=sumX1X2[1], maskCT=maskCT, maskTC=maskTC)
    
    ##############################
    #Print title rate/file
    name = samplename + ' | LOG rate@ ' + str('{0:.2}'.format(np.log(totalrate)/np.log(10)))
    #name = samplename+' \ '+PrintRate
    
    ##############################
    #Normalize the frequencies
    freq_difference = M_norm(freq_difference)
    #print(freq_difference)
    
    try: plot_heatmap(mutation_values=freq_difference, circle=True, showtxt=True, title=name)
    except: print('missed')
    
    if SAVEFIG:
        if maskCT: filename='noCT_'+filename
        if maskTC: filename='noTC_'+filename
            
        try: plt.savefig(str(directory + filename + '_' + str(sumX1X2[0]) + '_' 
                        + str(sumX1X2[1]) + '.svg'), transparent=True)
        except: pass
    #plt.savefig()
    plt.show()
    
    
    
def calculate_rate(freq_difference, Mnorm, x1=750, x2=1150, maskCT=False, maskTC=False):
    
    #we first determine how many A,T,C,Gs are in the data:
    M1     = Mnorm[:,x1:x2]
    counts = M_norm(np.sum(M1[[0,1,2,3],:]>.5, axis=1))
    ATCGcounts = np.repeat(counts ,4)
    
    #print(freq_difference)
    
    if maskCT: ATCGcounts[9]=0
    if maskTC: ATCGcounts[12]=0
        
    #print(ATCGcounts)
    
    #take the weighted average dot product
    return ATCGcounts.dot(freq_difference)
    
#this is sample code for running a single sample to obtain a mutation footprint 
if False: print("""
y=mutation_types(M1norm,x1=900,x2=1100)

y[9]=0
z = M_norm(y)
z[9]=0

print('{0:.2}'.format(np.log(sum(y))/np.log(10)))

plot_heatmap(mutation_values=z,circle=True,showtxt=True, title='ctrl')
     """)


def mutation_types(Mnorm,x1=750, x2=1150, STDEV=False):
    
    #input Mnorm, single matrix
    #output: A->:A, T, C, G, T->:A, T, C, G, C->:A, T, C, G, G->:A, T, C, G array
    
    lookup={'A':0,'T':1,'C':2,'G':3}
    Mmut  =np.zeros(16)
    Mstd  =np.zeros(16)
    
    for i in ['A','T','C','G']:
        
        letter = lookup[i]
        
        #select all of letter 'N'
        Msub = Mnorm[:4,x1:x2]
        xN = Msub[letter,:]>.51
        #take the cross-wise average of N->A, N->T, N->C, N->G
        y  = np.average(Msub[:,xN],axis=1) 
        z  = np.std(Msub[:,xN],axis=1)
        #print(y)
        
        #set identity base to zero!
        y[y>.51] = 0
        #this yields the following like array: [.999, .0001, .01, .004, 0.00]

        Mmut[letter*4:letter*4+4] = y
        Mstd[letter*4:letter*4+4] = z
    
    if STDEV:
        return Mmut, Mstd   
    
    return Mmut #16x1 array
    

def plot_heatmap(mutation_values=np.zeros([100,5]),circle=False, showtxt=False, title=''):
    #input format:
    """
    0 A->A
    1 A->T
    2 A->C
    3 A->G
    4 T -A
    5   -T
    6   -C
    7   -G
    8 C..
    9
    10
    11
    12
    13
    14...
    
    A->:A, T, C, G, T->:A, T, C, G, C->:A, T, C, G, G->:A, T, C, G
    
    """
    #initial values
    xy_values=np.array([[0, 0], [0, 1], [0, 2], [0, 3], [1, 0], [1, 1], [1, 2], [1, 3], [2, 0], [2, 1], 
        [2, 2], [2, 3], [3, 0], [3, 1], [3, 2], [3, 3]])*.2
    
    size_rect=.2
    size_radius=.098
    
    fig = plt.figure()
    ax = fig.add_subplot(111, aspect='equal')
    plt.axis('off')
    
    ax.text(.4, .9, title, verticalalignment='center', horizontalalignment='center')
    
    #coloring
    #norm=np.sum(mutation_values)
    #mutation_values/=(norm+01.00)
    
    #print(mutation_values)
    
    
    #mutation values
    
    
    norm = colors.Normalize(vmin=np.min(mutation_values), vmax=np.max(mutation_values))
    max_mv=np.max(mutation_values)
    
    for i, xy in enumerate(xy_values):
        
        mv=mutation_values[i]
        ##NEED TO SET COLOR VALUES HERE
        #Cvalue=[mutation_values[i]*2,1,mutation_values[i]*2]
        
        Cvalue=[1-norm(mv)*.5, 1-norm(mv)*.95-.05,1-norm(mv)*.95-.05]
    
        
        if circle == False:
            ax.add_patch(
                patches.Rectangle(xy, size_rect, size_rect, color=Cvalue)
            )
        
        else:
            try: r=np.sqrt(mutation_values[i]/max_mv)
            except: r=0
            ax.add_patch(
                patches.Circle(xy+.1, size_radius*r , color=Cvalue)
                #RADIUS is proporation to the mutation_value[i], scaled by max_mv ~ allowing the maximum circle to fill
                #a box... the scaler doesn't matter, as it stays constant in a given graph.
            )
            
        #ax.add_patch(patches.Rectangle(xy, size_rect, size_rect, fill=False, lw=.2))
        if showtxt: #and mv <= 0: 
            nump = '{0:.0%}'.format(max(mv,0))
            ax.text(xy[0]+.105,xy[1]+.1, nump ,verticalalignment='center', horizontalalignment='center')
    
    
    #ax.get_xaxis().set_visible(False)
    #ax.get_yaxis().set_visible(False)
    
    
#########################PLOTTING A TABLE FOR EACH MUTATION RATE
def table_output(Mnormlist, row_labels=[], subtract_zero=False, zero=[], x1=1170, x2=1270, rows=[]):
    
    #rows = row_labels
    if len(rows)==0: rows = ('wtT7-KO', 'pA154', 'APN1-KO', 'APN2-SHL', 'CIP1-KO', 'REV3-SHL', 'REV7-SHL', 
     '479/82/149int', 'g479', 'g480', 'g481', 'g482', 'g79/80/81', 'g63/66/69', 'wt-ctrl', 
     'wtT7-KO', 'pA152', 'pA153','pA158', 'pA159', 'polK', 'UDGmut', '79/82+148int', '79/82+681', '79/81+716')
    
    #are we going to execute the zero subtraction?? if so that will be the first matrix Mnormlist[0]
    normalize = subtract_zero
    if subtract_zero: baseline = np.round( mutation_types(zero, x1=x1, x2=x2) * 1E3, 3 )
    #############################################################
    #COMPUTE THE MUTATION TYPE FREQUENCIES and APPEND TO A MATRIX
    MUT = np.zeros((len(Mnormlist), 16))
    
    
    for i, norm in enumerate(Mnormlist):
        #calculate mutation types and round to 3rd decimnal place x10E3
        if normalize: MUT[i] =  np.round( mutation_types(norm, x1=x1, x2=x2) * 1E3, 3 ) - baseline
        else: MUT[i] =  np.round( mutation_types(norm, x1=x1, x2=x2) * 1E3, 3 )
        
        if normalize: MUT[i][MUT[i]<0] = 0
        #
        #print( str(i+1) + ': ' + str(MUT[i]))
    
    #####################
    #MAKE THE COLUMN PLOT
    normal = MUT/np.max(MUT+1E-8,axis=0)
    
    fig, axs =plt.subplots(1)
    collabel=("A:A", "A:T", "A:C", "A:G", "T:A", "T:T", "T:C", "T:G", "C:A", "C:T", 
          "C:C", "C:G", "G:A", "G:T", "G:C", "G:G")

    rowlabel= [str(i)+": "+rows[i] for i in range(len(MUT))]
    
    axs.axis('tight')
    axs.axis('off')
    the_table = axs.table(cellText=MUT,colLabels=collabel,rowLabels=rowlabel, 
                      loc='center', cellColours=plt.cm.summer(normal))
    the_table.scale(2,1)
    #axs.title("Mutation frequencies for samples x 10^-3")
    #plt.title('Mutation frequencies for samples x 10^-3')
    plt.show()

    
    
    
def CurvesDots(graph1, graph2, name1='0', name2='ctrl', baseStart='A', baseEND='T', window=20, Xrange=[0,1250], machine='none',  remove_out=True):
    
    #setting the range
    Xrange=Xrange
    
    #defining the dot area!
    area=.9
    
    if machine=='miseq': Ydic={'A':{'T':[1E-6,3E-4], 'C':[1E-6,1E-4], 'G':[1E-6,6E-4]},  
      'T':{'A':[1E-6,2E-4], 'C':[1E-6,6E-4], 'G':[1E-6,1E-4]}, 
      'C':{'A':[1E-6,1E-4], 'T':[1E-5,1E-1], 'G':[1E-6,1E-4]}, 
      'G':{'A':[1E-6,5E-4], 'T':[1E-6,2E-4], 'C':[1E-6,1E-4]} }
        
    elif machine=='hiseq': Ydic={'A':{'T':[1E-6,6E-4], 'C':[1E-6,6E-4], 'G':[1E-6,6E-4]},  
      'T':{'A':[1E-6,6E-4], 'C':[1E-6,6E-4], 'G':[1E-6,6E-4]}, 
      'C':{'A':[1E-6,6E-4], 'T':[1E-6,1E-1], 'G':[1E-6,6E-4]}, 
      'G':{'A':[1E-6,2E-3], 'T':[1E-6,6E-4], 'C':[1E-6,6E-4]} }
    
    else: Ydic={'A':{'T':[3E-5,1E-3], 'C':[3E-5,2E-3], 'G':[3E-5,1E-3]},  
      'T':{'A':[3E-5,1E-3], 'C':[3E-5,6E-4], 'G':[3E-5,2E-3]}, 
      'C':{'A':[3E-5,2E-3], 'T':[3E-5,1E-1], 'G':[3E-5,5E-4]}, 
      'G':{'A':[3E-5,2E-3], 'T':[3E-5,1E-3], 'C':[3E-5,6E-4]} }

    #print(baseStart, baseEnd)
    Yrange=Ydic[baseStart][baseEND]

    if baseStart=='C' and baseEND=='T': graphtype='log'
    else: graphtype='linear'

    ##############################################3
    lookup={'A':0,'T':1,'C':2,'G':3} #dictionar for bases
    MT=graph1 #pick out which data file to use 

    x_T=MT[lookup[baseStart]]>.51 #

    #plt.scatter(range(len(MT.T)), 1-np.max(MT.T[:,:],axis=1), s=area*2, color='green')
    plt.scatter(np.arange(-window, -window+len(MT.T)), MT.T[:,lookup[baseEND]]*x_T, s=area*3, color='blue')
    plt.scatter(np.arange(-window, -window+len(graph2.T)), graph2.T[:,lookup[baseEND]]*x_T[:len(graph2.T)], s=area*3, color='orange')

    plot_NtoN([graph1],baseStart=baseStart,baseEnd=baseEND,N=window,log=True, 
          lims=False, ASScolor=True, colors=['blue', 'blue'], linewidth=2.3, remove_out= remove_out)

    plot_NtoN([graph2],baseStart=baseStart,baseEnd=baseEND,N=window,log=True, 
          lims=False, ASScolor=True, colors=['orange', 'orange'], linewidth=2.3, remove_out= remove_out)

    plt.ylim(Yrange)
    plt.xlim(Xrange)
    plt.title(name1 + ' // ' + name2 + ' : ' + baseStart + ' -> ' + baseEND)
    #plt.yscale('log')
    plt.yscale(graphtype)
    return plt

def CurvesDots(graph1, graph2, name1='0', name2='ctrl', baseStart='A', baseEND='T', window=20, Xrange=[0,1250], machine='none',  remove_out=True, area=.9, lines=True):
    
    #setting the range
    Xrange=Xrange
    
    #defining the dot area!
    area=area
    
    if machine=='miseqALT': Ydic={'A':{'T':[1E-6,3E-4], 'C':[1E-6,1E-4], 'G':[1E-6,3E-4]},  
      'T':{'A':[1E-6,2E-4], 'C':[1E-6,4E-4], 'G':[1E-6,1E-4]}, 
      'C':{'A':[1E-6,1E-4], 'T':[1E-5,1E-1], 'G':[1E-6,1E-4]}, 
      'G':{'A':[1E-6,7E-4], 'T':[1E-6,2E-4], 'C':[1E-6,1E-4]} }
        
    if machine=='miseq': Ydic={'A':{'T':[1E-6,3E-4], 'C':[1E-6,1E-4], 'G':[1E-6,6E-4]},  
      'T':{'A':[1E-6,2E-4], 'C':[1E-6,6E-4], 'G':[1E-6,1E-4]}, 
      'C':{'A':[1E-6,1E-4], 'T':[1E-5,1E-1], 'G':[1E-6,1E-4]}, 
      'G':{'A':[1E-6,5E-4], 'T':[1E-6,2E-4], 'C':[1E-6,1E-4]} }
        
    elif machine=='hiseq': Ydic={'A':{'T':[1E-6,6E-4], 'C':[1E-6,6E-4], 'G':[1E-6,6E-4]},  
      'T':{'A':[1E-6,6E-4], 'C':[1E-6,6E-4], 'G':[1E-6,6E-4]}, 
      'C':{'A':[1E-6,6E-4], 'T':[1E-6,1E-1], 'G':[1E-6,6E-4]}, 
      'G':{'A':[1E-6,2E-3], 'T':[1E-6,6E-4], 'C':[1E-6,6E-4]} }
        
    elif machine=='dots': Ydic={'A':{'T':[.5E-6,5E-4], 'C':[.5E-6,5E-4], 'G':[.5E-6,5E-4]},  
      'T':{'A':[.5E-6,5E-4], 'C':[.5E-6,5E-4], 'G':[.5E-6,5E-4]}, 
      'C':{'A':[.5E-6,5E-4], 'T':[.5E-6,1E-1], 'G':[.5E-6,5E-4]}, 
      'G':{'A':[.5E-6,2E-3], 'T':[.5E-6,5E-4], 'C':[.5E-6,5E-4]} }
    
    else: Ydic={'A':{'T':[3E-5,1E-3], 'C':[3E-5,2E-3], 'G':[3E-5,1E-3]},  
      'T':{'A':[3E-5,1E-3], 'C':[3E-5,6E-4], 'G':[3E-5,2E-3]}, 
      'C':{'A':[3E-5,2E-3], 'T':[3E-5,1E-1], 'G':[3E-5,5E-4]}, 
      'G':{'A':[3E-5,2E-3], 'T':[3E-5,1E-3], 'C':[3E-5,6E-4]} }

    #print(baseStart, baseEnd)
    Yrange=Ydic[baseStart][baseEND]

    if baseStart=='C' and baseEND=='T': graphtype='log'
    else: graphtype='linear'

    ##############################################3
    lookup={'A':0,'T':1,'C':2,'G':3} #dictionar for bases
    MT=graph1 #pick out which data file to use 

    x_T=MT[lookup[baseStart]]>.51 #

    #plt.scatter(range(len(MT.T)), 1-np.max(MT.T[:,:],axis=1), s=area*2, color='green')
    plt.scatter(np.arange(-window, -window+len(MT.T)), MT.T[:,lookup[baseEND]]*x_T, s=area*3, color='blue')
    plt.scatter(np.arange(-window, -window+len(graph2.T)), graph2.T[:,lookup[baseEND]]*x_T[:len(graph2.T)], s=area*3, color='orange')

    if lines: plot_NtoN([graph1],baseStart=baseStart,baseEnd=baseEND,N=window,log=True, 
          lims=False, ASScolor=True, colors=['blue', 'blue'], linewidth=2.3, remove_out= remove_out)

    if lines: plot_NtoN([graph2],baseStart=baseStart,baseEnd=baseEND,N=window,log=True, 
          lims=False, ASScolor=True, colors=['orange', 'orange'], linewidth=2.3, remove_out= remove_out)

    plt.ylim(Yrange)
    plt.xlim(Xrange)
    plt.title(name1 + ' // ' + name2 + ' : ' + baseStart + ' -> ' + baseEND)
    #plt.yscale('log')
    plt.yscale(graphtype)
    return plt

def return_index(M, letter='A'):
    
    lookup={'A':0,'T':1,'C':2,'G':3}
    
    start = lookup[letter]
    xI = M[start,:] > .8
    
    return xI

def add_means(means, stdevs):
    mean = np.sum(means)
    std  = np.sum(stdevs**2)**(1/len(means))
    return mean, std

def divide_means(mean1, mean2, std1, std2):
    
    mean1, mean2 = np.array(mean1), np.array(mean2)
    
    mean1[mean1==0]=1E-6
    mean2[mean2==0]=1E-6
        
    mean = mean1 / mean2
    
    std  = mean * np.sqrt( (std1/mean1)**2 + (std2/mean2)**2 )
    return mean, std

def average_means(means, stdevs):
    #check if there are means to be averaged
    if len(means)<2: return means[0], stdevs[0]
    
    #average the means and proper stddev error prop
    mean = np.average(means)
    std  = np.sqrt( np.sum((means - mean)**2) / (len(means) - 1.0) )
    return mean, std

def ATCG_values(MnormA1, MnormA2, MnormB1, MnormB2, Xstart=800, Xend=1100, naming=''):

    #pull out mutation type information for replicate samples A1, 2 and B1, 2
    m1, m2 = mutation_types(MnormA1, x1=Xstart, x2=Xend), mutation_types(MnormA2, x1=Xstart, x2=Xend)
    m3, m4 = mutation_types(MnormB1, x1=Xstart, x2=Xend), mutation_types(MnormB2, x1=Xstart, x2=Xend)
    
    #generate and average and STDEV for the 2 replicates
    x1, y1 = np.average([m1, m2], axis=0), np.std([m1, m2], axis=0)
    x2, y2 = np.average([m3, m4], axis=0), np.std([m3, m4], axis=0)
    
    #divide the means to find a ratio of increase and error-prop standard deviations
    x, y = divide_means(x2, x1, y2, y1) 
    
    #print(x)
    
    #average means (A->N [1,2,3]) (T->N [4,6,7]) (C->N [8,9,11]) (C->T [9]) (C->!T [8, 11]) (G->N [12,13,14])
    lookup = {"A->N":[1,2,3], "T->N":[4,6,7], "C->!T":[8,11], "C->T":[9], "G->N":[12,13,14]}
    values, means, stds = {}, [], []
    
    for i in ["A->N", "T->N", "C->!T", "C->T", "G->N"]:
        z = average_means(x[lookup[i]], y[[lookup[i]]])
        
        means.append(z[0])
        stds.append(z[1])
        
        values[i] = z
        
    plt.barh(range(len(means))[::-1], means, .5, alpha=1, color='grey', 
                    xerr=stds, error_kw={'ecolor': '0.0'}, label='testing',log=True)
    
    plt.xlim([.1,100])
    plt.axvline(x=1, color='black', linestyle='--')
    plt.yticks(range(len(means))[::-1], ["A->N", "T->N", "C->!T", "C->T", "G->N"])
    plt.xlabel('fold increase in substitution rate +pT7/-pT7')
    plt.ylabel('substitution types')
    
    plt.title('ratio of subs rates +/-pT7 / 2-replicate / bp' + str(Xstart) + '-' + str(Xend) + '/' + naming)
        
    return values, means, stds, x, y



def runATCG(M1, M2, N1, N2, window=10, Xrange=[0,1250], savefig=False, machine='hiseq', remove_out=True, area=.9, lines=True):
    
    plt.figure(figsize=(25,25))
    k=1
    
    for i in ['G','C','T','A']:
        for j in ['A','T','C','G']:
            if i != j:
                #print(i,j)
                plt.subplot(4,4,k)
                fig = CurvesDots(M1, M2, name1=N1, name2=N2,
                    baseStart=j, baseEND=i, window=window, Xrange=Xrange, machine=machine, remove_out=remove_out, area=area, lines=lines)
                
                k+=1
            else: k+=1
    
    filename=N1+'__'+N2
    
    if savefig: plt.savefig(str('20180112_Analysis/' + filename + '.png'), transparent=False)
        #plt.savefig(str(directory + filename + '_' + str(sumX1X2[0]) + '_' 
                       # + str(sumX1X2[1]) + '.svg'), transparent=True)
            
            
def ATCG_values(MnormA1, MnormA2, MnormB1, MnormB2, Xstart=800, Xend=1100, naming=''):

    #pull out mutation type information for replicate samples A1, 2 and B1, 2
    m1, m2 = mutation_types(MnormA1, x1=Xstart, x2=Xend), mutation_types(MnormA2, x1=Xstart, x2=Xend)
    m3, m4 = mutation_types(MnormB1, x1=Xstart, x2=Xend), mutation_types(MnormB2, x1=Xstart, x2=Xend)
    
    #generate and average and STDEV for the 2 replicates
    x1, y1 = np.average([m1, m2], axis=0), np.std([m1, m2], axis=0)
    x2, y2 = np.average([m3, m4], axis=0), np.std([m3, m4], axis=0)
    
    #divide the means to find a ratio of increase and error-prop standard deviations
    x, y = divide_means(x2, x1, y2, y1) 
    
    #print(x)
    
    #average means (A->N [1,2,3]) (T->N [4,6,7]) (C->N [8,9,11]) (C->T [9]) (C->!T [8, 11]) (G->N [12,13,14])
    lookup = {"A->N":[1,2,3], "T->N":[4,6,7], "C->!T":[8,11], "C->T":[9], "G->N":[12,13,14]}
    values, means, stds = {}, [], []
    
    for i in ["A->N", "T->N", "C->!T", "C->T", "G->N"]:
        z = average_means(x[lookup[i]], y[[lookup[i]]])
        
        means.append(z[0])
        stds.append(z[1])
        
        values[i] = z
        
        
    plt.bar(range(len(means)), means, .5, alpha=1, color='grey', 
                    yerr=stds, error_kw={'ecolor': '0.0'}, label='testing',log=True)
    
    plt.axhline(y=1, color='black', linestyle='--')
    plt.xticks(range(len(means)), ["A->N", "T->N", "C->!T", "C->T", "G->N"])
    plt.ylabel('fold increase in substitution rate +pT7/-pT7')
    plt.xlabel('substitution types')
    
    plt.title('ratio of subs rates +/-pT7 / 2-replicate / bp' + str(Xstart) + '-' + str(Xend) + '/' + naming)
        
    return values, means, stds
        
