#!/usr/bin/python3 -u

##
## Author:			Matthew Tam
##						CARMA, University of Newcastle, Australia
## Created:			19th March 2014
## Last modified: 17 June, 2016 (to conform with PythProx)
##
##  Description:  Implementation of the Douglas--Rachford algorithm, and variants,
##                       for the protein reconstruction problem.
##
##

import numpy as np
#cimport numpy as np
import time as time
import os as os
import subprocess as sp
import proxtoolbox.Utilities.Procrustes as procrustes
import sys
import itertools
from shutil import which

#DTYPE = np.float64
#ctypedef np.float64_t DTYPE_t

# The DRprotein class is the main class for the protein reconstruction problem.
# Member variables:
#  randomSeed     = the initial random seed.
#  inputFile      = a PDB input file containing the protein from which the NMR experiment is to be simulated.
#  outputFile     = the file name for the reconstructed protein to be saved as a PDB file.
#  logFile        = file name for the log file.
#  GAP_EPS        = the algorithm will terminate when the "gap" vector is less than this
#  MAX_ITERATIONS = the maximum number of iterations befored by the algorithm
#  P1_EPS         = the projection onto the known bond length set will be computed to within this tolerence.
#  P2_RANK        = the rank which the protein should be embeded in (typically this should be set to 3).
#  useRAAR        = Should the iteration be regularized/relaxed?
#  betaRange      = a pair [start_beta,final_beta] of regularization/relaxation paramters, used only if useRAAR==True
#  typePC2        = a string indicating how the rank projection is performed. 
#  typePC2freq    = an integer which denoted how often the rank projection should be performed exactly, only applies
#                   only if typePC2!='default'
# Member functions:
#  solve          = solve the problem!
class DRprotein:  

    def __init__(self,problemName='DRprotein', startingSeed=False):
        #A label for which names the problem.
        self.problemName = problemName

        #Seed the random number generator.
        seed = startingSeed if startingSeed else np.random.random_integers(1)
        self.setRandomSeed(seed)

        self.inputFile   = ''
        self.outputFile = ''

        #Set the logFile and imageFile based on the problemName
        #self.logFile = 'LogFiles/' + self.problemName + '_'
        #self.imageFile = 'Images/' + self.problemName + '_'
        self.logFile = f'{self.problemName}_'
        self.imageFile = f'{self.problemName}_'
        fileSuffix = 0
        while(os.path.isfile(self.logFile + str(fileSuffix).zfill(3) + '.txt')):
            fileSuffix += 1
        self.logFile      =  self.logFile     + str(fileSuffix).zfill(3) + '.txt'
        self.imageFile = self.imageFile + str(fileSuffix).zfill(3) + '.jpg'
        self.fileSuffix = fileSuffix

        # Open the logfile for writting
        self.logFileHandle = open(self.logFile, "w")

        # Print & log initialisation info.
        self.drucken(" DRprotein ".center(80,"="))
        self.drucken(f'Run filename:           {os.path.basename(__file__)}')
        self.drucken('Current time:           ' + time.strftime("%Y-%m-%d %H:%M:%S"))
        self.drucken(f'Starting (random) seed: {str(self.randomSeed)}')
        self.drucken('')

        self.setGapEps()
        self.setMaxIterations()
        self.setP1Eps()
        self.setP2Rank()

        self.NOT_DATA = -100 # A placeholder use to indicator that a dataMatrix entry is unknown.

        self.useRAAR(False)
        self.setP2Type('default')

        if not which('jmol'):
            print(""" Warning: jmol not found in $PATH. Without Jmol, molecular images cannot be
viewed or saved. If you wish to view or saves images, please install Jmol from
http://jmol.sourceforge.net or via your distribution's package manager.
""")

    def __del__(self):
        #sys.stdout.close()
        self.logFileHandle.close()

    def drucken(self, message):
        print(message)
        self.logFileHandle.write(message + "\n")

    def readInputFile(self, inputFile):
        self.inputFile = inputFile
        fileType = os.path.splitext(inputFile)[-1]
        if fileType == ".pdb":
            self.readPDB()
        elif fileType == ".xyz":
            self.readXYZ()
        else:
            raise Exception(f'Input file type {fileType} not recognised.')

    def addOutputFile(self, outputFile):
        splitText = os.path.splitext(outputFile)
        self.outputFile = (
            f"{splitText[0]}_{str(self.fileSuffix).zfill(3)}{splitText[1]}"
        )

    def setRandomSeed(self, startingSeed=2):
        self.randomSeed = startingSeed
        np.random.seed(self.randomSeed)

    def setGapEps(self, gapEps=1e-3):
        self.GAP_EPS = gapEps

    def setMaxIterations(self, maxIterations=1000):
        self.MAX_ITERATIONS = maxIterations

    def setP1Eps(self, P1Eps=1e-12):
        self.P1_EPS = P1Eps

    def setP2Rank(self, rank=3):
        self.P2_RANK = rank

    def useRAAR(self,betaRange=False):
        if betaRange:
            self.useRAAR = True
            self.betaRange = betaRange
        else:
            self.useRAAR = False
            self.betaRange = [1, 1]

    def setP2Type(self, type='default', frequency=1):
        self.typePC2 = type
        self.typePC2freq = frequency

    # This function parses a *.pdb file and, from it, generates the coresponding EDM.
    # Input: self.inputFile         = a string specificing the *.pdb to be parsed.
    # Ouput: self.trueCoordinates   = a 3 by n (where n is the number of atom) matrix where the i-th row contains 
    #                                  the coordinates of the i-th atom. 
    #        self.trueMatrix        = the corresponding EDM
    def readPDB(self):
        self.drucken(" Loading PDB File ".center(80,"="))
        self.drucken(f"Input file:      {self.inputFile}")
        t0 = time.time()
        self.trueCoordinates = []
        with open(self.inputFile) as f:
            for line in f:
                if line.startswith('ATOM') or line.startswith('HETATM'):
                    coords = [line[30:37], line[38:45], line[46:54]]
                    self.trueCoordinates.append(np.array(list(map(float, coords))))
        numAtoms = len(self.trueCoordinates)
        self.trueMatrix = np.matrix( [ [ np.linalg.norm(coords1-coords2)**2 for coords1 in self.trueCoordinates ] for coords2 in self.trueCoordinates] )
        self.trueCoordinates = np.matrix(self.trueCoordinates)
        self.drucken(f"Number of atoms: {numAtoms}")
        self.drucken(f"Load time:       {str(time.time() - t0)} sec")
        self.drucken("")

    def readXYZ(self):
        self.drucken(" Loading XYZ File ".center(80,"="))
        self.drucken(f"Input file:      {self.inputFile}")
        t0 = time.time()
        self.trueCoordinates = []
        self.atomType = []
        lineNumber = 1
        with open(self.inputFile) as f:
            for line in f:
                if lineNumber > 2:# line.startswith('ATOM') or line.startswith('HETATM'):
                    coords = line.split()[1:4]
                    self.atomType.append(line.split()[0])
                    #coords = [line[30:37], line[38:45], line[46:54]]
                    self.trueCoordinates.append(np.array(map(float, coords)))
                lineNumber += 1
        numAtoms = len(self.trueCoordinates)
        self.trueMatrix = np.matrix( [ [ np.linalg.norm(coords1-coords2)**2 for coords1 in self.trueCoordinates ] for coords2 in self.trueCoordinates] )
        self.trueCoordinates = np.matrix(self.trueCoordinates)
        self.drucken(f"Number of atoms: {numAtoms}")
        self.drucken(f"Load time:       {str(time.time() - t0)} sec")
        self.drucken("")

#    def readXYZ2(self):
#        print " Loading XYZ File ".center(80,"=")
#        print "Input file:      " + self.inputFile
#        t0 = time.time()
#        self.trueCoordinates = []
#        atoms = []
#        lineNumber = 1        
#        with open(self.inputFile) as f:
#            for line in f:
#                if lineNumber > 2:# line.startswith('ATOM') or line.startswith('HETATM'):
#                    atom   = line.split()[0]
#                    coords = line.split()[1:3]
#                    #coords = [line[30:37], line[38:45], line[46:54]]
#                    self.trueCoordinates.append(np.array(map(float, coords)))
#                    atoms = append(atom)
#                lineNumber += 1
#        numAtoms = len(self.trueCoordinates)
#        self.trueMatrix = np.matrix( [ [ np.linalg.norm(coords1-coords2)**2 for coords1 in self.trueCoordinates ] for coords2 in self.trueCoordinates] )
#        self.trueCoordinates = np.matrix(self.trueCoordinates)
#        print "Number of atoms: " + str(numAtoms)
#        print "Load time:       " + str(time.time()-t0) + " sec"
#        print ''
        #Compute distances
        #Cation/Anion: N/O & C/H/N
#        self.dataMatrix = np.multiply( self.trueMatrix, (self.trueMatrix<=maxLength**2) ) + self.NOT_DATA * (self.trueMatrix>maxLength**2)        
#        nonzeroDistances = self.dataMatrix.size - self.dataMatrix.shape[0]
#        numKnownDistances = np.sum( self.dataMatrix !=  self.NOT_DATA ) - self.dataMatrix.shape[0]
#        print "Number of known non-zero inter-atomic distances:     " + str(numKnownDistances) +  ' of ' + str(nonzeroDistances)
 #       print "Percentage of known non-zero inter-atomic distances: " + str(float(numKnownDistances) / nonzeroDistances * 100) + ' %'
  #      print ''


    # This function simulates a NMR experiment in which all distances below maxLength are recorded.
    # Input:  self.trueMatrix  = the compelte EDM
    #         maxLength        = the maximum length which the NMR machine can record
    # Output: self.dataMatrix  = a partial EDM contain only distances less than maxLength
    def simulateNMR(self, maxLength, noiseLevel=0):
        self.drucken(' Simluating NMR experiment '.center(80,"="))
        if noiseLevel>0:
            theShape = self.trueMatrix.shape
            #noisyMatrix = self.trueMatrix +  np.matrix(np.random.uniform(-float(noiseLevel)/2,float(noiseLevel)/2,tuple(theShape)))
            noisyMatrix = self.trueMatrix +  np.matrix(np.random.normal(0,float(noiseLevel)/3,tuple(theShape)))
        else:
            noisyMatrix = self.trueMatrix
        self.dataMatrix = np.multiply( noisyMatrix, (noisyMatrix<=maxLength**2) ) + self.NOT_DATA * (noisyMatrix>maxLength**2)
        nonzeroDistances = self.dataMatrix.size - self.dataMatrix.shape[0]
        numKnownDistances = np.sum( self.dataMatrix !=  self.NOT_DATA ) - self.dataMatrix.shape[0]
        self.drucken(
            f"Maximum known bond length:                           {str(maxLength)}"
        )
        self.drucken(
            f"Number of known non-zero inter-atomic distances:     {str(numKnownDistances)} of {str(nonzeroDistances)}"
        )
        self.drucken(
            f"Percentage of known non-zero inter-atomic distances: {str(float(numKnownDistances) / nonzeroDistances * 100)} %"
        )
        self.drucken("")

    # This function all inter-atom distances within amino acids, from the PDB file. This is resonable since their
    # structures are known.
    # Input:
    # Output:
    def addAminoInfo(self):
        self.drucken(' Adding A Priori Amino Acid Information '.center(80,"="))
        t0 = time.time()
        seqResCount = []
        resSeqNum = -1
        with open(self.inputFile) as f:
            for line in f:
                if line.startswith('ATOM'):
                    thisResSeqNum = int(line[22:26])
                    thisResName = line[17:20]
                    if thisResSeqNum == resSeqNum:
                        seqResCount[-1] += 1
                    else:
                        seqResCount.append(1)
                        resSeqNum = thisResSeqNum
                elif line.startswith('HETATM'):
                    seqResCount.append(1)
        aminoInfoMatrix = np.matrix(np.ones(self.trueMatrix.shape)) * self.NOT_DATA

        #print len(seqResCount), np.sum( aminoInfoMatrix !=  self.NOT_DATA ), aminoInfoMatrix.shape
        #print seqResCount
        #print aminoInfoMatrix[0:5,0:5]


        for c in range(len(seqResCount)):
             start = sum(seqResCount[:c])
             end = start + seqResCount[c]
             #print start, end
             aminoInfoMatrix[start:end,start:end] = self.trueMatrix[start:end,start:end]

        #print len(seqResCount), np.sum( aminoInfoMatrix !=  self.NOT_DATA ), aminoInfoMatrix.shape
        #print seqResCount
        #print aminoInfoMatrix[6:13,6:13]


        #print np.linalg.norm(np.maximum(aminoInfoMatrix,self.dataMatrix)-self.dataMatrix)
        #print aminoInfoMatrix
        #print np.amax(aminoInfoMatrix)
        #self.dataMatix = np.matrix(np.maximum(np.array(self.dataMatrix),np.array(aminoInfoMatrix)))
        #print np.amax(self.dataMatrix)
        #self.dataMatrix = aminoInfoMatrix
        #print self.dataMatrix.size, self.dataMatrix.shape[0], np.sum( self.dataMatrix !=  self.NOT_DATA )
        #self.dataMatrix = np.maximum(self.dataMatrix,self.trueMatrix)
        self.dataMatrix = np.maximum(self.dataMatrix,aminoInfoMatrix)

        nonzeroDistances = self.dataMatrix.size - self.dataMatrix.shape[0]
        numKnownDistances = np.sum( self.dataMatrix !=  self.NOT_DATA ) - self.dataMatrix.shape[0]
        self.drucken(
            f"Number of known non-zero inter-atomic distances:     {str(numKnownDistances)} of {str(nonzeroDistances)}"
        )
        self.drucken(
            f"Percentage of known non-zero inter-atomic distances: {str(float(numKnownDistances) / nonzeroDistances * 100)} %"
        )

        #print self.dataMatrix[0:11,0:11]
        #print seqResCount
        #print sum(seqResCount), len(seqResCount)
        self.drucken("")

    # This function is similar to simulateNMR by return a partial EDM with the specified percentage of known entries.
    # Input:  desiredPercentage = the desire percentage of entries to known in the partial EDM
    # Output: 
    def addPercentOfDist(self, desiredPercentage=10):
        trueList = sorted(
            [
                self.trueMatrix[i, j]
                for i, j in list(
                    itertools.permutations(range(self.trueMatrix.shape[0]), 2)
                )
                if i < j
            ]
        )
        maxBondIndex = int(len(trueList)*desiredPercentage/100)-1
        maxBondLength = np.sqrt(trueList[maxBondIndex])
        self.simulateNMR(maxBondLength, noiseLevel=0)

    def addVanDerWaalInformation(self, threshold=1.0, vanDerWaalRadiiFile='InputFiles/vanDerWaalRadii.txt'): 
        self.drucken(' Adding A Priori van Waal Radii Information '.center(80,"="))
        vDW = np.loadtxt(vanDerWaalRadiiFile,dtype={'names':('atomType','radii'),'formats':('S2','f8')})
        vDWRadii = vDW['radii'].tolist()
        vDWAtomType = vDW['atomType'].tolist()
        xyzVDWRadii=[vDWRadii[vDWAtomType.index(atomType)] for atomType in self.atomType]
        radiiSq=np.matrix([[(r1+r2)**2 *  threshold**2 for r1 in xyzVDWRadii] for r2 in xyzVDWRadii])
        self.dataMatrix = np.multiply(self.trueMatrix,self.trueMatrix<=radiiSq) + self.NOT_DATA * (self.trueMatrix>radiiSq)
        nonzeroDistances = self.dataMatrix.size - self.dataMatrix.shape[0]
        numKnownDistances = np.sum( self.dataMatrix !=  self.NOT_DATA ) - self.dataMatrix.shape[0]
        self.drucken(
            f"Number of known non-zero inter-atomic distances:     {str(numKnownDistances)} of {str(nonzeroDistances)}"
        )
        self.drucken(
            f"Percentage of known non-zero inter-atomic distances: {str(float(numKnownDistances) / nonzeroDistances * 100)} %"
        )
        self.drucken("")        #print self.atomType

    # This function solves the loaded problem.
    # Input:  self.MAX_ITERATIONS = maximum number of iterations ot be performed.
    #         self.dataMatrix
    # Output: 
    def solve(self):
        # Print variable values from the input file.
        self.drucken(" Algorithm Parameters ".center(80,"="))
        memberVars = vars(self).keys()
        dontShow = ['trueCoordinates','trueMatrix','dataMatrix','inputFile','NOT_DATA', 'logFileHandle','fileSuffix','atomType']
        if self.typePC2 == 'default':
            dontShow.append('typePC2freq')
        if self.useRAAR == False:
            dontShow.append('betaRange')
        dontShow = set(dontShow)
        memberVars = list(set(memberVars)-dontShow)
        for memberVar in memberVars:
            self.drucken(f"{memberVar}:".ljust(30, " ") + str(vars(self)[memberVar]))
        self.drucken("")

        # Initiailize Variables (to take advantage of Cython we need to make typed copies)
        x0 = np.matrix( np.random.random_sample(size=self.dataMatrix.shape) )
        x0 = (x0 + x0.T) / 2

        xChange = self.NOT_DATA * np.ones(self.MAX_ITERATIONS)
        gap = self.NOT_DATA * np.ones(self.MAX_ITERATIONS)
        relGap = self.NOT_DATA * np.ones(self.MAX_ITERATIONS)
        dataMatrix = self.dataMatrix
        trueMatrix = self.trueMatrix
        P1_EPS = self.P1_EPS
        NOT_DATA = self.NOT_DATA
        GAP_EPS = self.GAP_EPS
        MAX_ITERATIONS = self.MAX_ITERATIONS
        P2_RANK = self.P2_RANK

        x   = x0
        P1x = P1(x0, dataMatrix, P1_EPS, NOT_DATA)
        P2R1x = np.empty(x0.shape)
        Tx    = np.empty(x0.shape)
        P1Tx  = np.empty(x0.shape)

        ## Variables used in computation of different heuristics. We compute them once here so we need not do it each time.
        dataMatrixLength = self.dataMatrix.shape[0]
        #cdef np.ndarray[DTYPE_t, ndim=2] y
        #cdef np.ndarray[DTYPE_t, ndim=2] U3 
        v = np.matrix( np.ones( (dataMatrixLength, 1)))
        v[-1] += np.sqrt(dataMatrixLength)
        houseHolderMatrix = np.eye(dataMatrixLength) - 2 * v * v.T / np.linalg.norm(v, 'fro')**2

        # Heuristics settings
        typePC2 = self.typePC2
        typePC2freq = self.typePC2freq

        ### Main Algorithm Loop ###
        self.drucken(" Running algorithm ".center(80,"="))

        t0 = time.time()
        iter = 0
        #while (iter==0) or ((iter < MAX_ITERATIONS) and (gap[iter-1] > GAP_EPS)):
        while (iter==0) or ((iter < MAX_ITERATIONS) and (relGap[iter-1] > GAP_EPS )):
            # Print progress so far.
            if iter==0:
                self.drucken(" Iterations ".center(10," ") + "Change".center(15," ") + "Gap".center(15," ") + "Rel. Gap".center(15," ") + "Error".center(15," "))
            elif (iter>0) and (np.mod( iter, 10 ) == 0):
              #self.drucken(str(iter).center(12," ") +  str(xChange[iter-1]).center(15," ") + str(gap[iter-1]).center(15," ")  + str(relGap[iter-1]).center(15," ") + str(np.linalg.norm(P1x-trueMatrix,'fro')).center(15," "))
                self.drucken( " " + str(iter).center(10, " ") + "{0:15e}".format(xChange[iter-1]) +  "{0:15e}".format(gap[iter-1]) +   "{0:15e}".format(relGap[iter-1]) + "{0:15e}".format(np.linalg.norm(P1x-trueMatrix,'fro')))

            # How should the projection onto C_2 be performed?
            if typePC2 in ['default', '']:
                P2R1x = P2(2*P1x-x, dataMatrixLength,  P2_RANK, houseHolderMatrix)
            elif typePC2 == 'lazyUpdate':
                #P2R1x = P2R1x
                if (iter<10) or np.mod(iter,typePC2freq) == 0:
                   P2R1x = P2(2*P1x-x, dataMatrixLength,  P2_RANK, houseHolderMatrix)
            elif typePC2 == 'lazyEigh':
                if (iter<10) or np.mod(iter,typePC2freq) == 0:
                   P2R1x, U3 = P2withU(2*P1x-x, dataMatrixLength,  P2_RANK, houseHolderMatrix)
                else:
                   P2R1x =  - P2R1x + PM(y, U3, houseHolderMatrix)
            elif typePC2 == 'adaptiveUpdate':
                if (iter<1000) or xChange[iter-1]<100:
                   P2R1x = P2(2*P1x-x, dataMatrixLength,  P2_RANK, houseHolderMatrix)
            else:
                raise Exception(f"typePC2 {typePC2} not found.")


            # Should Douglas--Rachford, or RAAR be used?
            if self.useRAAR:
                beta = 1
                Tx = beta * (x + P2R1x) + (1-2*beta) * P1x             
            else:
                Tx = x + P2R1x - P1x             

            P1Tx = P1(Tx, dataMatrix, P1_EPS, NOT_DATA)

            # When lazyEigh is used, we store the previous point.
            if self.typePC2 == 'lazyEigh':
                y = 2*P1Tx + P1x - x

            # Compute and append error statistics.
            xChange[iter] = np.linalg.norm( Tx-x, 'fro')
            gap[iter]     = np.linalg.norm( P2R1x-P1x, 'fro')
            relGap[iter] = gap[iter] / np.linalg.norm(P1x, 'fro')

            # Update.
            x = Tx
            P1x = P1Tx
            iter += 1

        # Collect final iterates and statistics.
        self.solveTime = time.time() - t0
        self.P1x   = P1x
        self.P2R1x = P2R1x
        self.xChange = xChange[xChange!=self.NOT_DATA]
        self.gap     = gap[gap!=self.NOT_DATA]
        self.relGap = relGap[relGap!=self.NOT_DATA]
        self.iter  = iter

        self.drucken("\n" + " Algorithm terminated ".center(80," ") + "\n")


    # This function generates a summary of statistics after solve() has been called.
    # Input:
    # Output:
    def printStats(self):
        self.drucken(' Summary statistics '.center(80,"="))
        self.drucken('Iterations:'.ljust(30)  + str(self.iter))
        self.drucken('Solve time:'.ljust(30) + str(self.solveTime) + ' sec')
        self.drucken('Change:'.ljust(30)     + str(self.xChange[self.iter-1]))
        self.drucken('Gap:'.ljust(30)        + str(self.gap[self.iter-1]))

        [fittedCoords, procrustesError ] = self.getCoordinates(True,True)
        fittedEDM = np.matrix([[np.linalg.norm(coords1-coords2)**2 for coords1 in fittedCoords] for coords2 in fittedCoords])
        EDMError = np.linalg.norm(fittedEDM-self.P1x,'fro')

        self.drucken('Procrustes (true) error:'.ljust(30) + str(procrustesError))
        self.drucken('Distance Matrix (true) error:'.ljust(30) + str(EDMError))
        self.drucken("")
        return [self.iter,self.solveTime,self.xChange[self.iter-1],self.gap[self.iter-1],procrustesError,EDMError]

    # This function generates specified dimensional atom coordiantes from the reconstructed EDM.
    # Requires: procrustes.py
    # Input:    self.P1x     = the final reconstruction
    #           self.P2_RANK = the dimension of the coordinates to be generated (default = 3)
    # Output:   
    def getCoordinates(self,fitCoords=True,procrustesError=True):
        #Generate the coordinates.
        n = self.P1x.shape[0]
        L = np.matrix(np.eye(n)) - (1. / n)
        tau = -0.5 * L * self.P1x * L
        U, s, V = np.linalg.svd(tau)
        coords = (U  * (np.diag(s)**0.5))[:,:self.P2_RANK]

        #Translate points to the positive orthant.
        #coords += np.matrix( -np.minimum(coords.min(axis=0), 0) )

        if fitCoords:
           d, coords, tform  = procrustes.procrustes(self.trueCoordinates, coords,scaling=False)
           d = np.linalg.norm(self.trueCoordinates-coords,'fro')

        return [coords, d] if procrustesError else coords

    # This function saves a *.pdb file containing the reconstructed molecule.
    # Input:  self.inputFile
    #         self.outputFile
    #         self.coord
    # Output:
    def savePDB(self):
        t0 = time.time()
        outputText = ""
        coords = self.getCoordinates(True,False)
        self.drucken( " Saving PDB File ".center(80,"="))
        with open(self.inputFile) as f:
            atomIndex = 0
            for line in f:
                if line.startswith('ATOM') or line.startswith('HETATM'):
                    #(x,y,z) = [line[30:37], line[38:45], line[46:54]]
                    thisLine = line[:30]
                    thisLine += "{0:8.3f}".format(coords[atomIndex, 0])
                    thisLine += "{0:8.3f}".format(coords[atomIndex, 1])
                    thisLine += "{0:8.3f}".format(coords[atomIndex, 2])
                    thisLine += line[54:]
                    outputText += thisLine
                    atomIndex += 1
                else:
                    outputText += line
        with open(self.outputFile, 'w') as f:
            f.write(outputText)
        self.drucken(f"Output file:     {str(self.outputFile)}")
        self.drucken(f"Save time:       {str(time.time() - t0)} sec")
        self.drucken("")


    def saveXYZ(self):
        t0 = time.time()
        outputText = ""
        coords = self.getCoordinates(True,False)
        #print coords
        self.drucken(" Saving XYZ File ".center(80,"="))
        lineNumber = 1
        with open(self.inputFile) as f:
            atomIndex = 0
            for line in f:
                if lineNumber > 2:
                    #(x,y,z) = [line[30:37], line[38:45], line[46:54]]
                    #print coords[atomIndex]
                    thisLine = line[:10]
                    thisLine += "{0:11.8f}".format(coords[atomIndex, 0])
                    thisLine += "   "
                    thisLine += "{0:11.8f}".format(coords[atomIndex, 1])
                    thisLine += "   "
                    thisLine += "{0:11.8f}".format(coords[atomIndex, 2])
                    thisLine += line[54:]
                    outputText += thisLine
                    atomIndex += 1
                else:
                    outputText += line
                lineNumber += 1
        with open(self.outputFile, 'w') as f:
            f.write(outputText)
        self.drucken(f"Output file:     {str(self.outputFile)}")
        self.drucken(f"Save time:       {str(time.time() - t0)} sec")
        self.drucken("")

    def saveJPG(self):
        self.drucken(" Saving JPG File ".center(80,"="))
        jmolCmd = [
            "jmol",
            "-n",
            self.outputFile,
            "-w",
            f"JPG:{self.imageFile}",
            "--exit",
        ]
        try:
            self.drucken(sp.check_output(jmolCmd))
        except sp.CalledProcessError as e:
            self.drucken(e.output.decode())
            self.drucken(f"Error while calling Jmol: {e}")
            self.drucken("")
    
    def show(self):
    	sp.call(["jmol", self.outputFile])

####################################################################################################################

# This function performs the projection onto "C_1" which is the data constraint.
# Input:  x           = the point to be projected
#         dataMatrix  = the partial matrix contain the known distances
#         P1_EPS      = the tolerance with which the distances are known
#         NOT_DATA    = a flag used to indicat which positions of dataMatrix are unknown
# Output: xx          = the projection of x onto "C_1"
#cdef np.ndarray[DTYPE_t, ndim=2] P1(np.ndarray[DTYPE_t, ndim=2] x, np.ndarray[DTYPE_t,  ndim=2] dataMatrix, DTYPE_t P1_EPS, int NOT_DATA):
#    cdef np.ndarray[DTYPE_t,  ndim=2] xx
def P1(x, dataMatrix, P1_EPS, NOT_DATA):
    xx = np.empty([x.shape[0], x.shape[1]])
    for j in range(x.shape[1]):
        for i in range(x.shape[0]):
            if dataMatrix[i, j] == NOT_DATA:
                xx[i, j] = max(x[i, j], 0)
            elif i == j:
                xx[i, j] = 0
            elif x[i, j] > dataMatrix[i, j] + P1_EPS:
                xx[i, j] = dataMatrix[i, j] + P1_EPS
            elif x[i, j] < max(dataMatrix[i, j] - P1_EPS, 0):
                xx[i, j] = max( dataMatrix[i, j] - P1_EPS,  0)
    return xx

# This function performs the projection onto "C_2" which is the EDM constraint.
# Input: x                 = the point to be projected
#        dataMatrixLength  = the number of atoms
#        P2_RANK           = the rank in which the EDM should be embeddable (default = 3)
#        houseHolderMatrix = the Householder matrix (computed elsewhere so that we don't have to each time)
# Output: the projection of x onto "C_2"
#cdef np.ndarray[DTYPE_t, ndim=2] P2(np.ndarray[DTYPE_t, ndim=2] x,  int dataMatrixLength,  int P2_RANK,  np.ndarray[DTYPE_t,  ndim=2] houseHolderMatrix):
#    cdef int n = dataMatrixLength
def P2(x,  dataMatrixLength,  P2_RANK,  houseHolderMatrix):
    Q = houseHolderMatrix
    n = dataMatrixLength

    #Do the rest
    QxQ = np.dot(np.dot(Q,x),Q)#Q * x * Q
#    cdef np.ndarray[DTYPE_t,  ndim=1] s
#    cdef np.ndarray[DTYPE_t,  ndim=2] U
    s,U = np.linalg.eigh( QxQ[:n-1,:n-1] ) #eigenvalues and eigenvectors of symmetric matrix

    s=np.minimum(s, 0) #ensure all eigenvalues are negative
    if P2_RANK>0:      #asumming ordered
        s[P2_RANK:]=0

    QxQ[:n-1, :n-1] = U * np.matrix(np.diag(s)) * U.T
    return Q * QxQ * Q

# 
#cdef P2withU(np.ndarray[DTYPE_t, ndim=2] x,  int dataMatrixLength,  int P2_RANK,  np.ndarray[DTYPE_t,  ndim=2] houseHolderMatrix):
#    cdef int n = dataMatrixLength
#    cdef np.ndarray[DTYPE_t, ndim=2] Q = houseHolderMatrix
#
#    #Do the rest
#    cdef np.ndarray[DTYPE_t, ndim=2]QxQ = np.dot(np.dot(Q,x),Q)#Q * x * Q
#    cdef np.ndarray[DTYPE_t,  ndim=1] s
#    cdef np.ndarray[DTYPE_t,  ndim=2] U
#    s,U = np.linalg.eigh( QxQ[:n-1,:n-1] ) #eigenvalues and eigenvectors of symmetric matrix
#
#    s=np.minimum(s, 0) #ensure all eigenvalues are negative
#    if P2_RANK>0: #asumming order ... not the case
#        s[P2_RANK:]=0
#
#    QxQ[:n-1, :n-1] = U * np.matrix(np.diag(s), dtype=DTYPE) * U.T
#    return Q * QxQ * Q, U[:,0:P2_RANK]
#
#cdef np.ndarray[DTYPE_t, ndim=2] PM(np.ndarray[DTYPE_t, ndim=2] y, np.ndarray[DTYPE_t, ndim=2] U3, np.ndarray[DTYPE_t, ndim=2] houseHolderMatrix):
#    QyQ = houseHolderMatrix * y * houseHolderMatrix
#    cdef np.ndarray[DTYPE_t, ndim=2] hatY
#    hatY = U3.T  * QyQ[:-1,:-1] * U3
#    hatY = np.diag(np.diag(hatY))
#    QyQ[:-1,:-1] = U3 * hatY * U3.T
#    return houseHolderMatrix * QyQ * houseHolderMatrix



####################################################################################################################

#def runExperiment():   
#    parameterFiles = ['1PTQ_parameters','1PTQ_parameters_lazyUpdate','1PTQ_parameters_lazyEigh','1PTQ_parameters_RAAR']
#    replications = 2
#    startingSeed = 2
#    for pF in parameterFiles:
#      stats = []
#      seed = startingSeed
#      for rep in xrange(replications):
#        problem = DRprotein(pF,seed)
#        problem.solve()
#        theseStats = problem.printStats()
#        stats.append(theseStats)
#        problem.savePDB()
#        del problem        
#        seed += 1
#        #[self.iter,self.solveTime,self.xChange[self.iter-1],self.gap[self.iter-1],procrustesError,EDMError]

def runExperiment():   
    proteins = ['1PTQ','1PHT','1HOE','1POA','1LFB','1AX8']
    suffix   = ['','_lazyUpdate','_lazyEigh','_RAAR']
#    parameterFiles = ['1PTQ_parameters','1PTQ_parameters_lazyUpdate','1PTQ_parameters_lazyEigh','1PTQ_parameters_RAAR']
    parameterFiles = [
        f'{pro}_parameters{suf}' for (pro, suf) in zip(proteins, suffix)
    ]
    replications = 5
    startingSeed = 2
    for pF in parameterFiles:
        stats = []
        seed = startingSeed
        for _ in range(replications):
            problem = DRprotein(pF,seed)
            problem.solve()
            theseStats = problem.printStats()
            stats.append(theseStats)
            problem.savePDB()
            del problem
            seed += 1
        #[self.iter,self.solveTime,self.xChange[self.iter-1],self.gap[self.iter-1],procrustesError,EDMError]       


def main():
    # Single demo.
    #parameterFile = '1PTQ_parameters'
    parameterFile = '1AX8_parameters'
    problem = DRprotein(parameterFile)
    problem.addAminoInfo()
#    problem.solve()
#    problem.savePDB()
#    problem.saveJPG()
    del problem

    #experiment demo.
    #runExperiment()

    #parameterFile = 'PAN_parameters'
    #problem = DRprotein(parameterFile)
    #problem.solve()
    #problem.printStats()
    #problem.saveXYZ()
    #del problem



if __name__ == '__main__':
  main()
