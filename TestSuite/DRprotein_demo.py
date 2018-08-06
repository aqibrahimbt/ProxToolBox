#Input, Output and Logfiles
import sys
sys.path.append('..')
problemName = '1PTQ'
inputFile   = '../InputData/Protein/1PTQ.pdb' #in InputFiles directory
outputFile  = '1PTQrecon.pdb'

#Parameters used by all variants of the algorithm
gapEps = 1e-5
maxIterations = 5000
P1Eps = 1e-12
P2Rank = 3

#Data simulation parameters
maxBondLength = 6

startingSeed = 2
replications = 1

from proxtoolbox.Problems import DRprotein

for rep in range(replications):
    problem = DRprotein(problemName,startingSeed+rep)

    problem.readInputFile(inputFile)
    problem.addOutputFile(outputFile)

    problem.setGapEps(gapEps)
    problem.setMaxIterations(maxIterations)
    problem.setP1Eps(P1Eps)
    problem.setP2Rank(P2Rank)

    problem.simulateNMR(maxBondLength)
    problem.addAminoInfo()

    problem.solve()
    problem.printStats()

    problem.savePDB()
    #problem.saveJPG()
    problem.show()

    del problem
