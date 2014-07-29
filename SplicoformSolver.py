"""
The MIT License (MIT)

Copyright (c) 2013 Oliver Serang

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
"""

import numpy as np
import itertools
import math
import random
import timeit
from scipy.signal import fftconvolve
import sys
import getopt

def normalized(arrayVec):
    return np.array(arrayVec) / sum(arrayVec)

class Peptide:
    def __init__(self, alpha, beta, pepScore):
        self.alpha = alpha
        self.beta = beta
        self.pepScore = pepScore
    def peptidePrior(self, count):
        return 1-(1-self.beta)*pow(1-self.alpha, count)
    def peptideLikelihood(self, count):
        prior = self.peptidePrior(count)
        # note: does not perform Fido's charge state prior removal
        return prior * self.pepScore + (1-prior) * (1-self.pepScore)
    def __repr__(self):
        return 'Peptide a=' + str(self.alpha) + " b=" + str(self.beta) + ' score=' + str(self.pepScore)

class Protein:
    def __init__(self, gamma, uniquePeptides):
        self.gamma = gamma
        self.uniquePeptides = uniquePeptides
        self.initUniquePeptideLikelihoodTable()
    def initUniquePeptideLikelihoodTable(self):
        result = []
        for intBool in (0,1):
            prod = 1.0
            for pep in self.uniquePeptides:
                prod *= pep.peptideLikelihood(intBool)
            result.append(prod)
        # normalize for greater precision
        self.uniquePeptideLikelihoodTable = normalized(result)
    def uniquePeptideLikelihood(self, intBool):
        return self.uniquePeptideLikelihoodTable[intBool]
    def prior(self, state):
        if state == 1:
            return self.gamma;
        return 1-self.gamma;
    def priorVector(self):
        return np.array( [self.prior(i) for i in (0,1)] )
    def likelihoodVector(self):
        return np.array( [self.uniquePeptideLikelihood(i) for i in (0,1)] )
    def __repr__(self):
        return 'Protein g=' + str(self.gamma)

class SplicoformGraph:
    def __init__(self, proteins, sharedEvidenceTable):
        self.proteins = proteins
        self.sharedEvidenceTable = sharedEvidenceTable
    def sharedLikelihoods(self, n):
        return self.sharedEvidenceTable[n]
    def allProteinPosteriors(self):
        result = []
        for protIndex in range(len(self.proteins)):
            result.append(self.posteriorForVariable(protIndex))
        return result
    def posteriorForVariable(self, protInd):
        raise Exception('Abstract: derived class needs to implement')
            
class ExponentialSplicoformGraph(SplicoformGraph):
    def joint(self, configuration):
        return self.prior(configuration) * self.uniqueLikelihoods(configuration) * self.sharedLikelihoods(sum(configuration))
    def prior(self, configuration):
        prod = 1.0
        for intBool, prot in zip(configuration, self.proteins):
            prod *= prot.prior(intBool)
        return prod
    def uniqueLikelihoods(self, configuration):
        prod = 1.0
        for protIndex, intBool in enumerate(configuration):
            prod *= self.proteins[protIndex].uniquePeptideLikelihood(intBool)
        return prod
    def allConfigurations(self):
        return list(itertools.product( *[ (0,1) for i in range(len(self.proteins)) ] ))
    def allConfigurationsAndJoints(self):
        result = []
        allConfigs = self.allConfigurations()
        for configuration in allConfigs:
            result.append( (configuration, self.joint(configuration)) )
        return result
    # overrides SplicoformGraph::allProteinPosteriors because it is
    # more efficient to compute all posterios in a single pass
    def allProteinPosteriors(self):
        total = 0.0
        totalWithProtein = np.array( [ 0.0 for x in range(len(self.proteins)) ] )
        for config, joint in self.allConfigurationsAndJoints():
            total += joint
            for proteinNumber in range(len(self.proteins)):
                if config[proteinNumber] == 1:
                    totalWithProtein[proteinNumber] += joint

        allPosteriors = totalWithProtein / total

        allPosteriorVectors = []
        for posterior in allPosteriors:
            allPosteriorVectors.append(np.array([1-posterior, posterior]))

        return allPosteriorVectors

class QuadraticCell:
    def __init__(self, edgeValues):
        self.edgeValues = edgeValues
        self.cumulativeFromLeft = None
        self.cumulativeFromRight = None
    def __str__(self):
        return 'E:' + str(self.edgeValues) +\
            ' ' + 'L:' + str(self.cumulativeFromLeft) + '\tR:' + str(self.cumulativeFromRight)
    def __repr__(self):
        return self.__str__()
 
class QuadraticTable:
    def __init__(self, nToSharedLikelihoods, proteins, edgeValuesForEachProtein):
        self.edgeValuesForEachProtein = edgeValuesForEachProtein

        # build the unique likelihood pyramid
        self.layerToCumulativeToCell = []
        for protInd in range(len(proteins)):
            edgeValues = self.edgeValuesForEachProtein[protInd]

            self.layerToCumulativeToCell.append( [] )
            # Layer i has i+1 cells
            for n in range(1 + protInd):
                self.layerToCumulativeToCell[protInd].append( QuadraticCell(edgeValues) )

        # add an n-sized layer for the nToSharedLikelihoods table
        self.layerToCumulativeToCell.append([])
        for n in range(len(proteins)+1):
            terminal_cell = QuadraticCell(None)
            terminal_cell.cumulativeFromRight = nToSharedLikelihoods[n]
            self.layerToCumulativeToCell[-1].append(terminal_cell)

        self.populateCumulativeFromRight()
        self.populateCumulativeFromLeft()

    def populateCumulativeFromLeft(self):
        # the initial node has no left hand side
        self.layerToCumulativeToCell[0][0].cumulativeFromLeft = 1.0

        # iterate up to the second-to-last element
        for protInd in list(range(len(self.layerToCumulativeToCell)))[:-1]:
            layer = self.layerToCumulativeToCell[protInd]
            nextLayer = self.layerToCumulativeToCell[protInd+1]
            # clear cumulative in next layer
            for cell in nextLayer:
                cell.cumulativeFromLeft = 0.0
            for count, cell in enumerate(layer):
                for countChange, weight in enumerate(self.edgeValuesForEachProtein[protInd]):
                    nextLayer[count+countChange].cumulativeFromLeft += weight * cell.cumulativeFromLeft

            # for greater precision, normalize the layer:
            totalAtNextLayer = sum([ cell.cumulativeFromLeft for cell in nextLayer ])
            for cell in nextLayer:
                cell.cumulativeFromLeft /= totalAtNextLayer

    def populateCumulativeFromRight(self):
        # work backwards using dynamic programming:
        # iterate in reversed order up to the second-to-last element
        for protInd in list(range(len(self.layerToCumulativeToCell)))[-2::-1]:
            layer = self.layerToCumulativeToCell[protInd]
            nextLayer = self.layerToCumulativeToCell[protInd+1]
            for count, cell in enumerate(layer):
                cell.cumulativeFromRight = sum([self.edgeValuesForEachProtein[protInd][countChange]*nextLayer[count+countChange].cumulativeFromRight for countChange, weight in enumerate(self.edgeValuesForEachProtein[protInd])])

            # for greater precision, normalize the layer:
            totalAtNextLayer = sum([ cell.cumulativeFromRight for cell in nextLayer ])
            for cell in nextLayer:
                cell.cumulativeFromRight /= totalAtNextLayer

    def posteriorForVariable(self, protIndex):
        layer = self.layerToCumulativeToCell[protIndex]
        nextLayer = self.layerToCumulativeToCell[protIndex+1]
        # for each node in the layer, compute the vector sum of all
        # edge vectors weighted by their cumulativeFromTheRight values
        numEdges = len(layer[0].edgeValues)
        edgeTotal = np.array( [0.0 for e in range(numEdges) ] )
        for cell in range(len(layer)):
            edgeTotal += np.array( [ layer[cell].cumulativeFromLeft*layer[cell].edgeValues[e]*nextLayer[cell+e].cumulativeFromRight for e in range(numEdges) ] )
        return (edgeTotal / sum(edgeTotal))

class QuadraticSplicoformGraph(SplicoformGraph):
    def __init__(self, proteins, sharedEvidenceTable):
        SplicoformGraph.__init__(self, proteins, sharedEvidenceTable)

        edgeValuesForEachProtein = []
        for protInd in range(len(self.proteins)):
            priors = np.array([ self.proteins[protInd].prior(i) for i in (0,1) ])
            likelihoods = np.array([ self.proteins[protInd].uniquePeptideLikelihood(i) for i in (0,1) ])

            # normalize for greater precision
            likelihoods = normalized(likelihoods)

            edgeValues = likelihoods * priors
            edgeValues = normalized(edgeValues)

            edgeValuesForEachProtein.append(edgeValues)

        # build the simple shared likelihood table
        nToSharedLikelihoods = []
        for n in range(len(self.proteins)+1):
            nToSharedLikelihoods.append(self.sharedLikelihoods(n))

        # for greater precision, normalize nToSharedLikelihoods
        nToSharedLikelihoods = normalized(np.array(nToSharedLikelihoods))
        self.quadraticTable = QuadraticTable(nToSharedLikelihoods, self.proteins, edgeValuesForEachProtein)

    def posteriorForVariable(self, protIndex):
        return self.quadraticTable.posteriorForVariable(protIndex)

class ConvolutionTreeNode:
    def __init__(self, jointAbove):
        # normalize for greater precision
        self.jointAbove = normalized(jointAbove)

        self.leftParent = None
        self.rightParent = None

        self.likelihoodBelow = None

    def __str__(self):
        result += "jointAbove: " + str(self.jointAbove) + "\n"
        result += "likelihoodBelow: " + str(self.likelihoodBelow) + "\n"

        return result

    @classmethod
    def createCountNode(cls, lhs, rhs):
        # pass messages down from lhs and rhs and create a new node

        ####### manual quadratic convolution of priors
        # prior = [0]*(len(lhs.prior)+len(rhs.prior)-1)
        # for i,valL in enumerate(lhs.prior):
        #     for j,valR in enumerate(rhs.prior):
        #         prior[i+j] += valL * valR
        ####### FFT version
#        prior = fftconvolve(lhs.prior, rhs.prior)

        ####### manual quadratic convolution of likelihoods above
        # likelihoodAbove = [0]*(len(lhs.prior)+len(rhs.prior)-1)
        # for i,valL in enumerate(lhs.likelihoodAbove*lhs.prior):
        #     for j,valR in enumerate(rhs.likelihoodAbove*rhs.prior):
        #         likelihoodAbove[i+j] += valL * valR / prior[i+j]

        ####### FFT version
        jointAbove = fftconvolve(lhs.jointAbove, rhs.jointAbove)
        result = cls(jointAbove)

        result.leftParent = lhs
        result.rightParent = rhs

        return result

    def messageUp(self, answerSize, otherJointVector):
        ####### manual quadratic convolution
        # note that this node is the N node
        # result = [0.0]*answerSize
        # for x in range(len(result)):
        #     for n in range(len(self.prior)):
        #         otherX = n-x
        #         if otherX not in range(len(otherJointVector)):
        #             continue
        #         result[x] += self.likelihoodBelow[n]*otherJointVector[otherX]

        ####### FFT
        startingPoint = len(otherJointVector)-1
        result = fftconvolve(otherJointVector[::-1], self.likelihoodBelow)[startingPoint:startingPoint+answerSize]

        return normalized(result)

    def messageUpLeft(self):
        return self.messageUp(len(self.leftParent.jointAbove), self.rightParent.jointAbove)
    def messageUpRight(self):
        return self.messageUp(len(self.rightParent.jointAbove), self.leftParent.jointAbove)

    # once all messages are received
    def posterior(self):
        return normalized(self.jointAbove * self.likelihoodBelow)

class ConvolutionTree:
    def __init__(self, nToSharedLikelihoods, proteins):
        self.nToSharedLikelihoods = nToSharedLikelihoods
        self.logLength = int(math.ceil(np.log2(float(len(proteins)))))
        self.allLayers = []
        self.buildFirstLayer(proteins)
        self.buildRemainingLayers()
        self.propagateBackward()

    def buildFirstLayer(self, proteins):
        # construct first layer (of proteins)
        layer = []
        for prot in proteins:
            protNode = ConvolutionTreeNode(prot.priorVector()*prot.likelihoodVector())
            layer.append( protNode )

        # pad with necessarily absent dummy variables so that the
        # number of variables is a power of 2; this is not the most
        # efficient method for this, but is simple for demonstration
        # purposes. because they are absent, they won't influence the
        # total sum, and thus Ds.
        for i in range(0, 2**self.logLength - len(proteins) ):
            # this protein cannot be present
            prot = Protein(0.0, [])
            layer.append( ConvolutionTreeNode(prot.priorVector()*prot.likelihoodVector()) )

        self.allLayers.append(layer)

    def buildRemainingLayers(self):
        # construct layers of count nodes
        for L in range(self.logLength):
            mostRecentLayer = self.allLayers[-1]
            layer = []
            for i in range(len(self.allLayers[0]) / (2**(L+1))):
                leftParent = mostRecentLayer[i*2]
                rightParent = mostRecentLayer[i*2+1]
                countNode = ConvolutionTreeNode.createCountNode(leftParent, rightParent)
                layer.append( countNode )

            # add connection to remaining nodes (when layer above is not a power of 2)
            self.allLayers.append(layer)

        # final node gets (Ds | N) multiplied into its likelihoodBelow
        finalNode = self.allLayers[-1][0]
        # normalize for greater precision
        finalNode.likelihoodBelow = normalized(self.nToSharedLikelihoods)

    def propagateBackward(self):
        # propagate backward, setting likelihoodBelow.
        # the loop has upper bound at logLength+1
        # because of the layer of proteins
        for L in range(1, self.logLength+1)[::-1]:
            layer = self.allLayers[L]

            for i in range(len(layer)):
                node = layer[i]

                leftParent = node.leftParent
                rightParent = node.rightParent

                leftParent.likelihoodBelow = node.messageUpLeft()
                rightParent.likelihoodBelow = node.messageUpRight()

        self.proteinLayer = self.allLayers[0]

    def posteriorForVariable(self, protInd):
        return self.proteinLayer[protInd].posterior()

class SubquadraticSplicoformGraph(SplicoformGraph):
    def __init__(self, proteins, sharedEvidenceTable):
        SplicoformGraph.__init__(self, proteins, sharedEvidenceTable)

        nToSharedLikelihoods = np.array([ self.sharedLikelihoods(n) for n in range(len(self.proteins)+1) ])

        self.convolutionTree = ConvolutionTree(nToSharedLikelihoods, self.proteins)

    def posteriorForVariable(self, protInd):
        return self.convolutionTree.posteriorForVariable(protInd)


def solvePowerSet():
    print 'Power-set (exponential)'
    sg = ExponentialSplicoformGraph(proteins, sharedEvidenceTable)
    allPosteriors = sg.allProteinPosteriors()
    if printResults:
        for protInd, posterior in enumerate(allPosteriors):
            print protInd, posterior
        
def solveQuadratic():
    print 'Dynamic programming (quadratic)'
    sg = QuadraticSplicoformGraph(proteins, sharedEvidenceTable)
    allPosteriors = sg.allProteinPosteriors()
    if printResults:
        for protInd, posterior in enumerate(allPosteriors):
            print protInd, posterior

def solveSubquadratic():
    print 'Convolution tree (subquadratic)'
    sg = SubquadraticSplicoformGraph(proteins, sharedEvidenceTable)
    allPosteriors = sg.allProteinPosteriors()
    if printResults:
        for protInd, posterior in enumerate(allPosteriors):
            print protInd, posterior

def printUsage():
    print """Usage: python SplicoformSolver.py --numberOfProteins = <int> [--printDescription] [--exponential] [--quadratic] [--subquadratic] [--printResults]"""

def printProblemDescription():
        print """***  ILLUSTRATION SPLICOFORM SOLVER ***
        Each protein has:
        \t-at least one unique peptide
        \t-a few shared peptides (shared by all splicoforms)

        Posterior probabilities for an arbitrary protein are computed using a generalized Fido model
        note that:
        \t-no pruning is performed (PSM scores are simulated so that pruning is inaccurate)
        \t-grouping is not performed (each protein has unique peptides, and cannot group)
        \t-the junction tree algorithm cannot improve this problem due to a clique in the moral graph
        """
def main(params):
    try:
        opts, args = getopt.getopt(params, '', ['printDescription', 'numberOfProteins=', 'exponential', 'quadratic', 'subquadratic', 'printResults'])

        # strip off the leading '--' and put into a dictionary
        optionsMap = dict( [ (flag[2:], value) for flag, value in opts ] )

        possibleOptions = ['printDescription', 'numberOfProteins', 'exponential', 'quadratic', 'subquadratic', 'printResults']
        if len(set(optionsMap) - set(possibleOptions)) > 0:
            raise Exception('Invalid options: ' + set(optionsMap) - set(possibleOptions))

    except getopt.GetoptError as error:
        printUsage()
        exit(1)

    if 'printDescription' in optionsMap:
        printProblemDescription()

    if 'numberOfProteins' not in optionsMap:
        raise(Exception('numberOfProteins must be specified'))
    n = int(optionsMap['numberOfProteins'])

    print 'Solving on ', n, 'proteins'
    # unfortunately, timeit doesn't accept parameters
    global proteins
    global sharedEvidenceTable
    global priorVectorOnProteins
    global printResults
    printResults = False
    # each protein has a random gamma, each peptide has a random
    # alpha and random beta
    proteins = [ Protein( random.uniform(0.1, 0.9), [ Peptide(random.uniform(0.1, 0.9), random.uniform(0.1, 0.9), random.uniform(0.1, 0.9)) for i in range(random.randint(1,1)) ]) for j in range(n) ]
    sharedEvidenceTable = normalized([ random.uniform(0.1, 0.9) for n in range(n+1) ])

    if 'printResults' in optionsMap:
        printResults = True

    if 'exponential' in optionsMap:
        t = timeit.timeit('solvePowerSet()', "from __main__ import solvePowerSet", number = 1)
        print 'Runtime (seconds):', t

    if 'quadratic' in optionsMap:
        t = timeit.timeit('solveQuadratic()', "from __main__ import solveQuadratic", number = 1)
        print 'Runtime (seconds):', t

    if 'subquadratic' in optionsMap:
        t = timeit.timeit('solveSubquadratic()', "from __main__ import solveSubquadratic", number = 1)
        print 'Runtime (seconds):', t

if __name__ == '__main__':
    main(sys.argv[1:])
