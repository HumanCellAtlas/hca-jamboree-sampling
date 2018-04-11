# Copyright © 2017 Anna Gilbert, Alexander Vargo, Umang Varma
# 
# This file is part of PicturedRocks.
#
# PicturedRocks is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# PicturedRocks is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with PicturedRocks.  If not, see <http://www.gnu.org/licenses/>.

import numpy as np
import scipy.spatial.distance as scipydist
from scipy.sparse.linalg import svds

# Currently need this for 1CS method (along with numpy)
#import cvxpy as cvx

# Helper functions for the 1 bit compressed sensing method

#### Correlations, as in SPA Stanardize
#
# This computes a vector. The $i$-th entry is the correlation (as defined in
# section 3.2.3 of Genzel's thesis) between feature vector corresponding to gene
# $i$ and the cluster data corresponding to the input clusters . 

# compute corrlation of columns of mat with vec Note that if a denominator is 0
# (from a standard deviation being 0), then the correlation will also be 0 (as
# vec = vecBar or the column of mat = the column of matBar) this is desired for
# this situation since this will only occur if vec = 0 or the column of mat = 0
def corrVec(vec, mat):
    matBar = mat.mean(axis=0)
    vecBar = vec.mean()
    
    r_num = np.sum( (mat-matBar)*(vec[:,None]-vecBar), axis=0) *1.0
    r_den = vec.shape[0]*np.std(vec)*np.std(mat,axis=0)

    if (len(list(r_den[r_den == 0])) != 0): 
        r_den[r_den == 0] = r_den[r_den == 0] + 0.000000001
            
    return r_num/(r_den)

# Generate a list of values of epsilon to test.
def genEpsList(nsteps=10, minEps=10.0**(-5), maxEps=10.0**(-3)):
    stepSize = (maxEps - minEps)/(nsteps*1.0)
        
    epsList = [maxEps - mult* stepSize for mult in range(nsteps+1)]
    return epsList

# soft threshold inVec by decreasing all entries by param
def softThreshold(inVec, param):
    
    outVec = np.copy(inVec)
    outVec[np.abs(inVec) < param] = 0.0
    outVec[inVec >= param] = inVec[inVec >= param] - param
    outVec[inVec <= -1.0*param] = inVec[inVec <= -1.0*param] + param
    
    return outVec

def ST(inVec, param):
    
    signs = np.sign(inVec)
    inVec = inVec - param * signs
    inVec[ np.invert(signs == np.sign(inVec)) ] = 0.0
    return inVec
    

def pca(Xin, dim=3):
    Xbar = Xin.mean(axis=0)
    Xcent = Xin - Xbar
    print("Computing Corvariance Matrix")
    Sigma = np.cov(Xcent.T)
    print("Performing SVDs")
    pcs = svds(Sigma, dim, return_singular_vectors="u")
    Xpca = Xcent.dot(pcs[0])
    return (Xcent, Sigma, pcs, Xpca)


class Rocks:
    # A Rocks object takes in the following parameters:
    # X: the gene expression matrix. We expect a numpy array of shape (N, P)
    #    containing data for N cells and P genes (note rows are cells and
    #    columns are genes.
    # y: cluster labels for N cells. We expect a numpy array of shape (N, 1) or
    #    (N,).
    # genes: names of genes. We expect an array of P strings, containing the
    #    names of various genes
    # verbose: verbosity level for debugging; defaults to 0.
    def __init__(self, X, y, genes=None, verbose=0):
        # self.X is the expression data
        # self.y is the cluster assignment
        # self.N is the number of cells
        # self.P is the number of genes
        # self.K is the number of clusters
        # self.genes contains the names of the genes
        self.verbose = verbose
        self.X = X
        self.y = y
        self.N, self.P = X.shape

        self.genes = genes
        assert (genes is None or len(genes) == self.P), \
                "genes must be an array of length P or None"

        self.Xcent, self.pcs, self.Xpca = (None, None, None)
        if y.shape == (self.N,):
            self.y = y.reshape((self.N, 1))
        assert self.y.shape == (self.N, 1), \
                "y should be a matrix of shape (N,1)"
        
        self.K = self.y.max() + 1
        assert np.array_equal(np.unique(self.y), range(self.K)), \
                "Cluster labels should be 0, 1, 2, ..., K -1"

        # Some extra elements that will be needed for OvA
        # they will be changed in the methods below
        self.cs_currY = self.y
        self.cs_currX = self.X

        # save the vectors of consts so that we don't need to deal with
        # manipulating the data over and over.  
        # Useful when testing multiple values of the parameter.
        self.cs_OvA = np.zeros((self.K, self.P))
        self.cs_scales = np.zeros((self.K, self.P))
        self.cs_generated = False
        
        self.clusterindices = {}
        for k in range(self.K):
            self.clusterindices[k] = np.nonzero(self.y == k)[0]
        nk = np.array([len(self.clusterindices[k]) for k in range(self.K)])
    
    def _debug(self, level, message):
        if self.verbose >= level:
            print(message, flush=True)

    # markers can be a list or 1-dim np array
    def markers_to_genes(self, markers):
        try:
            return [self.genes[a] for a in markers]
        except TypeError:
            raise ValueError("Gene names not specified. Set using object.genes")

    
    def normalize(self, totalexpr="median", log=True):
        cellsize = self.X.sum(axis=1).reshape((self.N,1))
        targetsize = np.median(cellsize) if totalexpr == "median" else totalexpr

        # avoid zero_divide issues
        cellsize[cellsize == 0] = 1

        self.X = (targetsize*self.X)/cellsize
        if log:
            self.X = np.log(self.X +1)
        self.Xcent, self.pcs, self.Xpca = (None, None, None)
    
    def pca(self, dims=3):
        self.Xcent, Sigma, self.pcs, self.Xpca = pca(self.X, dims)
        self.totalvariation = np.trace(Sigma)
        # Sigma is too big to be worth storing in memory

    
    def markers_mutualinfo(self, n, objective = "MRmR", pool = None):        
        import datetime
        X = np.log2(self.X+1).round().astype(int)
        if pool is None:
            pool = range(self.P)
        
        maxentry = max(X.max(), self.y.max())
        base = 10**int(round(np.log10(maxentry) + 1))

        def I(cx, cy):
            xy = np.array([cx, cy]).T
            xyconcat = xy.dot(np.array([base, 1]))
            values, counts = np.unique(xyconcat, return_counts=True)
            valcount = zip(list(values), list(counts))
            pxy = np.zeros((maxentry + 1, maxentry + 1))
            px = np.zeros(maxentry + 1)
            py = np.zeros(maxentry + 1)
            for value, count in valcount:
                xval = value // base
                yval = value % base
                pxy[xval, yval] += count
                px[xval] += count
                py[yval] += count
            with np.errstate(divide='ignore', invalid='ignore'):
                s = (pxy/np.expand_dims(px, axis=1))/np.expand_dims(py, axis=0) * n
                r = np.sum(np.nan_to_num((pxy/n) * np.log(s)))
            return r
        
        self._debug(1, "Computing I(x_i, y) values...")
        yflat = self.y.flatten()
        Ixiy = np.zeros(self.P)
        for i in pool:
            Ixiy[i] = I(X[:,i], yflat)
        
        Ixixj = np.zeros((self.P,self.P))
        IxixjChecked = np.zeros((self.P,self.P)).astype(bool)

        def getIxixj(i,j):
            if not IxixjChecked[i,j]:
                Ixixj[i, j] = I(X[:, i], X[:,j])
                Ixixj[j, i] = Ixixj[i, j]
                IxixjChecked[i,j] = IxixjChecked[j,i] = True
            return Ixixj[i, j]

        start = datetime.datetime.now()
        S = []
        Phi = 0
        self._debug(1, "Selecting candidate features...")
        assert objective in ["MR", "MRmMR", "MRmR"], "Invalid objective"
        for m in range(n):
            self._debug(1, "m = {}".format(m))
            maxDeltaPhi = float("-inf") # max increase in Phi
            argmaxDeltaPhi = -1 # index  i that corresponds to the value above
            for i in pool:
                if i in S:
                    continue
                DeltaPhi = Ixiy[i]
                if objective == "MR":
                    pass
                if objective == "MRmMR":
                    if m > 0:
                        DeltaPhi -= (1.0/m) * np.max([getIxixj(i, j) for j in S])
                elif objective == "MRmR":
                    if m > 0:
                        DeltaPhi -= (1.0/m) * np.sum([getIxixj(i, j) for j in S])
                if DeltaPhi > maxDeltaPhi:
                    maxDeltaPhi = DeltaPhi
                    argmaxDeltaPhi = i
                    self._debug(3,
                            "Feature provisionally selected m = {}, i = {}".\
                                    format(m, i))
            S.append(argmaxDeltaPhi)
            self._debug(2, "Features: {}".format(S))
        end = datetime.datetime.now()
        timedelta = end-start
        self._debug(1,
                "It took {:.2f} minutes to find the {} features via {}."\
                        .format(timedelta.total_seconds()/60, n, objective))
        return S


    # THE 1 BIT COMPRESSED SENSING METHOD (1CS, GENZEL)
    #
    # The methods below implement only one-vs-all (ova) multiclass
    # classification.  All-vs-all will be added soon
    #
    # TODO: Some functions below might be repeats of earlier functions

    # Get the indices of cells in cluster clustind
    # Note that the index that you input should be one MORE than the actual cluster 
    # you are looking for
    def clust2vec(self, clustind=1):
        returnvec = -1.0*np.ones(self.cs_currY.shape[0])
        returnvec[[i for i,x in enumerate(self.cs_currY) if x==clustind - 1]] \
                = 1.0
        return returnvec

    # Make the appropriate vector for 1CS methods
    def coeffs(self, clustind=1):
        """Creates the vector of coefficients used in the optimization problem.
        The solution is a soft-thresholding of this vector.  

        The method uses the data in self.cs_currX to get to the correct vector
        of coefficients, not the data in self.X

        :param clustind: The cluster that you would like to be separating.

        """
        return np.sum( self.clust2vec(clustind)[0:,np.newaxis]*self.cs_currX,
                axis=0 )

    def findXform(self, alpha=0.3, c=2.0, lambFac=1.0, clustOne=1, clustTwo=None):
        """Find the transform needed for scaling the data for the 1-bit CS methods.

        See Genzel (2015) for more information about the parameters.

        :param alpha: The CS alpha parameter
        :param lambFac: The CS lambda paramter
        :param clustOne: The first cluster that you want to consider
        :param clustTwo: [Optional] the second cluster to consider
        (for use with finding markers between two clusters; this functionality
        is not fully implemented)
        """
        # print("Finding transform with alpha = {}, c = {}, lambFac = {}, cluster
        # = {}".format(alpha, c, lambFac, clustOne), flush=True)
       
        # restrict the data if we are looking at two clusters
        if (clustTwo):
            self.cs_currY = np.concatenate( [self.y[self.y ==clustOne],
                self.y[self.y==clustTwo]] )
            self.cs_currX = np.concatenate( [self.X[
                np.squeeze(np.asarray(self.y==clustOne)) ],
                self.X[ np.squeeze(np.asarray(self.y==clustTwo)) ]] )
        
        # find vector of correlations
        rho = corrVec( self.clust2vec(clustOne), self.cs_currX )
        # print("Correlations: {}".format(rho[0:10]), flush=True)
        
        # find scaling vector
        sigma = np.std(self.cs_currX, axis=0)

        # the only time sigma = 0 for our situation is a gene that is never
        # expressed so this doesn't hurt anything.
        sigma[sigma==0] = 0.0000001
        alphaFac = alpha**(c * ( 1-np.abs(rho)))
        scaleFac = lambFac*alphaFac + (1 - alphaFac)/sigma
        
        # find center
        xbar = self.cs_currX.mean(axis=0)
        
        return xbar, scaleFac

    # Implement the optimization using the cvxpy methods
    #
    # When we call this, we need self.cs_currX to be exactly the data we want to
    # use in the optimization That is, it should be standardized and transformed
    # in every way that we want.  Also, lamb is the overal lambda constraint on
    # the 1-norm
    def findW(self, clustInd, lamb, tol=10.0**(-3)):
        consts = self.coeffs(clustInd)
        
        w = cvx.Variable(self.P) # self.P is the number of genes
        
        constraints = [ cvx.pnorm(w,1) <= lamb, cvx.norm(w,2) <= 1.0 ]
        objective = cvx.Maximize( cvx.sum_entries(
            cvx.mul_elemwise(consts[0:,np.newaxis],w) ))
        prob= cvx.Problem(objective,constraints)

        prob.solve(solver="SCS")
        
        # consistency checks
        if (np.sum(w.value) - lamb > tol):
            print("Warning for cluster {} - 1 norm condition failed: {}"\
                    .format(clustInd, np.sum(w.value)),flush=True)
        
        tmp = 0
        for i in range(self.cs_currX.shape[1]):
            tmp += w.value[i,0]**2
            
        if (tmp - 1 > tol):
            print("Warning for cluster {} - 2 norm condition failed: {}"\
                    .format(clustInd, tmp), flush=True)
        
        
        return w.value

    # The functions below are used when we transform the data back to a
    # consistent space and then find margins.
    #
    # Since we transform the data back to a consistent space, we don't need to
    # worry about keeping the transformed data around in order to find the
    # margins, so we can do everything on the full list of dubs

    # Input: one dubs vector (as a numGenes x 1 matrix), desired data in
    # consistent space (numCells x numGenes matrix)

    # TODO: Be careful which data you are using here.  Is self.cs_currX correct,
    # or do we need to input the data that we are using?  Answer: for now, we
    # are fine.  We transform the dubs back into the standard way of looking at
    # things.  So we just use self.cs_currX (which should be a centered version
    # of self.X for OVR),
    def dubs2margin(self, eps, dubs):
        
        epsDubs = np.copy(dubs)
        epsDubs[abs(epsDubs) < eps] = 0
        
        return np.squeeze(np.asarray(np.dot(self.cs_currX, epsDubs)))

    # dubs should be a list of (numGenes x 1) matrices, each matrix specifies an
    # optimal hyperplane
    #
    # eps is the cuttoff - we make all entries of the obtained dubs vector 0 if
    # they are smaller than eps
    #
    # desired data in consistent space (numCells x numGenes matrix)
    #
    # Output: all margins for a specific value of epsilon
    def genMarginsDubs(self, eps, dubs):
        
        margins = [self.dubs2margin(eps, w) for w in dubs]
        return margins

    # TODO: this is probably not the best way to accomplish this, since we are
    # really using a different clustering method to determine the accuracy of
    # the selected set of markers.  It should probably be removed.
    #
    # This checks for how many we get correct based on  the given cluster data
    # margins are the output from the above - all margins for a specific value
    # of epsilon.
    def numCorrect(self, margins):

        # some hacking to make this work when you are not cross-validating
        tmpY = self.cs_currY
        if (len(self.cs_currY.shape) < 2): tmpY = tmpY[:,None]

        final = np.argmax(np.array(margins), axis=0)
        # Debugging: some shape errors
        #print(tmpY.shape, flush=True)
        #print(final.shape, flush=True)
        # how many classified correctly
        return tmpY[final[:,np.newaxis]==tmpY].shape[0]

    # Run the method in a simple way
    # 
    # The following method runs the ovr method for simple cases.
    def runOvr(self, currLamb, lambFac, alpha, epsList=genEpsList()):

        numClusts = self.K
        
        self._debug(1, "Working on lamb = {}".format(currLamb))
        
        dubs = []
        margins = []

        for clust in range(1,numClusts+1):

            # Don't think that this is needed, but just being safe
            self.cs_currX = self.X
            cent,scale = self.findXform(alpha=alpha, c=2.0, lambFac=lambFac,
                    clustOne=clust)
            # print("Transform center and scale: {}, {}".format(cent, scale),
            # flush=True)

            # transform the data into the proper space for marker selection and
            # find the markers. Since this is OvR, we will currently have that
            # self.cs_currX = self.X.  Thus, there is some sketchy use of
            # self.cs_currX here. We reset self.cs_currX = self.X at the end to
            # get rid of the memory.
            self.cs_currX = (self.X - cent)*scale
            currDub = self.findW(clust, currLamb)

            # rescale the weight vector and save it
            dubs.append(np.squeeze(np.asarray(currDub))*scale)

            # reset currX to free up memory
            self.cs_currX = self.X
            
        
        # calculate all margins and how many cells we classify correctly for
        # each value in the epslist.  in OVR, our center stays the same (since
        # we are always looking at all of the data).
        # Note that we are putting the ORIGINAL DATA back into the
        # classification method
        # TODO: Should probably remove this step.  Just pick and eps and stick
        # with it.
        cent = self.X.mean(axis=0)
        self.cs_currX = self.X - cent
        margins = [ self.genMarginsDubs(eps, dubs) for eps in epsList ]

        ###
        #print("Margins:", flush=True)
        #for marg in margins:
        #        print("{}".format(marg[0:10]))

        final = [ self.numCorrect(marg) for marg in margins ]
        epsInd = np.argmax(np.array(final))

        bestEps = epsList[epsInd]
        currCorrect = final[epsInd]
        # assume monotonic - we always do better as we decrease epsilon
        # This is not really the case.
        change = (epsInd == len(final)-1)
        
        self._debug(1, "Summary: Number of correctly classified cells for"
                " all values of epsilon")
        self._debug(1, "Summary: {}".format(final))
        self._debug(1, "Summary: bestEps = {}, correctly classified cells = {}"\
                .format(bestEps, currCorrect))
        
        if (change):
            self._debug(1, "Warning: Optimal eps had change=True")
                
        return [dubs, np.array(margins[epsInd]), epsInd]

    
    # run the method and save a list of markers 
    # default behavior: alpha = 0 (no transform used on the data)
    # if you input a value for alpha, lambFac will default to 0
    def simpleCS(self, currLamb, writeOut = True, alpha = 0.0, lambFac = 0.0,
        epsList=genEpsList()):

        dubs,margins,epsInd = self.runOvr(currLamb=currLamb, lambFac=lambFac,
            alpha=alpha, epsList=epsList)

        self._debug(1,"Found markers")

        # the actual classification - note that we don't use this for anything
        # right now
        classes = np.argmax(margins, axis=0)
        currCorrect = self.y[ classes[:,np.newaxis]==self.y ].shape[0]

        # find the support genes by truncating according to the best value of
        # epsilion
        bestEps = epsList[epsInd]
        print("Testing: bestEps = {} with currCorrect = {}".format(bestEps,
            currCorrect), flush=True)
        for w in dubs:
            w[abs(w) < bestEps] = 0

        support_genes = [np.nonzero(w)[0] for w in dubs]

        # write the support genes to a data file
        if (writeOut):
            geneFile = "ovrGenes-lamb{}-lFac{}-alpha{}.dat".format(currLamb,
                    lambFac, alpha)
            gFile = open(geneFile, 'w')
            for genes in support_genes:
                for gene in genes:
                    gFile.write("{} ".format(gene))
                gFile.write("\n")
        
        # Return the classifcation and indices of the genes used in the
        # hyperplanes
        return [classes, support_genes]
    
    # Rescale inVec so that its 2-norm is less than param2 
    #        and its 1-norm is less than param1
    # Return the rescaled vector along with the norm (1 or 2)
    # that is "tightest"
    def rescale(self, inVec, param1, param2):
    
        tol = 1e-8
    
        whichNorm = 2
        outVec = inVec * param2 / np.linalg.norm(inVec)
    
        norm1 = np.sum(np.abs(outVec))
        if (norm1 > param1):
            whichNorm = 1
            outVec = outVec * param1/norm1
        
        # These are probably not needed and I have never seen them. 
        if (self.verbose >= 3):
            assert np.linalg.norm(outVec) <= param2 + tol, \
                "ERROR: 2-norm of rescaled vector exceeds desired parameter"
            assert np.sum(np.abs(outVec)) <= param1 + tol, \
                "ERROR: 1-norm of rescaled vector exceeds desired parameter"
    
        return outVec, whichNorm

    # Try soft thresholding with each coordinate to figure out how many
    # coordinates the optimization should actually take.  Assuming that we
    # start with the 2-norm (i.e. s >= 1 and consts != 0)
    #
    # We only report 1 norm for strict '>' in the above function; thus, we
    # never cut off too many indices (even if there is a point where the 2-norm
    # and 1 are equal
    #
    # absConsts should be the sorted norms of the consts vector
    # TODO: this shouldn't really be a class method...
    def findNormChanges(self, consts, absConsts, s=10):

        currNorm = 2
        prevNorm = 2
        switchInd = [] 
        numSwitch = 0
        ind = 1 # skip the largest index so that we point along an axis 

        # for simple tests, make sure that the index only changes once when
        # running on actual data just find the first place that the index
        # changes.
        while ind < consts.shape[0] and len(switchInd) == 0:

            val = absConsts[ind]
            curr, currNorm = self.rescale( softThreshold(consts, val), np.sqrt(s), 1.0 )
        
            if (currNorm != prevNorm):
                self._debug(1, "Switched from {} norm to {} norm at ind {}"\
                        .format(prevNorm, currNorm, ind))
                prevNorm = currNorm
                switchInd.append(ind)
                numSwitch = numSwitch + 1
        
            ind = ind + 1

        if (numSwitch > 1): self._debug(0, "Warning: norms switched {} times"\
                .format(numSwitch))

        return switchInd

    def findCSmarkers(self, currLamb, writeOut=True, alpha=0.0, lambFac=0.0):
        """Find markers via the 1-bit compressed sensing method following
        a "one vs all" (OvA) philosophy: for each cluster, find the markers
        that best separate the cluster from the rest of the data and return 
        the union of the markers.
        
        We use a soft-thresholding of a specific vector to get to a quick solution
        (instead of optimizing).  The first run will be slower, since we generate
        a transformed version of the data.  This information is saved, however, 
        so subsequent runs will be very fast.  This should allow you to find markers
        for many different values of the inputs quite quickly.

        :param currLamb: The sparsity parameter (controls how many features to select)
        :param writeOut: (bool) whether or not to write the selected markers to a file.
        :param alpha: the CS alpha paramter
        :param lambFac: the CS lambda paramter.
        """
        self._debug(1, "Working on lamb = {}".format(currLamb))

        marks = []

        # self.K = number of clusters
        for clust in range(1,self.K+1):

            self._debug(1, "Finding markers for cluster {}".format(clust))
            # Don't think that this is needed, but just being safe
            self.cs_currX = self.X

            cent,scale = self.findXform(alpha=alpha, c=2.0, lambFac=lambFac,
                clustOne=clust)

            self._debug(3, "Transform center and scale: {}, {}".\
                format(cent, scale))

            self.cs_currX = (self.X - cent)*scale
            consts = self.coeffs(clust)

            # flipud reverses the list
            # the sorting does not appear to take very long.  I could save this
            # information if it becomes the next bottleneck, however.
            absConsts = np.flipud(np.sort(np.abs(consts)))
            switchInd = self.findNormChanges(consts, absConsts, s=currLamb)
            print("Found all norm changes.", flush=True)
            if len(switchInd) > 1: self._debug(0, "More than one norm switch:\
                            assuming the first is the best")
            
            curr, currNorm = self.rescale(\
                softThreshold(np.asarray(consts).flatten(),\
                    absConsts[switchInd[0]]),\
                np.sqrt(currLamb), 1.0 )

            marks.append(curr.nonzero()[0])

        # delete final cs_currX
        self.cs_currX = self.X

        # write the support genes to a data file
        if (writeOut):
            geneFile = "ovrGenes-lamb{}-lFac{}-alpha{}.dat".format(currLamb,
                    lambFac, alpha)
            gFile = open(geneFile, 'w')
            for geneList in marks:
                for gene in geneList:
                    gFile.write("{} ".format(gene))
                gFile.write("\n")

        # return the list of markers
        return self.support2list(marks)

    # assuming OvA so that the center is just the center of the data
    # (and thus we are free to throw it out).
    def findCSmarkersForClust(self, currLamb, alpha=0.0, lambFac=0.0, clust=1):

        self._debug(1, "Finding markers for cluster {}".format(clust))

        consts = self._genOvaConsts(clust, alpha, lambFac)
        # flipud reverses the list
        # the sorting does not appear to take very long.  I could save this
        # information if it becomes the next bottleneck, however.
        absConsts = np.flipud(np.sort(np.abs(consts)))
        switchInd = self.findNormChanges(consts, absConsts, s=currLamb)
        print("Found all norm changes.", flush=True)
        if len(switchInd) > 1: self._debug(0, "More than one norm switch:\
                        assuming the first is the best")
        
        curr, currNorm = self.rescale(\
            softThreshold(np.asarray(consts).flatten(),\
                absConsts[switchInd[0]]),\
            np.sqrt(currLamb), 1.0 )

        return curr.nonzero()[0]

    def _genOvaConsts(self, clust, alpha, lambFac):
        """Transform the data into the proper space for marker selection
        and generate the correct consts vector in the transformed space.

        You only need to run this if you are transforming the data in a 
        specific way before finding the consts vector.
        """

        # don't do this if you don't have to.
        if self.cs_generated: 
            return self.cs_OvA[clust-1]

        # make sure to run this in order so that you don't short circuit
        if clust == self.K: self.cs_generated = True

        # Don't think that this is needed, but just being safe
        self.cs_currX = self.X

        cent,scale = self.findXform(alpha=alpha, c=2.0, lambFac=lambFac,
            clustOne=clust)

        self._debug(3, "Transform center and scale: {}, {}".\
            format(cent, scale))

        self.cs_currX = (self.X - cent)*scale
        consts = self.coeffs(clust)
        
        # save the vector so that you only do this once.
        self.cs_OvA[clust-1] = consts
        self.cs_scales[clust-1] = scale
        return consts


    # Makes a list out of a collection of support genes
    def support2list(self, sgenes):
        import itertools
        flattened = itertools.chain.from_iterable(sgenes)
        return list(set(flattened))

    def markers_CS(self, currLamb, writeOut=False, alpha = 0.0, lambFac = 0.0,
            epsList=genEpsList()):
        
        [classes, support_genes] =  self.simpleCS(currLamb, writeOut,
                alpha=alpha, lambFac=lambFac, epsList=epsList)
        return self.support2list(support_genes)

    # This finds `eps`, the offset that you need to soft threshold after you have
    # found the index when the norm constraint switches from the 2-norm to the
    # 1-norm.  (Note that, after the norm constraint switches from 2 to 1, the
    # entries are too large - that is, you have not soft thresholded enough at that
    # point).
    # 
    # * `xVec` should be the soft thresholded version of the solution that has the
    # support of the exact solution.  
    # * `kay` should be the number of nonzero entries in xVec 
    # * `ess` is the full sparsity parameter
    def findEps(self, xVec, kay, ess):
    
        assert xVec.nonzero()[0].shape[0] == kay,\
            "Error with inputs: xVec doesn't have kay nonzero entries"
    
        xNorm = np.sum(np.abs(xVec))
        c = xNorm*xNorm - ess * np.sum( np.abs(xVec) * np.abs(xVec) )
        b = -1.0 * xNorm * 2 * (kay-ess)
        a = kay*(kay - ess)
    
        # quadratic formula
        #return [(-1.0 * b - np.sqrt( b*b - 4*a*c )) / (2*a), (-1.0 * b + np.sqrt( b*b - 4*a*c )) / (2*a)]

        # use np.roots for more precision.
        return np.roots(np.array([a,b,c]))

    # curr should be the result of a rescale(softThreshold)
    def findWST(self, curr, ess):
        tol = 1e-4

        # number of nonzero entries
        kay = curr.nonzero()[0].shape[0]

        # np.roots returns array of length 2 even if root of multiplicity 2
        eps = self.findEps(curr, kay, ess)

        curr1 = ST(curr, eps[0])
        curr2 = ST(curr, eps[1])

        curr1 = curr1 / np.sum(np.abs(curr1)) * np.sqrt(ess) if np.sum(np.abs(curr1)) > 0 else 0.0
        curr2 = curr2 / np.sum(np.abs(curr2)) * np.sqrt(ess) if np.sum(np.abs(curr2)) > 0 else 0.0

        n1 = np.abs(1.0-np.linalg.norm(curr1))
        n2 = np.abs(1.0-np.linalg.norm(curr2))

        if (n1 < tol and n2 < tol): _debug(0, "Warning: two very similar\
                solutions for W found.  Arbitrarily choosing one.")

        return curr1 if n1 <= n2 else curr2

def pcafigure(celldata):
    import colorlover as cl
    import plotly.graph_objs as go
    if celldata.Xpca is None:
        celldata.pca(3)
    Xpca = celldata.Xpca
    clusterindices = celldata.clusterindices
    colscal = cl.scales['9']['qual']['Set1']
    plotdata = [go.Scatter3d(
            x=Xpca[inds,0],
            y=Xpca[inds,1],
            z=Xpca[inds,2],
            mode='markers',
            marker=dict(
                size=4,
                color=colscal[k % len(colscal)], # set color to an array/list
                #                                  of desired values
                opacity=1),
            name="Cluster {}".format(k),
            hoverinfo="name"
            )
            for k, inds in clusterindices.items()]

    layout = go.Layout(
        margin=dict(
            l=0,
            r=0,
            b=0,
            t=0
        )
    )
    return go.Figure(data=plotdata, layout=layout)
