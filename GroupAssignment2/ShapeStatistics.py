#Purpose: To implement a suite of 3D shape statistics and to use them for point
#cloud classification
#TODO: Fill in all of this code for group assignment 2
import sys
sys.path.append("S3DGLPy")
from Primitives3D import *
from PolyMesh import *

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm, colors
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from scipy.special import sph_harm


POINTCLOUD_CLASSES = ['biplane', 'desk_chair', 'dining_chair', 'fighter_jet', 'fish', 'flying_bird', 'guitar', 'handgun', 'head', 'helicopter', 'human', 'human_arms_out', 'potted_plant', 'race_car', 'sedan', 'shelves', 'ship', 'sword', 'table', 'vase']

NUM_PER_CLASS = 10

#########################################################
##                UTILITY FUNCTIONS                    ##
#########################################################

#Purpose: Export a sampled point cloud into the JS interactive point cloud viewer
#Inputs: Ps (3 x N array of points), Ns (3 x N array of estimated normals),
#filename: Output filename
def exportPointCloud(Ps, Ns, filename):
    N = Ps.shape[1]
    fout = open(filename, "w")
    fmtstr = "%g" + " %g"*5 + "\n"
    for i in range(N):
        fields = np.zeros(6)
        fields[0:3] = Ps[:, i]
        fields[3:] = Ns[:, i]
        fout.write(fmtstr%tuple(fields.flatten().tolist()))
    fout.close()

#Purpose: To sample a point cloud, center it on its centroid, and
#then scale all of the points so that the RMS distance to the origin is 1
def samplePointCloud(mesh, N):
    (Ps, Ns) = mesh.randomlySamplePoints(N)
    # Find centroid. A column matrix of x,y,z averages
    centroid = Ps.mean(1)[:,None] 
    # Center the points on the centroid with broadcasting
    Ps -= centroid
    # Calculate the scale factor
    distances = np.linalg.norm(Ps,axis=0)# Calculate the euclidean distance of all points
    summedSquared = np.sum(np.square(distances)) # Square all distances then sum them
    scaleFactor = math.sqrt(np.divide(N,summedSquared)) # Calculate the scaling factor from this
    # Scale the points
    Ps = np.multiply(Ps,scaleFactor)

    return (Ps, Ns)

#Purpose: To sample the unit sphere as evenly as possible.  The higher
#res is, the more samples are taken on the sphere (in an exponential 
#relationship with res).  By default, samples 66 points
def getSphereSamples(res = 2):
    m = getSphereMesh(1, res)
    return m.VPos.T

#Purpose: To compute PCA on a point cloud
#Inputs: X (3 x N array representing a point cloud)
def doPCA(X):
    ##TODO: Fill this in for a useful helper function
    eigs = np.array([1, 1, 1]) #Dummy value
    V = np.eye(3) #Dummy Value
    return (eigs, V)

# computes angle between numpy arrays a and b
def angleBetween(a, b):
        angle = np.arccos(a.dot(b) / (np.linalg.norm(a) * np.linalg.norm(b)))

#########################################################
##                SHAPE DESCRIPTORS                    ##
#########################################################

#Purpose: To compute a shape histogram, counting points
#distributed in concentric spherical shells centered at the origin
#Inputs: Ps (3 x N point cloud), Ns (3 x N array of normals) (not needed here
#but passed along for consistency)
#NShells (number of shells), RMax (maximum radius)
#Returns: hist (histogram of length NShells)
def getShapeHistogram(Ps, Ns, NShells, RMax):
	norms = np.linalg.norm(Ps,axis=0)
	return np.array(np.histogram(norms,NShells,(0,RMax))[0])
    
#Purpose: To create shape histogram with concentric spherical shells and
#sectors within each shell, sorted in decreasing order of number of points
#Inputs: Ps (3 x N point cloud), Ns (3 x N array of normals) (not needed here
#but passed along for consistency), NShells (number of shells), 
#RMax (maximum radius), SPoints: A 3 x S array of points sampled evenly on 
#the unit sphere (get these with the function "getSphereSamples")
def getShapeShellHistogram(Ps, Ns, NShells, RMax, SPoints):
    NSectors = SPoints.shape[1] #A number of sectors equal to the number of
    #points sampled on the sphere
    # First, dot product all points in Ps against all points in SPoints.
    # We're doing a matrix multiplication of Ps^T with SPoints. This will give an NxM matrix
    # Where the rows are the dot product magnitudes of each point of Ps with all the points in SPoints
    dots = np.dot(Ps.T, SPoints)
    # Then find the index of the point in SPoints that yields the largest dot product for each Ps point
    maximums = np.argmax(dots,axis=1).T
    #Create a 2D histogram that is NShells x NSectors
    NSectors = SPoints.shape[1]
    hist = np.zeros((NShells, 0)) 
    # Go through each sector and create a histogram
    for i in range(NSectors):
        sectorElems = Ps[:,maximums==i] # Select every element in the given sector
        sectorHistogram=np.array(np.histogram(np.linalg.norm(sectorElems,axis=0),NShells,(0,RMax))[0])[np.newaxis].T # Create a histogram for the sector in Column form
        hist=np.hstack((hist,sectorHistogram)) # Add the column to the histogram
    # Sort shells in descending order in rows
    hist=-np.sort(-hist,axis=1) 
    return hist.flatten() #Flatten the 2D histogram to a 1D array

#Purpose: To create shape histogram with concentric spherical shells and to 
#compute the PCA eigenvalues in each shell
#Inputs: Ps (3 x N point cloud), Ns (3 x N array of normals) (not needed here
#but passed along for consistency), NShells (number of shells), 
#RMax (maximum radius), sphereRes: An integer specifying points on the sphere
#to be used to cluster shells
def getShapeHistogramPCA(Ps, Ns, NShells, RMax):
    # Create our linear distance spacing for the histogram
    shellSpacing = np.linspace(0,RMax,num=NShells+1)
    # Create a distance array for all the points
    pointDistances = np.linalg.norm(Ps,axis=0)
    numPoints = Ps.shape[1]
    #Create a 2D histogram, with 3 eigenvalues for each shell
    hist = np.zeros((0, 3))
    # Compute the PCA values for the matrices of points by using pointDistances as a binary selector for points
    for i in range(1,shellSpacing.shape[0]):
        if(i<shellSpacing.shape[0]-1): # If we're not in the outer shell
            shellPoints = Ps[:,np.logical_and(pointDistances>shellSpacing[i-1],pointDistances<shellSpacing[i])]
        else: # If we're in the outer shell, choose all points in that shell and greater. This was a design choice for our system
            shellPoints = Ps[:,pointDistances>shellSpacing[i]]
        # Now, get PCA on points
        D = shellPoints.dot(shellPoints.T)
        eigs = np.sort(np.real(np.linalg.eigvals(D)))
        hist = np.vstack((hist,[eigs])) # Put the PCA eigvenvalues in the histogram
    return hist.flatten() #Flatten the 2D histogram to a 1D array

#Purpose: To create shape histogram of the pairwise Euclidean distances between
#randomly sampled points in the point cloud
#Inputs: Ps (3 x N point cloud), Ns (3 x N array of normals) (not needed here
#but passed along for consistency), DMax (Maximum distance to consider), 
#NBins (number of histogram bins), NSamples (number of pairs of points sample
#to compute distances)
def getD2Histogram(Ps, Ns, DMax, NBins, NSamples):
    numPoints = Ps.shape[1]
    r1 = np.random.random_integers(0, numPoints-1, NSamples)
    r2 = np.random.random_integers(0, numPoints-1, NSamples)
    distances = []
    for i in range(NSamples):
        d = np.linalg.norm(Ps[:, r1[i]] - Ps[:, r2[i]])
        if(d > DMax):
            continue
        distances.append(d)
        
    hist = np.histogram(distances, NBins)[0]
    return hist

#Purpose: To create shape histogram of the angles between randomly sampled
#triples of points
#Inputs: Ps (3 x N point cloud), Ns (3 x N array of normals) (not needed here
#but passed along for consistency), NBins (number of histogram bins), 
#NSamples (number of triples of points sample to compute angles)
def getA3Histogram(Ps, Ns, NBins, NSamples):
    numPoints = Ps.shape[1]
    r1 = np.random.random_integers(0, numPoints-1, NSamples)
    r2 = np.random.random_integers(0, numPoints-1, NSamples)
    r3 = np.random.random_integers(0, numPoints-1, NSamples)

    angles = []
    for i in range(NSamples):
        ba = Ps[:, r1[i]] - Ps[:, r2[i]]
        bc = Ps[:, r3[i]] - Ps[:, r2[i]]
        baNorm = np.linalg.norm(ba)
        bcNorm = np.linalg.norm(bc)
        if baNorm == 0 or bcNorm == 0:
            angles.append(0)
            continue
        cosTheta = ba.dot(bc) / (baNorm * bcNorm)
        if(cosTheta*cosTheta >= 1):
            angles.append(0)
            continue
        angle = np.arccos(cosTheta)
        angles.append(angle)
        
    hist = np.histogram(angles, NBins)[0]
    return hist

#Purpose: To create the Extended Gaussian Image by binning normals to
#sphere directions after rotating the point cloud to align with its principal axes
#Inputs: Ps (3 x N point cloud) (use to compute PCA), Ns (3 x N array of normals), 
#SPoints: A 3 x S array of points sampled evenly on the unit sphere used to 
#bin the normals
#def getEGIHistogram(Ps, Ns, SPoints):
def getEGIHistogram(Ps, Ns, SPoints):
    A = (Ps).dot(Ps.T)
    w, v = np.linalg.eig(A)
    
    maxEIndex = np.argmax(w)
    maxEVector = v[maxEIndex]
    
    #sorts the eigenvectors smallest to largest
    idx = w.argsort()[::1]   
    w = w[idx]
    v = v[:,idx]
    
    B = np.array([[1,0,0],[0,1,0],[0,0,1]])
    
    # need to align the v axis to the coordinate axis, but there are 4 right-handed
    # orientations : need to check all four of them ? maybe later
    # R dot A = B
    # R = B dot A^-1 = B dot A.T (because A and B are orthonormal)
    # note that this assumes our vectors are columns, and currently they're in rows
    R = B.dot(v.T)
    
    # now lets rotate all of the normals by this rotation matrix
    rotatedNs = R.dot(Ns)

    dots = np.dot(rotatedNs.T, SPoints)
    # Then find the index of the point in SPoints that yields the largest dot product for each Ps point
    maximums = np.argmax(dots,axis=1).T
    #Create a 2D histogram that is 1 x NSectors
    NSectors = SPoints.shape[1]
    hist = np.zeros((SPoints.shape[1])) 
    # Go through each sector and create a histogram
    for i in range(NSectors):
        sectorElems = Ps[:,maximums==i] # Select every element in the given sector
        hist[i]=sectorElems.shape[1] # Add the # elems to the histogram
    return hist

#Purpose: To create an image which stores the amalgamation of rotating
#a bunch of planes around the largest principal axis of a point cloud and 
#projecting the points on the minor axes onto the image.
#Inputs: Ps (3 x N point cloud), Ns (3 x N array of normals, not needed here),
#NAngles: The number of angles between 0 and 2*pi through which to rotate
#the plane, Extent: The extent of each axis, Dim: The number of pixels along
#each minor axis
def getSpinImage(Ps, Ns, NAngles, Extent, Dim):
    A = (Ps).dot(Ps.T)
    w, v = np.linalg.eig(A)
    
    maxEIndex = np.argmax(w)
    maxEVector = v[maxEIndex]
    
    idx = w.argsort()[::1]   
    w = w[idx]
    v = v[:,idx]
    
    #print w
    #print v
        
    # need to align the v axis to the coordinate axis
    # R dot A = B
    # R = B dot A^-1 = B dot A.T (because A and B are orthonormal)
    B = np.array([[1,0,0],[0,1,0],[0,0,1]]) # aligns largest eigenvector to z axis
    R = B.dot(v.T)
    
    RotatedPs = R.dot(Ps)
    
    A = (RotatedPs).dot(RotatedPs.T)
    w, v = np.linalg.eig(A)
    
    maxEIndex = np.argmax(w)
    maxEVector = v[maxEIndex]
    
    idx = w.argsort()[::1]   
    w = w[idx]
    v = v[:,idx]
    
    #print w
    #print v
    
    theta = 2*3.14159654/NAngles
    c = np.cos(theta)
    s = np.sin(theta)
    Rz = np.array([[c,-s,0],[s,c,0],[0,0,1]])
    
    hist = np.zeros((Dim, Dim))
    for i in range(NAngles):
        c = np.cos(theta*i)
        s = np.sin(theta*i)
        Rz = np.array([[c,-s,0],[s,c,0],[0,0,1]])
        r = Rz.dot(RotatedPs)
        h = np.histogram2d(r[0,:], r[1,:], Dim, [[-Extent, Extent], [-Extent, Extent]])[0]
        hist = hist + h
        
    #plt.imshow(hist)
    #plt.show()
    return hist.flatten()


#Purpose: To create a histogram of spherical harmonic magnitudes in concentric spheres
#Inputs: Ps (3 x N point cloud), Ns (3 x N array of normals, not used here), RMax: maximum 
# radius, NHarmonics: the number of spherical harmonics, NSpheres: the number of concentric 
# spheres to take
def getSphericalHarmonicMagnitudes(Ps, Ns, RMax, NHarmonics, NSpheres):
    #m = PolyMesh()
    #m.loadFile("models_off/biplane0.off") #Load a mesh
    #(Ps, Ns) = samplePointCloud(m, 20000) #Sample 20,000 points and associated normals

    res = 2 # this can (should be?) changed
    
    SPoints = getSphereSamples(res)
    B = SPoints.shape[1]
    
    # point n = SPoints[:,n]
    Bs = np.zeros((B, 2))
    for i in range(0, B):
        x = SPoints[0, i]
        y = SPoints[1, i]
        z = SPoints[2, i]
        # construct the Bs matrix
        if(x == 0 and y>0):
            theta = np.pi / 2
        elif (x == 0):
            theta = -np.pi/2
        else:
            theta = np.arctan(y/x)
            
        if(x < 0):
            theta = theta + np.pi
        
        phi = np.arccos(z/np.sqrt(x*x + y*y + z*z))
        Bs[i] = np.array([theta, phi])
    
    # Calculate the spherical harmonics matrix
    F = np.zeros((NHarmonics, B))
    for m in range(0, NHarmonics):
        # the paper ignores "l" (the degree) by summing up degrees from -m to m
        F0 = np.absolute(sph_harm(np.abs(-m), m, Bs[:,0] , Bs[:,1]))
        for l in range(-m+1, m+1):
            F0 += np.absolute(sph_harm(np.abs(l), m, Bs[:,0] , Bs[:,1]))
        F[m] = F0
    
    # now we compute h (a BxN matrix, B = number of sampled sphere points, N = number of shells)(see getShapeShellHistogram for detailed comments)
    dots = np.dot(Ps.T, SPoints)
    maximums = np.argmax(dots,axis=1).T
    h = np.zeros((NSpheres, 0)) 
    for i in range(B):
        sectorElems = Ps[:,maximums==i] # Select every element in the given sector
        sectorHistogram = np.array(np.histogram(np.linalg.norm(sectorElems,axis=0), NSpheres, (0,RMax))[0])[np.newaxis].T # Create a histogram for the sector in Column form
        h=np.hstack((h,sectorHistogram)) # Add the column to the histogram
    
    h = h.T
    
    # finally, we compute the spherical harmonic coefficient matrix
    H = F.dot(h)

    return H
    
#Purpose: Utility function for wrapping around the statistics functions.
#Inputs: PointClouds (a python list of N point clouds), Normals (a python
#list of the N corresponding normals), histFunction (a function
#handle for one of the above functions), *args (addditional arguments
#that the descriptor function needs)
#Returns: AllHists (A KxN matrix of all descriptors, where K is the length
#of each descriptor)
def makeAllHistograms(PointClouds, Normals, histFunction, *args):
    N = len(PointClouds)
    #Call on first mesh to figure out the dimensions of the histogram
    print "Computing histogram 1 of %i..."%(N)
    h0 = histFunction(PointClouds[0], Normals[0], *args)
    K = h0.size
    AllHists = np.zeros((K, N))
    AllHists[:, 0] = h0
    for i in range(1, N):
        print "Computing histogram %i of %i..."%(i+1, N)
        AllHists[:, i] = histFunction(PointClouds[i], Normals[i], *args)
    return AllHists

#########################################################
##              HISTOGRAM COMPARISONS                  ##
#########################################################

#Normalizes each histogram so that sum(h[k]) = 1
#AllHists: (K x N matrix of histograms, where K is the length
#of each histogram and N is the number of point clouds)
#Returns: hists (a K x N matrix of the normalized histograms)
def normalizeHists(AllHists):
	sums = np.sum(AllHists,axis=0)
	hists = 1.0 * AllHists / sums
	return hists

#Purpose: To compute the euclidean distance between a set
#of histograms
#Inputs: AllHists (K x N matrix of histograms, where K is the length
#of each histogram and N is the number of point clouds)
#Returns: D (An N x N matrix, where the ij entry is the Euclidean
#distance between the histogram for point cloud i and point cloud j)
def compareHistsEuclidean(AllHists):
    hists = normalizeHists(AllHists);
    aa = np.sum(np.multiply(hists,hists),0);
    bb = np.sum(np.multiply(hists,hists),0);
    ab = np.dot(hists.T,hists);
    
    D2 = (aa - (2*ab).T).T + bb;
    D = np.sqrt(D2)
    return D

#Purpose: To compute the cosine distance between a set
#of histograms
#Inputs: AllHists (K x N matrix of histograms, where K is the length
#of each histogram and N is the number of point clouds)
#Returns: D (An N x N matrix, where the ij entry is the cosine
#distance between the histogram for point cloud i and point cloud j)
def compareHistsCosine(AllHists):
    hists = normalizeHists(AllHists);
    N = hists.shape[1]
    numerator = np.dot(hists.T,hists);    
    
    norms = np.linalg.norm(hists,axis=0)
    norms = norms.reshape(1,N)
    denominator = np.dot(norms.T,norms)
    
    D = np.divide(numerator,denominator)
    return 1-D

#Purpose: To compute the chi square distance between a set
#of histograms
#Inputs: AllHists (K x N matrix of histograms, where K is the length
#of each histogram and N is the number of point clouds)
#Returns: D (An N x N matrix, where the ij entry is the chi squared
#distance between the histogram for point cloud i and point cloud j)
def compareHistsChiSquared(AllHists):
    hists = normalizeHists(AllHists);
    sub = hists[:, None, :].T - hists.T # used this resource: http://stackoverflow.com/questions/32473635/how-do-i-calculate-all-pairs-of-vector-differences-in-numpy
    sub2 = np.square(sub)
    add = hists[:, None, :].T + hists.T
    div = np.divide(1.0*sub2,add)
    D = 0.5*np.sum(div,2)
    return D

#Purpose: To compute the 1D Earth mover's distance between a set
#of histograms (note that this only makes sense for 1D histograms)
#Inputs: AllHists (K x N matrix of histograms, where K is the length
#of each histogram and N is the number of point clouds)
#Returns: D (An N x N matrix, where the ij entry is the earth mover's
#distance between the histogram for point cloud i and point cloud j)
def compareHistsEMD1D(AllHists):
    hists = normalizeHists(AllHists);
    histsC = np.cumsum(hists,axis=0)
    sub = histsC[:, None, :].T - histsC.T
    absSub = np.abs(sub)
    D = np.sum(absSub,2)
    return D


#########################################################
##              CLASSIFICATION CONTEST                 ##
#########################################################

#Purpose: To implement your own custom distance matrix between all point
#clouds for the point cloud clasification contest
#Inputs: PointClouds, an array of point cloud matrices, Normals: an array
#of normal matrices
#Returns: D: A N x N matrix of distances between point clouds based
#on your metric, where Dij is the distnace between point cloud i and point cloud j
def getMyShapeDistances(PointClouds, Normals):
    #TODO: Finish this
    #This is just an example, but you should experiment to find which features
    #work the best, and possibly come up with a weighted combination of 
    #different features
    HistsD2 = makeAllHistograms(PointClouds, Normals, getD2Histogram, 3.0, 30, 100000)
    DEuc = compareHistsEuclidean(HistsD2)
    return DEuc

#########################################################
##                     EVALUATION                      ##
#########################################################

#Purpose: To return an average precision recall graph for a collection of
#shapes given the similarity scores of all pairs of histograms.
#Inputs: D (An N x N matrix, where the ij entry is the earth mover's distance
#between the histogram for point cloud i and point cloud j).  It is assumed
#that the point clouds are presented in contiguous chunks of classes, and that
#there are "NPerClass" point clouds per each class (for the dataset provided
#there are 10 per class so that's the default argument).  So the program should
#return a precision recall graph that has 9 elements
#Returns PR, an (NPerClass-1) length array of average precision values for all 
#recalls
def getPrecisionRecall(D, NPerClass = 10):
    N = D.shape[0]
    PR = np.zeros(NPerClass-1)
    ind = np.argsort(D)
    
    rowIndex = 0
    for row in ind:
        print "row index:",rowIndex
        numFound = 0
        numSearched = 0
        position = 0

        while position < N:
            # skip the value if it is being compared to itself
            if rowIndex == row[position]:
                position += 1
                continue
            
            numSearched += 1
            
            originalPosition = row[position]
            if (rowIndex / NPerClass) == (originalPosition / NPerClass):
                numFound += 1
                toAdd = 1.0 * numFound / numSearched
                PR[numFound-1] += toAdd
            position += 1
        
        print "\tPR:",PR
        
        rowIndex += 1
    
    print PR
    
    # at the end, divide PR by rowIndex (number of rows)
    PR = 1.0 * PR / rowIndex
    
    return PR
    
def teehee():
    #return np.array([[1,3,5],[2,4,6]])
    return np.array([[1,3],[2,5],[3,7]])

def runDistanceMetricsExperiments():
    SPoints = getSphereSamples(2)
    #HistsEGI = makeAllHistograms(PointClouds, Normals, getEGIHistogram, SPoints)
    #HistsSpin = makeAllHistograms(PointClouds, Normals, getSpinImage, 100, 2, 40)
    #HistsShellSector = makeAllHistograms(PointClouds, Normals, getShapeShellHistogram, 10, 2, SPoints)
    #HistsShell = makeAllHistograms(PointClouds, Normals, getShapeHistogram, 10, 2)
    HistsA3 = makeAllHistograms(PointClouds, Normals, getA3Histogram, 30, 100000)
    #HistsD2 = makeAllHistograms(PointClouds, Normals, getD2Histogram, 3.0, 30, 100000)

    Hists = HistsA3

    DSpin1 = compareHistsEuclidean(Hists)
    DSpin2 = compareHistsCosine(Hists)
    DSpin3 = compareHistsChiSquared(Hists)
    DSpin4 = compareHistsEMD1D(Hists)

    PRSpin1 = getPrecisionRecall(DSpin1)
    PRSpin2 = getPrecisionRecall(DSpin2)
    PRSpin3 = getPrecisionRecall(DSpin3)
    PRSpin4 = getPrecisionRecall(DSpin4)
 
    recalls = np.linspace(1.0/9.0, 1.0, 9)
    plt.plot(recalls, PRSpin1, 'c', label='Euclidean')
    plt.hold(True)    
    plt.plot(recalls, PRSpin2, 'k', label='Cosine')
    plt.plot(recalls, PRSpin3, 'r', label='ChiSquared')
    plt.plot(recalls, PRSpin4, 'b', label='EMD1D')
    plt.xlabel('Recall')
    plt.ylabel('Precision')
    plt.legend()
    plt.show()

def runExperiments():
    SPoints = getSphereSamples(2)
    HistsShell = makeAllHistograms(PointClouds, Normals, getShapeHistogram, 10, 2)
    HistsShellSector= makeAllHistograms(PointClouds, Normals, getShapeShellHistogram, 10, 2, SPoints)
    HistsShellPCA = makeAllHistograms(PointClouds, Normals, getShapeHistogramPCA, 10, 2)
    HistsEGI = makeAllHistograms(PointClouds, Normals, getEGIHistogram, SPoints)
    HistsA3 = makeAllHistograms(PointClouds, Normals, getA3Histogram, 30, 100000)
    HistsD2 = makeAllHistograms(PointClouds, Normals, getD2Histogram, 3.0, 30, 100000)
    HistsSpin = makeAllHistograms(PointClouds, Normals, getSpinImage, 100, 2, 40)

    DS = compareHistsEuclidean(HistsShell)
    DSS = compareHistsEuclidean(HistsShellSector)
    DSPCA = compareHistsEuclidean(HistsShellPCA)
    DEGI = compareHistsEuclidean(HistsEGI)
    DA3 = compareHistsEuclidean(HistsA3)
    DD2 = compareHistsEuclidean(HistsD2)
    DSpin = compareHistsEuclidean(HistsSpin)

    PRS = getPrecisionRecall(DS)
    PRSS = getPrecisionRecall(DSS)
    PRSPCA = getPrecisionRecall(DSPCA)
    PREGI = getPrecisionRecall(DEGI)
    PRA3 = getPrecisionRecall(DA3)
    PRD2 = getPrecisionRecall(DD2)
    PRSpin = getPrecisionRecall(DSpin)
 
    recalls = np.linspace(1.0/9.0, 1.0, 9)
    plt.plot(recalls, PRS, 'c', label='Shell')
    plt.hold(True)
    plt.plot(recalls, PRSS, 'g', label='ShellSector')
    plt.plot(recalls, PRSPCA, 'y', label='ShellPCA')
    plt.plot(recalls, PREGI, 'm', label='EGI')    
    plt.plot(recalls, PRA3, 'k', label='A3')
    plt.plot(recalls, PRD2, 'r', label='D2')
    plt.plot(recalls, PRSpin, 'b', label='Spin')
    plt.xlabel('Recall')
    plt.ylabel('Precision')
    plt.legend()
    plt.show()

def runShellNumberExperiments():
    HistsShell = makeAllHistograms(PointClouds, Normals, getShapeHistogram, 1, 2)
    DS = compareHistsEuclidean(HistsShell)
    PRS = getPrecisionRecall(DS)
    recalls = np.linspace(1.0/9.0, 1.0, 9)
    plt.plot(recalls, PRS, 'c', label='Shells = 1')
    plt.hold(True)
    
    HistsShell = makeAllHistograms(PointClouds, Normals, getShapeHistogram, 2, 2)
    DS = compareHistsEuclidean(HistsShell)
    PRS = getPrecisionRecall(DS)
    plt.plot(recalls, PRS, 'g', label='Shells = 2')

    
    HistsShell = makeAllHistograms(PointClouds, Normals, getShapeHistogram, 4, 2)
    DS = compareHistsEuclidean(HistsShell)
    PRS = getPrecisionRecall(DS)
    plt.plot(recalls, PRS, 'y', label='Shells = 4')
    
    HistsShell = makeAllHistograms(PointClouds, Normals, getShapeHistogram, 8, 2)
    DS = compareHistsEuclidean(HistsShell)
    PRS = getPrecisionRecall(DS)
    plt.plot(recalls, PRS, 'm', label='Shells = 8')
    
    HistsShell = makeAllHistograms(PointClouds, Normals, getShapeHistogram, 16, 2)
    DS = compareHistsEuclidean(HistsShell)
    PRS = getPrecisionRecall(DS)
    plt.plot(recalls, PRS, 'k', label='Shells = 16')
    
    HistsShell = makeAllHistograms(PointClouds, Normals, getShapeHistogram, 32, 2)
    DS = compareHistsEuclidean(HistsShell)
    PRS = getPrecisionRecall(DS)
    plt.plot(recalls, PRS, 'r', label='Shells = 32')
    
    HistsShell = makeAllHistograms(PointClouds, Normals, getShapeHistogram, 64, 2)
    DS = compareHistsEuclidean(HistsShell)
    PRS = getPrecisionRecall(DS)
    plt.plot(recalls, PRS, 'b', label='Shells = 64')
   
    plt.xlabel('Recall')
    plt.ylabel('Precision')
    plt.legend()
    plt.show()
    
def runD2SampleExperiments():
    HistsD2 = makeAllHistograms(PointClouds, Normals, getD2Histogram, 3.0, 30, 512)
    DD2 = compareHistsEuclidean(HistsD2)
    PRD2 = getPrecisionRecall(DD2)
    recalls = np.linspace(1.0/9.0, 1.0, 9)
    plt.plot(recalls, PRD2, 'c', label='Samples = 512')
    plt.hold(True)
    
    HistsD2 = makeAllHistograms(PointClouds, Normals, getD2Histogram, 3.0, 30, 1024)
    DD2 = compareHistsEuclidean(HistsD2)
    PRD2 = getPrecisionRecall(DD2)
    plt.plot(recalls, PRD2, 'g', label='Samples = 1024')
     
    HistsD2 = makeAllHistograms(PointClouds, Normals, getD2Histogram, 3.0, 30, 2048)
    DD2 = compareHistsEuclidean(HistsD2)
    PRD2 = getPrecisionRecall(DD2)
    plt.plot(recalls, PRD2, 'y', label='Samples = 2048')
    
    HistsD2 = makeAllHistograms(PointClouds, Normals, getD2Histogram, 3.0, 30, 4096)
    DD2 = compareHistsEuclidean(HistsD2)
    PRD2 = getPrecisionRecall(DD2)
    plt.plot(recalls, PRD2, 'm', label='Samples = 4096')
    
    HistsD2 = makeAllHistograms(PointClouds, Normals, getD2Histogram, 3.0, 30, 8192)
    DD2 = compareHistsEuclidean(HistsD2)
    PRD2 = getPrecisionRecall(DD2)
    plt.plot(recalls, PRD2, 'k', label='Samples = 8192')
    
    HistsD2 = makeAllHistograms(PointClouds, Normals, getD2Histogram, 3.0, 30, 16384)
    DD2 = compareHistsEuclidean(HistsD2)
    PRD2 = getPrecisionRecall(DD2)
    plt.plot(recalls, PRD2, 'r', label='Samples = 16384')

    HistsD2 = makeAllHistograms(PointClouds, Normals, getD2Histogram, 3.0, 30, 32768)
    DD2 = compareHistsEuclidean(HistsD2)
    PRD2 = getPrecisionRecall(DD2)
    plt.plot(recalls, PRD2, 'b', label='Samples = 32768')
    
    plt.xlabel('Recall')
    plt.ylabel('Precision')
    plt.legend()
    plt.show()
    
def runEGIExperiments():
    SPoints = getSphereSamples(1)
    HistsEGI = makeAllHistograms(PointClouds, Normals, getEGIHistogram, SPoints)
    DEGI = compareHistsEuclidean(HistsEGI)
    PREGI = getPrecisionRecall(DEGI)
    recalls = np.linspace(1.0/9.0, 1.0, 9)
    plt.plot(recalls, PREGI, 'c', label='Resolution=1')
    plt.hold(True)
    
    SPoints = getSphereSamples(2)
    HistsEGI = makeAllHistograms(PointClouds, Normals, getEGIHistogram, SPoints)
    DEGI = compareHistsEuclidean(HistsEGI)
    PREGI = getPrecisionRecall(DEGI)
    plt.plot(recalls, PREGI, 'g', label='Resolution=2')

    SPoints = getSphereSamples(3)
    HistsEGI = makeAllHistograms(PointClouds, Normals, getEGIHistogram, SPoints)
    DEGI = compareHistsEuclidean(HistsEGI)
    PREGI = getPrecisionRecall(DEGI)
    plt.plot(recalls, PREGI, 'y', label='Resolution=3')

    SPoints = getSphereSamples(4)
    HistsEGI = makeAllHistograms(PointClouds, Normals, getEGIHistogram, SPoints)
    DEGI = compareHistsEuclidean(HistsEGI)
    PREGI = getPrecisionRecall(DEGI)
    plt.plot(recalls, PREGI, 'm', label='Resolution=4')

    SPoints = getSphereSamples(5)
    HistsEGI = makeAllHistograms(PointClouds, Normals, getEGIHistogram, SPoints)
    DEGI = compareHistsEuclidean(HistsEGI)
    PREGI = getPrecisionRecall(DEGI)
    plt.plot(recalls, PREGI, 'k', label='Resolution=5')

    plt.xlabel('Recall')
    plt.ylabel('Precision')
    plt.legend()
    plt.show()
    
def runRandomComparisonExperiments():
    HistsShell = makeAllHistograms(PointClouds, Normals, getShapeHistogram, 10, 2)
    DS = compareHistsEuclidean(HistsShell)
    PRS = getPrecisionRecall(DS)
 
    HistsRandom = np.random.random(HistsShell.shape)
    DR = compareHistsEuclidean(HistsRandom)
    PRR = getPrecisionRecall(DR)
 
    recalls = np.linspace(1.0/9.0, 1.0, 9)
    plt.plot(recalls, PRS, 'c', label='Shell')
    plt.hold(True)
    plt.plot(recalls, PRR, 'g', label='Random')

    plt.xlabel('Recall')
    plt.ylabel('Precision')
    plt.legend()
    plt.show()
    
#########################################################
##                     MAIN TESTS                      ##
#########################################################


if __name__ == '__main__':  
    NRandSamples = 10000 #You can tweak this number
    np.random.seed(100) #For repeatable results randomly sampling
    #Load in and sample all meshes
    
    numClasses = len(POINTCLOUD_CLASSES)
    #numClasses = 2 ## NC - comment this in if you don't want to load all 20 classes
    
    PointClouds = []
    Normals = []
    for i in range(numClasses):
        print "LOADING CLASS %i of %i..."%(i+1, numClasses)
        PCClass = []
        for j in range(NUM_PER_CLASS):
            m = PolyMesh()
            filename = "models_off/%s%i.off"%(POINTCLOUD_CLASSES[i], j)
            print "Loading ", filename
            m.loadOffFileExternal(filename)
            (Ps, Ns) = samplePointCloud(m, NRandSamples)
            PointClouds.append(Ps)
            Normals.append(Ps)

    #runRandomComparisonExperiments()
    #runExperiments()
    runDistanceMetricsExperiments()
    
    #TODO: Finish this, run experiments.  Also in the above code, you might
    #just want to load one point cloud and test your histograms on that first
    #so you don't have to wait for all point clouds to load when making
    #minor tweaks
