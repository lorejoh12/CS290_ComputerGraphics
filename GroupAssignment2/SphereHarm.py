import sys
sys.path.append("S3DGLPy")
from Primitives3D import *
from PolyMesh import *
import numpy as np
from scipy.special import sph_harm
import json

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

#Purpose: To create shape histogram with concentric spherical shells and
#sectors within each shell, sorted in decreasing order of number of points
#Inputs: Ps (3 x N point cloud), Ns (3 x N array of normals) (not needed here
#but passed along for consistency), NShells (number of shells), 
#RMax (maximum radius), SPoints: A 3 x S array of points sampled evenly on
#RMin (minimum radius), 
#the unit sphere (get these with the function "getSphereSamples")
def getShapeShellHistogram(Ps, Ns, NShells,RMin, RMax, SPoints):
    # Remove any points below RMin from the origin
    pointDistances = np.linalg.norm(Ps, axis=0)
    filteredElems = Ps[:,pointDistances>=RMin]

    NSectors = SPoints.shape[1] #A number of sectors equal to the number of
    #points sampled on the sphere
    # First, dot product all points in Ps against all points in SPoints.
    # We're doing a matrix multiplication of Ps^T with SPoints. This will give an NxM matrix
    # Where the rows are the dot product magnitudes of each point of Ps with all the points in SPoints
    dots = np.dot(filteredElems.T, SPoints)
    # Then find the index of the point in SPoints that yields the largest dot product for each Ps point
    maximums = np.argmax(dots,axis=1).T
    #Create a 2D histogram that is NShells x NSectors
    NSectors = SPoints.shape[1]
    hist = np.zeros((NShells, 0)) 
    # Go through each sector and create a histogram
    for i in range(NSectors):
        sectorElems = filteredElems[:,maximums==i] # Select every element in the given sector
        sectorHistogram=np.array(np.histogram(np.linalg.norm(sectorElems,axis=0),NShells,(RMin,RMax))[0])[np.newaxis].T # Create a histogram for the sector in Column form
        hist=np.hstack((hist,sectorHistogram)) # Add the column to the histogram
      
    return hist

#Purpose: To create a shape shell histogram using differential sampling of concentric spherical shells.
# meant to differentiate the sampling volume difference for further away sectors
# Note: Rmax is automatically calculated by this algorithm
#Inputs: Ps (3 x N point cloud), Ns (3 x N array of normals) (not needed here
#but passed along for consistency), NShells (number of shells), 
#RMax (maximum radius), SPointRange: A range of spoint sampling rates to be spread linearly
#amongst the shells in RMax. i.e. if 4 numbers are given in the range, the first will be used for 0-RMax/4
#Returns an array of tuples of (innerShellRadius, outerShellRadius, SValue, ShapeShellHistogram, ShellWidth)

def getDifferentialShapeShellHistogram(Ps, Ns, NShells, SPointRange):
    pointDistances = np.linalg.norm(Ps, axis=0)
    Rmax = np.amax(pointDistances,) 
    # Figure out how many different sampling rates we're dealing with
    numSampleIncrements = len(SPointRange)
    distanceIncrement = 1.0*Rmax/numSampleIncrements
    shellIncrement = 1.0*NShells/numSampleIncrements
    if (math.ceil(shellIncrement)!=math.floor(shellIncrement)):
        print "Warning: NShells should be integer divisible by the number of sample increments"
    shellIncrement=math.floor(shellIncrement)
    histogram = []
    # Run the shape shell histogram on each shell
    for i in range(numSampleIncrements):
        RMin=i*distanceIncrement
        RMax=(i+1)*distanceIncrement
        sampleHistogram=getShapeShellHistogram(Ps,Ns,shellIncrement,RMin, RMax, getSphereSamples(SPointRange[i]))
        histogram.append((RMin,RMax,SPointRange[i],sampleHistogram, distanceIncrement/shellIncrement))
    return histogram


#Returns the squared distance between two 3D points
def sqDistBetweenTwoPoints(p1, p2):
    a = (p2[0] - p1[0])**2
    b = (p2[1] - p1[1])**2
    c = (p2[2] - p1[2])**2 
    return a + b + c;

#Returns the angle between two unit vectors
def getAngleBetweenTwoUnitVectors(v1, v2):
    return np.arccos(v1.dot(v2))

#Returns approximately the smallest angle between any two points
#Assumes that SPoints is evenly sampled for a sphere
def getMinimumAngleAmongPoints(SPoints):
    firstPoint = SPoints[:, 0]
    minDist = 5 # assume 2 is the maximum distance btwn any two points (b/c unit sphere)
    # go through each point and find the point closest to the first point
    for i in range(1, SPoints.shape[1]): #1 is intentional!
        point = SPoints[:, i]
        sqDist = sqDistBetweenTwoPoints(firstPoint,point)
        if sqDist < minDist:
            minDist = sqDist
            neighborPoint = point
    # now that we have the first point and the closest point to it, return the angle between them
    return getAngleBetweenTwoUnitVectors(firstPoint, neighborPoint)
    
#Returns a point cloud from a 2D histogram
#The point cloud is generated with random sampling
# numSamples is the overall number of sectors to use in the sampling
# hist is the histogram in the format ((Rmin,Rmax,SValue, Histogram)) for a special shell sampled with differential sampling
def getSamplePointsFrom2DHistogram(numSamples, histogram):
    numSections = len(histogram)
    pointsPerSection = numSamples/numSections
    # initialize Ps
    Ps = np.zeros((3,numSamples))
    sampleNum=0
    for j in range(numSections):
        shellHist = histogram[j]
        RMin = shellHist[0]
        RMax = shellHist[1]
        SVal = shellHist[2]
        hist = shellHist[3]
        shellWidth = shellHist[4]
        SPoints = getSphereSamples(SVal)
        # this is the angle between the two nearest sampled points
        # which is the same as the angle of the sampled cone around the vectors
        angle = getMinimumAngleAmongPoints(SPoints)
        # get flattened matrix of probability from the densities
        flatProbabilities = (1.0*hist/np.sum(hist)).flatten()
        # get random samples of indexes of sectors
        randomSectorIndexes = np.random.choice(len(flatProbabilities), pointsPerSection, p=flatProbabilities)
        histWidth = hist.shape[1]
        
        for i in range(0, len(randomSectorIndexes)):
            index = randomSectorIndexes[i]
            # calculate the INDEXES of the shell and sector
            shellIndex = index / histWidth # index of shell in hist
            sectorIndex = index % histWidth  # index of sector in hist
            
            # Radius range, phi range, and theta range, and randomly sample in them to determine new point location
            # shell is zero indexed
            radiusMin = RMin+shellIndex*shellWidth
            radiusMax = RMin+(shellIndex+1)*shellWidth

            # phi, theta location somewhere around the sphere sampled vector
            sector = SPoints[:, sectorIndex]
            [_theta, _phi] = convertRectangularToPolar(sector[0], sector[1], sector[2])
            phiMin = _phi - angle/2.0
            phiMax = _phi + angle/2.0
            thetaMin = _theta - angle/2.0
            thetaMax = _theta + angle/2.0
            
            # randomly sample these
            randomRadius = (radiusMax - radiusMin) * np.random.random_sample() + radiusMin
            randomPhi = (phiMax - phiMin) * np.random.random_sample() + phiMin
            randomTheta = (thetaMax - thetaMin) * np.random.random_sample() + thetaMin
            
            # convert returns [x, y, z]
            ret = convertPolarToRectangular(randomRadius, randomPhi, randomTheta)
            Ps[0,sampleNum] = ret[0]
            Ps[1,sampleNum] = ret[1]
            Ps[2,sampleNum] = ret[2]
            sampleNum+=1
        
    return Ps
        
def convertRectangularToPolar(x, y, z):
    if(x == 0 and y>0):
        theta = np.pi / 2
    elif (x == 0):
        theta = -np.pi/2
    else:
        theta = np.arctan(y/x)

    if(x < 0):
        theta = theta + np.pi

    phi = np.arccos(z/np.sqrt(x*x + y*y + z*z))
    
    return [theta,phi]

def convertPolarToRectangular(r, phi, theta):
    x = r * np.cos(theta) * np.sin(phi)
    y = r * np.sin(theta) * np.sin(phi)
    z = r * np.cos(phi)
    return [x, y, z]

def computeUnraveledFMatrix(NHarmonics, res, magnitude = True):
    
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
    numCols = NHarmonics*NHarmonics + 2*NHarmonics + 1
    F = np.zeros((numCols, B))+0j
    i = 0
    for m in range(0, NHarmonics):
        for l in range(-m, m+1):
            if(magnitude):
                F0 = np.absolute(sph_harm(l, m, Bs[:,0] , Bs[:,1]))
            else:
                F0 = sph_harm(l, m, Bs[:,0] , Bs[:,1])
            F[i, :] = F0
            i = i+1

    return F
 
 
#Purpose: To create a histogram of spherical harmonic magnitudes in concentric spheres. RMax is calculated based on input data
#Inputs: Ps (3 x N point cloud), Ns (3 x N array of normals, not used here), 
# NHarmonics: the number of spherical harmonics, 
# WindowSize: The size of each shell's window, HopSize: the hop between windows
def getSphericalHarmonicMagnitudes(Ps, Ns, NHarmonics, NSpheres, SPointRange):
    # Calculate RMax based on the input shape
    pointDistances = np.linalg.norm(Ps, axis=0)
    Rmax = np.amax(pointDistances,) 
    # Figure out how many different sector sample rates we have
    numSampleIncrements = len(SPointRange)
    distanceIncrement = 1.0*Rmax/numSampleIncrements
    shellIncrement = 1.0*NSpheres/numSampleIncrements
    if (math.ceil(shellIncrement)!=math.floor(shellIncrement)):
        print "Warning: NShells should be integer divisible by the number of sample increments"
    shellIncrement=math.floor(shellIncrement)
    allH = []
    allHApprox = []
    # Go through all the shells and calculate the spherical harmonic on each.
    for j in range(numSampleIncrements):
        res=SPointRange[j]
        RMin=j*distanceIncrement
        RMax=(j+1)*distanceIncrement
        F = computeUnraveledFMatrix(NHarmonics, res, False)
        SPoints = getSphereSamples(res)
        B = SPoints.shape[1]
        # Filter out the points below the min radius
        filteredElems = Ps[:,pointDistances>=RMin]
        # now we compute h (a BxN matrix, B = number of sampled sphere points, N = number of shells)(see getShapeShellHistogram for detailed comments)
        dots = np.dot(filteredElems.T, SPoints)
        maximums = np.argmax(dots,axis=1).T
        h = np.zeros((shellIncrement, 0)) 
        print "computing h"
        for i in range(B):
            sectorElems = filteredElems[:,maximums==i] # Select every element in the given sector
            sectorHistogram = np.array(np.histogram(np.linalg.norm(sectorElems,axis=0), shellIncrement, (RMin,RMax))[0])[np.newaxis].T # Create a histogram for the sector in Column form
            h=np.hstack((h,sectorHistogram)) # Add the column to the histogram
        
        h = h.T

        H = F.dot(h)

        h_approx = np.absolute(F.T.dot(H))
                
        h_approx = np.floor(np.sum(h)*10.0*h_approx/np.sum(h_approx))/10.0
        h_approx = np.sum(h)*1.0*h_approx/np.sum(h_approx)
        h_approx = h_approx.T
        h=h.T
        
        # Add the shell's data to the overall data output
        allH.append((RMin,RMax,res,h,distanceIncrement/shellIncrement))
        allHApprox.append((RMin,RMax,res,h_approx,distanceIncrement/shellIncrement))



    # return H.flatten()
    return [allHApprox, allH]

    
# Calculate the size of a histogram. This is used for extimating space savings
def calculateHistogramSize(histogram):
    size=0
    for i in range(len(histogram)):
        size+=histogram[i][3].shape[0]*histogram[i][3].shape[1]
    return size


#########################################################
##                     MAIN TESTS                      ##
#########################################################

if __name__ == '__main__':  
    np.random.seed(100) #For repeatable results randomly sampling
    # Read in the args
    args = sys.argv
    if len(args)<7:
        print "Need more args... Usage: "
        print "python SphereHarm.py InputFilename.off OutputFilename.pts NumSamplePoints"
        print "NumHarmonics NumShells ShellSampleAmount"
        print "Where ShellSampleAmount is a JSON array of all the sphere sampling values"
        exit(1)
    filename = args[1]
    output_name = args[2]
    approx_name = output_name.replace(".pts","_approx.pts")
    original_name = output_name.replace(".pts","_original.pts")
    n_rand_samples = int(args[3])
    num_harmonics = int(args[4])
    num_shells = int(args[5])
    shell_samples = json.loads(args[6])
    #Load in the mesh we care about
    m = PolyMesh()
    m.loadOffFileExternal(filename)
    (Ps, Ns) = samplePointCloud(m, n_rand_samples)
    # Get the spherical harmonic magnitudes
    [h_approx,h] = getSphericalHarmonicMagnitudes(Ps, Ns, num_harmonics, num_shells, shell_samples)
    new_approx_point_cloud = getSamplePointsFrom2DHistogram(n_rand_samples, h_approx)
    new_point_cloud = getSamplePointsFrom2DHistogram(n_rand_samples, h)
    # Now we'll export this point cloud. Normals will be completely wrong because we never calculated them
    exportPointCloud(new_point_cloud,Ns,output_name)
    exportPointCloud(new_approx_point_cloud,Ns,approx_name)
    exportPointCloud(Ps,Ns,original_name)
