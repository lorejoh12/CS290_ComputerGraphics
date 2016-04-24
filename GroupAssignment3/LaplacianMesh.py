import sys
sys.path.append("S3DGLPy")
from PolyMesh import *
from Primitives3D import *
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import lsqr, cg, eigsh
import matplotlib.pyplot as plt
import scipy.io as sio


##############################################################
##                  Laplacian Mesh Editing                  ##
##############################################################

#Purpose: To return a sparse matrix representing a Laplacian matrix with
#the graph Laplacian (D - A) in the upper square part and anchors as the
#lower rows
#Inputs: mesh (polygon mesh object), anchorsIdx (indices of the anchor points)
#Returns: L (An (N+K) x N sparse matrix, where N is the number of vertices
#and K is the number of anchors)
WEIGHT = 1

def getLaplacianMatrixUmbrella(mesh, anchorsIdx):
    N = len(mesh.vertices)
    K = len(anchorsIdx)
    I = []
    J = []
    V = []
    for i in range(0, N):
        vertex = mesh.vertices[i]
        # append the D matrix
        I.append(vertex.ID)
        J.append(vertex.ID)
        V.append(len(vertex.getVertexNeighbors()))
        for j in range(len(vertex.getVertexNeighbors())):
            neighbor = vertex.getVertexNeighbors()[j]
            # append the -A matrix
            I.append(neighbor.ID)
            J.append(vertex.ID)
            V.append(-1)
    
    # append the anchors
    row = N
    for index in anchorsIdx:
        I.append(row)
        J.append(index)
        V.append(WEIGHT)
        row+=1      

    L = sparse.coo_matrix((V, (I, J)), shape=(N+K, N)).tocsr()
    return L

#Purpose: To return a sparse matrix representing a laplacian matrix with
#cotangent weights in the upper square part and anchors as the lower rows
#Inputs: mesh (polygon mesh object), anchorsIdx (indices of the anchor points)
#Returns: L (An (N+K) x N sparse matrix, where N is the number of vertices
#and K is the number of anchors)
def getLaplacianMatrixCotangent(mesh, anchorsIdx):
    I = []
    J = []
    V = []
    N = len(mesh.vertices)
    K = len(anchorsIdx)
    
    for vi in mesh.vertices:
        iiEntry = 0 # sum (subtract) up as we go through the j's
        for neighbor in vi.getVertexNeighbors():
            edge = getEdgeInCommon(vi, neighbor)
            cotangentA = 0 # assume 0 if the face is null (TODO: is this OK?)
            cotangentB = 0 # assume 0 if the face is null (TODO: is this OK?)
            if edge.f1:
                vertices = edge.f1.getVertices()
                for vk in vertices:
                    if(vk != vi and vk != neighbor):
                        # vk is the third vertex
                        cotangentA = computeCotAngle(mesh.VPos[vi.ID],mesh.VPos[neighbor.ID],mesh.VPos[vk.ID])
            if edge.f2:
                vertices = edge.f2.getVertices()
                for vk in vertices:
                    if(vk != vi and vk != neighbor):
                        # vk is the third vertex
                        cotangentB = computeCotAngle(mesh.VPos[vi.ID],mesh.VPos[neighbor.ID],mesh.VPos[vk.ID])

            ijEntry = -0.5*(cotangentB + cotangentA)
            iiEntry -= ijEntry
            # set the ij entry
            I.append(vi.ID)
            J.append(neighbor.ID)
            V.append(ijEntry)

        # set the ii entry
        I.append(vi.ID)
        J.append(vi.ID)
        V.append(iiEntry)

    # append the anchors
    row = len(mesh.vertices)
    for index in anchorsIdx:
        I.append(row)
        J.append(index)
        V.append(WEIGHT)
        row+=1

    L = sparse.coo_matrix((V, (I, J)), shape=(N+K, N)).tocsr()
    return L

# computes the cotangent of the angle between CA and CB
# assume A = 1 x 3 matrix of [Ax, Bx, Cx]
def computeCotAngle(A, B, C):
    CA = np.subtract(A, C)
    CB = np.subtract(B, C)
    dotProd = np.dot(CA, CB)
    crossProdMag = np.linalg.norm(np.cross(CA, CB))
    return dotProd/crossProdMag

#Purpose: Given a mesh, to perform Laplacian mesh editing by solving the system
#of delta coordinates and anchors in the least squared sense
#Inputs: mesh (polygon mesh object), anchors (a K x 3 numpy array of anchor
#coordinates), anchorsIdx (a parallel array of the indices of the anchors)
#Returns: Nothing (should update mesh.VPos)
def solveLaplacianMesh(mesh, anchors, anchorsIdx):
    # Solve for the Laplacian and delta matrix.
    N = len(mesh.vertices)
    K = len(anchorsIdx)
    laplacian_matrix = getLaplacianMatrixCotangent(mesh, anchorsIdx)
    delta = np.array(laplacian_matrix.dot(mesh.VPos))
    # Now update the anchors in the delta matrix
    for i in range(K):
        delta[N+i,:]=anchors[i,:].T*WEIGHT
    # Update vpos with the new vertex locations
    x = lsqr(laplacian_matrix,delta[:,0])
    y = lsqr(laplacian_matrix,delta[:,1])
    z = lsqr(laplacian_matrix,delta[:,2])
    mesh.VPos[:,0] = x[0]
    mesh.VPos[:,1] = y[0]
    mesh.VPos[:,2] = z[0]

#Purpose: Given a few RGB colors on a mesh, smoothly interpolate those colors
#by using their values as anchors and 
#Inputs: mesh (polygon mesh object), anchors (a K x 3 numpy array of anchor
#coordinates), anchorsIdx (a parallel array of the indices of the anchors)
#Returns: Nothing (should update mesh.VPos)
def smoothColors(mesh, colors, colorsIdx):
    N = len(mesh.vertices)
    K = len(colorsIdx)
    laplacian_matrix = getLaplacianMatrixCotangent(mesh, colorsIdx)
    delta = np.zeros((N+K, 3))
    # Now update the anchors in the delta matrix
    for i in range(K):
        delta[N+i,:]=colors[i,:].T*WEIGHT
    # Update vpos with the new vertex locations
    x = lsqr(laplacian_matrix,delta[:,0])
    y = lsqr(laplacian_matrix,delta[:,1])
    z = lsqr(laplacian_matrix,delta[:,2])
    _colors = np.zeros((N, 3))
    _colors[:,0] = x[0]
    _colors[:,1] = y[0]
    _colors[:,2] = z[0]
    
    return _colors

#Purpose: Given a mesh, to smooth it by subtracting off the delta coordinates
#from each vertex, normalized by the degree of that vertex
#Inputs: mesh (polygon mesh object)
#Returns: Nothing (should update mesh.VPos)
def doLaplacianSmooth(mesh):
    # Solve for the Laplacian and delta matrix.
    N = len(mesh.vertices)
    anchorsIdx = []
    L = getLaplacianMatrixCotangent(mesh, anchorsIdx)
    (I,J,V) = sparse.find(L)
    V2 = V
    L_N = L
    for i in range(0, len(I)):
        V2[i] = V[i] / L[I[i], I[i]]
    L_N = sparse.coo_matrix((V2, (I, J)), shape=(N, N)).tocsr()
    mesh.VPos = mesh.VPos - (np.array(L_N.dot(mesh.VPos)))
    

#Purpose: Given a mesh, to sharpen it by adding back the delta coordinates
#from each vertex, normalized by the degree of that vertex
#Inputs: mesh (polygon mesh object)
#Returns: Nothing (should update mesh.VPos)
def doLaplacianSharpen(mesh):
    # Solve for the Laplacian and delta matrix.
    N = len(mesh.vertices)
    anchorsIdx = []
    L = getLaplacianMatrixCotangent(mesh, anchorsIdx)
    (I,J,V) = sparse.find(L)
    V2 = V
    L_N = L
    for i in range(0, len(I)):
        V2[i] = V[i] / L[I[i], I[i]]
    L_N = sparse.coo_matrix((V2, (I, J)), shape=(N, N)).tocsr()
    mesh.VPos = mesh.VPos + (np.array(L_N.dot(mesh.VPos)))

#Purpose: Given a mesh and a set of anchors, to simulate a minimal surface
#by replacing the rows of the laplacian matrix with the anchors, setting
#those "delta coordinates" to the anchor values, and setting the rest of the
#delta coordinates to zero
#Inputs: mesh (polygon mesh object), anchors (a K x 3 numpy array of anchor
#coordinates), anchorsIdx (a parallel array of the indices of the anchors)
#Returns: Nothing (should update mesh.VPos)
def makeMinimalSurface(mesh, anchors, anchorsIdx):
    N = len(mesh.vertices)
    # we don't want any anchors at the bottom of the matrix
    temp = []
    laplacian_matrix = getLaplacianMatrixCotangent(mesh, temp)
    # update the anchors in the L matrix
    (I,J,V) = sparse.find(laplacian_matrix)
    for i in range(0, len(I)):
        if I[i] in anchorsIdx:
            if I[i] == J[i]:
                V[i] = 1
            else:
                V[i] = 0
                
    laplacian_matrix = sparse.coo_matrix((V, (I, J)), shape=(N, N)).tocsr()
        
    delta = np.zeros((N, 3))
    # Now update the anchors in the delta matrix
    for i in range(0, len(anchorsIdx)):
        delta[anchorsIdx[i],:]=anchors[i,:].T
    # Update vpos with the new vertex locations
    x = lsqr(laplacian_matrix,delta[:,0])
    y = lsqr(laplacian_matrix,delta[:,1])
    z = lsqr(laplacian_matrix,delta[:,2])
    mesh.VPos[:,0] = x[0]
    mesh.VPos[:,1] = y[0]
    mesh.VPos[:,2] = z[0]

##############################################################
##        Spectral Representations / Heat Flow              ##
##############################################################

#Purpose: Given a mesh, to compute first K eigenvectors of its Laplacian
#and the corresponding eigenvalues
#Inputs: mesh (polygon mesh object), K (number of eigenvalues/eigenvectors)
#Returns: (eigvalues, eigvectors): a tuple of the eigenvalues and eigenvectors
def getLaplacianSpectrum(mesh, K):
    temp = []
    L = getLaplacianMatrixUmbrella(mesh,temp)
    (eigvalues, eigvectors) = eigsh(L.asfptype(), K, which='LM', sigma = 0)
    return (eigvalues, eigvectors)

#Purpose: Given a mesh, to use the first K eigenvectors of its Laplacian
#to perform a lowpass filtering
#Inputs: mesh (polygon mesh object), K (number of eigenvalues/eigenvectors)
#Returns: Nothing (should update mesh.VPos)
def doLowpassFiltering(mesh, K):
    (eigvalues, eigvectors) = getLaplacianSpectrum(mesh, K)
    mesh.VPos = (eigvectors.dot(eigvectors.T)).dot(mesh.VPos)
    
#Purpose: Given a mesh, to simulate heat flow by projecting initial conditions
#onto the eigenvectors of the Laplacian matrix, and then to sum up the heat
#flow of each eigenvector after it's decayed after an amount of time t
#Inputs: mesh (polygon mesh object), eigvalues (K eigenvalues), 
#eigvectors (an NxK matrix of eigenvectors computed by your laplacian spectrum
#code), t (the time to simulate), initialVertices (indices of the verticies
#that have an initial amount of heat), heatValue (the value to put at each of
#the initial vertices at the beginning of time
#Returns: heat (a length N array of heat values on the mesh)
def getHeat(mesh, eigvalues, eigvectors, t, initialVertices, heatValue = 100.0):
    N = mesh.VPos.shape[0]
    heat = np.zeros(N) #Dummy value
    return heat #TODO: Finish this

#Purpose: Given a mesh, to approximate its curvature at some measurement scale
#by recording the amount of heat that stays at each vertex after a unit impulse
#of heat is applied.  This is called the "Heat Kernel Signature" (HKS)
#Inputs: mesh (polygon mesh object), K (number of eigenvalues/eigenvectors to use)
#t (the time scale at which to compute the HKS)
#Returns: hks (a length N array of the HKS values)
def getHKS(mesh, K, t):
    N = mesh.VPos.shape[0]
    hks = np.zeros(N) #Dummy value
    return hks #TODO: Finish this

##############################################################
##                Parameterization/Texturing               ##
##############################################################

#Purpose: Given 4 vertex indices on a quadrilateral, to anchor them to the 
#square and flatten the rest of the mesh inside of that square
#Inputs: mesh (polygon mesh object), quadIdxs (a length 4 array of indices
#into the mesh of the four points that are to be anchored, in CCW order)
#Returns: nothing (update mesh.VPos)
def doFlattening(mesh, quadIdxs):
    # we assume the points are given as in the picture, where the shared edge is between vertices 1 and 3
    if(not len(quadIdxs) == 4):
        print '  Error, points do not conform to standard'
    N = len(mesh.vertices)
    # we don't want any anchors at the bottom of the matrix
    temp = []
    laplacian_matrix = getLaplacianMatrixUmbrella(mesh, temp)
    # update the anchors in the L matrix
    (I,J,V) = sparse.find(laplacian_matrix)
    
    if mesh.vertices[quadIdxs[2]] in mesh.vertices[quadIdxs[0]].getVertexNeighbors():
        idx = [0, 2]
        print '  0 and 2'
    elif mesh.vertices[quadIdxs[1]] in mesh.vertices[quadIdxs[3]].getVertexNeighbors():
        idx = [1, 3]
        print '  1 and 3'
    else:
        print '  Error, points do not conform to standard'

    for i in range(0, len(I)):
        #overwrite anchor rows in matrix
        if I[i] in quadIdxs:
            if I[i] == J[i]:
                V[i] = 1
            else:
                V[i] = 0
        # if I[i] == quadIdxs[idx[0]]:
            # if J[i] == quadIdxs[idx[1]]:
                # V[i] = 0
            # if J[i] == I[i]:
                # V[i] = V[i] - 1
        # if I[i] == quadIdxs[idx[1]]:
            # if J[i] == quadIdxs[idx[0]]:
                # V[i] = 0
            # if J[i] == I[i]:
                # V[i] = V[i] - 1
                
    laplacian_matrix = sparse.coo_matrix((V, (I, J)), shape=(N, N)).tocsr()
        
    delta = np.zeros((N, 3))
    # Now update the anchors in the delta matrix
    delta[quadIdxs[0],:]=np.array([0, 0, 0])
    delta[quadIdxs[1],:]=np.array([0, 1, 0])
    delta[quadIdxs[2],:]=np.array([1, 1, 0])
    delta[quadIdxs[3],:]=np.array([1, 0, 0])

    # Update vpos with the new vertex locations
    x = lsqr(laplacian_matrix,delta[:,0])
    y = lsqr(laplacian_matrix,delta[:,1])
    z = lsqr(laplacian_matrix,delta[:,2])
    mesh.VPos[:,0] = x[0]
    mesh.VPos[:,1] = y[0]
    mesh.VPos[:,2] = z[0]

#Purpose: Given 4 vertex indices on a quadrilateral, to anchor them to the 
#square and flatten the rest of the mesh inside of that square.  Then, to 
#return these to be used as texture coordinates
#Inputs: mesh (polygon mesh object), quadIdxs (a length 4 array of indices
#into the mesh of the four points that are to be anchored, in CCW order)
#Returns: U (an N x 2 matrix of texture coordinates)
def getTexCoords(mesh, quadIdxs):
    # we assume the points are given as in the picture, where the shared edge is between vertices 1 and 3
    if(not len(quadIdxs) == 4):
        print '  Error, points do not conform to standard'
    N = len(mesh.vertices)
    # we don't want any anchors at the bottom of the matrix
    temp = []
    laplacian_matrix = getLaplacianMatrixUmbrella(mesh, temp)
    # update the anchors in the L matrix
    (I,J,V) = sparse.find(laplacian_matrix)
    
    if mesh.vertices[quadIdxs[2]] in mesh.vertices[quadIdxs[0]].getVertexNeighbors():
        idx = [0, 2]
        print '  0 and 2'
    elif mesh.vertices[quadIdxs[1]] in mesh.vertices[quadIdxs[3]].getVertexNeighbors():
        idx = [1, 3]
        print '  1 and 3'
    else:
        print '  Error, points do not conform to standard'

    for i in range(0, len(I)):
        #overwrite anchor rows in matrix
        if I[i] in quadIdxs:
            if I[i] == J[i]:
                V[i] = 1
            else:
                V[i] = 0
        # if I[i] == quadIdxs[idx[0]]:
            # if J[i] == quadIdxs[idx[1]]:
                # V[i] = 0
            # if J[i] == I[i]:
                # V[i] = V[i] - 1
        # if I[i] == quadIdxs[idx[1]]:
            # if J[i] == quadIdxs[idx[0]]:
                # V[i] = 0
            # if J[i] == I[i]:
                # V[i] = V[i] - 1
                
    laplacian_matrix = sparse.coo_matrix((V, (I, J)), shape=(N, N)).tocsr()
        
    delta = np.zeros((N, 3))
    # Now update the anchors in the delta matrix
    delta[quadIdxs[0],:]=np.array([0, 0, 0])
    delta[quadIdxs[1],:]=np.array([0, 1, 0])
    delta[quadIdxs[2],:]=np.array([1, 1, 0])
    delta[quadIdxs[3],:]=np.array([1, 0, 0])

    # Update vpos with the new vertex locations
    x = lsqr(laplacian_matrix,delta[:,0])
    y = lsqr(laplacian_matrix,delta[:,1])
    z = lsqr(laplacian_matrix,delta[:,2])
    U = np.zeros((N, 2))
    U[:,0] = x[0]
    U[:,1] = y[0]
    return U

if __name__ == '__main__':
    print "TODO"
