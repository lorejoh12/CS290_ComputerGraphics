//Purpose: A file that holds the code that students fill in

//Given a ray described by an initial point P0 and a direction V both in
//world coordinates, check to see 
//if it intersects the polygon described by "vertices," an array of vec3
//values describing the location of the polygon vertices in its child frame.
//mvMatrix is a matrix describing how to transform "vertices" into world coordinates
//which you will have to do to get the correct intersection in world coordinates.
//Be sure to compute the plane normal only after you have transformed the points,
//and be sure to only compute intersections which are inside of the polygon
//(you can assume that all polygons are convex and use the area method)
function rayIntersectPolygon(P0, V, vertices, mvMatrix) {
    // If we have fewer than 3 points, we can't do this. Return null
  	if(vertices.length<3){
      return null;
    }
  
    //Step 1: Make a new array of vec3s which holds "vertices" transformed 
    //to world coordinates (hint: vec3 has a function "transformMat4" which is useful)
  	var worldVertices = [];
	for(var i=0; i<vertices.length; i++){
      var tempVertex = vec3.create();
      vec3.transformMat4(tempVertex,vertices[i],mvMatrix);
      worldVertices.push(tempVertex);
    }
    //Step 2: Compute the plane normal of the plane spanned by the transformed vertices
  
  	// Create 2 arbitrary vectors between points on plane. 
  	var vec1 = vec3.create();
  	var vec2 = vec3.create();
  	vec3.subtract(vec1, worldVertices[0], worldVertices[1]);
  	vec3.subtract(vec2, worldVertices[1], worldVertices[2]);
    // Cross them to get our normal
  	var normal = vec3.create();
  	vec3.cross(normal,vec1,vec2);
    //Step 3: Perform ray intersect plane
    
  	//First make sure V isn't parallel to the normal.  If it is, return null
  	var vDotn = vec3.dot(V, normal);
  	if(vDotn == 0){
     	return null; 
    }
  
  	// Calculate the distance along the ray where it intersects with the point at the normal base
  	var qMinusP = vec3.create();
  	vec3.subtract(qMinusP,worldVertices[0], P0);
  	var t = vec3.dot(qMinusP, normal)/vDotn;
  	
  	// If t is negative, we're on the wrong side of the ray so we don't intersect
  	if(t<0){
      return null;
    }

  	//Get our point of intersection with the plane
	var intersectionPoint = vec3.create();
  	vec3.scale(intersectionPoint,V,t);
  	vec3.add(intersectionPoint, P0, intersectionPoint);
    
    //Step 4: Check to see if the intersection point is inside of the transformed polygon
    //You can assume that the polygon is convex.  If you use the area test, you can
    //allow for some wiggle room in the two areas you're comparing (e.g. absolute difference
    //not exceeding 1e-4)
 
  	if(!isInPolygonFromVertices(intersectionPoint, worldVertices)){
      return null;
    }
  	
    //Step 5: Return the intersection point if it exists or null if it's outside
    //of the polygon or if the ray is perpendicular to the plane normal (no intersection)
    
    return {t:t, P:intersectionPoint}; //These are dummy values, but you should return 
    //both an intersection point and a parameter t.  The parameter t will be used to sort
    //intersections in order of occurrence to figure out which one happened first
}

function isInPolygonFromVertices(point, vertices) {
	var numVertices = vertices.length;

	// find area of triangles made from the point to the rest of the face
	var pointTrianglesArea = computeAreaOfTriangle(point, vertices[numVertices-1], vertices[0]);
	for (var i = 0; i < numVertices - 1; i++) {
		pointTrianglesArea += computeAreaOfTriangle(point, vertices[i], vertices[i+1]);
	}
  
	// find area of the whole face
	var faceArea = 0;
	var chosenPoint = vertices[numVertices-1];
	for (var i = 0; i < numVertices - 2; i++) {
		 faceArea += computeAreaOfTriangle(chosenPoint, vertices[i], vertices[i+1]);
	}
	
	// compare the two areas
	var areaThreshold = 0.00001;
	return Math.abs(pointTrianglesArea - faceArea) <= areaThreshold;
}

// Returns the area of the triangle formed by points p1, p2, p3
// e.g. p1 = [1, 2.6, 3.3]
function computeAreaOfTriangle(p1, p2, p3) {
	var a = getDistBtwnPoints(p1,p2);
	var b = getDistBtwnPoints(p2,p3);
	var c = getDistBtwnPoints(p1,p3);
  
	var p = (a + b + c)/2.0;
	var area = Math.sqrt(p*(p-a)*(p-b)*(p-c));
	
	return area;
}

// Returns the distance between two 3D points
// e.g. p1 = [1, 2.6, 3.3]
function getDistBtwnPoints(p1, p2) {
  var a2 = Math.pow(p1[0]-p2[0],2);
  var b2 = Math.pow(p1[1]-p2[1],2);
  var c2 = Math.pow(p1[2]-p2[2],2);
  
  return Math.sqrt(a2 + b2 + c2);
}


function addImageSourcesFunctions(scene) {
    //Setup all of the functions that students fill in that operate directly
    //on the scene
    
    //Purpose: A recursive function provided which helps to compute intersections of rays
    //with all faces in the scene, taking into consideration the scene graph structure
    //Inputs: P0 (vec3): Ray starting point, V (vec3): ray direction
    //node (object): node in scene tree to process, 
    //mvMatrix (mat4): Matrix to put geometry in this node into world coordinates
    //excludeFace: Pointer to face object to be excluded (don't intersect with
    //the face that this point lies on)
    //Returns: null if no intersection,
    //{tmin:minimum t along ray, PMin(vec3): corresponding point, faceMin:Pointer to mesh face hit first}
    
    //NOTE: Calling this function with node = scene and an identity matrix for mvMatrix
    //will start the recursion at the top of the scene tree in world coordinates
    scene.rayIntersectFaces = function(P0, V, node, mvMatrix, excludeFace, pruneByBoundingBox) {
      	if (pruneByBoundingBox) scene.computeBoundingBoxes();
        var tmin = Infinity;//The parameter along the ray of the nearest intersection
        var PMin = null;//The point of intersection corresponding to the nearest interesection
        var faceMin = null;//The face object corresponding to the nearest intersection
        if (node === null) {
            return null;
        }
      
      	if (pruneByBoundingBox && !scene.rayIntersectsBoundingBox(P0, V, node, mat4.create())) {
            return null;
        }
      
        if ('mesh' in node) { //Make sure it's not just a dummy transformation node
            var mesh = node.mesh;
            for (var f = 0; f < mesh.faces.length; f++) {
                if (mesh.faces[f] == excludeFace) {
                    continue;//Don't re-intersect with the face this point lies on
                }
                //Intersect the ray with this polygon
                var res = rayIntersectPolygon(P0, V, mesh.faces[f].getVerticesPos(), mvMatrix);
                if (!(res === null) && (res.t < tmin)) {
                    tmin = res.t;
                    PMin = res.P;
                    faceMin = mesh.faces[f];
                }
            }
        }
        
        if ('children' in node) {
            //Recursively check the meshes of the children to make sure the ray
            //doesn't intersect any of them first
            for (var i = 0; i < node.children.length; i++) {
                var nextmvMatrix = mat4.create();
                //Multiply on the right by the next transformation of the child
                //node
                mat4.mul(nextmvMatrix, mvMatrix, node.children[i].transform);
                //Recursively intersect with the child node
                var cres = scene.rayIntersectFaces(P0, V, node.children[i], nextmvMatrix, excludeFace, pruneByBoundingBox);
                if (!(cres === null) && (cres.tmin < tmin)) {
                    tmin = cres.tmin;
                    PMin = cres.PMin;
                    faceMin = cres.faceMin;
                }
            }
        }
        if (PMin === null) {
            return null;
        }
        return {tmin:tmin, PMin:PMin, faceMin:faceMin};
    }
    
    // Returns a boolean of whether the point lies within the convex polygon face
    // e.g. point = [1, 2.6, 3.3]. 
    scene.isInPolygonFromFace = function(point, face, worldTransform) {
      	var vertices = scene.getWorldCoordinateVertices(face, worldTransform);
      	return isInPolygonFromVertices(point, vertices);
    }
    
    // Returns the world coordinates of the vertices of a face given a world transform
    scene.getWorldCoordinateVertices = function(face, worldTransform) {
      	var rawVertices = face.getVerticesPos();
      	var vertices = []
        for (var i = 0; i < rawVertices.length; i++) {
          	var temp = vec3.create();
			vec3.transformMat4(temp, rawVertices[i], worldTransform);
          	vertices.push(temp);
        }
      	return vertices;
    }
    
    scene.rayIntersectsBoundingBox = function(P0, V, node, mvMatrix) {
        // check if the ray intersects the bounding box
        var box = node.boundingBox;
      	var v1 = [box.minX, box.minY, box.minZ];
      	var v2 = [box.minX, box.minY, box.maxZ];
      	var v3 = [box.minX, box.maxY, box.minZ];
      	var v4 = [box.minX, box.maxY, box.maxZ];
      	var v5 = [box.maxX, box.minY, box.minZ];
      	var v6 = [box.maxX, box.minY, box.maxZ];
      	var v7 = [box.maxX, box.maxY, box.minZ];
      	var v8 = [box.maxX, box.maxY, box.maxZ];
      	
      	// face 1: minX
      	var vertices1 = [v1,v2,v4,v3];
      	// face 2: maxX
      	var vertices2 = [v5,v6,v8,v7];
      	// face 3: minY
      	var vertices3 = [v1,v5,v6,v2];
      	// face 4: maxY
      	var vertices4 = [v3,v7,v8,v4];
      	// face 5: minZ
      	var vertices5 = [v1,v3,v7,v5];
      	// face 6: maxZ
      	var vertices6 = [v2,v6,v8,v4];
      	
      	if (!rayIntersectPolygon(P0, V, vertices1, mvMatrix) &&
	      	!rayIntersectPolygon(P0, V, vertices2, mvMatrix) && 
   	   		!rayIntersectPolygon(P0, V, vertices3, mvMatrix) && 
      		!rayIntersectPolygon(P0, V, vertices4, mvMatrix) && 
      		!rayIntersectPolygon(P0, V, vertices5, mvMatrix) && 
      		!rayIntersectPolygon(P0, V, vertices6, mvMatrix)
           ) 
        {
          return false;
        }
      
      	return true;
      	//sorry Duvall
    }
    
    scene.computeBoundingBoxes = function() {
       if (!scene.boundingBox) {
            var identity = mat4.create();
            scene.findBoundingBox(scene, identity);
        }
  	}

   scene.findBoundingBox = function(node) {
        if (node === null) {
            return null;
        }
     
        var maxX = Number.MIN_VALUE;
        var minX = Number.MAX_VALUE;
        var maxY = Number.MIN_VALUE;
        var minY = Number.MAX_VALUE;
        var maxZ = Number.MIN_VALUE;
        var minZ = Number.MAX_VALUE;
     
     	// compute max and min for the node itself
        if ('mesh' in node) { //Make sure it's not just a dummy transformation node          
            var mesh = node.mesh;
            for (var v = 0; v < mesh.vertices.length; v++) {
                var temp = mesh.vertices[v].pos;
              	var vertexPos = vec3.create();
              	vec3.transformMat4(vertexPos, temp, node.worldTransform);
              
              	if (vertexPos[0] < minX) minX = vertexPos[0];
                if (vertexPos[0] > maxX) maxX = vertexPos[0];
              	if (vertexPos[1] < minY) minY = vertexPos[1];
              	if (vertexPos[1] > maxY) maxY = vertexPos[1];
              	if (vertexPos[2] < minZ) minZ = vertexPos[2];
              	if (vertexPos[2] > maxZ) maxZ = vertexPos[2];
            }
        }
        
     	// compute max and min for all of the node's children
        if ('children' in node) {
            //Recursively check the meshes of the children to compute the min and max
            for (var i = 0; i < node.children.length; i++) {
                //Recursively intersect with the child node
                var box = scene.findBoundingBox(node.children[i]);
              
              	if (box.minX < minX) minX = box.minX;
              	if (box.maxX > maxX) maxX = box.maxX;
              	if (box.minY < minY) minY = box.minY;
              	if (box.maxY > maxY) maxY = box.maxY;
              	if (box.minZ < minZ) minZ = box.minZ;
              	if (box.maxZ > maxZ) maxZ = box.maxZ;
            }
        }
     	
     	node.boundingBox = {minX:minX, maxX:maxX, minY:minY, maxY:maxY, minZ:minZ, maxZ:maxZ};
     
        return node.boundingBox;
    }
    
    scene.computeWorldTransforms = function(node) {
    	if (node == null) return;
    	if('children' in node){
			for(var k=0;k<node.children.length;k++){
				node.children[k].worldTransform = mat4.create();
				//Multiply on the right by the next transformation of the child node
                mat4.mul(node.children[k].worldTransform, node.worldTransform, node.children[k].transform);
                scene.computeWorldTransforms(node.children[k]);
            }
        }  
    }
    
    //Purpose: Fill in the array scene.imsources[] with a bunch of source
    //objects.  It's up to you what you put in the source objects, but at
    //the very least each object needs a field "pos" describing its position
    //in world coordinates so that the renderer knows where to draw it
    //You will certainly also need to save along pointers from an image source
    //to its parent so that when you trace paths back you know where to aim
    //Recursion is highly recommended here, since you'll be making images of 
    //images of images (etc...) reflecting across polygon faces.
    
    //Inputs: order (int) : The maximum number of bounces to take
    scene.computeImageSources = function(order) {
      	scene.maxOrder = order;
        scene.source.order = 0;//Store an order field to figure out how many 
        //bounces a particular image represents
        scene.source.rcoeff = 1.0;//Keep track of the reflection coefficient of the node that
        //gave rise to this source
        scene.source.parent = null;//Keep track of the image source's parent
        scene.source.genFace = null;//Keep track of the mesh face that generated this image
        //Remember not to reflect an image across the face that just generated it, 
        //or you'll get its parent image.  This information can also be used later
        //when tracing back paths
      	scene.worldTransform = mat4.create(); // create Identity matrix for initial scene transform
      	
      	scene.computeWorldTransforms(scene);
      	
        scene.imsources = [scene.source];
      	// For all orders               	
      	
      	for (var currentOrder = 1; currentOrder <= order; currentOrder++) {
          	// Go through all the sources
          	var sourceIndex = 0;
          	while(sourceIndex<scene.imsources.length) {
              var source = scene.imsources[sourceIndex];
              sourceIndex++;
              if (source.order == currentOrder - 1) {
                var stack = [];
		      	stack.push(scene);
              	// go through every source, and for each source, find every face in the scene graph
            	// check to make sure that the face we're on is not the genFace
              	// if it's not, then find the reflected point and add it to the 
              
              	// pull a node off of the stack
              	while(stack.length>0){
                  	var node = stack.pop(); //get our node
                  	//Throw its children into a pile
                  	if('children' in node){
                    	for(var k=0;k<node.children.length;k++){
                          stack.push(node.children[k]);
                        }
                    }                  	
                  	//Now do some serious polygonal processing shit
                  	if ('mesh' in node) {
                      	var mesh = node.mesh;
                        for (var i = 0; i < mesh.faces.length; i++) {
                          	var face = mesh.faces[i];
                            if (face == source.genFace) {
                                continue;// don't re-reflect across the old genFace
                            }
                          	else {
                              	var normal = face.getNormal(); // Normal of the face
                              	var tempCentroid = face.getCentroid(); // Arbitrary point on the face
                              	var p = source.pos;                              	
                              
                              	// we need to make sure we're all in the world reference
                              
                              	var q = vec3.create();
                              
                              	// transform by the worldTransform to get the world location of the centroid point
                              	vec3.transformMat4(q, tempCentroid, node.worldTransform);

                              	//newSource = p - 2*(p-q)dotn * n;
                              	var newSource = vec3.create();
                                var pMinusQ = vec3.create();
                              	var t1 = vec3.create();
                              	vec3.subtract(pMinusQ, p, q);
                                vec3.scale(t1, normal, vec3.dot(pMinusQ, normal)*2);
                              	vec3.subtract(newSource, p, t1);
                              
                              	// find the point of reflection 
                              	// (the point on the plane of the face that is the midpoint of the source and the new source)
                              	var t2 = vec3.create();
                              	var pointOfReflection = vec3.create();
                                vec3.scale(t2, normal, vec3.dot(pMinusQ, normal));
                              	vec3.subtract(pointOfReflection, p, t2);

                              	// create the mirror image source and push it to the list of sources
                              	var mirrorImage = {pos:newSource, order:currentOrder, parent:source, genFace:face, rcoeff:node.rcoeff};
                                scene.imsources.push(mirrorImage);
                              	
                            }
                        }
                    }
                  }
                }
          	}
        }        
    }    
    
    //Purpose: Based on the extracted image sources, trace back paths from the
    //receiver to the source, checking to make sure there are no occlusions
    //along the way.  Remember, you're always starting by tracing a path from
    //the receiver to the image, and then from the intersection point with
    //that image's corresponding face to the image's parent, and so on
    //all the way until you get back to the original source.
    
    //Fill in the array scene.paths, where each element of the array is itself
    //an array of objects describing vertices along the path, starting
    //with the receiver and ending with the source.  Each object in each path
    //array should contain a field "pos" which describes the position, as well
    //as an element "rcoeff" which stores the reflection coefficient at that
    //part of the path, which will be used to compute decays in "computeInpulseResponse()"
    //Don't forget the direct path from source to receiver!
    scene.extractPaths = function(pruneByBoundingBox) {
        scene.paths = [];
      	var path = [scene.receiver];
      	// Check all possible orders of paths
      	for (var i=0; i<=scene.maxOrder;i++){
      		scene.checkForSourcePath(scene.receiver.pos,null,i,path, pruneByBoundingBox);
        }

        //TODO: Finish this. Extract the rest of the paths by backtracing from
        //the image sources you calculated.  Return an array of arrays in
        //scene.paths.  Recursion is highly recommended
        //Each path should start at the receiver and end at the source
        //(or vice versa), so scene.receiver should be the first element 
        //and scene.source should be the last element of every array in 
        //scene.paths
    }
    
    // Recursive function for path checking.
    // Point is the point to start checking from, face is the current source's face to ignore when doing the intersection, order is the order of source we are examining
    // and path is the current path we are traversing
    // Automatically adds the path to the scene.paths variable.
    scene.checkForSourcePath = function(point, face, order, path, pruneByBoundingBox){
      	// If we've found our source, this is also great. Let's return the path
        if(point == scene.source.pos){
          scene.paths.push(path);
        }
      	//If our order is less than 0, we're done. This point is worthless to us
        if(order<0){
          return;
        }
		
      	// Go through all the image sources
      	for (var i = 0; i<scene.imsources.length;i++){
          	var source = scene.imsources[i];
          	// If the image source is of the order that we're looking for
          	if (source.order == order) {
          		// Get the unit direction from point to the source
                var dirVector = vec3.create();
              	vec3.subtract(dirVector, source.pos, point);
              	vec3.normalize(dirVector, dirVector);
              
                //If we have at least 2 points in the path (receiver+something else), check the dot product versus the normal to make sure same
              	// Otherwise, this path is not the path you're looking for, continue
              	var THRESHOLD = 0.0004
				if(face != null && path.length >=2) {
					// Get the vector from the old source to the current point
					var oldCurrentDir = vec3.create();
					vec3.subtract(oldCurrentDir, path[path.length-2].pos, point);
					vec3.normalize(oldCurrentDir, oldCurrentDir);
					
					var dot1 = vec3.dot(dirVector, face.getNormal());
					var dot2 = vec3.dot(oldCurrentDir, face.getNormal());
					if (path.length>=2 && Math.abs(dot1-dot2) > THRESHOLD){
					  continue;
					}
				}
          
          		// Pass into scene.rayIntersectFaces to get the points of intersection. Be sure to ignore the face the point is currently on
				
				//Returns: null if no intersection,
				//{tmin:minimum t along ray, PMin(vec3): corresponding point, faceMin:Pointer to mesh face hit first}
          		var intersect = scene.rayIntersectFaces(point, dirVector, scene, mat4.create(), face, pruneByBoundingBox); 
              
          		// If the closest intersection is the source's parent face (or if the node is the source), we're good. Duplicate the path array, add the point of intersection to the path,
    			// subtract 1 from order, and call the recursive method on the point of intersection with the source's parent face.
				var hitSource = (intersect == null || getDistBtwnPoints(point, scene.source.pos) < getDistBtwnPoints(point, intersect.PMin)) && (source == scene.source);
				if (hitSource) {
					var newPath = path.slice(0);
					newPath.push({pos:scene.source.pos});
					scene.checkForSourcePath(scene.source.pos, null, order-1, newPath, pruneByBoundingBox);
				}
              	else if (intersect != null && intersect.faceMin == source.genFace) {
					var newPath = path.slice(0);
					newPath.push({pos:intersect.PMin});
					scene.checkForSourcePath(intersect.PMin, intersect.faceMin, order-1, newPath, pruneByBoundingBox);
				}
            }
        }
    }
    
    
    //Inputs: Fs: Sampling rate (samples per second)
    scene.computeImpulseResponse = function(Fs) {
        var SVel = 340;//Sound travels at 340 meters/second
        //TODO: Finish this.  Be sure to scale each bounce by 1/(1+r^p), 
        //where r is the length of the line segment of that bounce in meters
        //and p is some integer less than 1 (make it smaller if you want the 
        //paths to attenuate less and to be more echo-y as they propagate)
        //Also be sure to scale by the reflection coefficient of each material
        //bounce (you should have stored this in extractPaths() if you followed
        //those directions).  Use some form of interpolation to spread an impulse
        //which doesn't fall directly in a bin to nearby bins
        //Save the result into the array scene.impulseResp[]
    }
}
