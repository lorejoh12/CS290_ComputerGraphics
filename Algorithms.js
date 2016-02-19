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
    //TODO: Fill this in
    
    //Step 1: Make a new array of vec3s which holds "vertices" transformed 
    //to world coordinates (hint: vec3 has a function "transformMat4" which is useful)

    
    //Step 2: Compute the plane normal of the plane spanned by the transformed vertices
    
    //Step 3: Perform ray intersect plane
    
    
    //Step 4: Check to see if the intersection point is inside of the transformed polygon
    //You can assume that the polygon is convex.  If you use the area test, you can
    //allow for some wiggle room in the two areas you're comparing (e.g. absolute difference
    //not exceeding 1e-4)
    
    
    //Step 5: Return the intersection point if it exists or null if it's outside
    //of the polygon or if the ray is perpendicular to the plane normal (no intersection)
    
    return {t:1e9, P:vec3.fromValues(0, 0, 0)}; //These are dummy values, but you should return 
    //both an intersection point and a parameter t.  The parameter t will be used to sort
    //intersections in order of occurrence to figure out which one happened first
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
    scene.rayIntersectFaces = function(P0, V, node, mvMatrix, excludeFace) {
        var tmin = Infinity;//The parameter along the ray of the nearest intersection
        var PMin = null;//The point of intersection corresponding to the nearest interesection
        var faceMin = null;//The face object corresponding to the nearest intersection
        if (node === null) {
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
                var cres = scene.rayIntersectFaces(P0, V, node.children[i], nextmvMatrix, excludeFace);
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
    
    scene.isInPolygon = function(point,face) {
      	// TODO NC
      	// find area of triangles made from the point to the rest of the face
      
      	// find area of the whole face
      	scene.computeAreaOfTriangle(point,point,point);
    }
    
    scene.computeAreaOfTriangle = function(p1, p2, p3) {
      	// TODO NC
      	var a = 1;
      	var b = 2;
      	var c = 3;
      	var p = (a + b + c)/2.0;
      	var area = Math.sqrt(p*(p-a)*(p-b)*(p-c));
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
        scene.source.order = 0;//Store an order field to figure out how many 
        //bounces a particular image represents
        scene.source.rcoeff = 1.0;//Keep track of the reflection coefficient of the node that
        //gave rise to this source
        scene.source.parent = null;//Keep track of the image source's parent
        scene.source.genFace = null;//Keep track of the mesh face that generated this image
        //Remember not to reflect an image across the face that just generated it, 
        //or you'll get its parent image.  This information can also be used later
        //when tracing back paths
        scene.imsources = [scene.source];
      
      	for (var i = 1; i <= order; i++) {
          	var sourceIndex = 0;
          	while(sourceIndex<scene.imsources.length) {
              //var
              if (s.order == order - 1) {
                var stack = [];
		      	stack.push(scene);
              	// go through every source, and for each source, find every face in the scene graph
            	// check to make sure that the face we're on is not the genFace
              	// if it's not, then find the reflected point and add it to the 
              
              	// pull a node off of the stack
              	while(stack.length>0){
                  	var node = stack.pop(); //get our node
                  	console.log(node); // TODO delete... testing only
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
                              	var q = face.getCentroid(); // Arbitrary point on the face
                              	var p = source.pos;
                              	
                              	// we need to make sure we're all in the world reference
                              	
                              	//var newSource = p - 2*(p-q)dotn * n;
                              	var newSource = vec3.create();
                                var t = vec3.create();
                              	var t2 = vec3.create();
                              	vec3.subtract(t, p, q);
                                vec3.scale(t2, normal, vec3.dot(t, normal)*2);
                              	vec3.subtract(newSource, p, t2);
                              
                              
                              
                              	var point = source.pos; //temp!!! this point should be the point that the original point is reflected over
                              	if (scene.isInPolygon(point,face)) {
                                  
                                  	// TODO add the reflection to the list of points
                                }
								
                                
                                // "face" is a face we want to reflect over; it's made up of vertices & edges
                              	// 1. find normal vector of face
                              	// 2. find a vector from the point to the face
                              	// 3. do the dot product equation to reflect over the polygon
                              	// 4. check to make sure that the spot on the face reflected over is inside the polygon
                              	
                              
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
    scene.extractPaths = function() {
        scene.paths = [];
        
        //TODO: Finish this. Extract the rest of the paths by backtracing from
        //the image sources you calculated.  Return an array of arrays in
        //scene.paths.  Recursion is highly recommended
        //Each path should start at the receiver and end at the source
        //(or vice versa), so scene.receiver should be the first element 
        //and scene.source should be the last element of every array in 
        //scene.paths
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
