
load "utilities.m";
 
QQ := Rationals();
interval := Polytope([[-1],[1]]);
polygon4 := Polytope([[-1,0],[0,1],[1,0],[0,-1]]);
polygon5 := Polytope([[-1,-1],[-1,0],[0,1],[1,0],[0,-1]]);
polygon6 := Polytope([[1,0],[1,1],[0,1],[-1,0],[-1,-1],[0,-1]]);
polygon7 := Polytope([[-1,-1],[1,-1],[1,0],[0,1],[-1,1]]);
polygon8 := Polytope([[-1,-1],[1,-1],[-1,1],[1,1]]);
polygon9 := Polytope([[-1,-1],[2,-1],[-1,2]]);

P := polygon6*polygon6;

// Cache all the dual face information.
DualFaces := AssociativeArray();
for i in [1,2] do
    for face in Faces(P,i) do
        DualFaces[face] := DualFace(P,face);
    end for;
end for;

// The ordering of the vertices of edges of Polar(P) orients each edge. Fixing a clockwise orientation
// of each 2d face of Polar(P) we obtain a list of relative orientations. Note that we form OrderedVertices
// in case repeated calls to Vertices(edge) does not consistently return vertices of 'edge' in the same order.
RelativeOrientations := AssociativeArray();
OrderedVertices := AssociativeArray();
for edge in Edges(Polar(P)) do
    OrderedVertices[edge] := [Vertices(edge)[1],Vertices(edge)[2]];
end for;
for face in Faces(Polar(P),2) do
    for edge in [edge_in_face : edge_in_face in Edges(face)] do
        RelativeOrientations[[edge,face]] := GetOrientation(edge,face,OrderedVertices[edge]);
    end for;
end for;
 
// A sequence of all all the Minkowski decompositions possible on each face.
all_decompositions := [[[D : D in M | Dimension(D) gt 0] : M in MinkowskiDecomposition((NumberOfPoints(DualFaces[F])-1)*F)] : F in Faces(P,2)];
sizes := [#M : M in all_decompositions];
 
// The Minkowski decompositions for each face are stored in index_vector.
index_vector := [1 : M in [1..#all_decompositions]];
 
// A fixed list of vertices for polygon5.
hex_vertices := [[1,0],[1,1],[0,1],[-1,0],[-1,-1],[0,-1]]; 
 
// Test sets of Minkowski decompositions. Each entry consists of a pair of subsets of [1,2,3,4,5,6].
interesting_values := [
[[],[]], [[1],[]], [[1,2],[]], [[1,2,3],[]], [[1,2,3,4],[]], [[1,2,3,4,5],[]],[[1,2,3,4,5,6],[]],
[[1],[1]], [[1,2],[1]], [[1,2,3],[1]], [[1,2,3,4],[1]], [[1,2,3,4,5],[1]],[[1,2,3,4,5,6],[1]],
[[1,2],[1,2]], [[1,2,3],[1,2]], [[1,2,3,4],[1,2]], [[1,2,3,4,5],[1,2]],[[1,2,3,4,5,6],[1,2]],
[[1,2,3],[1,2,3]], [[1,2,3,4],[1,2,3]], [[1,2,3,4,5],[1,2,3]],[[1,2,3,4,5,6],[1,2,3]],
[[1,2,3,4],[1,2,3,4]], [[1,2,3,4,5],[1,2,3,4]], [[1,2,3,4,5,6],[1,2,3,4]],
[[1,2,3,4,5],[1,2,3,4,5]], [[1,2,3,4,5,6],[1,2,3,4,5]],
[[1,2,3,4,5,6],[1,2,3,4,5,6]],
];

results := [1 : i in [1..#interesting_values]];

for choice in interesting_values do
    index_vector := [1 : M in [1..#all_decompositions]];
    fail := false;
 	
    for val_a in choice[1] do
	   index := Index(Faces(P,2),polygon6*Polytope([hex_vertices[val_a]]));
        index_vector[index] := 2;
    end for;
    for val_b in choice[2] do
	   index := Index(Faces(P,2),Polytope([hex_vertices[val_b]])*polygon6);
   	    index_vector[index] := 2;
    end for;
 
    // The sequence of Minkowski decompositions determined by index_vector.
    decomps := [all_decompositions[i][index_vector[i]] : i in [1..#all_decompositions]];
 
    // Solutions of the linear conditions imposed on the slopes of the PL function are realised as 
    // the kernel of the matrix 'adjacency'.
    adjacency := Matrix(QQ,#Edges(P),0,[]);
 
    for i in [1..#Faces(P,2)] do
 
        // Store the face indexed by i
        face := Faces(P,2)[i];
 
        // adjacency has a block decomposition. Each block corresponds to a 2D face F of P.
        // Columns of each block are indexed by summands of F in the chosen decomposition.
        tall_block := ZeroMatrix(QQ,#Edges(P),#decomps[i]);
 
        // k indexes columns of tall_block, while j indexes rows.
        for k in [1..#decomps[i]] do
            for j in [1..#Edges(P)] do
 
                // Check two conditions:
                // 1) whether the edge E of P indexed by is contained in the edges of the face F indexed by i,
                // 2) if so, whether E is an edge of F corresponding to the Minkowski factor indexed by i and k.
                if Edges(P)[j] in Edges(face) and IsEdgeInducedByFactor(Edges(P)[j],decomps[i][k],face)  then
                    tall_block[j,k] := RelativeOrientations[ [DualFaces[face], DualFaces[Edges(P)[j]] ] ];
                end if;
            end for;
        end for;
 
        adjacency := HorizontalJoin(adjacency,tall_block);
    end for;
 
    // We test the kernel of adjacency. Each vector in the kernel of adjacency corresponds to a tuple of slopes
    // for a PL function. To construct a strictly convex function we need to be able to choose different slopes
    // for the various Minkowski factors of a single 2D face of P.
 
    // 'position' records the location in each vector in the kernel corresponding to a 2D face F.
    // 'positions' stores these locations for later use.
    position := 1;
    positions := [1 : i in Faces(P,2)];
 
    for i in [1..#Faces(P,2)] do
        // 'zeroes' records the pairs of Minkowski factors which have the same slope for every member of a basis
        // of the kernel of adjacency so far tested. We start assuming no pair of slopes can be seperated.
        zeroes := {1..#decomps[i]*(#decomps[i]-1)/2};
 
        // store the indices of the factors of F in 'indices'
        indices := [position..position+#decomps[i]-1];
 
        // Iterate over elements in a basis of the kernel of 'adjacency'.
        for b in Basis(Kernel(Transpose(adjacency))) do
 
            // Form the differences of slopes corresponding to various Minkowski factors of F. 
            seq := [Eltseq(b)[j] - Eltseq(b)[k] : j in indices, k in indices | j lt k];
 
            // Reduce the list of inseperable pairs of slopes using seq.
            zeroes := zeroes meet {j : j in [1..#seq] | seq[j] eq 0};
        end for;
 
        // If there is at least one pair of Minkowski factors for which slopes must coincide, record failure as a zero.
        if #zeroes ge 1 then
            fail := true;
        end if;
 
        positions[i] := position;
        position := position + #decomps[i];
    end for;
 
    // If this choice of Minkowski decompositions was successful then print the configuration of faces.
    if not fail then
        print "Success:", choice;

    // Uncomment the following to print the Minkowski decompositions of hexagonal faces. 
	// for j in [1..5] do
	//	decomps[Index(Faces(P,2),polygon6*Polytope([pent_vertices[j]]))];
	//end for;
    end if;
end for;
