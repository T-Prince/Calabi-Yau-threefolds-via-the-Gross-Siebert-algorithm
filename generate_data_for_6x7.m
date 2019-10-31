
load "utilities.m";
 
QQ := Rationals();
interval := Polytope([[-1],[1]]);
polygon4 := Polytope([[-1,0],[0,1],[1,0],[0,-1]]);
polygon5 := Polytope([[-1,-1],[-1,0],[0,1],[1,0],[0,-1]]);
polygon6 := Polytope([[1,0],[1,1],[0,1],[-1,0],[-1,-1],[0,-1]]);
polygon7 := Polytope([[-1,-1],[1,-1],[1,0],[0,1],[-1,1]]);
polygon8 := Polytope([[-1,-1],[1,-1],[-1,1],[1,1]]);
polygon9 := Polytope([[-1,-1],[2,-1],[-1,2]]);
 
// P := polygon6*polygon6;
P := polygon6*polygon7;
 
// Cache all the dual face information.
DualFaces := AssociativeArray();
for i in [1,2] do
    for face in Faces(P,i) do
        DualFaces[face] := DualFace(P,face);
    end for;
end for;
 
// Cache all the dual face information.
PolarDualFaces := AssociativeArray();
Q := Polar(P);
for i in [1,2] do
    for face in Faces(Q,i) do
        PolarDualFaces[face] := DualFace(Q,face);
    end for;
end for;

// 'bad_triangles' contains the indices of all edges E of P dual to triangles of Polar(P) such that either:
// 1) Dual(E) is 2 x standard triangle and E has length 1 or,
// 2) Dual(E) is 1 x standard triangle and E has length 2.
bad_triangles := [Index(Edges(P),E) : E in Edges(P) | IsBadTriangle(DualFaces[E],E)];

 
// The ordering of the vertices of edges of Polar(P) orients each edge. Fixing a clockwise orientation
// of each 2d face of Polar(P) we obtain a list of relative orientations.
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

pent_vertices := [[-1,-1],[1,-1],[1,0],[0,1],[-1,1]];
 
interesting_values := [
    [],[1],[2],[1,2],[1,3],[3,4],[1,2,3],[2,3,4],[1,2,3,4],[2,3,4,5],[1,2,3,4,5] 
];

results := [1 : i in [1..#interesting_values]];
 
for choice in interesting_values do
    index_vector := [1 : M in [1..#all_decompositions]];
    fail := false;
 
    for val_a in choice do
    index_vector[Index(Faces(P,2),polygon6*Polytope([pent_vertices[val_a]]))] := 2;
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
 
    for i in bad_triangles do
 
        // 'F' stores the bad_triangle in Polar(P) itself.
        F := DualFaces[Edges(P)[i]];
 
    	// We test six different possible quadrilateral domains of linearity.
        tests := {1..6};

	   // the pair 'slopes' stores the slopes of a PL function on segments neighbouring a vertex
	   // of the given bad triangle.
        slope := [0,0];
 
        // Iterate over elements in a basis of the kernel of 'adjacency'. These store slopes of possible PL functions.
        for b in Basis(Kernel(Transpose(adjacency))) do
            for vert in Vertices(F) do
 
                // Obtain the edges neighbouring 'vert', and the edge opposite 'vert'.
		        // These are edges of F, and hence of Polar(P).
                neighbours := [edge : edge in Edges(F) | vert in edge];
                omitted_edge := Explode([edge : edge in Edges(F) | vert notin Vertices(edge)]);
 
                // 'omitted_position' holds the column index of the edge opposite 'vert'.
                omitted_position := positions[Index(Faces(P,2), PolarDualFaces[omitted_edge])];
                
                // We determine the slope of the PL function along each segment containing 'vert'.
                for j in [1,2] do
 
                    // 'local_position' holds the column index of the jth edge containing 'vert'.
                    local_position := positions[Index(Faces(P,2), PolarDualFaces[neighbours[j]])];
 
                    // If vert is the first vertex of this edge, then the associated slope is the minimal
                    // of the two associated to segments along the jth edge containing 'vert'.
                    if Index(OrderedVertices[neighbours[j]],vert) eq 1 then
 
                        // Set the slope to be the minimum of the two, then adjust the sign if necessary.
                        // This sign change indicates that, ordering the triangle clockwise, whether we
                        // are travelling into or out from the given vertex.
                        slope[j] := RelativeOrientations[ [neighbours[j],F] ]*
                            Min([Eltseq(b)[local_position],Eltseq(b)[local_position + 1]]);
 
                    else
                        
                        // Set the slope to be the maxmum of the two, again adjusting the sign if necessary.
                        slope[j] := RelativeOrientations[ [neighbours[j],F] ]*
                            Max([Eltseq(b)[local_position],Eltseq(b)[local_position + 1]]); 
 
                    end if;
                end for;
 
                // We now need to check whether the sum of these slopes match the slope corresponding to
                // either of the two segments of the edge of 'F' oppositive 'vert'.
                if slope[1] + slope[2] + RelativeOrientations[[omitted_edge,F]]*
                            Min([Eltseq(b)[omitted_position], Eltseq(b)[omitted_position + 1]]) ne 0 then
 
                    zeroes := zeroes diff {2*Index(Vertices(F),vert)-1};
 
                end if;
                if slope[1] + slope[2] + RelativeOrientations[[omitted_edge,F]]*
                            Max([Eltseq(b)[omitted_position], Eltseq(b)[omitted_position + 1]]) ne 0 then
 
                    zeroes := zeroes diff {2*Index(Vertices(F),vert)};
 
                end if;
            end for;
        end for;
         
        // If there is at least one pair of Minkowski factors for which slopes must coincide, record failure as a zero.
        if #zeroes ge 1 then
            fail := true;
        end if;
    end for;
 
    // If this choice of Minkowski decompositions was successful then print the configuration of faces.
    if not fail then
        print "Success:", choice, #choice;
    end if;
end for;
