
// Given an edge e and Minkowski factor f of a polygon poly, determine whether e determines an edge of f.
function IsEdgeInducedByFactor(e,f,poly)

	// We consider two cases, depending on the dimension of f. If f is a segment
	// we check whether e is parallel to f.
	if Dimension(f) eq 1 then
			E := Explode(Edges(f));
			if AreProportional(Vertices(E)[1] - Vertices(E)[2], Vertices(e)[1] - Vertices(e)[2]) or
			   AreProportional(Vertices(E)[1] - Vertices(E)[2], Vertices(e)[2] - Vertices(e)[1]) then
				return true;
			end if;
	end if;

	// If f is two dimensional, then check whether the normal vector of e in poly is a ray generator
	// of the normal fan of f.
	if Dimension(f) eq 2 then
		if Explode(Rays(NormalCone(poly,e))) in Rays(NormalFan(f)) then
			return true;
		end if;
	end if;

	return false;
end function;

function GetOrientation(edge,face,vs)

	// Put 'face' in a 2D lattice. This implicitly fixes an orientation of 'face'.
	flat_face,emb,shift := PolyhedronInSublattice(face);
	flat_edge := Preimage(emb,edge-shift);

	// We construct an interior point of 'flat_face' which we treat as the origin
	// to construct an orientation of its boundary
	interior_point := (Vertices(flat_face)[1] + Vertices(flat_face)[2] + Vertices(flat_face)[3])/3;

	// Return the sign of the determinant of the 2x2 matrix formed by the vertices of 'edge'
	// after projection to the minimal sublattice containing 'face'.
	
	ordered_vertices :=[Preimage(emb,vs[i]-shift) - interior_point : i in [1,2]];
	return Sign(Determinant(Matrix(Rationals(),2,2,[Eltseq(ordered_vertices[i]) : i in [1,2]])));
end function;

// Return the dual face of a reflexive polytope.
function DualFace(P,face)
	nc := NormalCone(P,face);
	return Polytope(RGenerators(nc));
end function;


function IsBadRectangle(face,dual_face)

	// Check the face is a quadrilateral.
	if #Edges(face) ne 4 then
		return false;
	end if;

	// Check the length of dual edge is equal to 1.
	if NumberOfPoints(dual_face) gt 2 then
		return false;
	end if;

	// Test whether either pair of opposite edges both have length one.
	// Equivalently, we ask whether every vertex is contained in an edge of length one.
	if #(&join[SequenceToSet(Vertices(edge)) : edge in Edges(face) | NumberOfPoints(edge) eq 2 ]) eq 4 then
		return true;
	end if;

	return false;
end function;

function IsBadTriangle(face,dual_face)

	// Check the face is a triangle.
	if #Edges(face) ne 3 then
		return false;
	end if;

	// Check the first case: length two dual edge, standard triangle face.
	if NumberOfPoints(dual_face) eq 3 and NumberOfPoints(face) eq 3 then
		return true;
	end if;

	// Check the second case: length one dual edge, dilated triangluar face.
	// First check the polygon is a triangle with 6 points, is hollow, and has an edge of length 2. 
	if NumberOfPoints(dual_face) eq 2 and NumberOfPoints(face) eq 6 and NumberOfPoints(Edges(face)[1]) eq 3 and NumberOfInteriorPoints(face) eq 0 then
		return true;
	end if;

	return false;
end function;

// Start here

QQ := Rationals();
interval := Polytope([[-1],[1]]);
polygon4 := Polytope([[-1,0],[0,1],[1,0],[0,-1]]);
polygon5 := Polytope([[-1,-1],[-1,0],[0,1],[1,0],[0,-1]]);
polygon6 := Polytope([[1,0],[1,1],[0,1],[-1,0],[-1,-1],[0,-1]]);
polygon7 := Polytope([[-1,-1],[1,-1],[1,0],[0,1],[-1,1]]);
polygon8 := Polytope([[-1,-1],[1,-1],[-1,1],[1,1]]);
polygon9 := Polytope([[-1,-1],[2,-1],[-1,2]]);

// P := polygon6*polygon6;
P := polygon6*polygon8;

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

// Find all faces which fall into special cases.
bad_rectangles := [Index(Edges(P),E) : E in Edges(P) | IsBadRectangle(DualFaces[E],E)];
bad_triangles := [Index(Edges(P),E) : E in Edges(P) | IsBadTriangle(DualFaces[E],E)];

// Check I've found all the bad triangles.
#bad_triangles;

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

// A fixed list of vertices for polygon6.
tri_vertices := [[-1,-1],[-1,2],[2,-1]];
wide_quad_vertices := [[-1,-1],[1,-1],[1,1],[-1,1]];
hex_vertices := [[1,0],[1,1],[0,1],[-1,0],[-1,-1],[0,-1]];
pent_vertices := [[-1,-1],[-1,0],[0,1],[1,0],[0,-1]];
quad_vertices := [[-1,0],[0,1],[1,0],[0,-1]];


// Test sets of Minkowski decompositions.
//interesting_values := [
//[[1],[]]
//];

//interesting_values := [
//[],[1],[1,2],[1,2,3],[1,2,3,4],[1,2,3,4,5]
//];

interesting_values := [
[],[1],[1,2],[1,3],[1,2,3],[1,2,3,4]
];

results := [1 : i in [1..#interesting_values]];

for choice in interesting_values do
	index_vector := [1 : M in [1..#all_decompositions]];
	fail := false;

	for val_a in choice do
	index_vector[Index(Faces(P,2),polygon6*Polytope([wide_quad_vertices[val_a]]))] := 2;
	end for;
	//for val_b in choice[2] do
	//	index_vector[Index(Faces(P,2),Polytope([hex_vertices[val_b]])*polygon6)] := 2;
	//end for;

	//index_vector;

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



	for i in bad_rectangles do

		// 'F' stores the bad_ractangle itself.
		F := DualFaces[Edges(P)[i]];
	
		opposite_pair := [Edges(F)[1],Edges(F)[
			Explode([j : j in [1..4] | #SequenceToSet(Vertices(Edges(F)[j]) cat Vertices(Edges(F)[1])) eq 4 ])] ];

		if [NumberOfPoints(opposite_pair[j]) : j in [1,2]] eq [2,2] then
			edge_pair := [edge : edge in Edges(F) | edge notin opposite_pair];;
		else 
			edge_pair := opposite_pair;
		end if;

		// We now find the indices of faces dual to the edges in edge_pair.
		dual_face_indices := [Index(Faces(P,2), DualFaces[edge]) : edge in edge_pair];
		local_positions := [positions[j] : j in dual_face_indices];

		zeroes := {1..#decomps[dual_face_indices[1]]*#decomps[dual_face_indices[2]]};

		// store the indices of the factors of F in 'indices'
		indices_1 := [local_positions[1]..local_positions[1]+#decomps[dual_face_indices[1]]-1];
		indices_2 := [local_positions[2]..local_positions[2]+#decomps[dual_face_indices[2]]-1];
		local_orientations := [RelativeOrientations[[edge,face]]: edge in edge_pair];

		// Iterate over elements in a basis of the kernel of 'adjacency'.
		for b in Basis(Kernel(Transpose(adjacency))) do

			// Form the differences of slopes corresponding to various Minkowski factors. 
			seq := [local_orientations[1]*Eltseq(b)[j] + local_orientations[2]*Eltseq(b)[k] : j in indices_1, k in indices_2];

			// Reduce the list of inseperable pairs of slopes using seq.
			zeroes := zeroes meet {j : j in [1..#seq] | seq[j] eq 0};
		end for;

		// If there is at least one pair of Minkowski factors for which slopes must coincide, record failure as a zero.
		if #zeroes ge 1 then
			fail := true;
		end if;
	end for;

	for i in bad_triangles do

		// 'F' stores the bad_triangle itself.
		F := DualFaces[Edges(P)[i]];

		zeroes := {1..6};
		slope := [0,0];

		// Iterate over elements in a basis of the kernel of 'adjacency'.
		for b in Basis(Kernel(Transpose(adjacency))) do
			for vert in Vertices(F) do

				// Obtain the edges neighbouring 'vert', and the edge opposite 'vert'.
				neighbours := [edge : edge in Edges(F) | vert in Vertices(F)];
				omitted_edge := Explode([edge : edge in Edges(F) | vert notin Vertices(edge)]);

				// 'omitted_position' holds the column index of the edge opposite 'vert'.
				omitted_position := positions[Index(Faces(P,2), PolarDualFaces[omitted_edge])];

				for j in [1,2] do

					// 'local_position' holds the column index of the jth edge containing 'vert'.
					local_position := positions[Index(Faces(P,2), PolarDualFaces[neighbours[j]])];

					// If vert is the first vertex of this edge, then the associated slope is the minimal
					// of the two associated to segments along the jth edge containing 'vert'.
					if Index(OrderedVertices[neighbours[j]],vert) eq 1 then

						slope[j] := RelativeOrientations[ [neighbours[j],F] ]*
							Min([Eltseq(b)[local_position],Eltseq(b)[local_position + 1]]);

					else

						slope[j] := RelativeOrientations[ [neighbours[j],F] ]*
							Max([Eltseq(b)[local_position],Eltseq(b)[local_position + 1]]); 

					end if;
				end for;

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
		print "Success:", choice;//#choice[1], #choice[2];
	end if;
end for;