
// Given an edge e and Minkowski factor f of a polygon poly, determine whether e determines an edge of f.
function IsEdgeInducedByFactor(e,f,poly)

	// We consider two cases, depending on the dimension of f. If f is a segment
	// we check whether e is parallel to f.
	if Dimension(f) eq 1 then
		return (AreProportional(Vertices(f)[1] - Vertices(f)[2], Vertices(e)[1] - Vertices(e)[2]) or
			   AreProportional(Vertices(f)[1] - Vertices(f)[2], Vertices(e)[2] - Vertices(e)[1]));
	end if;

	// If f is two dimensional, then check whether the normal vector of e in poly is a ray generator
	// of the normal fan of f.
	if Dimension(f) eq 2 then
		return (NormalCone(poly,e) in NormalFan(f));
	end if;

	return false;
end function;

// Takes as input an edge, a face containing the edge and an ordered list
// containing the vertices of the edge. The latter input fixes an orientation
// of the edge.
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
	return Polytope(Rays(nc));
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