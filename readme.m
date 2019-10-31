These files contain Magma source code used to generate the data displayed in Tables 3-6 in the article titled
'Smoothing Calabi-Yau toric hypersurfaces using the Gross-Siebert algorithm' by T. Prince.

File 'generate_data_for_mxn.m' contains source code to generate regular pairs (P_{m,n},D), where P_{m,n} is the
product of reflexive polygons of normalised volume m and n respectively which admit a standard decomposition.

-------------------------------------

Notes for 'generate_data_for_4x6.m'

Input: a disjoint pair A and B of subsets of [1,2,3,4].

Dilating a hexagonal face of P_{4,6} by a factor of two, the resulting polygon admits Minkowski decompositions into
1) six line segments,
2) three line segments and a pair of triangles, or,
3) four triangles.

Vertices of the integral square with vertices
[1,0],[0,1],[-1,0],[0,-1]
are numbered from 1 to 4. Hexagonal faces of P_{4,6} are in bijection with the vertices of this integral square.
Vertices in A are assigned the Minkowski decomposition (2).
Vertices in B are assigned the Minkowski decomposition (3).
All other vertices are assigned the first Minkowski decomposition (into 6 line segments).

-------------------------------------

Notes for 'generate_data_for_5x6.m'

Input:  disjoint pair A and B of subsets of [1,2,3,4,5].

Hexagonal faces deocompose admit Minkowski decompositions into
1) a triple of line segments, or,
2) a pair of triangles.

Dilating such a face by a factor of two the resulting polygon admits Minkowski decompositions into
1) six line segments,
2) three line segments and a pair of triangles, or,
3) four triangles.

Vertices of the integral pentagon with vertices
[1,0],[0,1],[-1,0],[-1,-1],[0,-1]
are numbered from 1 to 5. Hexagonal faces of P_{5,6} are in bijection with the vertices of this integral pentagon.
Vertices in A are assigned the Minkowski decomposition (2) on each list.
Vertices in B are assigned the Minkowski decomposition (3) (into 4 triangles).
All other vertices are assigned the first Minkowski decomposition (into 3 or 6 line segments, 
depending on the length of the edge of Polar(P_{5,6}) dual to the corresponding 2d face of P.)
Note that only the values 1 and 2 may appear in the second sequence.

Output: a list of such pairs for which the corresponding standard decomposition D is regular.

-------------------------------------

Notes for 'generate_data_for_6x6.m'

Input:  a pair A and B of subsets of [1,2,3,4,5,6].

Hexagonal faces deocompose admit Minkowski decompositions into
1) a triple of line segments, or,
2) a pair of triangles.

Vertices of the integral hexagon P_6 with vertices
[1,0],[1,1],[0,1],[-1,0],[-1,-1],[0,-1]
are numbered from 1 to 6. There are twelve Hexagonal faces of P_{6,6}, which are in bijection
with the vertices of this integral pentagon, together with an additional binary choice.

2-faces of the form (a x P_6), for a in A are assigned the Minkowski decomposition (2) on each list.
2-faces of the form (P_6 x a), for a in A are also assigned the Minkowski decomposition (2) on each list.
All other vertices are assigned the first Minkowski decomposition (into 3 line segments).

Output: a list of such pairs for which the corresponding standard decomposition D is regular.

-------------------------------------

Notes for 'generate_data_for_6x7.m'

Input: a subset A of [1,2,3,4,5].

Hexagonal faces deocompose admit Minkowski decompositions into
1) a triple of line segments, or,
2) a pair of triangles.

Vertices of the integral pentagon with vertices
[-1,-1],[1,-1],[1,0],[0,1],[-1,1]
are numbered from 1 to 5. Hexagonal faces of P_{5,6} are in bijection with the vertices of this integral pentagon.
Vertices in A are assigned the Minkowski decomposition (2) on each list.
All other vertices are assigned the first Minkowski decomposition.

Output: a list of such pairs for which the corresponding standard decomposition D is regular.

