function pOp = point(symOp, r)
% determines point group symmetry at a given position
%
% pOp = SWSYM.POINT(symOp, r)
%
% The function determines point group symmetry in an arbitrary position in
% the unit cell in any space group. Returns all the generators of the point
% group.
%
% Input:
%
% symOp         Symmetry operators of the space group stored in a matrix
%               with dimensions of [3 4 nOp].
% r             Position in the unit cell, dimensions are [3 1].
%
% Output:
%
% pOp           Point group operators, dimensions are [3 3 npOp], these
%               operators act on the relative atomic positions (they are in
%               the lattice coordinate system). To convert them to
%               Cartesian coordinate system, use:
%                   R = A*pOp(:,:,ii)*inv(A)
%               Where A is a 3x3 matrix, containing the basis vectors of
%               the lattice as column vectors.
%
% See also SWSYM.GENERATOR, SWSYM.OPERATOR, SWSYM.POSITION.
%

if nargin == 0
    help swsym.point
    return
end

[~,~,info] = swsym.position(symOp,r);
% point group operators are the ones that does NOT move the atom
pOp = symOp(:,1:3,~info.ismoved{1});

end