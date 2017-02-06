function varargout = newcell(obj,bvect, bshift)
% changes lattice vectors while keeping atoms
%
% {T} = NEWCELL(obj, bvect, {bshift})
%
% The function defines new unit cell using the 3 vectors contained in
% bvect. The three vectors in lattice units define a parallelepiped. This
% will be the new unit cell. The atoms from the original unit cell will
% fill the new unit cell. Also the magnetic structure and bond and single
% ion property definitions will be erased from the structure.
%
% Input:
%
% obj       spinw class object.
% bvect     Defines the new lattice vectors in the original lattice
%           coordinate system. Cell with the following elements
%           {v1 v2 v3}.
% bshift    Vector defines a shift of the position of the unit cell.
%           Optional.
%
% Output:
%
% T     is a transformation matrix that converts Q points (in reciprocal
%       lattice units) from the old reciprocal lattice to the new
%       reciprocal lattice as follows:
%           Qrlu_new = T * Qrlu_old,
%       where the dimensions of the Q vectors are [1 3].
%
% Example:
%
% tri = spinw;
% tri.genlattice('lat_const',[3 3 5],'angled',[90 90 120])
% tri.addatom('r',[0 0 0])
% tri.newcell({[1 0 0] [1 2 0] [0 0 1]})
% plot(tri)
%
% The example show how to convert a triangular lattice into orthorhombic
% lattice vectors and plots the new unit cell.
%
% See also SPINW.GENLATTICE, SPINW.GENCOUPLING, SPINW.NOSYM.
%

if nargin <= 1
    help spinw.newcell
    return
end

%warning('spinw:newcell:PossibleBug','There might be an error, if there are atoms at the faces of the original cell!')

if ~iscell(bvect) || numel(bvect)~=3
    error('spinw:newcell:WrongInput','Input has to be cell type with 3 vectors inside!');
end

% shift
if nargin == 2
    bshift = [0;0;0];
else
    bshift = bshift(:);
end

% here 3 coordinate systems are used:
% - xyz real space, Cartesian
% - original unit cell a,b and c vectors
% - new unit cell a',b' and c' vectors


% transformation matrix from the new lattice units into the original
% lattice
% v_orig = Tn_o*v_new
Tn_o = [bvect{1}(:) bvect{2}(:) bvect{3}(:)];

% transformation from the original lattice into xyz real space
% xyz = To_xyz * v_orig
To_xyz = obj.basisvector;

% the new basis vectors
basisvector2 = To_xyz * Tn_o;

% new lattice parameters
obj.lattice.lat_const = sqrt(sum(basisvector2.^2,1));
% new angles
bnorm = bsxfun(@rdivide,basisvector2,obj.lattice.lat_const);
obj.lattice.angle = acos(sum([bnorm(:,2).*bnorm(:,3) bnorm(:,1).*bnorm(:,3) bnorm(:,1).*bnorm(:,2)],1));

% coordinates of the corners of the new coordinate system
pp = [zeros(3,1) Tn_o Tn_o(:,1)+Tn_o(:,3) Tn_o(:,2)+Tn_o(:,3) Tn_o(:,1)+Tn_o(:,2) sum(Tn_o,2)];

% number of cells needed for the extension
nExt  = ceil(max(pp,[],2) - min(pp,[],2))'+2;
%obj.mag_str.N_ext = int32(nExt);


% generated atoms
atomList   = obj.atom;
% original number of atoms in the unit cell
nAtom0 = numel(atomList.idx);
atomList.S = obj.unit_cell.S(atomList.idx);
atomList = sw_extendlattice(nExt,atomList);
rExt   = bsxfun(@plus,bsxfun(@times,atomList.RRext,nExt'),(min(pp,[],2)-1));
idxExt = atomList.idxext;



rExt = bsxfun(@plus,rExt,bshift);
% atomic positions in the new unit cell
rNew = inv(Tn_o)*rExt; %#ok<MINV>

epsilon = 10*eps;
% cut atoms outside of the unit cell
%idxCut = any((rNew<0) | (rNew>=(1-eps)),1);
idxCut = any((rNew<-epsilon) | (rNew>(1-epsilon)),1);
rNew(:,idxCut) = [];
idxExt(idxCut) = [];
atomList.Sext(idxCut)   = [];

% atoms are close to the face or origin --> put them exactly
rNew(rNew<epsilon) = 0;

% new lattice simmetry is no-symmetry
obj.lattice.sym     = zeros(3,4,0);
obj.lattice.label   = 'P 0';
% no symmetry operations
obj.sym = false;

% new unit cell defined
obj.unit_cell.r     = rNew;
obj.unit_cell.S     = atomList.Sext;

obj.unit_cell.ff  = obj.unit_cell.ff(:,:,idxExt);

fNames = fieldnames(obj.unit_cell);
fNames = setdiff(fNames,{'r' 'S' 'ff'});

for ii = 1:numel(fNames)
    obj.unit_cell.(fNames{ii}) = obj.unit_cell.(fNames{ii})(:,idxExt);
end

% reset the magnetic structure
% obj.mag_str.F     = zeros(3,0,0);
% obj.mag_str.k     = zeros(3,0);
% obj.mag_str.N_ext = int32([1 1 1]);

% reset the magnetic structure and the bonds
obj.initfield({'coupling' 'mag_str'});

% correct for formula units
obj.unit.nformula = obj.unit.nformula*numel(obj.unit_cell.S)/nAtom0;

% transformation from the original reciprocal lattice into the new
% reciprocal lattice
if nargout>0
    varargout{1} = Tn_o;
end

end