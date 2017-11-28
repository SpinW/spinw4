function setmatrix(obj, varargin)
% sets exchange tensor values
% 
% ### Syntax
% 
% `setmatrix(obj,Name,Value)`
% 
% ### Description
% 
% `setmatrix(obj,Name,Value)` sets the value of a selected matrix based on
% symmetry analysis.
% 
% ### Examples
% 
% This example will set 'J1' coupling to the 6th symmetry allowed matrix,
% with prefactor 0.235.
% ```
% setmatrix(crystal,'label','J1','pref',{[6 0.235]})
% ```
% This will set 'J2' to antiferromagnetic Heisenberg exchange, with value
% of 1.25 meV.
% ```
% setmatrix(crystal,'label','J2','pref',{1.25})
% ```
%
% ### Input Arguments
% 
% `obj`
% : [spinw] object.
% 
% ### Name-Value Pair Arguments
% 
% One of the below options has to be given:
% 
% `'mat'`
% : Label or index of the matrix that is already assigned to
%   a bond, anisotropy or g-tensor.
% 
% `'bond'`
% : Index of the bond in `spinw.coupling.idx`, e.g. 1 for first neighbor.
% 
% `'subIdx'`
% : Selects a certain bond within the symmetry equivalent bonds, within
%   default value is 1.
% 
% `'aniso'`
% : Label or index of the magnetic atom that has a single ion
%   anisotropy matrix is assigned, e.g. `'Cr3'` to select anisotropy on
%   atoms with this label.
% 
% `'gtensor'`
% : Label or index of the magnetic atom that has a g-tensor is 
%   assigned.
% 
% Optional inputs:
% 
% `'pref'`
% : Defines prefactors as a vector for the symmetry allowed
%           components in a row vector with $n_{symMat}$ number of elements. Alternatively, if only
%           a few of the symmetry allowed matrices have non-zero
%           prefactors, use:
%   ```
%   {[6 0.1 5 0.25]}
%   ```
%   This means, the 6th symmetry allowed matrix have prefactor 0.1,
%           the 5th symmetry allowed matrix have prefactor 0.25. Since
%           Heisenberg isotropic couplings are always allowed, a cell with
%           a single element will create a Heisenberg coupling, example:
%   ```
%   {0.1}
%   ```
%   This is identical to `obj.matrix.mat = eye(3)*0.1`.
%           For DM interactions (antisymmetric coupling matrices), use
%           three element vector in the cell:
%   ```
%   {[D1 D2 D3]}
%   ```
%   In this case, these will be the prefactors of the 3
%           antisymmetric symmetry allowed matrices. In case no crystal
%           symmetry is defined, these will define directly the components
%           of the  DM interaction in the $xyz$ coordinate system. Be
%           carefull with the sign of the DM interaction, it depends on the
%           order of the two interacting atoms! Default value is `{1}`.
%           For anisotropy matrices antisymmetric matrices are not allowed.
% 
% ### Output Arguments
% 
% The selected `obj.matrix.mat` will contain the new value.
% 
% ### See Also
% 
% [spinw] \| [spinw.gencoupling] \| [spinw.getmatrix]
% 
% *[DM]: Dzyaloshinskii-Moriya
%

if nargin == 1
    help spinw.setmatrix
    return
end

[aMat, param] = obj.getmatrix(varargin{:},'fid',0);

if isempty(param.pref)
    % Identity matrix (Heisenberg coupling) is always allowed by symmetry!
    if param.g_idx
        % g-tensor default is 2 for spin only magnetic moment
        aMat = 2*eye(3);
    else
        % identity matrix
        aMat = eye(3);
    end
end

if param.matIdx > 0
    obj.matrix.mat(:,:,param.matIdx) = aMat;
else
    error('spinw:setmatrix:WrongInput','It is not possible to unambiguously select a matrix from the input options!')
end

end