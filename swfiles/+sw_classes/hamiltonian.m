classdef hamiltonian < sw_classes.hamiltonian_base
    % Class to calculate the (bilinear) Hamiltonian: sum_ij Si x Jij x Sj
    % after eq (25) and (26) of Toth & Lake
    properties(SetAccess=private)
        % (Properties defined in sw_classes.hamiltonian_base)
        % User input quantities:
        %dR         % The difference position vector between atom i and j
        % Deduced quantities:
        %nCoupling  % Number of couplings
        %nMagExt    % Number of magnetic atoms in the extended supercell
        %subblocks  % A vector of matrix elements (nBlock x 1 cell array, each block is nCoupling x 1 vector)
        %subidx     % A (nBlock*nCoupling) x 2 list of indices to convert .subblocks to Hamiltonian matrix.
        %diagblocks % A (2*nMagExt x 2*nMagExt) matrix of the -C hkl-independent diagonal blocks
    end
    methods
        function self = hamiltonian(JJ, dR, atom1, atom2, u, v, S_mag)
            % Calculates the (bilinear) Hamiltonian: sum_ij Si x Jij x Sj
            % Syntax: ham = sw_classes.hamiltonia(JJ, dR, atom1, atom2, u, v, S_mag)
            % Inputs:
            %   JJ    % The magnetic couplings (3 x 3 x nCoupling array)
            %   dR    % The difference position vector between atom i and j
            %   atom1 % The atom1 indices (1 x nCoupling vector)
            %   atom2 % The atom2 indices (1 x nCoupling vector)
            %   u     % The u vectors (zed in original code, 3 x nMagExt)
            %   v     % The v vectors (eta in original code, 3 x nMagExt)
            %   S_mag % The magnetic moment magnitude (1 x nMagExt vector)

            % Don't know why we need to do this, but otherwise need to take a conj later
            u = conj(u);  % (or the equations don't match Toth & Lake...)

            self.dR = dR;
            self.nMagExt = max([atom1; atom2]);

            % Input u and v are in base order, reorder them according to input
            % indices and expand dimension to match JJ (implicit expansion doesn't work)
            uT_i = repmat(permute(u(:,atom1),[1 3 2]),[1 3 1]);
            u_j = repmat(permute(u(:,atom2),[3 1 2]),[3 1 1]);
            vT_i = repmat(permute(v(:,atom1),[1 3 2]),[1 3 1]);
            v_j = repmat(permute(v(:,atom2),[3 1 2]),[3 1 1]);

            % The prefactors
            S_mag = S_mag(:);
            SiSj = sqrt(S_mag(atom1).*S_mag(atom2));

            % Computes the submatrices as vectors and their indices
            % The Hamiltonian matrix is (2*nMagExt) x (2*nMagExt)
            % with the basis eq(24) formed first of all the creation
            % then all the anihilation operators.
            % The indices atom1 and atom2 cover only one set (max(atom1,atom2)==nMagExt).
            nMagExt = self.nMagExt;

            % The upper left block A^ij in eq (26) of Toth & Lake:
            AD0 =  SiSj.*squeeze(sum(sum(uT_i.*JJ.*conj(u_j),2),1));
            idxA1 = [atom1 atom2];            % For upper left block
            idxD1 = [atom1 atom2] + nMagExt;  % For lower right block

            % The upper right block B^ij in eq (26) of Toth & Lake:
            BC0 =  SiSj.*squeeze(sum(sum(uT_i.*JJ.*u_j,2),1));
            idxB = [atom1 (atom2 + nMagExt)]; % For upper right block

            self.subblocks = {AD0 2*BC0 conj(AD0)};
            self.subidx = [idxA1; idxB; idxD1];

            % The diagonal elements C^ij in eq (26) of Toth & Lake:
            AD = squeeze(sum(sum(vT_i.*JJ.*v_j,2),1));
            A20 = -S_mag(atom2).*AD;  % This is -C in A(k)-C in upper left of eq(25)
            D20 = -S_mag(atom1).*AD;  % This is -C in conj(A(-k))-C in lower right of eq(25)
            idxA2 = [atom1 atom1];            % For upper left block
            idxD2 = [atom2 atom2] + nMagExt;  % For lower right block

            self.diagblocks = accumarray([idxA2; idxD2], 2*[A20; D20], [1 1]*2*nMagExt);
        end
    end
end
