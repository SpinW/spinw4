classdef biquadratic_hamiltonian < sw_classes.hamiltonian_base
    % Class to calculate the biquadratic Hamiltonian sum_ij Bij * (Si.Sj)^2
    methods
        function self = biquadratic_hamiltonian(JJ, dR, atom1, atom2, u, v, S_mag)
            % Calculates the biquadratic Hamiltonian
            % Syntax: bq_ham = sw_classes.biquadratic_hamiltonia(JJ, dR, atom1, atom2, u, v, S_mag)
            % Inputs:
            %   JJ     % The biquadratic couplings (1 x nCoupling vector)
            %   dR     % The difference position vector between atom i and j
            %   atom1  % The atom1 indices (1 x nCoupling vector)
            %   atom2  % The atom2 indices (1 x nCoupling vector)
            %   u      % The u vectors (zed in original code, 3 x nMagExt)
            %   v      % The v vectors (eta in original code, 3 x nMagExt)
            %   S_mag  % The magnetic moment magnitude (1 x nMagExt vector)

            self.dR = dR;
            self.nMagExt = max([atom1 atom2]);
            JJ = transpose(JJ);
            u = u';
            v = transpose(v);

            % In this case the coupling constants JJ are scalars
            % (In the bilinear Hamiltonian they are 3x3 matrices)
            % So we just compute the uv factors converting from
            % the spins in the rotating to boson operators here
            bqM = sum(v(atom1,:) .* v(atom2,:), 2);
            bqN = sum(v(atom1,:) .* u(atom2,:), 2);
            bqO = sum(u(atom1,:) .* u(atom2,:), 2);
            bqP = sum(conj(u(atom1,:)) .* u(atom2,:), 2);
            bqQ = sum(u(atom1,:) .* v(atom2,:), 2);

            % The prefactors
            Si = S_mag(atom1);
            Sj = S_mag(atom2);
            nMagExt = self.nMagExt;

            % Now compute the blocks of the Hamiltonian matrix
            bqA0 = (Si.*Sj).^(3/2) .* (bqM.*conj(bqP) + bqQ.*conj(bqN)) .* JJ;
            bqB0 = (Si.*Sj).^(3/2) .* (bqM.*bqO + bqQ.*bqN) .* JJ;
            % Creates the serial indices for every matrix element in ham matrix.
            % Aij(k) matrix elements (b^+ b)
            idxbqA  = [atom1' atom2'];
            % b b^+ elements
            idxbqA2 = [atom1' atom2'] + nMagExt;
            % Bij(k) matrix elements (b^+ b^+)
            idxbqB  = [atom1' (atom2' + nMagExt)];

            self.subblocks = {bqA0 conj(bqA0) 2*bqB0};
            self.subidx = [idxbqA; idxbqA2; idxbqB];

            % The hkl-independent diagonal blocks
            bqC  = Si.*Sj.^2 .* (conj(bqQ).*bqQ - 2*bqM.^2) .* JJ;
            bqD  = Si.*Sj.^2 .* (bqQ).^2 .* JJ;
            idxbqC  = [atom1' atom1'];
            idxbqC2 = [atom1' atom1'] + nMagExt;
            idxbqD  = [atom1' (atom1' + nMagExt)];

            self.diagblocks = accumarray([idxbqC; idxbqC2; idxbqD], [bqC; bqC; 2*bqD], [1 1]*2*nMagExt);
        end
    end
end
