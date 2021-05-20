classdef (Abstract) hamiltonian_base < handle
    % An abstract base class for Hamiltonian calculations.
    % It only defines the ham_k() method to calculate the hkl-dependent part of the Hamiltonian.
    % The constructor is abstract and must be defined concretely by derived classes.
    % The constructor must also set the attributes defined here.

    properties(SetAccess=protected)
        dR         % The difference position vector between atom i and j
        nMagExt    % Number of magnetic atoms in the extended supercell
        subblocks  % A vector of matrix elements (nBlock x 1 cell array, each block is nCoupling x 1 vector)
        subidx     % A (nBlock*nCoupling) x 2 list of indices to convert .subblocks to Hamiltonian matrix.
        diagblocks % A (2*nMagExt x 2*nMagExt) matrix of the -C hkl-independent diagonal blocks
    end
    
    methods
        function ham = ham_k(self, hkl)
            % Returns the Fourier transformed Hamiltonian matrix
            %
            % ### Syntax
            % `hamk = hamiltonian.ham_k(hkl)
            %
            % ### Input Arguments
            %
            % `hkl`
            % : a 3 x nHkl matrix of hkl q-vectors
            %
            % ### Output Arguments
            %
            % `hamk`
            % : a N x N x nHkl matrix of transformed Hamiltonian matrices

            ExpF = exp(1i * inner(self.dR, hkl));

            % Multiplies subblocks by the Fourier phase factor
            % with implicit dimension expansion:
            % Each element has dimension: nCoupling x 1; ExF is nCoupling x nHkl
            ABCD = cellfun(@(c) c .* ExpF, self.subblocks, 'UniformOutput', false);
            % Append the subblocks and indices into one
            ABCD = cat(1, ABCD{:});
            % Expand the indices to include the hkl indices
            nHkl = size(hkl,2);
            idxAll = extend_indices(self.subidx, nHkl);

            % Compute the square Hamiltonian matrix from the matrix element vector and indices
            ham = accumarray(idxAll,ABCD(:),[2*self.nMagExt 2*self.nMagExt nHkl]);
            % Subtract the C^ij hkl-independent diagonal blocks
            ham = ham + repmat(self.diagblocks,[1 1 nHkl]);
        end
    end
end

function out = inner(A, B)
    % Inner product with dimension expansion.
    % A must be M x N_A and B must be M x N_B
    % N_A and N_B do not have to be the same (M must be).
    % Output is N_A x N_B array
    assert(size(A,1) == size(B,1));
    out = squeeze(sum(A .* permute(B, [1 3 2]), 1));
end

function out = extend_indices(idx0, nHkl)
    nCoupling = size(idx0, 1);
    idx3 = repmat(1:nHkl, [nCoupling 1]);
    out  = [repmat(idx0,[nHkl 1]) idx3(:)];
    out  = out(:,[2 1 3]);
end
