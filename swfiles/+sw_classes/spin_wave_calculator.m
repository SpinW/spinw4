classdef spin_wave_calculator < handle
    properties
        spinWaveObject     % spinw class parent object
        parameters         % a struct of parameters
        magnetic_structure % a sw_classes.magnetic_structure object
        qvectors           % a sw_classes.qvectors object
        hamiltonian        % a sw_classes.hamiltonian object
        bq_hamiltonian     % a sw_classes.biquadratic_hamiltonian object
    end

    properties(Access=private)
        gtensor % a (3 x 3 x nMagExt) matrix
        eta     % The v vector in eq (9) of Toth & Lake
        zed     % The u vector in eq (9) of Toth & Lake
        SI      % A struct with single-ion (SI) properties
        RR      % A 3 x nMagExt matrix with positions of the spins/atoms
        orthWarn = false   % Flag if orthogonal warning is triggered
        posdefWarn = false % Flag if positive-definite warning triggered
    end

    methods(Access=public)
        function self = spin_wave_calculator(spinWaveObject, varargin)
            self.spinWaveObject = spinWaveObject;
            self.spinwave_parse_input(varargin{:});
            self.parseMagneticStructure();
            self.prepareHamiltonian();
        end
        function spectra = calculateSpinWave(self, hkl)
            self.qvectors = sw_classes.qvectors(hkl);
            nChunk = self.computeNumChunk();
            self.qvectors.nChunk = nChunk;
            nTwins = self.checkTwins();
            omega_twin = cell(1,nTwins);
            Sab_twin = cell(1,nTwins);
            Hsave_twin = cell(1,nTwins);
            Vsave_twin = cell(1,nTwins);
            Sabp_twin = cell(1,nTwins);
            omegap_twin = cell(1,nTwins);
            if self.magnetic_structure.incomm
                k_incomm_vec = [-self.magnetic_structure.km', zeros(3,1), self.magnetic_structure.km'];
                calc_type = 'INCOMMENSURATE';
            else
                k_incomm_vec = zeros(3,1);
                calc_type = 'COMMENSURATE';
            end
            fprintf0(self.parameters.fid,['Calculating %s spin wave spectra '...
                '(nMagExt = %d, nHkl = %d, nTwin = %d)...\n'], calc_type, ...
                self.magnetic_structure.nMagExt, self.qvectors.nHkl, nTwins);
            % Gets the transformation matrices for the hkl for each twin
            [~, rotQ] = self.spinWaveObject.twinq([0;0;0]);
            for iTwin = 1:nTwins
                rotC = self.spinWaveObject.twin.rotc(:,:,iTwin);
                % Incommensurate loop
                sz_incomm = size(k_incomm_vec,2);
                Hsave = cell(1, sz_incomm);
                Vsave = cell(1, sz_incomm);
                omega_cell = cell(1, sz_incomm);
                Sab_cell = cell(1, sz_incomm);
                incomm_idx = 0;
                for k_incomm = k_incomm_vec
                    Sab = cell(1, nChunk);
                    incomm_idx = incomm_idx+1;
                    for ii = 1:nChunk
                        hklIdxMEM = self.qvectors.getIdx(ii);
                        % Get chunk and apply user defined transformation in qmat (see also spinw.newcell)
                        hklChunk = self.spinWaveObject.unit.qmat * self.qvectors.getChunk(ii);
                        % If twin, then rotate the current hkl set by its rotation matrix
                        if ~self.parameters.notwin
                            hklChunk = (hklChunk' * rotQ(:,:,iTwin))';
                        end
                        ham = self.calculateHamiltonian(hklChunk, rotC, k_incomm);
                        if self.parameters.saveH
                            Hsave{incomm_idx}(:,:,hklIdxMEM) = ham;
                        end
                        [V, omega_cell{incomm_idx}(:,hklIdxMEM)] = self.diagonaliseHamiltonian(ham);
                        if self.parameters.saveV
                            Vsave{incomm_idx}(:,:,hklIdxMEM) = V;
                        end
                        Sab{ii} = self.spinspincorrel(V, hklChunk);
                        if self.magnetic_structure.incomm && self.parameters.saveSabp && incomm_idx == 2
                            Sabp_twin{iTwin}(:,:,:,hklIdxMEM) = Sab{ii};
                        end
                    end
                    Sab_cell{incomm_idx} = cat(4, Sab{:});
                    if self.magnetic_structure.incomm
                        % Transforms from the rotating frame to the lab frame according to eq (40) of Toth and Lake (2015).
                        if incomm_idx == 1
                            Sab_cell{1} = mmat(Sab_cell{1}, self.magnetic_structure.R1);
                        elseif incomm_idx == 2
                            Sab_cell{2} = mmat(Sab_cell{2}, self.magnetic_structure.R2);
                        elseif incomm_idx == 3
                            Sab_cell{3} = mmat(Sab_cell{3}, conj(self.magnetic_structure.R1));
                        end
                    end
                end
                if self.magnetic_structure.incomm && self.parameters.saveSabp
                    omegap_twin{iTwin} = omega_cell{2};
                end
                omega = cat(1, omega_cell{:});
                Sab = cat(3, Sab_cell{:});
                if self.parameters.sortMode
                    % sort the spin wave modes
                    [omega, Sab] = sortmode(omega,reshape(Sab,9,size(Sab,3),[]));
                    Sab = reshape(Sab,3,3,size(Sab,2),[]);
                end
                if ~self.parameters.notwin
                    omega_twin{iTwin} = omega;
                    Sab_twin{iTwin} = Sab;
                    if self.parameters.saveV
                        Vsave_twin{iTwin} = cat(3, Vsave{:});
                    end
                    if self.parameters.saveH
                        Hsave_twin{iTwin} = cat(3, Hsave{:});
                    end
                    if ~(isdiag(rotC) && all(abs(diag(rotC) - 1) < self.parameters.tol))  % Not identity
                        % Rotates the spin-spin correlation tensor using the twin rotation matrix
                        sSab = size(Sab_twin{iTwin});
                        % TODO: consider using a simple loop of 'reshape->arrayfun->reshape' operation
                        Sab = reshape(Sab_twin{iTwin}, 3, 3, []);   % Originally 3x3xNModexNQ now 3x3x(NMode*NQ)
                        Sab = arrayfun(@(idx)(rotC*Sab(:,:,idx)*(rotC')), 1:size(Sab,3), 'UniformOutput', false);
                        % Sab is now a cell array of 3x3 matrices
                        Sab_twin{iTwin} = reshape(cat(3, Sab{:}), sSab);  % Converts back to normal matrix
                    end
                end
            end
            if self.posdefWarn && ~self.parameters.fitmode
                warning('spinw:spinwave:NonPosDefHamiltonian',['To make the Hamiltonian '...
                    'positive definite, a small omega_tol value was added to its diagonal!'])
            end
            % issue eigorth warning
            if self.orthWarn
                warning('spinw:spinwave:NoOrth','Eigenvectors of defective eigenvalues cannot be orthogonalised at some q-point!');
            end
            % If number of formula units are given per cell normalize to formula unit
            if self.spinWaveObject.unit.nformula > 0
                Sab = Sab / double(self.spinWaveObject.unit.nformula);
            end
            % Creates output structure with the calculated values.
            if self.parameters.notwin
                spectra.omega    = omega;
                spectra.Sab      = Sab;
            else
                spectra.omega    = omega_twin;
                spectra.Sab      = Sab_twin;
            end
            spectra.hkl      = self.qvectors.hkl;
            spectra.hklA     = 2*pi*((self.spinWaveObject.unit.qmat*spectra.hkl)' / self.spinWaveObject.basisvector)';
            spectra.helical  = false; % TODO: sort out incommensurate and helical
            spectra.norm     = false;
            spectra.nformula = double(self.spinWaveObject.unit.nformula);
            % Save different intermediate results.
            if self.parameters.saveV
                if self.parameters.notwin
                    spectra.V = cat(3, Vsave{:});
                else
                    spectra.V = Vsave_twin;
                end
            end
            if self.parameters.saveH
                if self.parameters.notwin
                    spectra.H = cat(3, Hsave{:});
                else
                    spectra.H = cat(3, Hsave_twin{:});
                end
            end
            if self.parameters.saveSabp
                spectra.Sabp = cat(4, Sabp_twin{:});
                spectra.omegap = cat(2, omega_twin{:});
            end
            % save the important parameters
            spectra.param.sortMode  = self.parameters.sortMode;
            spectra.param.tol       = self.parameters.tol;
            spectra.param.omega_tol = self.parameters.omega_tol;
            spectra.param.hermit    = self.parameters.hermit;
            spectra.title           = self.parameters.title;
            spectra.gtensor         = self.parameters.gtensor;
            if ~self.parameters.fitmode
                spectra.dateend = datestr(now);
                spectra.obj = copy(self.spinWaveObject);
            end
        end
    end


    methods(Access=private)
        function nTwin = checkTwins(self)
            % checks for twins
            if self.parameters.notwin
                nTwin = 1;
            else
                nTwin = size(self.spinWaveObject.twin.vol,2);
                % if there is only one twin and it is the identity set param.notwin true
                rotc1 = self.spinWaveObject.twin.rotc(:,:,1) - eye(3);
                if (nTwin == 1) && norm(rotc1(:))==0
                    self.parameters.notwin = true;
                end
            end
        end

        function spinwave_parse_input(self, varargin)
            pref = swpref();
            title0 = 'Numerical LSWT spectrum';

            inpForm.fname  = {'fitmode' 'notwin' 'sortMode' 'optmem' 'tol' 'hermit'};
            inpForm.defval = {false     false    true       0        1e-4  true    };
            inpForm.size   = {[1 1]     [1 1]    [1 1]      [1 1]    [1 1] [1 1]   };

            inpForm.fname  = [inpForm.fname  {'omega_tol' 'saveSabp' 'saveV' 'saveH'}];
            inpForm.defval = [inpForm.defval {1e-5        false      false   false  }];
            inpForm.size   = [inpForm.size   {[1 1]       [1 1]      [1 1]   [1 1]  }];

            inpForm.fname  = [inpForm.fname  {'formfact' 'formfactfun' 'title' 'gtensor'}];
            inpForm.defval = [inpForm.defval {false       @sw_mff      title0  false    }];
            inpForm.size   = [inpForm.size   {[1 -1]      [1 1]        [1 -2]  [1 1]    }];

            inpForm.fname  = [inpForm.fname  {'cmplxBase' 'tid' 'fid' }];
            inpForm.defval = [inpForm.defval {false       -1    -1    }];
            inpForm.size   = [inpForm.size   {[1 1]       [1 1] [1 1] }];

            self.parameters = sw_readparam(inpForm, varargin{:});
            %if ~self.parameters.fitmode
            %    % save the time of the beginning of the calculation
            %    self.spectra.datestart = datestr(now);
            %end
            if self.parameters.fitmode
                self.parameters.sortMode = false;
                self.parameters.tid = 0;
            end
            if self.parameters.tid == -1
                self.parameters.tid = pref.tid;
            end
            if self.parameters.fid == -1
                self.parameters.fid = pref.fid;
            end
            pref = swpref;
            % use mex file by default?
            self.parameters.useMex = pref.usemex;
        end

        function parseMagneticStructure(self)
            % generate magnetic structure in the rotating notation
            self.magnetic_structure = sw_classes.magnetic_structure(self.spinWaveObject.magstr(), self.parameters.tol);
            if self.parameters.cmplxBase
                % The coordinate system is fixed by the complex magnetisation vectors:
                % e1 = imag(M), e3 = real(M), e2 = cross(e3,e1)
                % TODO: put this into the magnetic_structure class when it is refactored.
                F0  = self.spinWaveObject.mag_str.F;
                RF0 = vecnorm(real(F0));
                IF0 = vecnorm(imag(F0));
                % e3 = real(M)
                e3  = real(F0)./repmat(RF0,[3 1]);
                % e1 = imag(M) perpendicular to e3
                e1  = imag(F0)./repmat(IF0,[3 1]);
                e1 = e1 - (sum(e1.*e3).*e3);
                e1  = e1./vecnorm(e1);
                % e2 = cross(e3,e1)
                e2  = cross(e3,e1);
            else
                [e1, e2, e3] = self.magnetic_structure.get_local_basis();
            end
            % assign complex vectors that define the rotating coordinate system on
            % every magnetic atom
            self.zed = e1 + 1i*e2;
            self.eta = e3;
        end

        function nChunk = computeNumChunk(self)
            nMagExt = self.magnetic_structure.nMagExt;
            nHkl = self.qvectors.nHkl;
            if self.parameters.optmem == 0
                freeMem = sw_freemem;
                if freeMem > 0
                    nChunk = ceil(nMagExt^2*nHkl*6912/freeMem*2);
                else
                    nChunk = 1;
                    if ~self.parameters.fitmode
                        warning('spinw:spinwave:FreeMemSize','The size of the free memory is unkown, no memory optimisation!');
                    end
                end
            else
                nChunk = self.parameters.optmem;
            end
            if nHkl < nChunk
                fprintf0(self.parameters.fid,['Memory allocation is not optimal, nMagExt is'...
                    ' too large compared to the free memory!\n']);
                nChunk = nHkl;
            elseif nChunk > 1
                fprintf0(self.parameters.fid,['To optimise memory allocation, Q is cut'...
                    ' into %d pieces!\n'],nChunk);
            end
        end

        function prepareHamiltonian(self)
            % Create the interaction matrix and atomic positions in the extended
            % magnetic unit cell.
            nMagExt = self.magnetic_structure.nMagExt;
            S_mag = self.magnetic_structure.S_mag';
            [SS, self.SI, self.RR] = self.spinWaveObject.intmatrix('fitmode',true,'conjugate',true);

            % add the dipolar interactions to SS.all
            SS.all = [SS.all SS.dip];

            bq = SS.all(15,:)==1;  % SS is an array of Js. Each row corresponds to each unique property.
            if (any(bq))
                if (self.magnetic_structure.incomm)
                    error('spinw:spinwave:Biquadratic','Biquadratic exchange can be only calculated for k=0 structures!');
                end
                bqdR    = SS.all(1:3, bq);
                bqAtom1 = SS.all(4, bq);
                bqAtom2 = SS.all(5, bq);
                bqJJ    = SS.all(6, bq);
                self.bq_hamiltonian = sw_classes.biquadratic_hamiltonian(bqJJ, bqdR, bqAtom1, bqAtom2, self.zed, self.eta, S_mag);
                % Remove the biquadratic interactions for subsequent calculations
                SS.all = SS.all(1:14,SS.all(15,:)==0);
            end

            dR    = [SS.all(1:3,:) zeros(3,nMagExt)];
            atom1 = transpose([SS.all(4,:)   1:nMagExt]);
            atom2 = transpose([SS.all(5,:)   1:nMagExt]);

            % magnetic couplings, 3x3xnJ; Jij in eq (26) of Toth & Lake.
            JJ = cat(3, reshape(SS.all(6:14,:),3,3,[]), self.SI.aniso);

            if self.magnetic_structure.incomm
                % transform JJ due to the incommensurate wavevector
                [~, K] = sw_rot(self.magnetic_structure.n, self.magnetic_structure.km*dR*2*pi);
                % multiply JJ with K matrices for every interaction
                % and symmetrising JJ for the rotating basis
                JJ = (mmat(JJ,K)+mmat(K,JJ))/2;
            end

            if self.parameters.gtensor
                self.gtensor = self.SI.g;
                if self.magnetic_structure.incomm
                    self.gtensor = self.magnetic_structure.transform_frame(self.gtensor);
                end
            end

            % Calculates the hkl-independent part of the Hamiltonian in a separate class
            self.hamiltonian = sw_classes.hamiltonian(JJ, dR, atom1, atom2, self.zed, self.eta, S_mag);
        end

        function ham = calculateHamiltonian(self, hkl, rotC, k_incomm)
            nExt = self.magnetic_structure.N_ext;

            % calculate all magnetic form factors
            if self.parameters.formfact
                obj = self.spinWaveObject;
                % Angstrom^-1 units for Q
                hklA0 = 2*pi*(hkl'/obj.basisvector)';
                % store form factor per Q point for each atom in the magnetic supercell
                self.parameters.FF = repmat(self.parameters.formfactfun(permute(obj.unit_cell.ff(1,:,obj.matom.idx),[3 2 1]),hklA0),[prod(nExt) 1]);
            end

            % Converts wavevector list into the extended unit cell
            hklExt  = 2*pi*(hkl.*nExt' + k_incomm);

            % Calculates the hkl-dependent Hamiltonian
            ham = self.hamiltonian.ham_k(hklExt);

            if ~isempty(self.bq_hamiltonian)
                ham = ham + self.bq_hamiltonian.ham_k(hklExt);
            end
            if any(self.SI.field)
                % Calculates the contribution of the magnetic field (Zeeman term) to the Hamiltonian
                field = self.SI.field*rotC*self.spinWaveObject.unit.muB;
                MF = repmat(field * permute(mmat(self.SI.g,permute(self.eta,[1 3 2])),[1 3 2]),[1 2]);
                ham = ham + diag(MF);
            end

            ham = (ham + conj(permute(ham,[2 1 3])))/2;
        end

        function [V, omega] = diagonaliseHamiltonian(self, ham)
            nMagExt = self.magnetic_structure.nMagExt;
            if self.parameters.hermit
                % All the matrix calculations are according to Colpa's paper
                % J.H.P. Colpa, Physica 93A (1978) 327-353
                [V, omega, self.posdefWarn] = spinwave_hermit(ham, self.parameters, nMagExt);

            else
                % All the matrix calculations are according to White's paper
                % R.M. White, et al., Physical Review 139, A450?A454 (1965)

                % diagonal of the boson commutator matrix
                gCommd = [ones(nMagExt,1); -ones(nMagExt,1)];
                % boson commutator matrix
                gComm  = diag(gCommd);
                %gd = diag(g);

                gham = mmat(gComm,ham);

                [V, omega, orthWarn1] = eigorth(gham, self.parameters.omega_tol, self.parameters.useMex);

                self.orthWarn = orthWarn1 || self.orthWarn;

                for ii = 1:size(omega, 2)
                    % multiplication with g removed to get negative and positive
                    % energies as well
                    %omega(:,end+1) = D(:,ii); %#ok<AGROW>
                    M              = diag(gComm*V(:,:,ii)'*gComm*V(:,:,ii));
                    V(:,:,ii)      = V(:,:,ii)*diag(sqrt(1./M));
                end
            end
        end

        function Sab = spinspincorrel(self, V, hkl)
            nMagExt = self.magnetic_structure.nMagExt;
            S0 = self.magnetic_structure.S_mag;
            nExt = self.magnetic_structure.N_ext;
            hklExt0MEM = 2*pi*hkl.*nExt';
            nHklMEM = size(V,3);

            % Calculates correlation functions.
            % V right
            VExtR = repmat(permute(V  ,[4 5 1 2 3]),[3 3 1 1 1]);
            % V left: conjugate transpose of V
            VExtL = conj(permute(VExtR,[1 2 4 3 5]));

            % Introduces the exp(-ikR) exponential factor.
            ExpF =  exp(-1i*sum(repmat(permute(hklExt0MEM,[1 3 2]),[1 nMagExt 1]).*repmat(self.RR,[1 1 nHklMEM]),1));
            % Includes the sqrt(Si/2) prefactor.
            ExpF = ExpF.*repmat(sqrt(S0/2),[1 1 nHklMEM]);

            ExpFL = repmat(permute(ExpF,[1 4 5 2 3]),[3 3 2*nMagExt 2]);
            % conj transpose of ExpFL
            ExpFR = conj(permute(ExpFL,[1 2 4 3 5]));

            zeda = repmat(permute([self.zed conj(self.zed)],[1 3 4 2]),[1 3 2*nMagExt 1 nHklMEM]);
            % conj transpose of zeda
            zedb = conj(permute(zeda,[2 1 4 3 5]));

            % calculate magnetic structure factor using the hklExt0 Q-values
            % since the S(Q+/-k,omega) correlation functions also belong to the
            % F(Q)^2 form factor
            if self.parameters.formfact
                % include the form factor in the z^alpha, z^beta matrices
                zeda = zeda.*repmat(permute(self.parameters.FF,[3 4 5 1 2]),[3 3 2*nMagExt 2 1]);
                zedb = zedb.*repmat(permute(self.parameters.FF,[3 4 1 5 2]),[3 3 2 2*nMagExt 1]);
            end

            if self.parameters.gtensor
                % include the g-tensor
                zeda = mmat(repmat(permute(self.gtensor,[1 2 4 3]),[1 1 1 2]),zeda);
                zedb = mmat(zedb,repmat(self.gtensor,[1 1 2]));
            end
            % Dynamical structure factor from S^alpha^beta(k) correlation function.
            % Sab(alpha,beta,iMode,iHkl), size: 3 x 3 x 2*nMagExt x nHkl.
            % Normalizes the intensity to single unit cell.
            %Sab = cat(4,Sab,squeeze(sum(zeda.*ExpFL.*VExtL,4)).*squeeze(sum(zedb.*ExpFR.*VExtR,3))/prod(nExt));
            Sab = squeeze(sum(zeda.*ExpFL.*VExtL,4)) .* squeeze(sum(zedb.*ExpFR.*VExtR,3)) / prod(nExt);

            if self.magnetic_structure.incomm && self.magnetic_structure.helical
                % integrating out the arbitrary initial phase of the helix
                Sab = self.magnetic_structure.transform_frame(Sab);
            end
        end

    end
end


function [V, omega, warn1] = spinwave_hermit(ham, param, nMagExt)
    warn1 = false;
    nHklMEM = size(ham, 3);

    % diagonal of the boson commutator matrix
    gCommd = [ones(nMagExt,1); -ones(nMagExt,1)];
    % boson commutator matrix
    gComm  = diag(gCommd);
    %gd = diag(g);

    % basis functions of the magnon modes
    V = zeros(2*nMagExt,2*nMagExt,nHklMEM);
    omega = zeros(2*nMagExt,0);

    if param.useMex && nHklMEM>1
        % use mex files to speed up the calculation
        % mex file will return an error if the matrix is not positive definite.
        [K2, invK] = chol_omp(ham,'Colpa','tol',param.omega_tol);
        [V, omega] = eig_omp(K2,'sort','descend');
        % the inverse of the para-unitary transformation V
        for ii = 1:nHklMEM
            V(:,:,ii) = V(:,:,ii)*diag(sqrt(gCommd.*omega(:,ii)));
        end
        V = sw_mtimesx(invK,V);
    else
        for ii = 1:nHklMEM
            [K, posDef]  = chol(ham(:,:,ii));
            if posDef > 0
                try
                    % get tolerance from smallest negative eigenvalue
                    tol0 = eig(ham(:,:,ii));
                    tol0 = sort(real(tol0));
                    tol0 = abs(tol0(1));
                    % TODO determine the right tolerance value
                    tol0 = tol0*sqrt(nMagExt*2)*4;
                    if tol0>param.omega_tol
                        error('spinw:spinwave:NonPosDefHamiltonian','Very baaaad!');
                    end
                    try
                        K = chol(ham(:,:,ii)+eye(2*nMagExt)*tol0);
                    catch
                        K = chol(ham(:,:,ii)+eye(2*nMagExt)*param.omega_tol);
                    end
                    warn1 = true;
                catch PD
                    if param.tid == 2
                        % close timer window
                        sw_timeit(100,2,param.tid);
                    end
                    error('spinw:spinwave:NonPosDefHamiltonian',...
                        ['Hamiltonian matrix is not positive definite, probably'...
                        ' the magnetic structure is wrong! For approximate'...
                        ' diagonalization try the param.hermit=false option']);
                end
            end

            K2 = K*gComm*K';
            K2 = 1/2*(K2+K2');
            % Hermitian K2 will give orthogonal eigenvectors
            [U, D] = eig(K2);
            D      = diag(D);

            % sort modes accordign to the real part of the energy
            [~, idx] = sort(real(D),'descend');
            U = U(:,idx);
            % omega dispersion
            omega(:,end+1) = D(idx); %#ok<AGROW>

            % the inverse of the para-unitary transformation V
            V(:,:,ii) = inv(K)*U*diag(sqrt(gCommd.*omega(:,end))); %#ok<MINV>
        end
    end
end
