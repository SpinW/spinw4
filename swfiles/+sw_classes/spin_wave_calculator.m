classdef spin_wave_calculator < handle
    properties
        spinWaveObject
        parameters
        spectra
        magnonEnergies
        magnetic_structure
        qvectors
        hamiltonian  % this could be its own class which contains its own calculation methon for the eigenvectors and eigenvalues (magnonEnergies)
    end

    properties(Access=private)
        eta
        zed
        hamVars
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
            % Gets the transformation matrices for the hkl for each twin
            [~, rotQ] = self.spinWaveObject.twinq([0;0;0]);
            for iTwin = 1:nTwins
                rotC = self.spinWaveObject.twin.rotc(:,:,iTwin);
                omega = [];
                Sab = cell(1, nChunk);
                Hsave = [];
                Vsave = [];
                for ii = 1:nChunk
                    hklIdxMEM = self.qvectors.getIdx(ii);
                    hklChunk = self.qvectors.getChunk(ii);
                    % If twin, then rotate the current hkl set by its rotation matrix
                    if ~self.parameters.notwin
                        hklChunk = (hklChunk' * rotQ(:,:,iTwin))';
                    end
                    [ham, hklExt] = self.calculateHamiltonian(hklChunk, rotC);
                    if self.parameters.saveH
                        Hsave(:,:,hklIdxMEM) = ham; %#ok<AGROW>
                    end
                    [V, omega(:,hklIdxMEM)] = self.diagonaliseHamiltonian(ham);  %#ok<AGROW>
                    if self.parameters.saveV
                        Vsave(:,:,hklIdxMEM) = V; %#ok<AGROW>
                    end
                    Sab{ii} = self.spinspincorrel(V, hklExt);
                end
                Sab = cat(4, Sab{:});
                if self.parameters.sortMode
                    % sort the spin wave modes
                    [omega, Sab] = sortmode(omega,reshape(Sab,9,size(Sab,3),[]));
                    Sab = reshape(Sab,3,3,size(Sab,2),[]);
                end
                if ~self.parameters.notwin
                    omega_twin{iTwin} = omega;
                    Sab_twin{iTwin} = Sab;
                    if self.parameters.saveV
                        Vsave_twin{iTwin} = Vsave;
                    end
                    if self.parameters.saveH
                        Hsave_twin{iTwin} = Hsave;
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
                    spectra.V = Vsave;
                else
                    spectra.V = Vsave_twin;
                end
            end
            if self.parameters.saveH
                if self.parameters.notwin
                    spectra.H = Hsave;
                else
                    spectra.H = cat(3, Hsave_twin{:});
                end
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
            if ~self.parameters.fitmode
                % save the time of the beginning of the calculation
                self.spectra.datestart = datestr(now);
            end
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
            S_mag = self.magnetic_structure.S_mag;
            [SS, SI, RR] = self.spinWaveObject.intmatrix('fitmode',true,'conjugate',true);

            % add the dipolar interactions to SS.all
            SS.all = [SS.all SS.dip];

            %fprintf0(fid,['Calculating COMMENSURATE spin wave spectra '...
            %    '(nMagExt = %d, nHkl = %d, nTwin = %d)...\n'],nMagExt, nHkl0, nTwin);

            dR    = [SS.all(1:3,:) zeros(3,nMagExt)];
            atom1 = [SS.all(4,:)   1:nMagExt];
            atom2 = [SS.all(5,:)   1:nMagExt];
            % magnetic couplings, 3x3xnJ
            JJ = cat(3,reshape(SS.all(6:14,:),3,3,[]),SI.aniso);

            nCoupling = size(JJ,3);

            zedL = repmat(permute(self.zed(:,atom1),[1 3 2]),[1 3 1]);
            zedR = repmat(permute(self.zed(:,atom2),[3 1 2]),[3 1 1]);

            etaL = repmat(permute(self.eta(:,atom1),[1 3 2]),[1 3 1]);
            etaR = repmat(permute(self.eta(:,atom2),[3 1 2]),[3 1 1]);

            SiSj = sqrt(S_mag(atom1).*S_mag(atom2));

            % Creates temporary values for calculating matrix elements.
            AD  =  shiftdim(sum(sum(etaL.*JJ.*etaR,2),1),1);
            A20 = -S_mag(atom2).*AD;
            D20 = -S_mag(atom1).*AD;
            BC0 =  SiSj.*shiftdim(sum(sum(zedL.*JJ.*     zedR ,2),1),1);
            AD0 =  SiSj.*shiftdim(sum(sum(zedL.*JJ.*conj(zedR),2),1),1);

            % Creates the serial indices for every matrix element in ham matrix.
            idxA1 = [atom1'         atom2'         ];
            idxA2 = [atom1'         atom1'         ];
            idxB  = [atom1'         atom2'+nMagExt ];
            % transpose of idxB
            %idxC  = [atom2'+nMagExt atom1'         ]; % SP1
            idxD1 = idxA1+nMagExt;
            idxD2 = [atom2'+nMagExt atom2'+nMagExt ];
            idxMF = [(1:2*nMagExt)' (1:2*nMagExt)' ];

            % Variables required for other parts of the calculations
            self.hamVars = struct('SI', SI, 'RR', RR, 'dR', dR, ...
                'AD0', AD0, 'BC0', BC0, 'A20', A20, 'D20', D20, ...
                'idxA1', idxA1, 'idxB', idxB, 'idxD1', idxD1, ...
                'idxA2', idxA2, 'idxD2', idxD2, 'idxMF', idxMF, ...
                'nCoupling', nCoupling);
        end

        function [ham, hklExt] = calculateHamiltonian(self, hkl, rotC)
            nExt = self.magnetic_structure.N_ext;
            nMagExt = self.magnetic_structure.nMagExt;
            % Copies the variables in hamVars to this workspace
            assign_vars(self.hamVars);

            % Transform the momentum values to the new lattice coordinate system
            hkl = self.spinWaveObject.unit.qmat*hkl;

            % calculate all magnetic form factors
            if self.parameters.formfact
                obj = self.spinWaveObject;
                % Angstrom^-1 units for Q
                hklA0 = 2*pi*(hkl'/obj.basisvector)';
                % store form factor per Q point for each atom in the magnetic supercell
                % TODO check prod(nExt)? instead of nExt
                %FF = repmat(param.formfactfun(permute(obj.unit_cell.ff(1,:,obj.matom.idx),[3 2 1]),hklA0),[1 nExt]);
                self.parameters.FF = repmat(self.parameters.formfactfun(permute(obj.unit_cell.ff(1,:,obj.matom.idx),[3 2 1]),hklA0),[prod(nExt) 1]);
            end

            % Calculates momentum transfer in A^-1 units.
            %hklA = 2*pi*(hkl'/self.spinWaveObject.basisvector)';

            % number of Q points
            nHkl0 = size(hkl,2);

            %hkl0   = hkl;
            nHkl   = nHkl0;

            % Converts wavevector list into the extended unit cell
            % hklExt  = bsxfun(@times,hklExt,nExt')*2*pi;
            hklExt  = 2*pi*hkl.*nExt';
            % q values without the +/-k_m value
            %hklExt0 = hklExt;

            % Creates the matrix of exponential factors nCoupling x nHkl size.
            % Extends dR into 3 x 3 x nCoupling x nHkl
            %     ExpF = exp(1i*permute(sum(repmat(dR,[1 1 nHklMEM]).*repmat(...
            %         permute(hklExtMEM,[1 3 2]),[1 nCoupling 1]),1),[2 3 1]))';
            ExpF = exp(1i*permute(sum(bsxfun(@times,dR,permute(hklExt,[1 3 2])),1),[2 3 1]))';

            % Creates the matrix elements containing zed.
            A1 = bsxfun(@times,     AD0 ,ExpF);
            B  = bsxfun(@times,     BC0 ,ExpF);
            D1 = bsxfun(@times,conj(AD0),ExpF);

            % Store all indices
            % SP1: speedup for creating the matrix elements
            %idxAll = [idxA1; idxB; idxC; idxD1]; % SP1
            idxAll   = [idxA1; idxB; idxD1];
            % Store all matrix elements
            %ABCD   = [A1     B     conj(B)  D1]; % SP1
            ABCD   = [A1     2*B      D1];

            % Stores the matrix elements in ham.
            %idx3   = repmat(1:nHklMEM,[4*nCoupling 1]); % SP1
            idx3   = repmat(1:nHkl,[3*nCoupling 1]);
            idxAll = [repmat(idxAll,[nHkl 1]) idx3(:)];
            idxAll = idxAll(:,[2 1 3]);

            ABCD   = ABCD';

            % quadratic form of the boson Hamiltonian stored as a square matrix
            ham = accumarray(idxAll,ABCD(:),[2*nMagExt 2*nMagExt nHkl]);

            ham = ham + repmat(accumarray([idxA2; idxD2],2*[A20 D20],[1 1]*2*nMagExt),[1 1 nHkl]);

            if any(SI.field)
                % Calculates the contribution of the magnetic field (Zeeman term) to the Hamiltonian
                MF = repmat(SI.field*rotC*self.spinWaveObject.unit.muB * permute(mmat(SI.g,permute(self.eta,[1 3 2])),[1 3 2]),[1 2]);
                ham = ham + repmat(accumarray(idxMF,MF,[1 1]*2*nMagExt),[1 1 nHkl]);
            end

            ham = (ham + conj(permute(ham,[2 1 3])))/2;
        end

        function [V, omega] = diagonaliseHamiltonian(self, ham)
            nMagExt = self.magnetic_structure.nMagExt;
            if self.parameters.hermit
                % All the matrix calculations are according to Colpa's paper
                % J.H.P. Colpa, Physica 93A (1978) 327-353
                [V, omega] = spinwave_hermit(ham, self.parameters, nMagExt);

            else
                % All the matrix calculations are according to White's paper
                % R.M. White, et al., Physical Review 139, A450?A454 (1965)

                % diagonal of the boson commutator matrix
                gCommd = [ones(nMagExt,1); -ones(nMagExt,1)];
                % boson commutator matrix
                gComm  = diag(gCommd);
                %gd = diag(g);

                gham = mmat(gComm,ham);

                [V, omega, orthWarn] = eigorth(gham, self.parameters.omega_tol, self.parameters.useMex);

                %orthWarn0 = orthWarn || orthWarn0;

                for ii = 1:size(omega, 2)
                    % multiplication with g removed to get negative and positive
                    % energies as well
                    %omega(:,end+1) = D(:,ii); %#ok<AGROW>
                    M              = diag(gComm*V(:,:,ii)'*gComm*V(:,:,ii));
                    V(:,:,ii)      = V(:,:,ii)*diag(sqrt(1./M));
                end
            end
        end

        function Sab = spinspincorrel(self, V, hklExt0MEM)
            nMagExt = self.magnetic_structure.nMagExt;
            S0 = self.magnetic_structure.S_mag;
            nExt = self.magnetic_structure.N_ext;
            nHklMEM = size(V,3);
            % Copies the variables in hamVars to this workspace
            assign_vars(self.hamVars);

            % Calculates correlation functions.
            % V right
            VExtR = repmat(permute(V  ,[4 5 1 2 3]),[3 3 1 1 1]);
            % V left: conjugate transpose of V
            VExtL = conj(permute(VExtR,[1 2 4 3 5]));

            % Introduces the exp(-ikR) exponential factor.
            ExpF =  exp(-1i*sum(repmat(permute(hklExt0MEM,[1 3 2]),[1 nMagExt 1]).*repmat(RR,[1 1 nHklMEM]),1));
            % Includes the sqrt(Si/2) prefactor.
            ExpF = ExpF.*repmat(sqrt(S0/2),[1 1 nHklMEM]);

            ExpFL =      repmat(permute(ExpF,[1 4 5 2 3]),[3 3 2*nMagExt 2]);
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
                gtensor = SI.g;
                % include the g-tensor
                zeda = mmat(repmat(permute(gtensor,[1 2 4 3]),[1 1 1 2]),zeda);
                zedb = mmat(zedb,repmat(gtensor,[1 1 2]));
            end
            % Dynamical structure factor from S^alpha^beta(k) correlation function.
            % Sab(alpha,beta,iMode,iHkl), size: 3 x 3 x 2*nMagExt x nHkl.
            % Normalizes the intensity to single unit cell.
            %Sab = cat(4,Sab,squeeze(sum(zeda.*ExpFL.*VExtL,4)).*squeeze(sum(zedb.*ExpFR.*VExtR,3))/prod(nExt));
            Sab = squeeze(sum(zeda.*ExpFL.*VExtL,4)) .* squeeze(sum(zedb.*ExpFR.*VExtR,3)) / prod(nExt);
        end

        function incommensurateCase(self)
            error('Not implemented');
        end

        function twinCase(self)
            error('Not implemented');
        end

    end
end


function [V, omega] = spinwave_hermit(ham, param, nMagExt)
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
        % V = bsxfun(@times,invK,V);
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

function assign_vars(vars)
    % Loads save variables to workspace
    vnames = fieldnames(vars);
    for ii = 1:numel(vnames)
        assignin('caller', vnames{ii}, vars.(vnames{ii}));
    end
end
