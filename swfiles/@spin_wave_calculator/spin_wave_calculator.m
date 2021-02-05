classdef spin_wave_calculator
    properties
        spinWaveObject
        parameters
        spectra
        magnonEnergies
        hamiltonian  % this could be its own class which contains its own calculation methon for the eigenvectors and eigenvalues (magnonEnergies)
    end

    % properties (Constant)
    %     numberOfSlices
    % end

    methods(Access=public)
        function self = SpinWaveCalculator(spinWaveObject, varargin)
            self.spinWaveObject = spinWaveObject;
            self.setup(varargin{:}) % varargin{:} is a MATLAB cell array
            self.parseMagneticStructure
            self.prepareHamiltonian
        end

        function calculateSpinWave(self, hkl)
            self.prepareHamiltonian_hkl(hkl)
            self.doCalculateSpinWave
            self.prepateOutputs
        end
    end

    methods(Access=private)
        function parseMagneticStructure(self)
            % generate magnetic structure in the rotating notation
            magStr = self.spinWaveObject.magstr;

            % Calculates parameters eta and zed.
            if isempty(magStr.S)
                error('spinw:spinwave:NoMagneticStr','No magnetic structure defined in self.spinWaveObject!');
            end

            Si = magStr.S;
            S_mag = vecnorm(Si);
            % normal to rotation of the magnetic moments
            n  = magStr.n;
            nMagExt = size(Si,2);


            % size of the extended magnetic unit cell
            nExt = magStr.N_ext;
            % magnetic ordering wavevector in the extended magnetic unit cell
            km = magStr.k.*nExt;

            % whether the structure is incommensurate
            self.incomm = any(abs(km-round(km)) > param.tol);

            % Check for 2*km
            tol = self.parameters.tol*2;
            self.helical =  sum(abs(mod(abs(2*km)+tol,1)-tol).^2) > tol;

            if self.parameters.cmplxBase
                % The coordinate system is fixed by the complex magnetisation vectors:
                % e1 = imag(M), e3 = real(M), e2 = cross(e3,e1)
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
                % Local (e1,e2,e3) coordinate system fixed to the moments,
                % e3||Si,ata
                % e2 = Si x [1,0,0], if Si || [1,0,0] --> e2 = [0,0,1]
                % e1 = e2 x e3
                %
                % e3 || Si
                e3 = Si./S_mag;
                % e2 = Si x [1,0,0], if Si || [1,0,0] --> e2 = [0,0,1]
                e2  = [zeros(1,nMagExt); e3(3,:); -e3(2,:)];
                e2(3,~any(abs(e2)>1e-10)) = 1;
                e2  = e2./vecnorm(e2);
                % e1 = e2 x e3
                e1  = cross(e2,e3);
            end

            % assign complex vectors that define the rotating coordinate system on
            % every magnetic atom
            self.zed = e1 + 1i*e2;
            self.eta = e3;
            % self.nMagExt = nMagExt;
        end

        function prepareHamiltonian_hkl(self, hkl)

            % rotC = eye(3);

            % Transform the momentum values to the new lattice coordinate system
            hkl = self.spinWaveObject.unit.qmat*hkl;

            % Calculates momentum transfer in A^-1 units.
            hklA = 2*pi*(hkl'/self.spinWaveObject.basisvector)';

            % number of Q points
            nHkl0 = size(hkl,2);

            nTwin = 1;

            hkl0   = hkl;
            nHkl   = nHkl0;

            % Converts wavevector list into the extended unit cell
            % hklExt  = bsxfun(@times,hklExt,nExt')*2*pi;
            hklExt  = 2*pi*hkl.*nExt';
            % q values without the +/-k_m value
            hklExt0 = hklExt;
        end

        function nSlice = computeNumSlice(self)
            nMagExt = self.spinWaveObject.mag_str.nExt;
            if param.optmem == 0
                freeMem = sw_freemem;
                if freeMem > 0
                    nSlice = ceil(nMagExt^2*nHkl*6912/freeMem*2);
                else
                    nSlice = 1;
                    if ~param.fitmode
                        warning('spinw:spinwave:FreeMemSize','The size of the free memory is unkown, no memory optimisation!');
                    end
                end
            else
                nSlice = param.optmem;
            end

            if nHkl < nSlice
                fprintf0(fid,['Memory allocation is not optimal, nMagExt is'...
                    ' too large compared to the free memory!\n']);
                nSlice = nHkl;
            elseif nSlice > 1
                fprintf0(fid,['To optimise memory allocation, Q is cut'...
                    ' into %d pieces!\n'],nSlice);
            end

        end


        function prepareHamiltonian(self)

            % Create the interaction matrix and atomic positions in the extended
            % magnetic unit cell.
            [SS, SI, RR] = self.spinWaveObject.intmatrix('fitmode',true,'conjugate',true);

            % add the dipolar interactions to SS.all
            SS.all = [SS.all SS.dip];

            fprintf0(fid,['Calculating COMMENSURATE spin wave spectra '...
                '(nMagExt = %d, nHkl = %d, nTwin = %d)...\n'],nMagExt, nHkl0, nTwin);

            dR    = [SS.all(1:3,:) zeros(3,nMagExt)];
            atom1 = [SS.all(4,:)   1:nMagExt];
            atom2 = [SS.all(5,:)   1:nMagExt];
            % magnetic couplings, 3x3xnJ
            JJ = cat(3,reshape(SS.all(6:14,:),3,3,[]),SI.aniso);

            nCoupling = size(JJ,3);

            zedL = repmat(permute(zed(:,atom1),[1 3 2]),[1 3 1]);
            zedR = repmat(permute(zed(:,atom2),[3 1 2]),[3 1 1]);

            etaL = repmat(permute(eta(:,atom1),[1 3 2]),[1 3 1]);
            etaR = repmat(permute(eta(:,atom2),[3 1 2]),[3 1 1]);

            SiSj = sqrt(S_mag(atom1).*S_mag(atom2));

            % Creates temporary values for calculating matrix elements.
            AD  =  shiftdim(sum(sum(etaL.*JJ.*etaR,2),1),1);
            A20 = -S_mag(atom2).*AD;
            D20 = -S_mag(atom1).*AD;
            BC0 =  SiSj.*shiftdim(sum(sum(zedL.*JJ.*     zedR ,2),1),1);
            AD0 =  SiSj.*shiftdim(sum(sum(zedL.*JJ.*conj(zedR),2),1),1);

            % Calculates the contribution of the magnetic field (Zeeman term) to the Hamiltonian
            % MF = repmat(SI.field*rotC*self.spinWaveObject.unit.muB * permute(mmat(SI.g,permute(eta,[1 3 2])),[1 3 2]),[1 2]);

            % Creates the serial indices for every matrix element in ham matrix.
            idxA1 = [atom1'         atom2'         ];
            idxA2 = [atom1'         atom1'         ];
            idxB  = [atom1'         atom2'+nMagExt ];
            % transpose of idxB
            %idxC  = [atom2'+nMagExt atom1'         ]; % SP1
            idxD1 = idxA1+nMagExt;
            idxD2 = [atom2'+nMagExt atom2'+nMagExt ];
            idxMF = [(1:2*nMagExt)' (1:2*nMagExt)' ];

        end

        function spinwave_parse_input(self, varargin)

            pref = swpref;

            % use mex file by default?
            self.useMex = pref.usemex;

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
        end

        function incommensurateCase(self)
        end

        function twinCase(self)
        end

        function doCalculateSpinWave(self)
        end

        function prepareOutputs(self)
        end

    end
end
