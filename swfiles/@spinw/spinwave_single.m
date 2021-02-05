function spectra = spinwave_single(obj, hkl, param, rotC)
% calculates spin correlation function using linear spin wave theory

fid = param.fid;
pref = swpref;
useMex = pref.usemex;

if nargin < 4
    rotC = eye(3);
end

% save warning of eigorth
orthWarn0 = false;

% save warning for singular matrix
singWarn0 = warning('off','MATLAB:nearlySingularMatrix');

% generate magnetic structure in the rotating notation
magStr = obj.magstr;

% size of the extended magnetic unit cell
nExt = magStr.N_ext;
% magnetic ordering wavevector in the extended magnetic unit cell
km = magStr.k.*nExt;

% Transform the momentum values to the new lattice coordinate system
hkl = obj.unit.qmat*hkl;

% Calculates momentum transfer in A^-1 units.
hklA = 2*pi*(hkl'/obj.basisvector)';

% Check for 2*km
tol = param.tol*2;
helical =  sum(abs(mod(abs(2*km)+tol,1)-tol).^2) > tol;

% number of Q points
nHkl0 = size(hkl,2);

nTwin = 1;

hkl0   = hkl;
hklExt = hkl;
nHkl   = nHkl0;

% Create the interaction matrix and atomic positions in the extended
% magnetic unit cell.
[SS, SI, RR] = obj.intmatrix('fitmode',true,'conjugate',true);

% add the dipolar interactions to SS.all
SS.all = [SS.all SS.dip];

% Converts wavevctor list into the extended unit cell
hklExt  = bsxfun(@times,hklExt,nExt')*2*pi;
% q values without the +/-k_m value
hklExt0 = bsxfun(@times,hkl0,nExt')*2*pi;

% Calculates parameters eta and zed.
if isempty(magStr.S)
    error('spinw:spinwave:NoMagneticStr','No magnetic structure defined in obj!');
end

M0 = magStr.S;
S0 = sqrt(sum(M0.^2,1));
% normal to rotation of the magnetic moments
n  = magStr.n;
nMagExt = size(M0,2);

fprintf0(fid,['Calculating COMMENSURATE spin wave spectra '...
        '(nMagExt = %d, nHkl = %d, nTwin = %d)...\n'],nMagExt, nHkl0, nTwin);

if param.cmplxBase
    % The coordinate system is fixed by the complex magnetisation vectors:
    % e1 = imag(M), e3 = real(M), e2 = cross(e3,e1)
    F0  = obj.mag_str.F;
    RF0 = sqrt(sum(real(F0).^2,1));
    IF0 = sqrt(sum(imag(F0).^2,1));
    % e3 = real(M)
    e3  = real(F0)./repmat(RF0,[3 1]);
    % e1 = imag(M) perpendicular to e3
    e1  = imag(F0)./repmat(IF0,[3 1]);
    e1  = e1-bsxfun(@times,sum(e1.*e3,1),e3);
    e1  = e1./repmat(sqrt(sum(e1.^2,1)),[3 1]);
    % e2 = cross(e3,e1)
    e2  = cross(e3,e1);
else
    % Local (e1,e2,e3) coordinate system fixed to the moments,
    % e3||Si,ata
    % e2 = Si x [1,0,0], if Si || [1,0,0] --> e2 = [0,0,1]
    % e1 = e2 x e3
    %
    % e3 || Si
    e3 = bsxfun(@rdivide,M0,S0);
    % e2 = Si x [1,0,0], if Si || [1,0,0] --> e2 = [0,0,1]
    e2  = [zeros(1,nMagExt); e3(3,:); -e3(2,:)];
    e2(3,~any(abs(e2)>1e-10)) = 1;
    e2  = bsxfun(@rdivide,e2,sqrt(sum(e2.^2,1)));
    % e1 = e2 x e3
    e1  = cross(e2,e3);
end
% assign complex vectors that define the rotating coordinate system on
% every magnetic atom
zed = e1 + 1i*e2;
eta = e3;

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

SiSj = sqrt(S0(atom1).*S0(atom2));

% Creates temporary values for calculating matrix elements.
AD  =  shiftdim(sum(sum(etaL.*JJ.*etaR,2),1),1);
A20 = -S0(atom2).*AD;
D20 = -S0(atom1).*AD;
BC0 =  SiSj.*shiftdim(sum(sum(zedL.*JJ.*     zedR ,2),1),1);
AD0 =  SiSj.*shiftdim(sum(sum(zedL.*JJ.*conj(zedR),2),1),1);

% Calculates the contribution of the magnetic field (Zeeman term) to the Hamiltonian
MF = repmat(SI.field*rotC*obj.unit.muB * permute(mmat(SI.g,permute(eta,[1 3 2])),[1 3 2]),[1 2]);

% Creates the serial indices for every matrix element in ham matrix.
idxA1 = [atom1'         atom2'         ];
idxA2 = [atom1'         atom1'         ];
idxB  = [atom1'         atom2'+nMagExt ];
% transpose of idxB
%idxC  = [atom2'+nMagExt atom1'         ]; % SP1
idxD1 = idxA1+nMagExt;
idxD2 = [atom2'+nMagExt atom2'+nMagExt ];
idxMF = [(1:2*nMagExt)' (1:2*nMagExt)' ];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MEMORY MANAGEMENT LOOP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

% message for magnetic form factor calculation
yesNo = {'No' 'The'};
fprintf0(fid,[yesNo{param.formfact+1} ' magnetic form factor is'...
    ' included in the calculated structure factor.\n']);
% message for g-tensor calculation
fprintf0(fid,[yesNo{param.gtensor+1} ' g-tensor is included in the '...
    'calculated structure factor.\n']);

hklIdx = [floor(((1:nSlice)-1)/nSlice*nHkl)+1 nHkl+1];

% Empty omega dispersion of all spin wave modes, size: 2*nMagExt x nHkl.
omega = zeros(2*nMagExt,0);

% empty Sab
Sab = zeros(3,3,2*nMagExt,0);

% Empty matrices to save different intermediate results for further
% analysis: Hamiltonian, eigenvectors, dynamical structure factor in the
% rotating frame
if param.saveV
    Vsave = zeros(2*nMagExt,2*nMagExt,nHkl);
end
if param.saveH
    Hsave = zeros(2*nMagExt,2*nMagExt,nHkl);
end

sw_timeit(0,1,param.tid,'Spin wave spectrum calculation');

warn1 = false;

% calculate all magnetic form factors
if param.formfact
    spectra.formfact = true;
    % Angstrom^-1 units for Q
    hklA0 = 2*pi*(hkl0'/obj.basisvector)';
    % store form factor per Q point for each atom in the magnetic supercell
    % TODO check prod(nExt)? instead of nExt
    %FF = repmat(param.formfactfun(permute(obj.unit_cell.ff(1,:,obj.matom.idx),[3 2 1]),hklA0),[1 nExt]);
    param.FF = repmat(param.formfactfun(permute(obj.unit_cell.ff(1,:,obj.matom.idx),[3 2 1]),hklA0),[prod(nExt) 1]);
else
    spectra.formfact = false;
end

for jj = 1:nSlice
    % q indices selected for every chunk
    hklIdxMEM  = hklIdx(jj):(hklIdx(jj+1)-1);
    % q values contatining the k_m vector
    hklExtMEM  = hklExt(:,hklIdxMEM);
    % q values without the +/-k_m vector
    hklExt0MEM = hklExt0(:,hklIdxMEM);
    nHklMEM = size(hklExtMEM,2);

    % Creates the matrix of exponential factors nCoupling x nHkl size.
    % Extends dR into 3 x 3 x nCoupling x nHkl
    %     ExpF = exp(1i*permute(sum(repmat(dR,[1 1 nHklMEM]).*repmat(...
    %         permute(hklExtMEM,[1 3 2]),[1 nCoupling 1]),1),[2 3 1]))';
    ExpF = exp(1i*permute(sum(bsxfun(@times,dR,permute(hklExtMEM,[1 3 2])),1),[2 3 1]))';

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
    idx3   = repmat(1:nHklMEM,[3*nCoupling 1]);
    idxAll = [repmat(idxAll,[nHklMEM 1]) idx3(:)];
    idxAll = idxAll(:,[2 1 3]);

    ABCD   = ABCD';


    % quadratic form of the boson Hamiltonian stored as a square matrix
    ham = accumarray(idxAll,ABCD(:),[2*nMagExt 2*nMagExt nHklMEM]);

    ham = ham + repmat(accumarray([idxA2; idxD2],2*[A20 D20],[1 1]*2*nMagExt),[1 1 nHklMEM]);

    if any(SI.field)
        ham = ham + repmat(accumarray(idxMF,MF,[1 1]*2*nMagExt),[1 1 nHklMEM]);
    end

    ham = (ham + conj(permute(ham,[2 1 3])))/2;

    if param.hermit
        % All the matrix calculations are according to Colpa's paper
        % J.H.P. Colpa, Physica 93A (1978) 327-353
        [V, omega(:,hklIdxMEM)] = spinwave_hermit(ham, param, useMex, nMagExt);

    else
        % All the matrix calculations are according to White's paper
        % R.M. White, et al., Physical Review 139, A450?A454 (1965)

        gham = mmat(gComm,ham);

        [V, D, orthWarn] = eigorth(gham,param.omega_tol,useMex);

        orthWarn0 = orthWarn || orthWarn0;

        for ii = 1:nHklMEM
            % multiplication with g removed to get negative and positive
            % energies as well
            omega(:,end+1) = D(:,ii); %#ok<AGROW>
            M              = diag(gComm*V(:,:,ii)'*gComm*V(:,:,ii));
            V(:,:,ii)      = V(:,:,ii)*diag(sqrt(1./M));
        end
    end

    if param.saveV
        Vsave(:,:,hklIdxMEM) = V;
    end
    if param.saveH
        Hsave(:,:,hklIdxMEM) = ham;
    end

    Sab = cat(4, Sab, spinwave_spinspincorrel(V, RR, param, hklExt0MEM, hklIdxMEM, nHklMEM, nMagExt, zed, SI, nExt, S0));

    sw_timeit(jj/nSlice*100,0,param.tid);
end

if param.sortMode
    % sort the spin wave modes
    [omega, Sab] = sortmode(omega,reshape(Sab,9,size(Sab,3),[]));
    Sab          = reshape(Sab,3,3,size(Sab,2),[]);
end

[~,singWarn] = lastwarn;
% restore warning for singular matrix
warning(singWarn0.state,'MATLAB:nearlySingularMatrix');

% If number of formula units are given per cell normalize to formula
% unit
if obj.unit.nformula > 0
    Sab = Sab/double(obj.unit.nformula);
end

sw_timeit(100,2,param.tid);

fprintf0(fid,'Calculation finished.\n');

if warn1 && ~param.fitmode
    warning('spinw:spinwave:NonPosDefHamiltonian',['To make the Hamiltonian '...
        'positive definite, a small omega_tol value was added to its diagonal!'])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END MEMORY MANAGEMENT LOOP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

helical = false;

% Creates output structure with the calculated values.
spectra.omega    = omega;
spectra.Sab      = Sab;
spectra.hkl      = obj.unit.qmat\hkl(:,1:nHkl0);
spectra.hklA     = hklA;
spectra.helical  = helical;
spectra.norm     = false;
spectra.nformula = double(obj.unit.nformula);

% Save different intermediate results.
if param.saveV
    spectra.V = Vsave;
end
if param.saveH
    spectra.H = Hsave;
end

% save the important parameters
spectra.param.sortMode  = param.sortMode;
spectra.param.tol       = param.tol;
spectra.param.omega_tol = param.omega_tol;
spectra.param.hermit    = param.hermit;
spectra.title           = param.title;
spectra.gtensor         = param.gtensor;

if ~param.fitmode
    spectra.dateend = datestr(now);
    spectra.obj = copy(obj);
end

if ~param.gtensor && any(obj.single_ion.g)
    warning('spinw:spinwave:NonZerogTensor',['The SpinW model defines a '...
        'g-tensor that is not included in the calculation. Anisotropic '...
        'g-tensor values cannot be applied afterwards as they change relative'...
        'spin wave intensities!'])
end

% issue eigorth warning
if orthWarn0
    warning('spinw:spinwave:NoOrth','Eigenvectors of defective eigenvalues cannot be orthogonalised at some q-point!');
end

if strcmp(singWarn,'MATLAB:nearlySingularMatrix')
    lineLink = 'line 846';
    if feature('HotLinks')
        lineLink = ['<a href="matlab:opentoline([''' sw_rootdir 'swfiles' filesep '@spinw' filesep 'spinwave.m''' '],846,0)">' lineLink '</a>'];
    end
    warning('spinw:spinwave:nearlySingularMatrix',['Matrix is close '...
        'to singular or badly scaled. Results may be inaccurate.\n> In spinw/spinwave (' lineLink ')']);
    %fprintf(repmat('\b',[1 30]));
end

end

function [V, omega] = spinwave_hermit(ham, param, useMex, nMagExt)
    nHklMEM = size(ham, 3);

    % diagonal of the boson commutator matrix
    gCommd = [ones(nMagExt,1); -ones(nMagExt,1)];
    % boson commutator matrix
    gComm  = diag(gCommd);
    %gd = diag(g);

    % basis functions of the magnon modes
    V = zeros(2*nMagExt,2*nMagExt,nHklMEM);
    omega = zeros(2*nMagExt,0);

    if useMex && nHklMEM>1
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

function Sab = spinwave_spinspincorrel(V, RR, param, hklExt0MEM, hklIdxMEM, nHklMEM, nMagExt, zed, SI, nExt, S0)
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

    zeda = repmat(permute([zed conj(zed)],[1 3 4 2]),[1 3 2*nMagExt 1 nHklMEM]);
    % conj transpose of zeda
    zedb = conj(permute(zeda,[2 1 4 3 5]));

    % calculate magnetic structure factor using the hklExt0 Q-values
    % since the S(Q+/-k,omega) correlation functions also belong to the
    % F(Q)^2 form factor

    if param.formfact
        % include the form factor in the z^alpha, z^beta matrices
        zeda = zeda.*repmat(permute(param.FF(:,hklIdxMEM),[3 4 5 1 2]),[3 3 2*nMagExt 2 1]);
        zedb = zedb.*repmat(permute(param.FF(:,hklIdxMEM),[3 4 1 5 2]),[3 3 2 2*nMagExt 1]);
    end

    if param.gtensor
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

