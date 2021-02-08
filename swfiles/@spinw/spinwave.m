function spectra = spinwave(obj, hkl, varargin)
% calculates spin correlation function using linear spin wave theory
%
% ### Syntax
%
% `spectra = spinwave(obj,Q)`
%
% `spectra = spinwave(___,Name,Value)`
%
% ### Description
%
% `spinwave(obj,Q,Name,Value)` calculates spin wave dispersion and
% spin-spin correlation function at the reciprocal space points $Q$. The
% function can solve any single-k magnetic structure exactly and any
% multi-k magnetic structure appoximately and quadratic spinw-spin
% interactions as well as single ion anisotropy and magnetic field.
% Biquadratic exchange interactions are also implemented, however only for
% $k_m=0$ magnetic structures.
%
% If the magnetic ordering wavevector is non-integer, the dispersion is
% calculated using a coordinate system rotating from unit cell to unit
% cell. In this case the spin Hamiltonian has to fulfill this extra
% rotational symmetry which is not checked programatically.
%
% Some of the code of the function can run faster if mex files are used. To
% switch on mex files, use the `swpref.setpref('usemex',true)` command. For
% details see the [sw_mex] and [swpref.setpref] functions.
%
% ### Examples
%
% To calculate and plot the spin wave dispersion of the
% triangular lattice antiferromagnet ($S=1$, $J=1$) along the $(h,h,0)$
% direction in reciprocal space we create the built in triangular lattice
% model using `sw_model`.
%
% ```
% >>tri = sw_model('triAF',1)
% >>spec = tri.spinwave({[0 0 0] [1 1 0]})
% >>sw_plotspec(spec)
% >>snapnow
% ```
%
% ### Input Arguments
%
% `obj`
% : [spinw] object.
%
% `Q`
% : Defines the $Q$ points where the spectra is calculated, in reciprocal
%   lattice units, size is $[3\times n_{Q}]$. $Q$ can be also defined by
%   several linear scan in reciprocal space. In this case `Q` is cell type,
%   where each element of the cell defines a point in $Q$ space. Linear scans
%   are assumed between consecutive points. Also the number of $Q$ points can
%   be specified as a last element, it is 100 by defaults.
%
%   For example to define a scan along $(h,0,0)$ from $h=0$ to $h=1$ using
%   200 $Q$ points the following input should be used:
%   ```
%   Q = {[0 0 0] [1 0 0]  200}
%   ```
%
%   For symbolic calculation at a general reciprocal space point use `sym`
%   type input.
%
%   For example to calculate the spectrum along $(h,0,0)$ use:
%   ```
%   Q = [sym('h') 0 0]
%   ```
%   To calculate spectrum at a specific $Q$ point symbolically, e.g. at
%   $(0,1,0)$ use:
%   ```
%   Q = sym([0 1 0])
%   ```
%
% ### Name-Value Pair Arguments
%
% `'formfact'`
% : If true, the magnetic form factor is included in the spin-spin
%   correlation function calculation. The form factor coefficients are
%   stored in `obj.unit_cell.ff(1,:,atomIndex)`. Default value is `false`.
%
% `'formfactfun'`
% : Function that calculates the magnetic form factor for given $Q$ value.
%   value. Default value is `@sw_mff`, that uses a tabulated coefficients
%   for the form factor calculation. For anisotropic form factors a user
%   defined function can be written that has the following header:
%   ```
%   F = formfactfun(atomLabel,Q)
%   ```
%   where the parameters are:
%   * `F`           row vector containing the form factor for every input
%                   $Q$ value
%   * `atomLabel`   string, label of the selected magnetic atom
%   * `Q`           matrix with dimensions of $[3\times n_Q]$, where each
%                   column contains a $Q$ vector in $\\ang^{-1}$ units.
%
% `'gtensor'`
% : If true, the g-tensor will be included in the spin-spin correlation
%   function. Including anisotropic g-tensor or different
%   g-tensor for different ions is only possible here. Including a simple
%   isotropic g-tensor is possible afterwards using the [sw_instrument]
%   function.
%
% `'fitmode'`
% : If `true`, function is optimized for multiple consecutive calls (e.g.
%   the output spectrum won't contain the copy of `obj`), default is
%   `false`.
%
% `'notwin'`
% : If `true`, the spectra of the twins won't be calculated. Default is
% `false`.
%
% `'sortMode'`
% : If `true`, the spin wave modes will be sorted by continuity. Default is
%   `true`.
%
% `'optmem'`
% : Parameter to optimise memory usage. The list of Q values will be cut
%   into `optmem` number of pieces and will be calculated piece by piece to
%   decrease peak memory usage. Default value is 0, when the number
%   of slices are determined automatically from the available free memory.
%
% `'tol'`
% : Tolerance of the incommensurability of the magnetic ordering wavevector.
%   Deviations from integer values of the ordering wavevector smaller than
%   the tolerance are considered to be commensurate. Default value is
%   $10^{-4}$.
%
% `'omega_tol'`
% : Tolerance on the energy difference of degenerate modes when
%   diagonalising the quadratic form, default value is $10^{-5}$.
%
% `'hermit'`
% : Method for matrix diagonalization with the following logical values:
%
%   * `true`    using Colpa's method (for details see [J.H.P. Colpa, Physica 93A (1978) 327](http://www.sciencedirect.com/science/article/pii/0378437178901607)),
%               the dynamical matrix is converted into another Hermitian
%               matrix, that will give the real eigenvalues.
%   * `false`   using the standard method (for details see [R.M. White, PR 139 (1965) A450](https://journals.aps.org/pr/abstract/10.1103/PhysRev.139.A450))
%               the non-Hermitian $\mathcal{g}\times \mathcal{H}$ matrix
%               will be diagonalised, which is computationally less
%               efficient. Default value is `true`.
%
% {{note Always use Colpa's method, except when imaginary eigenvalues are
%   expected. In this case only White's method work. The solution in this
%   case is wrong, however by examining the eigenvalues it can give a hint
%   where the problem is.}}
%
% `'saveH'`
% : If true, the quadratic form of the Hamiltonian is also saved in the
%   output. Be carefull, it can take up lots of memory. Default value is
%   `false`.
%
% `'saveV'`
% : If true, the matrices that transform the normal magnon modes into the
%   magnon modes localized on the spins are also saved into the output. Be
%   carefull, it can take up lots of memory. Default value is `false`.
%
% `'saveSabp'`
% : If true, the dynamical structure factor in the rotating frame
%   $S'(k,\omega)$ is saved. Default value is `false`.
%
% `'title'`
% : Gives a title string to the simulation that is saved in the output.
%
% `'fid'`
% : Defines whether to provide text output. The default value is determined
%   by the `fid` preference stored in [swpref]. The possible values are:
%   * `0`   No text output is generated.
%   * `1`   Text output in the MATLAB Command Window.
%   * `fid` File ID provided by the `fopen` command, the output is written
%           into the opened file stream.
%
% `'tid'`
% : Determines if the elapsed and required time for the calculation is
%   displayed. The default value is determined by the `tid` preference
%   stored in [swpref]. The following values are allowed (for more details
%   see [sw_timeit]):
%   * `0` No timing is executed.
%   * `1` Display the timing in the Command Window.
%   * `2` Show the timing in a separat pup-up window.
%
% ### Output Arguments
%
% `spectra`
% : structure, with the following fields:
%   * `omega`   Calculated spin wave dispersion with dimensions of
%               $[n_{mode}\times n_{Q}]$.
%   * `Sab`     Dynamical structure factor with dimensins of
%               $[3\times 3\times n_{mode}\times n_{Q}]$. Each
%               `(:,:,i,j)` submatrix contains the 9 correlation functions
%               $S^{xx}$, $S^{xy}$, $S^{xz}$, etc. If given, magnetic form
%               factor is included. Intensity is in \\hbar units, normalized
%               to the crystallographic unit cell.
%   * `H`       Quadratic form of the Hamiltonian. Only saved if `saveH` is
%               true.
%   * `V`       Transformation matrix from the normal magnon modes to the
%               magnons localized on spins using the following:
%               $x_i = \sum_j V_{ij} \times x_j'$
%               Only saved if `saveV` is true.
%   * `Sabp`    Dynamical structure factor in the rotating frame,
%               dimensions are $[3\times 3\times n_{mode}\times n_{Q}]$,
%               but the number of modes are equal to twice the number of
%               magnetic atoms.
%   * `formfact`  Cell containing the labels of the magnetic ions if form
%               factor in included in the spin-spin correlation function.
%   * `cmplxBase` The local coordinate system on each magnetic moment is
%               defined by the complex magnetic moments:
%               $\begin{align}  e_1 &= \Im(\hat{M})\\
%                               e_3 &= Re(\hat{M})\\
%                               e_2 &= e_3\times e_1
%               \end{align}$
%
%   * `hkl`     Contains the input $Q$ values, dimensions are $[3\times n_{Q}]$.
%   * `hklA`    Same $Q$ values, but in $\\ang^{-1}$ unit, in the
%               lab coordinate system, dimensins are $[3\times n_{Q}]$.
%   * `incomm`  Logical value, tells whether the calculated spectra is
%               incommensurate or not.
%   * `obj`     The copy (clone) of the input `obj`, see [spinw.copy].
%
% The number of magnetic modes (labeled by `nMode`) for commensurate
% structures is double the number of magnetic atoms in the magnetic cell.
% For incommensurate structures this number is tripled due to the
% appearance of the $(Q\pm k_m)$ Fourier components in the correlation
% functions. For every $Q$ points in the following order:
% $(Q-k_m,Q,Q+k_m)$.
%
% If several twins exist in the sample, `omega` and `Sab` are packaged into
% a cell, that contains $n_{twin}$ number of matrices.
%
% ### See Also
%
% [spinw] \| [spinw.spinwavesym] \| [sw_mex] \| [spinw.powspec] \| [sortmode]
%

% calculate symbolic spectrum if obj is in symbolic mode
if obj.symbolic
    if nargin > 1
        spectra = spinwave_symbolic(obj, hkl, varargin{:});
    else
        spectra = spinwave_symbolic(obj, [], varargin{:});
    end
    return
end

% help when executed without argument
if nargin==1
    swhelp('spinw.spinwave');
    spectra = [];
    return
end

sw_calculator = sw_classes.spin_wave_calculator(obj, varargin{:});
spectra = sw_calculator.calculateSpinWave(hkl);

end

function spectra = spinwave_symbolic(obj, hkl, varargin)
    if obj.symbolic
        if numel(hkl) == 3
            hkl = sym(hkl);
        end

        if ~isa(hkl,'sym')
            inpForm.fname  = {'fitmode'};
            inpForm.defval = {false    };
            inpForm.size   = {[1 1]    };
            param0 = sw_readparam(inpForm, varargin{:});

            if ~param0.fitmode
                warning('spinw:spinwave:MissingInput','No hkl value was given, spin wave spectrum for general Q (h,k,l) will be calculated!');
            end
            spectra = obj.spinwavesym(varargin{:});
        else
            spectra = obj.spinwavesym(varargin{:},'hkl',hkl);
        end
        return
    end
end

% parse input arguments
[param, useMex] = spinwave_parse_input(varargin);

% checks for twins
if param.notwin
    nTwin = 1;
else
    nTwin = size(obj.twin.vol,2);
    % if there is only one twin and it is the identity set param.notwin true
    rotc1 = obj.twin.rotc(:,:,1) - eye(3);
    if (nTwin == 1) && norm(rotc1(:))==0
        param.notwin = true;
    end
end

hkl    = cell2mat(hkl);
hkl0   = cell2mat(hkl0);
hklExt = cell2mat(hklExt);

% determines a twin index for every q point
twinIdx = repmat(1:nTwin,[nHkl 1]);
twinIdx = twinIdx(:);

% Create the interaction matrix and atomic positions in the extended
% magnetic unit cell.
[SS, SI, RR] = obj.intmatrix('fitmode',true,'conjugate',true);

% add the dipolar interactions to SS.all
SS.all = [SS.all SS.dip];

% is there any biquadratic exchange
bq = SS.all(15,:)==1;

% Biquadratic exchange only supported for commensurate structures
if incomm && any(bq)
    error('spinw:spinwave:Biquadratic','Biquadratic exchange can be only calculated for k=0 structures!');
end

if any(bq)
    % Separate the biquadratic couplings
    % Just use the SS.bq matrix produced by intmatrix(), it won't contain
    % the transpose matrices (not necessary for biquadratic exchange)
    % TODO check whether to keep the transposed matrices to be sure
    SS.bq = SS.all(1:6,bq);
    % Keep only the quadratic exchange couplings
    SS.all = SS.all(1:14,SS.all(15,:)==0);
end

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

if incomm
    fprintf0(fid,['Calculating INCOMMENSURATE spin wave spectra '...
        '(nMagExt = %d, nHkl = %d, nTwin = %d)...\n'],nMagExt, nHkl0, nTwin);
if param.notwin
    spectra = obj.spinwave_single(hkl, param);
else
    hkl  = obj.twinq(hkl);
    for ii = 1:nTwin
        % Need to pass rotC to spinwave_single to calculate
        rotC = obj.twin.rotc(:,:,ii);
        spec(ii) = obj.spinwave_single(hkl{ii}, param, rotC);
        if ii == 1
            spectra = spec(1);
            spectra.omega = {spectra.omega};
            spectra.Sab = {spectra.Sab};
            if isfield(spectra, 'V')
                spectra.V = {spectra.V};
            end
        else
            spectra = append_twin(spectra, spec(ii));
        end
        if sum(sum(abs(rotC - eye(3)))) > 1e-10
            % Rotates the spin-spin correlation tensor using the twin rotation matrix
            sSab = size(spectra.Sab{ii});
            Sab = reshape(spectra.Sab{ii}, 3, 3, []);   % Originally 3x3xNModexNQ now 3x3x(NMode*NQ)
            Sab = arrayfun(@(idx)(rotC*Sab(:,:,idx)*(rotC')), 1:size(Sab,3), 'UniformOutput', false);
            % Sab is now a cell array of 3x3 matrices
            spectra.Sab{ii} = reshape(cat(3, Sab{:}), sSab);  % Converts back to normal matrix
        end
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [param, useMex] = spinwave_parse_input(arg_in)

    pref = swpref;

    % use mex file by default?
    useMex = pref.usemex;

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

    param = sw_readparam(inpForm, arg_in{:});

    if ~param.fitmode
        % save the time of the beginning of the calculation
        spectra.datestart = datestr(now);
    end

    if param.fitmode
        param.sortMode = false;
        param.tid = 0;
    end

    if param.tid == -1
        param.tid = pref.tid;
    end

    if param.fid == -1
        param.fid = pref.fid;
    end
    % Dynamical structure factor from S^alpha^beta(k) correlation function.
    % Sab(alpha,beta,iMode,iHkl), size: 3 x 3 x 2*nMagExt x nHkl.
    % Normalizes the intensity to single unit cell.
    Sab = cat(4,Sab,squeeze(sum(zeda.*ExpFL.*VExtL,4)).*squeeze(sum(zedb.*ExpFR.*VExtR,3))/prod(nExt));

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

end

if ~param.notwin
    % Rotate the calculated correlation function into the twin coordinate
    % system using rotC
    SabAll = cell(1,nTwin);
    for ii = 1:nTwin
        % select the ii-th twin from the Q points
        idx    = (1:nHkl0) + (ii-1)*nHkl0;
        % select correlation function of twin ii
        SabT   = Sab(:,:,:,idx);
        % size of the correlation function matrix
        sSabT  = size(SabT);
        % convert the matrix into cell of 3x3 matrices
        SabT   = reshape(SabT,3,3,[]);
        % select the rotation matrix of twin ii
        rotC   = obj.twin.rotc(:,:,ii);
        % rotate correlation function using arrayfun
        SabRot = arrayfun(@(idx)(rotC*SabT(:,:,idx)*(rotC')),1:size(SabT,3),'UniformOutput',false);
        SabRot = cat(3,SabRot{:});
        % resize back the correlation matrix
        SabAll{ii} = reshape(SabRot,sSabT);
    end
    Sab = SabAll;

    if nTwin == 1
        Sab = Sab{1};
    else
        omega = mat2cell(omega,size(omega,1),repmat(nHkl0,[1 nTwin]));
    end

function spectra = spinwave_symbolic(obj, hkl, varargin)
    if obj.symbolic
        if numel(hkl) == 3
            hkl = sym(hkl);
        end

        if ~isa(hkl,'sym')
            inpForm.fname  = {'fitmode'};
            inpForm.defval = {false    };
            inpForm.size   = {[1 1]    };
            param0 = sw_readparam(inpForm, varargin{:});

            if ~param0.fitmode
                warning('spinw:spinwave:MissingInput','No hkl value was given, spin wave spectrum for general Q (h,k,l) will be calculated!');
            end
            spectra = obj.spinwavesym(varargin{:});
        else
            spectra = obj.spinwavesym(varargin{:},'hkl',hkl);
        end
        return
    end
end

function spectra = append_twin(spectra, spec_new)
    spectra.omega = [spectra.omega {spec_new.omega}];
    spectra.Sab = [spectra.Sab {spec_new.Sab}];
    if isfield(spectra, 'V');
        spectra.V = [spectra.V {spec_new.V}];
    end
    if isfield(spectra, 'H');
        spectra.H = cat(3, spectra.H, spec_new.H);
    end
    warning('spinw:spinwave:nearlySingularMatrix',['Matrix is close '...
        'to singular or badly scaled. Results may be inaccurate.\n> In spinw/spinwave (' lineLink ')']);
    %fprintf(repmat('\b',[1 30]));
end

end
