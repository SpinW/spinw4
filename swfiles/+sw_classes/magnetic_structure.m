classdef magnetic_structure < handle
    % Class to handle information required by the spin wave calculations from the magnetic structure
    properties
        S      % Magnetic moments, 3xN array; N=number of spins
        k      % Magnetic propagation vector, 3x1 vector
        n      % Normal vector to spin plane, 3x1 vector
        N_ext  % Supercell size, 3x1 vector
        tol    % Tolerance factor for incommensurate determination
    end

    properties(SetAccess=private)
        S_mag    % Magnitude of spins, 1xN vector
        nMagExt  % Number of spins in supercell, scalar
        km       % Magnetic propagation vector in the supercell, 3x1 vector
        incomm   % Flag to indicate if structure in incommensurate
        helical  % Flag to indicate if structure is helical
        nx       % Rodriguez's rotation matrix to convert from rotating to lab frame
        R1       % R1 matrix in eq (39) of Toth and Lake (2015)
        R2       % R2 matrix in eq (39) of Toth and Lake (2015)
        constructed = false
    end

    methods
        % magstr should be output of spinw.magstr
        % TODO: This class should be refactored to do what spinw/genmagstr
        %       and spin/magstr does at present
        function self = magnetic_structure(magStr, tol)
            % Class to encapsulate the magnetic structure
            %
            % ### Syntax
            % `mag_str_obj = magnetic_structure(magstr, tol)
            %
            % ### Description
            % This class is an implementation detail and is subject to change.
            % Please use the genmagstr() method of the spinw class to define the magnetic structure
            %
            % ### Input arguments
            %
            % `magstr`
            % : A structure with fields S, k and n defining a spinw single-k magnetic structure
            %
            % `tol`
            % : Scalar tolerance to determine if a propagation vector is incommensurate or not (default: 1e-5)
            %
            % ### Output arguments
            %
            % `mag_str_obj`
            % : A `magnetic_structure` object
            if nargin < 2
                tol = 1e-5;
            end
            if isstruct(magStr)
                if isempty(magStr.S)
                    error('spinw:spinwave:NoMagneticStr','No magnetic structure defined in self.spinWaveObject!');
                end
                self.S = magStr.S;
                self.k = magStr.k;
                self.n = magStr.n;
                self.N_ext = magStr.N_ext;
            else 
                error('non-struct constructor not implemented yet');
            end
            self.tol = tol;
            self.constructed = true;
            self.set_dependent_properties();
        end

        function set.S(self, S)
            self.S = S;
            self.set_dependent_properties();
        end
        
        function set.k(self, k)
            self.k = k;
            self.set_dependent_properties();
        end
        
        function set.n(self, n)
            self.n = n;
            self.set_dependent_properties();
        end
        
        function set.N_ext(self, N)
            self.N_ext = N;
            self.set_dependent_properties();
        end
        
        function [e1, e2, e3] = get_local_basis(self)
            % Returns the local (e1,e2,e3) coordinate system fixed to the moments,
            %
            % ### Syntax
            % `[e1, e2, e3] = magnetic_structure.get_local_basis()
            %
            % ### Output Arguments
            %
            % `e1`, `e2`, `e3`
            % : 3-vectors defining the local coordinate systems with:
            %   e3||Si,ata
            %   e2 = Si x [1,0,0], if Si || [1,0,0] --> e2 = [0,0,1]
            %   e1 = e2 x e3
 
            % e3 || Si
            e3 = self.S ./ self.S_mag;
            % e2 = Si x [1,0,0], if Si || [1,0,0] --> e2 = [0,0,1]
            e2  = [zeros(1, self.nMagExt); e3(3,:); -e3(2,:)];
            e2(3,~any(abs(e2)>1e-10)) = 1;
            e2  = e2./vecnorm(e2);
            % e1 = e2 x e3
            e1  = cross(e2,e3);
        end
        
        function out = transform_frame(self, in)
            % Applies the transformation from the lab to the rotating frame
            %
            % ### Syntax
            % `out_vec = magnetic_structure.transform_frame(in_vec)`
            %
            % ### Input Arguments
            %
            % `in_vec`
            % : a 3-vector input to be transformed
            %
            % ### Output Arguments
            %
            % `out_vec`
            % : the transformed 3-vector
            m1 = eye(3);
            out = 1/2*in - 1/2*mmat(mmat(self.nx,in),self.nx) + ...
                  1/2*mmat(mmat(self.R2-m1,in),self.R2) + ...  % R2 also called nxn in old code
                  1/2*mmat(mmat(self.R2,in),2*self.R2-m1);
        end
    end

    methods(Access=private)
        % Calculates some dependent properties
        function set_dependent_properties(self)
            if ~self.constructed
                return;
            end
            self.S_mag = vecnorm(self.S);     % Magnitude of spins
            self.nMagExt = size(self.S, 2);   % Number of spins - name should probably be changed
            self.km = self.k .* self.N_ext;   % The magnetic propagation vector
            self.incomm = any(abs(self.km-round(self.km)) > self.tol);   % Whether structure is incommensurate
            km = self.km;
            tol = self.tol * 2;
            self.helical = sum(abs(mod(abs(2*km)+tol,1)-tol).^2) > tol;  % Whether incomm wavevector is multiple of 2
            self.nx  = [0 -self.n(3) self.n(2);
                        self.n(3) 0 -self.n(1);
                       -self.n(2) self.n(1) 0];
            self.R2 = self.n'*self.n;  % Previously also called nxn in code
            self.R1 = 1/2*(eye(3) - self.R2 - 1i*self.nx);
        end
    end
end
