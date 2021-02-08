classdef magnetic_structure < handle
    properties
        S
        k
        n
        N_ext
        tol
    end

    properties(SetAccess=private)
        S_mag
        nMagExt
        km
        incomm
        helical
        constructed = false
    end

    methods
        % magstr should be output of spinw.magstr
        % TODO: This class should be refactored to do what spinw/genmagstr
        %       and spin/magstr does at present
        function self = magnetic_structure(magStr, tol)
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
            % Local (e1,e2,e3) coordinate system fixed to the moments,
            % e3||Si,ata
            % e2 = Si x [1,0,0], if Si || [1,0,0] --> e2 = [0,0,1]
            % e1 = e2 x e3
            %
            % e3 || Si
            e3 = self.S ./ self.S_mag;
            % e2 = Si x [1,0,0], if Si || [1,0,0] --> e2 = [0,0,1]
            e2  = [zeros(1, self.nMagExt); e3(3,:); -e3(2,:)];
            e2(3,~any(abs(e2)>1e-10)) = 1;
            e2  = e2./vecnorm(e2);
            % e1 = e2 x e3
            e1  = cross(e2,e3);
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
        end
    end
end
