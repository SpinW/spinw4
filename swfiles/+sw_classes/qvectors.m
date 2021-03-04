classdef qvectors < handle
    % Class for describing qvectors, which are the hkl coordinates in
    % three-dimensional space that describes the orientation of one unit cell
    % of the crystal lattice with respect to the origin
    properties
        hkl
        nChunk   % number of times to process (by taking into account the avalable memory)
    end
    properties(SetAccess=private)
        nHkl     % number of hkl points to calculate for (scalar)
        hklIdx   % (2 x nChunk) array with start and end indices of hkl points for a chunk
    end
    methods
        function self = qvectors(hkl)
            self.hkl = hkl;
        end
        function hklIdxMEM = getIdx(self, jj)
            hklIdxMEM = self.hklIdx(jj):(self.hklIdx(jj+1)-1);
        end
        function hklChunk = getChunk(self, jj)
            hklChunk = self.hkl(:, self.hklIdx(jj):(self.hklIdx(jj+1)-1));
        end
        function set.hkl(self, hkl)
            if iscell(hkl)
                % Converts a list of q-segments to a 3xN array
                self.hkl = sw_qscan(hkl);
            else
                self.hkl = hkl;
            end
            self.nHkl = size(self.hkl, 2);
        end
        function set.nChunk(self, nChunk)
            self.nChunk = nChunk;
            self.hklIdx = [floor(((1:nChunk)-1)/nChunk*self.nHkl)+1 self.nHkl+1];
        end
    end
end
