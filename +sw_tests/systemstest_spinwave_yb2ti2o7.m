classdef systemstest_spinwave_yb2ti2o7 < matlab.unittest.TestCase

    properties
        generate_reference_data = false;
        reference_data = [];
        reference_data_dir = fullfile('.', 'test_data');
        reference_data_file = 'systemstest_spinwave_yb2ti2o7.mat';
        relToll = 0.01;
        absToll = 1e-6;
        yto = [];
    end

    methods (TestClassSetup)
        function get_reference_data(testCase)
            fname = fullfile(testCase.reference_data_dir, testCase.reference_data_file);
            if ~exist(testCase.reference_data_dir, 'dir')
                mkdir(testCase.reference_data_dir);
            end
            if ~exist(fname, 'file') || testCase.generate_reference_data
                testCase.generate_reference_data = true;
                tmp = []; save(fname, 'tmp');
                warning('Generating reference data');
            else
                testCase.reference_data = load(fname);
            end
        end
    end

    methods (TestClassTeardown)
        function save_reference_data(testCase)
            if testCase.generate_reference_data
                fname = fullfile(testCase.reference_data_dir, testCase.reference_data_file);
                ref_dat = load(fname);
                ref_dat = rmfield(ref_dat, 'tmp');
                save(fname, '-struct', 'ref_dat');
            end
        end
    end

    properties (TestParameter)
        B = {2 5};
        Q = {{[-0.5 -0.5 -0.5] [2 2 2]} {[1 1 -2] [1 1 1.5]} {[2 2 -2] [2 2 1.5]} {[-0.5 -0.5 0] [2.5 2.5 0]} {[0 0 1] [2.3 2.3 1]}};
    end

    methods (TestMethodSetup)
        function prepareForRun(testCase)
            % From Tutorial 20, Yb2Ti2O7, based on PRX 1, 021002 (2011)
            % First set up the crystal structure
            symStr = '-z, y+3/4, x+3/4; z+3/4, -y, x+3/4; z+3/4, y+3/4, -x; y+3/4, x+3/4, -z; x+3/4, -z, y+3/4; -z, x+3/4, y+3/4';
            yto = spinw;
            a = 10.0307;
            yto.genlattice('lat_const',[a a a],'angled',[90 90 90],'spgr',symStr,'label','F d -3 m Z')
            yto.addatom('label','Yb3+','r',[1/2 1/2 1/2],'S',1/2)
            % We generate the list of bonds.
            yto.gencoupling
            % We create two 3x3 matrix, one for the first neighbor anisotropic exchange
            % and one for the anisotropic g-tensor. And assign them appropriately.
            yto.addmatrix('label', 'J1', 'color', [255 0 0], 'value', 1)
            yto.addmatrix('label', 'g0', 'color', [0 0 255], 'value', -0.84*ones(3)+4.32*eye(3));
            yto.addcoupling('mat', 'J1', 'bond', 1)
            yto.addg('g0')
            % Sets the correct values for the matrix elements of J1
            J1 = -0.09; J2 = -0.22; J3 = -0.29; J4 = 0.01;
            yto.setmatrix('mat','J1','pref',[J1 J3 J2 -J4]);
            testCase.yto = yto;
        end
    end

    methods (Test)
        function test_yto(testCase, B, Q)
            n = [1 -1 0];
            % set magnetic field
            testCase.yto.field(n/norm(n)*B);
            % create fully polarised magnetic structure along the field direction
            testCase.yto.genmagstr('S',n','mode','helical');
            % find best structure using steepest descendend
            testCase.yto.optmagsteep;
            ytoSpec = testCase.yto.spinwave([Q {200}],'gtensor',true);
            ytoSpec = sw_neutron(ytoSpec);
            % bin the spectrum in energy
            ytoSpec = sw_egrid(ytoSpec,'Evect',linspace(0,2,500),'component','Sperp');
            %figure; sw_plotspec(ytoSpec,'axLim',[0 0.5],'mode',3,'dE',0.09,'colorbar',false,'legend',false); title(''); caxis([0 60]); colormap(jet);
            import matlab.unittest.constraints.IsEqualTo
            import matlab.unittest.constraints.RelativeTolerance
            import matlab.unittest.constraints.AbsoluteTolerance
            theseBounds = RelativeTolerance(testCase.relToll) | AbsoluteTolerance(testCase.absToll);
            fieldname = ['d' reshape(dec2hex(get_hash({B Q})),1,[])];
            if testCase.generate_reference_data
                data.input = struct(testCase.yto);
                data.spec = {ytoSpec.omega ytoSpec.Sab ytoSpec.swConv ytoSpec.swInt};
                filename = fullfile(testCase.reference_data_dir, testCase.reference_data_file);
                tmpstr.(fieldname) = data;
                save(filename, '-append', '-struct', 'tmpstr');
            else
                ref_data = testCase.reference_data.(fieldname);
                test_data.input = struct(testCase.yto);
                test_data.spec = {ytoSpec.omega ytoSpec.Sab ytoSpec.swConv ytoSpec.swInt};
                testCase.verifyThat(test_data, IsEqualTo(ref_data, 'Within', theseBounds));
            end
        end
    end

end


function out = get_hash(obj)
    % Calculates a hash for an object or struct using undocumented built-ins
    % Based on DataHash (https://uk.mathworks.com/matlabcentral/fileexchange/31272-datahash)
    Engine = java.security.MessageDigest.getInstance('MD5');
    Engine.update(getByteStreamFromArray(obj));
    out = typecast(Engine.digest, 'uint8');
end
