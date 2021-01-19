classdef systemstest_spinwave < matlab.unittest.TestCase
    % Base class for all systems test of spinwave.m based on tutorials

    properties
        generate_reference_data = false;
        reference_data = [];
        reference_data_dir = fullfile('.', 'test_data');
        relToll = 0.01;
        absToll = 1e-6;
        swobj = [];
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

    methods (Static)
        function out = get_hash(obj)
            % Calculates a hash for an object or struct using undocumented built-ins
            % Based on DataHash (https://uk.mathworks.com/matlabcentral/fileexchange/31272-datahash)
            Engine = java.security.MessageDigest.getInstance('MD5');
            Engine.update(getByteStreamFromArray(obj));
            out = typecast(Engine.digest, 'uint8');
        end
    end

    methods
        function generate_or_verify(testCase, spec, pars)
            import matlab.unittest.constraints.IsEqualTo
            import matlab.unittest.constraints.RelativeTolerance
            import matlab.unittest.constraints.AbsoluteTolerance
            theseBounds = RelativeTolerance(testCase.relToll) | AbsoluteTolerance(testCase.absToll);
            fieldname = ['d' reshape(dec2hex(testCase.get_hash(pars)),1,[])];
            if testCase.generate_reference_data
                data.input = struct(testCase.swobj);
                data.spec = {spec.omega spec.Sab};
                if isfield(spec, 'swConv'); data.spec = [data.spec {spec.swConv}]; end
                if isfield(spec, 'swInt');  data.spec = [data.spec {spec.swInt}];  end
                filename = fullfile(testCase.reference_data_dir, testCase.reference_data_file);
                tmpstr.(fieldname) = data;
                save(filename, '-append', '-struct', 'tmpstr');
            else
                ref_data = testCase.reference_data.(fieldname);
                test_data.input = struct(testCase.swobj);
                test_data.spec = {spec.omega spec.Sab};
                if isfield(spec, 'swConv'); test_data.spec = [test_data.spec {spec.swConv}]; end
                if isfield(spec, 'swInt');  test_data.spec = [test_data.spec {spec.swInt}];  end
                testCase.verifyThat(test_data, IsEqualTo(ref_data, 'Within', theseBounds));
            end
        end
    end

end


