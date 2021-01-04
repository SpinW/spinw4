function result = run_tests(out_dir)
    if ~exist('spinw', 'file') && exist('swfiles', 'file')
        addpath(genpath('swfiles'));
        addpath(genpath('external'));
        addpath('dat_files');
    end
    if nargin == 0
        out_dir = fullfile(pwd, 'test_output');
        if ~exist(out_dir, 'dir')
            mkdir(out_dir);
        end
    end

    import matlab.unittest.TestSuite
    import matlab.unittest.TestRunner
    import matlab.unittest.plugins.CodeCoveragePlugin
    import matlab.unittest.plugins.codecoverage.CoberturaFormat

    reportFormat = CoberturaFormat(fullfile(out_dir, 'spinw_coverage.xml'));
    coverage_plugin = CodeCoveragePlugin.forFolder('swfiles', 'Producing', reportFormat);
    %coverage_plugin = CodeCoveragePlugin.forFile('swfiles/@spinw/spinwave.m', 'Producing', reportFormat);

    suite = TestSuite.fromPackage('sw_tests');
    runner = TestRunner.withTextOutput;
    runner.addPlugin(coverage_plugin);
    result = runner.run(suite);
end
