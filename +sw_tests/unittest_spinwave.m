classdef unittest_spinwave < matlab.mock.TestCase
    % Runs through unit test for spinwave.m

    properties
        swobj = [];
        relToll = 0.01;
        absToll = 1e-6;
    end

    methods (TestClassSetup)
        function setup_spinw_model(testCase)
            % Just create a very simple FM 1D chain model
            testCase.swobj = spinw;
            testCase.swobj.genlattice('lat_const', [3 8 8], 'angled', [90 90 90]);
            testCase.swobj.addatom('r', [0 0 0],'S', 1, 'label', 'MNi2', 'color', 'blue');
            testCase.swobj.gencoupling('maxDistance', 7);
            testCase.swobj.addmatrix('value', -eye(3), 'label', 'Ja', 'color', 'green');
            testCase.swobj.addcoupling('mat', 'Ja', 'bond', 1);
            testCase.swobj.genmagstr('mode', 'direct', 'k', [0 0 0], 'n', [1 0 0], 'S', [0; 1; 0]);
        end
    end

    methods (Test)
        function test_noInput(testCase)
            % Tests that if call spinwave with no input, it calls the help
            % First mock the help call
            help_function = sw_tests.mock_function('swhelp');
            testCase.swobj.spinwave();
            testCase.assertEqual(help_function.n_calls, 1);
            testCase.assertEqual(help_function.arguments, {{'spinw.spinwave'}});
        end
        function test_symbollic(testCase)
            % Test that spinw.spinwavesym() is called if spinw.symbolic==true and spinw.spinwave() is called
            % Mock the spinwavesym and symbolic methods
            [mocksw, bh] = testCase.createMock(?spinw, 'MockedMethods', {'symbolic', 'spinwavesym'});
            % Make sure that spinw thinks it is symbolic
            testCase.assignOutputsWhen(withAnyInputs(bh.symbolic), true);
            % Make sure that spinwavesym outputs a sample spectra struct
            hkl = [1 2; 3 4; 5 6];
            out_spec = struct('hkl', hkl, 'swConv', [1; 2]);
            testCase.assignOutputsWhen(withAnyInputs(bh.spinwavesym), out_spec);
            % Call spinwave
            spec = mocksw.spinwave(hkl);
            testCase.assertCalled(withExactInputs(bh.spinwavesym()));
            testCase.assertEqual(spec, out_spec);
        end
        function test_incommensurate(testCase)
            % Tests that incommensurate calculation is ok
            hkl = {[0 0 0] [0 1 0] [1 0 0] 50};
            commensurate_spec = testCase.swobj.spinwave(hkl);
            testCase.swobj.genmagstr('mode', 'direct', 'k', [0.123 0 0], 'n', [1 0 0], 'S', [0; 1; 0]);
            % Forcing incomm struct means exchange parameters don't agree, so we need 'hermit' false
            incomm_spec = testCase.swobj.spinwave(hkl, 'hermit', false);
            testCase.assertEqual(size(incomm_spec.omega, 1), size(commensurate_spec.omega, 1) * 3);
            % Rests the structure for the following tests
            testCase.swobj.genmagstr('mode', 'direct', 'k', [0 0 0], 'n', [1 0 0], 'S', [0; 1; 0]);
        end
        function test_twin(testCase)
            % Tests that setting twins gives correct outputs (cell arrays)
            testCase.swobj.addtwin('axis', [0 0 1], 'phid', [60 120], 'vol', [1 1]);
            hkl = [1 2; 3 4; 5 6];
            out_spec = testCase.swobj.spinwave(hkl);
            testCase.assertTrue(iscell(out_spec.omega))
            testCase.assertEqual(numel(out_spec.omega), 3);  % Two _additional_ twins
            % Recalculate without twins for each set of hkl's and compare
            [~, rotQ] = testCase.swobj.twinq([0;0;0]);
            testCase.swobj.twin = struct('vol', 1, 'rotc', eye(3));
            for ii = 1:3
                spec_single = testCase.swobj.spinwave((hkl' * rotQ(:,:,ii))');
                testCase.assertEqual(spec_single.omega, out_spec.omega{ii});
            end
        end
        function test_formfact(testCase)
            % Tests that the form factor calculation is applied correctly
            hkl = {[0 0 0] [10 0 0] 100};
            % Runs calculation with/without formfactor
            spec_no_ff = sw_neutron(testCase.swobj.spinwave(hkl, 'formfact', false));
            spec_ff = sw_neutron(testCase.swobj.spinwave(hkl, 'formfact', true));
            import matlab.unittest.constraints.IsEqualTo
            import matlab.unittest.constraints.RelativeTolerance
            import matlab.unittest.constraints.AbsoluteTolerance
            theseBounds = RelativeTolerance(testCase.relToll) | AbsoluteTolerance(testCase.absToll);
            % The form factor is calculated using sw_mff, and the scaling is F(Q)^2 not F(Q).
            implied_ff = spec_ff.Sperp ./ spec_no_ff.Sperp;
            ff = sw_mff(testCase.swobj.unit_cell.label{1}, spec_ff.hklA);
            testCase.assertThat(implied_ff(1,:), IsEqualTo(ff.^2, 'Within', theseBounds));
        end
        function test_hermit(testCase)
            % Tests that the 'hermit' option to switch to a non-hermitian calculation works
            % First make the model non-Hermitian by adding a large axial SIA perpendicular to the spins
            testCase.swobj.addmatrix('label', 'K', 'value', diag([-1 0 0]));
            testCase.swobj.addaniso('K');
            hkl = {[0 0 0] [0 1 0] [1 0 0] 50};
            % Check that calling it with 'hermit' on gives an error
            testCase.assertError(@() testCase.swobj.spinwave(hkl, 'hermit', true), ?MException);
            % Now check that there are imaginary eigenvalues/energies in the output with 'hermit' off
            spec = testCase.swobj.spinwave(hkl, 'hermit', false);
            testCase.assertGreaterThan(sum(abs(imag(spec.omega(:)))), 0);
        end
    end

end
