classdef TestForce < matlab.unittest.TestCase

	% All of these tests will be for AbstractCellSimulation,
	% but they will be done using CellGrowing
   
	methods (Test)

		function TestBruteForceNodeNode(testCase)

			f = BruteForceNodeNode(1, 0, 1, 1.5);

			n1 = Node(0.5,0,0);
			n2 = Node(0,0,0);

			f.ApplyForce(n1,n2);

			testCase.verifyGreaterThan(n1.force(1), 0);
			testCase.verifyLessThan(n2.force(1), 0);

		end

	end

end
