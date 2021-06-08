classdef TestSimulation < matlab.unittest.TestCase

	% All of these tests will be for AbstractCellSimulation,
	% but they will be done using CellGrowing
   
	methods (Test)

		function TestTorqueConversion(testCase)

			t = TestSim();

			n1 = Node(1,2,3);
			n2 = Node(2,2,3);
			n3 = Node(3,3,3);
			e1 = Edge(n1,n2);
			e2 = Edge(n2,n3);
			e3 = Edge(n3,n1);
			s = Surface(e1,e2,e3);

			t.AddNodesToSimulation([n1,n2,n3]);
			t.AddEdgesToSimulation([e1,e2,e3]);
			t.AddSurfacesToSimulation(s);

			u = s.GetUnitNormal();
			F = 2*u;
			A = [1.1,2.2,3];
			
			s.ApplyForce(F,A);

			t.ConvertTorquesToForces();

			testCase.verifyEqual(s.torque, [0,0,0]);

			rD = [2, 7/3, 3];

			rDA = A - rD;

			torque = cross(rDA, F);

			testCase.verifyEqual(F, n1.force);
			testCase.verifyEqual(F, n2.force);
			testCase.verifyEqual(F, n3.force);


		end

	end

end