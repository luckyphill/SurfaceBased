classdef TestFace < matlab.unittest.TestCase

	% All of these tests will be for AbstractCellSimulation,
	% but they will be done using CellGrowing
   
	methods (Test)

		function TestProperties(testCase)

			% Test the basic features

			n1 = Node(1,1,0);
			n2 = Node(2,2,0);
			n3 = Node(2,1,0);
			e1 = Edge(n1,n2);
			e2 = Edge(n2,n3);
			e3 = Edge(n3,n1);
			s = Face(e1,e2,e3);

			n1.eta = 1;
			n2.eta = 1;
			n3.eta = 1;

			testCase.verifyEqual(s.GetTotalDrag(), 3);
			testCase.verifyEqual(s.GetCentreOfDrag(), [5/3,4/3,0]);

			% At this point the nodes have been put in their correct
			% order in the surface. I'm not sure if they will necessarily
			% be in the order given above, but we'll try

			[rD1, rD2, rD3] = s.GetNodeCoordinatesRelative();
			testCase.verifyEqual(rD1, [-2/3,-1/3,0], 'RelTol', 1e-10);
			testCase.verifyEqual(rD2, [1/3,2/3,0], 'RelTol', 1e-10);
			testCase.verifyEqual(rD3, [1/3,-1/3,0], 'RelTol', 1e-10);



			ID = s.GetDragMatrix();

			IDcalc = [2/3, -1/3, 0;
						-1/3, 2/3, 0;
							0, 0, 4/3];

			testCase.verifyEqual(ID, IDcalc, 'RelTol', 1e-10);

			u = s.GetUnitNormal();

			% This specific direction because of the way nodes are
			% stored in the surface.
			testCase.verifyEqual(u, [0,0,-1]);

			[u,v,w] = s.GetCoordinateSystem();

			testCase.verifyEqual(u, [-2,1,0]/sqrt(5), 'RelTol', 1e-10);
			testCase.verifyEqual(v, [1,2,0]/sqrt(5), 'RelTol', 1e-10);
			testCase.verifyEqual(w, [0,0,-1], 'RelTol', 1e-10);


			[r1, r2, r3] = s.GetNodeCoordinatesLocal();

			testCase.verifyEqual(r1, [1,-4/3,0]/sqrt(5), 'RelTol', 1e-10);
			testCase.verifyEqual(r2, [0,sqrt(5)/3,0], 'RelTol', 1e-10);
			testCase.verifyEqual(r3, [-1,-1/3,0]/sqrt(5), 'RelTol', 1e-10);


			% Don't really need to test IDz at this point because
			% its just IDx+IDy
			[IDx,IDy,IDxy,IDz] = s.GetDragMatrixLocal();

			testCase.verifyEqual(IDx, 42/45 , 'RelTol', 1e-10);
			testCase.verifyEqual(IDy, 2/5, 'RelTol', 1e-10);
			testCase.verifyEqual(IDxy, -1/5, 'RelTol', 1e-10);


			% Identical to rD for given etas
			testCase.verifyEqual(s.GetMidPoint(), [5/3,4/3,0]);



			% Retest the whole lot for a differently oriented surface




		end

		function TestApplyForce(testCase)

			% Need to test that by applying a force
			% at certain points, the force on the nodes
			% and the torque is updated accordingly

			n1 = Node(1,2,3);
			n2 = Node(2,2,3);
			n3 = Node(3,3,3);
			e1 = Edge(n1,n2);
			e2 = Edge(n2,n3);
			e3 = Edge(n3,n1);
			s = Face(e1,e2,e3);

			u = s.GetUnitNormal();
			A = [1.1,2.2,3];
			F = 2*u;
			s.ApplyForce(F, A);

			testCase.verifyEqual(F, n1.force);
			testCase.verifyEqual(F, n2.force);
			testCase.verifyEqual(F, n3.force);

			rD = [2, 7/3, 3];

			rDA = A - rD;

			torque = cross(rDA, F);

			testCase.verifyEqual(s.torque, torque);

		end

	end

end