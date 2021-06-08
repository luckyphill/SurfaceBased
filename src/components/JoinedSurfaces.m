classdef  JoinedSurfaces < AbstractSimulation


	properties

		% None yet

	end


	methods


		function obj = JoinedSurfaces()


			obj.stochasticJiggle = false;
			% A single node to provide force
			n = Node(1.5,1.3,0.3);

			% Construct a surface
			n1 = Node(1,1,0);
			n2 = Node(1.5,2,0);
			n3 = Node(2,1,0);
			n4 = Node(1.5,1,1);
			e1 = Edge(n1,n2);
			e2 = Edge(n2,n3);
			e3 = Edge(n3,n1);
			e4 = Edge(n1,n4);
			e5 = Edge(n4,n3);
			e6 = Edge(n4,n2);

			s1 = Surface(e1,e2,e3);
			s2 = Surface(e3,e4,e5);
			s3 = Surface(e2,e5,e6);

			obj.AddNodesToSimulation([n,n1,n2,n3,n4]);
			obj.AddEdgesToSimulation([e1,e2,e3,e4,e5,e6]);
			obj.AddSurfacesToSimulation([s1,s2]);

			springRate = 1;
			dAsym = 0;
			dSep = 1;
			dLim = 2;

			b = BruteForceNodeSurface(springRate, dAsym, dSep, dLim);

			obj.AddTissueBasedForce(b);


		end



	end


end