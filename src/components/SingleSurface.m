classdef  SingleSurface < AbstractSimulation


	properties

		% None yet

	end


	methods


		function obj = SingleSurface()


			obj.stochasticJiggle = false;
			% A single node to provide force
			n = Node(1.5,1.1,0.5);

			% Construct a surface
			n1 = Node(1,1,0);
			n2 = Node(2,2,0);
			n3 = Node(2,1,0);
			e1 = Edge(n1,n2);
			e2 = Edge(n2,n3);
			e3 = Edge(n3,n1);
			s = Surface(e1,e2,e3);

			obj.AddNodesToSimulation([n,n1,n2,n3]);
			obj.AddEdgesToSimulation([e1,e2,e3]);
			obj.AddSurfacesToSimulation(s);

			springRate = 1;
			dAsym = 0;
			dSep = 1;
			dLim = 2;

			b = BruteForceNodeSurface(springRate, dAsym, dSep, dLim);

			obj.AddTissueBasedForce(b);


		end



	end


end