classdef  SphereSurface < AbstractSimulation


	properties

		% None yet

	end


	methods


		function obj = SphereSurface()


			s = SphereShell([0,0,0]);

			obj.AddNodesToSimulation(s.nodeList);
			obj.AddEdgesToSimulation(s.edgeList);
			obj.AddSurfacesToSimulation(s.surfList);

			springRate = 1;
			dAsym = 0;
			dSep = 1;
			dLim = 2;

			b = BruteForceNodeSurface(springRate, dAsym, dSep, dLim);

			obj.AddTissueBasedForce(b);


		end



	end


end