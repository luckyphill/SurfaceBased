classdef  SphereSurface < AbstractSimulation


	properties

		% None yet

	end


	methods


		function obj = SphereSurface()


			s1 = SphereShell([1.1,0,0]);

			obj.AddNodesToSimulation(s1.nodeList);
			obj.AddEdgesToSimulation(s1.edgeList);
			obj.AddFacesToSimulation(s1.faceList);
			obj.AddCellsToSimulation(s1);

			s2 = SphereShell([-1.1,0,0]);

			obj.AddNodesToSimulation(s2.nodeList);
			obj.AddEdgesToSimulation(s2.edgeList);
			obj.AddFacesToSimulation(s2.faceList);
			obj.AddCellsToSimulation(s2);

			for i = 1:length(obj.edgeList)

				e = obj.edgeList(i);

				e.naturalLength = e.GetLength();

			end

			springRate = 1;
			dAsym = -0.2;
			dSep = 0.1;
			dLim = 0.2;

			obj.AddTissueBasedForce(  BruteForceNodeFace(springRate, dAsym, dSep, dLim)  );
			obj.AddEdgeBasedForce( EdgeLinearSpringForce(1) );
			obj.AddCellBasedForce( ConstantInternalPressure(.3) );


		end



	end


end