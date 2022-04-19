classdef  SphereSurface < AbstractSimulation


	properties

		pathName

	end


	methods


		function obj = SphereSurface()


			obj.dt = 0.001;

			springRate = 1;
			dAsym = -0.2;
			dSep = 0.1;
			dLim = 0.2;

			updatePeriod = 0.1; % 0.1 hours
			updateFrequency = floor(updatePeriod / obj.dt);

			s1 = SphereShell([1.1,0,0]);

			for i = 1:length(s1.nodeList)
				s1.nodeList(i).NewUpdateFrequency(updateFrequency);
				s1.nodeList(i).NewNodeNeighbourhoodRadius(dLim);
				s1.nodeList(i).NewFaceNeighbourhoodRadius(dLim);
			end

			for i = 1:length(s1.edgeList)
				s1.edgeList(i).NewUpdateFrequency(updateFrequency);
				s1.edgeList(i).NewEdgeNeighbourhoodRadius(dLim);
			end

			obj.AddNodesToSimulation(s1.nodeList);
			obj.AddEdgesToSimulation(s1.edgeList);
			obj.AddFacesToSimulation(s1.faceList);
			obj.AddCellsToSimulation(s1);

			s2 = SphereShell([-1.1,0,0]);

			for i = 1:length(s2.nodeList)
				s2.nodeList(i).NewUpdateFrequency(updateFrequency);
				s2.nodeList(i).NewNodeNeighbourhoodRadius(dLim);
				s2.nodeList(i).NewFaceNeighbourhoodRadius(dLim);
			end

			for i = 1:length(s2.edgeList)
				s2.edgeList(i).NewUpdateFrequency(updateFrequency);
				s2.edgeList(i).NewEdgeNeighbourhoodRadius(dLim);
			end

			obj.AddNodesToSimulation(s2.nodeList);
			obj.AddEdgesToSimulation(s2.edgeList);
			obj.AddFacesToSimulation(s2.faceList);
			obj.AddCellsToSimulation(s2);

			for i = 1:length(obj.edgeList)

				e = obj.edgeList(i);

				e.naturalLength = e.GetLength();

			end

			% The 0 is no attraction force
			obj.AddTissueBasedForce(  BruteForceSphere(0, springRate, dAsym, dSep, dLim)  );
			obj.AddEdgeBasedForce( EdgeLinearSpringForce(1) );
			obj.AddCellBasedForce( ConstantInternalPressure(.3) );

			obj.pathName = sprintf('SphereSurface/');
			obj.AddSimulationData(SpatialState());

			obj.AddDataWriter(WriteSpatialState(10, obj.pathName));


		end



	end


end