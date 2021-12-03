classdef  Tumour3D < AbstractSimulation


	properties

		% None yet

	end


	methods


		function obj = Tumour3D(radius, p, g, f, nns, nms, mep, seed)

			% Input variables:
			% radius: the initial radius of the membrane
			% p: The pause phase duration
			% g: The growth phase duration
			% f: The contact inhibition volume fraction
			% nns: node-to-node spring rate
			% nms: node-to-membrane spring rate
			% mep: The membrane edge parameter
			% seed: RNG seed

			% Each cell is represented by a single node
			% The membrane is represented by a sphere shell
			% made up of faces

			% The radius of a NodeCell when free floating in space
			freeRadius = 0.25;

			s = SphereShell([0,0,0],radius);

			obj.AddNodesToSimulation(s.nodeList);
			obj.AddEdgesToSimulation(s.edgeList);
			obj.AddFacesToSimulation(s.faceList);
			obj.AddCellsToSimulation(s);

			for i = 1:length(obj.edgeList)

				e = obj.edgeList(i);

				e.naturalLength = e.GetLength();

			end

			n = Node(0,0,0);
			ccm = NodeCellCycleSlowGrowth(p, g, f, freeRadius, obj.dt);
			c =  NodeCell(n, ccm);


			c.AddCellData(NodeCellVolume(obj, freeRadius));

			obj.AddCellsToSimulation(c);
			obj.AddNodesToSimulation(c.nodeList);

			%BruteForceNodeFace( spring rate, force asymptote, natural separation, limit of interactions) 
			obj.AddTissueBasedForce(  BruteForceNodeFace(nms, 0,     freeRadius, 1.05 * freeRadius)  );
			obj.AddTissueBasedForce(  BruteForceNodeNode(nns, 0, 2 * freeRadius,   3 * freeRadius)  );
			obj.AddEdgeBasedForce( EdgeLinearSpringForce(mep) );
			% obj.AddCellBasedForce( ConstantInternalPressure(.3) );

			pathName = sprintf('Tumour3D/r%gp%gg%gf%gnns%gnms%gmep%g_seed%g/',radius, p, g, f, nns, nms, mep, seed);
			obj.AddSimulationData(SpatialState());
			obj.AddDataWriter(WriteSpatialState(20, pathName));


		end

		function RunToConfluence(obj, tmin, tmax)

			% This function runs the simulation until all cells are stopped
			% by contact inhibition.
			% tmin: The simulation must run for this time before stopping is permitted
			% tmax: A catch all time limit in case the stopping condition is never reached

			% Use tmin > 0 if starting from a very small number of cells (say, 1)
			% to give the simulation a chance to get past conditions where all cells
			% stop at the same time briefly, before starting again


			obj.AddStoppingCondition(ConfluentStoppingCondition(tmin));

			obj.RunToTime(tmax);

		end



	end


end