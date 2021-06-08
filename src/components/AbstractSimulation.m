classdef (Abstract) AbstractSimulation < matlab.mixin.SetGet
	% A parent class that contains all the functions for running a simulation
	% The child/concrete class will only need a constructor that assembles the cells

	properties

		seed

		dt = 0.005  % Time step size
		t = 0		% Current time
		step = 0	% step count
        
		nodeList
		nextNodeId = 1

		edgeList
		nextEdgeId = 1

		surfList
		nextSurfId = 1

		cellList
		nextCellId = 1

		stochasticJiggle = true
		epsilon = 0.0001; % The size of the jiggle force

		% cellBasedForces AbstractCellBasedForce
		% edgeBasedForces AbstractEdgeBasedForce
		% surfBasedForces AbstractSurfBasedForce
		% neighbourhoodBasedForces AbstractNeighbourhoodBasedForce
		tissueBasedForces AbstractTissueBasedForce

		% stoppingConditions AbstractStoppingCondition

		stopped = false

		% tissueLevelKillers AbstractTissueLevelCellKiller
		% cellKillers AbstractCellKiller

		% simulationModifiers AbstractSimulationModifier

		% A collection of objects that store data over multiple time steps
		% with also the potential to write to file
		% dataStores AbstractDataStore

		% dataWriters AbstractDataWriter

		% A collection objects for calculating data about the simulation
		% stored in a map container so each type of data can be given a
		% meaingful name
		simData = containers.Map;

		% partition SpacePartition;

		usingPartition = false;

		writeToFile = true;
		
	end

	methods

		function SetRNGSeed(obj, seed)

			obj.seed = seed;
			rng(seed);

		end

		function NextTimeStep(obj)
			
			% First do the things to advance the simulation
			obj.Movement();

			% Need to implement these properly
			% obj.Expansion();

			% obj.Contraction();

			% Make any specific changes, i.e. boundary conditions
			% obj.Modify();

			obj.AdvanceAge();

			% obj.CollateData();

			% if obj.IsStoppingConditionMet()
			% 	obj.stopped = true;
			% end

		end

		function NTimeSteps(obj, n)
			% Advances a set number of time steps
			
			for i = 1:n
				% Do all the calculations
				obj.NextTimeStep();

				if mod(obj.step, 1000) == 0
					fprintf('Time = %3.3fhr\n',obj.t);
				end

				if obj.stopped
					fprintf('Stopping condition met at t=%3.3f\n', obj.t);
					break;
				end

			end
			
		end

		function RunToTime(obj, t)

			% Given a time, run the simulation until we reach said time
			if t > obj.t
				n = ceil((t-obj.t) / obj.dt);
				NTimeSteps(obj, n);
			end

		end

		function Movement(obj)

			% Function related to movement

			% obj.GenerateEdgeBasedForces();
			% obj.GenerateSurfaceBasedForces();
			% obj.GenerateCellBasedForces();
			obj.GenerateTissueBasedForces();

			% if obj.usingPartition
			% 	obj.GenerateNeighbourhoodBasedForces();
			% end

			% obj.GenerateForcesInTissues();
			% obj.GenerateNeighbourhoodBasedForces();
			% obj.GenerateOrganBasedForces();

			obj.ConvertTorquesToForcesGlobal();

			obj.MakeNodesMove();

		end

		%**********************************************************
		% Not checked in this prototype so far
		%**********************************************************
		function Expansion(obj)

			% Cells or other tissue components divide or
			% increase their number here

			newNodes 	= Node.empty();
			newEdges 	= Edge.empty();
			
			removedNodes = Node.empty();
			removedEdges = Edge.empty();
			

			for i = 1:length(obj.tissueList)
				
				[newNodeList, newEdgeList, removedNodeList, removedEdgeList] = obj.tissueList(i).TissueExpansion(obj);
				
				newNodes = [newNodes, newNodeList];
				newEdges = [newEdges, newEdgeList];
				
				removedNodes = [removedNodes, removedNodeList];
				removedEdges = [removedEdges, removedEdgeList];

			end

			obj.AddComponents(newNodes, newEdges);
			obj.DeleteComponents(removeNodes, removeEdges);

		end

		%**********************************************************
		% Not checked in this prototype so far
		%**********************************************************
		function Contraction(obj)

			% Cells or other tissue components are
			% removed from the simulation here

			newNodes 	= Node.empty();
			newEdges = Edge.empty();
			
			removedNodes 	= Node.empty();
			removedEdges = Edge.empty();
			

			for i = 1:length(obj.tissueList)
				
				[newNodeList, newEdgeList, removedNodeList, removedEdgeList] = obj.tissueList(i).TissueContraction(obj);
				
				newNodes = [newNodes, newNodeList];
				newEdges = [newEdges, newEdgeList];
				
				removedNodes = [removedNodes, removedNodeList];
				removedEdges = [removedEdges, removedEdgeList];

			end

			obj.AddComponents(newNodes, newEdges);
			obj.DeleteComponents(removeNodes, removeEdges);

		end

		%**********************************************************
		% Not checked in this prototype so far
		%**********************************************************
		function AdvanceAge(obj)

			obj.t = obj.t + obj.dt;
			obj.step = obj.step + 1;

			% for i = 1:length(obj.tissueList)
			% 	obj.tissueList(i).AdvanceAge(obj.dt);
			% end

		end

		function CollateData(obj)

			% All the processing to track and write data
			obj.StoreData();
			
			if obj.writeToFile
				obj.WriteData();
			end

		end

		function GenerateTissueBasedForces(obj)
			
			for i = 1:length(obj.tissueBasedForces)
				obj.tissueBasedForces(i).AddTissueBasedForces(obj);
			end

		end

		function GenerateCellBasedForces(obj)
			
			for i = 1:length(obj.cellBasedForces)
				obj.cellBasedForces(i).AddCellBasedForces(obj.cellList);
			end

		end

		function GenerateEdgeBasedForces(obj)

			for i = 1:length(obj.edgeBasedForces)
				obj.edgeBasedForces(i).AddEdgeBasedForces(obj.edgeList);
			end

		end

		function GenerateNeighbourhoodBasedForces(obj)
			
			if isempty(obj.partition)
				error('ACS:NoBoxes','Space partition required for NeighbourhoodForces, but none set');
			end

			for i = 1:length(obj.neighbourhoodBasedForces)
				obj.neighbourhoodBasedForces(i).AddNeighbourhoodBasedForces(obj.nodeList, obj.partition);
			end

		end

		function ConvertTorquesToForcesLocal(obj)

			% Loop through the edges and surfaces and
			% convert the torques to equivalent forces

			for i = 1:length(obj.surfList)

				% For each surface, use the rigid body properties to
				% find the equivalent force that can be applied to the
				% nodes around the boundary

				% Transform the torque into the coordinate system of the surface

				s = obj.surfList(i);

				t = s.torque;

				[u,v,w] = s.GetCoordinateSystem();

				tD = [dot(u,t),dot(v,t),0];

				if abs(dot(w,t)) > 1e-6
					error('Torque about the surface normal should be zero at this point')
				end

				[IDx,IDy,IDxy,IDz] = s.GetDragMatrixLocal();

				% We solve the components of the angular velocity in the body
				% system of coordinates
				omegaDx = (tD(1) * IDy + tD(2) * IDxy) / (IDx*IDy - IDxy^2);
				omegaDy = (tD(2) * IDx + tD(1) * IDxy) / (IDx*IDy - IDxy^2);

				% We now have the angular velocity vector in the body system
				% which also defines the rotation axis
				omegaD = [omegaDx,omegaDy,0];
				R = omegaD/norm(omegaD);

				% The change in angle is now
				a = obj.dt * norm(omegaD);

				% For a rotation about an arbitrary axis, the transformation matrix
				% is quite hefty, but reduced somewhat due to our choice of axes
				Rot = [cos(a) + R(1) * R(1) * (1 - cos(a)), R(1) * R(2) * (1 - cos(a)),  R(2) * sin(a); 
						R(1) * R(2) * (1 - cos(a)), cos(a) + R(2) * R(2) * (1 - cos(a)),  -R(1) * sin(a);
							-R(2) * sin(a), R(1) * sin(a), cos(a)];

				
				% Get the nodes in the local coordinates
				[r1, r2, r3] = s.GetNodeCoordinatesLocal();

				% Rotate the nodes in the local coordinates
				r1new = Rot * r1';
				r2new = Rot * r2';
				r3new = Rot * r3';

				% We now need to transform the local node coordinates into
				% the global coordinates
				T = [u;v;w];
				rD1new = T\r1new;
				rD2new = T\r2new;
				rD3new = T\r3new;

				% And finally shift back to correct origin
				rD = s.GetCentreOfDrag();

				rD1new = rD1new' + rD;
				rD2new = rD2new' + rD;
				rD3new = rD3new' + rD;

				% This gives us the point where the nodes will end up
				% after rotation. We can't rotate directly, so convert it
				% to an equivalent force.

				dr1 = rD1new - s.Node1.pos;
				dr2 = rD2new - s.Node2.pos;
				dr3 = rD3new - s.Node3.pos;

				Fe1 = s.Node1.eta * dr1 / obj.dt;
				Fe2 = s.Node2.eta * dr2 / obj.dt;
				Fe3 = s.Node3.eta * dr3 / obj.dt;

				% Finally, add these forces to the nodes directly

				s.Node1.ApplyForce(Fe1);
				s.Node2.ApplyForce(Fe2);
				s.Node3.ApplyForce(Fe3);

				% and reset the torque

				s.prevTorque = t;
				s.torque = [0,0,0];

			end

		end

		function ConvertTorquesToForcesGlobal(obj)

			% Does the convesion, except keeps everything
			% in the global axis, not the local one


			for i = 1:length(obj.surfList)

				% For each surface, use the rigid body properties to
				% find the equivalent force that can be applied to the
				% nodes around the boundary

				% Transform the torque into the coordinate system of the surface

				s = obj.surfList(i);

				t = s.torque;

				% At some point the thing fails if there is zero torque
				if sum(t) ~= 0

					ID = s.GetDragMatrix();

					omega = ID\t';

					R = omega / norm(omega);

					a = obj.dt * norm(omega);

					% The rotation matrix for a rotation of a about an axis R
					Rot = [cos(a) + R(1)^2 * (1 - cos(a)),  R(1) * R(2) * (1 - cos(a)) - R(3) * sin(a),  R(1) * R(3) * (1 - cos(a)) + R(2) * sin(a); 
							R(2) * R(1) * (1 - cos(a)) + R(3) * sin(a),  cos(a) + R(2)^2 * (1 - cos(a)),  R(2) * R(3) * (1 - cos(a)) - R(1) * sin(a); 
							R(3) * R(1) * (1 - cos(a)) - R(2) * sin(a),  R(3) * R(2) * (1 - cos(a)) + R(1) * sin(a),  cos(a) + R(3)^2 * (1 - cos(a))];


					% Get the nodes positions relative to the centre of drag
					[rD1, rD2, rD3] = s.GetNodeCoordinatesRelative();

					% Rotate the nodes about the axis, and translate
					rD = s.GetCentreOfDrag();
					rD1new = Rot * rD1' + rD';
					rD2new = Rot * rD2' + rD';
					rD3new = Rot * rD3' + rD';


					% This gives us the point where the nodes will end up
					% after rotation. We can't rotate directly, so convert it
					% to an equivalent force.

					dr1 = rD1new' - s.Node1.pos;
					dr2 = rD2new' - s.Node2.pos;
					dr3 = rD3new' - s.Node3.pos;

					Fe1 = s.Node1.eta * dr1 / obj.dt;
					Fe2 = s.Node2.eta * dr2 / obj.dt;
					Fe3 = s.Node3.eta * dr3 / obj.dt;

					% Finally, add these forces to the nodes directly

					s.Node1.ApplyForce(Fe1);
					s.Node2.ApplyForce(Fe2);
					s.Node3.ApplyForce(Fe3);

					% and reset the torque

					s.prevTorque = t;
					s.torque = [0,0,0];

				end

			end

		end

		function MakeNodesMove(obj)

			for i = 1:length(obj.nodeList)
				
				n = obj.nodeList(i);

				eta = n.eta;
				force = n.force;
				if obj.stochasticJiggle
					% Add in a tiny amount of stochasticity to the force calculation
					% to nudge it out of unstable equilibria

					% Make a random direction vector
					v = [rand-0.5,rand-0.5,rand-0.5];
					v = v/norm(v);

					% Add the random vector, and make sure it is orders of magnitude
					% smaller than the actual force
					force = force + v * obj.epsilon;

				end

				newPosition = n.pos + obj.dt/eta * force;

				n.MoveNode(newPosition);

				if obj.usingPartition
					obj.partition.UpdateBoxForNode(n);
				end

			end

		end

		
		%**********************************************************
		% Not checked in this prototype so far
		%**********************************************************
		function AdjustNodePosition(obj, n, newPos)

			% Only used by modifiers. Do not use to
			% progress the simulation

			% This will move a node to a given position regardless
			% of forces, but after all force movement has happened

			% Previous position and previous force are not modified

			% Make sure the node and edges are in the correct partition
			n.AdjustPosition(newPos);
			if obj.usingPartition
				obj.partition.UpdateBoxForNodeAdjusted(n);
			end

		end


		%**********************************************************
		% Not checked in this prototype so far
		%**********************************************************
		function MakeCellsDivide(obj)

			% Call the divide process, and update the lists
			newCells 	= AbstractCell.empty();
			newEdges = Edge.empty();
			newNodes 	= Node.empty();
			for i = 1:length(obj.cellList)
				c = obj.cellList(i);
				if c.IsReadyToDivide()
					[newCellList, newNodeList, newEdgeList] = c.Divide();
					newCells = [newCells, newCellList];
					newEdges = [newEdges, newEdgeList];
					newNodes = [newNodes, newNodeList];
				end
			end

			obj.AddNewCells(newCells, newEdges, newNodes);

		end


		%**********************************************************
		% Not checked in this prototype so far
		%**********************************************************
		function AddNewCells(obj, newCells, newEdges, newNodes)
			% When a cell divides, need to make sure the new cell object
			% as well as the new edges and nodes are correctly added to
			% their respective lists and partition if relevant

			for i = 1:length(newNodes)

				n = newNodes(i);
				n.id = obj.GetNextNodeId();
				if obj.usingPartition
					obj.partition.PutNodeInBox(n);
				end

			end

			for i = 1:length(newEdges)

				e = newEdges(i);
				e.id = obj.GetNextEdgeId();
				if obj.usingPartition && ~e.internal
					obj.partition.PutEdgeInBoxes(e);
				end

			end

			for i = 1:length(newCells)
				
				nc = newCells(i);
				nc.id = obj.GetNextCellId();

				if obj.usingPartition

					% Not necessary, as external edges never become internal edges
					% even considering cell death. As long as the internal flag is set
					% correctly during division, there is no issue
					% % If the cell type is joined, then we need to make sure the
					% % internal edge is labelled as such, and that the edge
					% % is removed from the partition.

					% if strcmp(class(nc), 'SquareCellJoined')

					% 	if ~nc.edgeRight.internal
					% 		nc.edgeRight.internal = true;
					% 		obj.partition.RemoveEdgeFromPartition(nc.edgeRight);
					% 	end

					% end

					% When a division occurs, the nodes and edges of the sister cell
					% (which was also the parent cell before division), may
					% have been modified to have a different node. This screws
					% with the space partition, so we have to fix it
					oc = nc.sisterCell;

					% Repair modified edges goes first because that adjusts nodes
					% in the function
					for j = 1:length(oc.edgeList)
						e = oc.edgeList(j);
						
						if e.modifiedInDivision
							obj.partition.RepairModifiedEdge(e);
						end

					end

					for j = 1:length(oc.nodeList)
						n = oc.nodeList(j);

						if n.nodeAdjusted
							obj.partition.UpdateBoxForNodeAdjusted(n);
						end

					end

				end

			end


			obj.cellList = [obj.cellList, newCells];

			obj.edgeList = [obj.edgeList, newEdges];

			obj.nodeList = [obj.nodeList, newNodes];

		end

		function MakeCellsAge(obj)

			for i = 1:length(obj.cellList)
				obj.cellList(i).AgeCell(obj.dt);
			end

		end

		function AddCellBasedForce(obj, f)

			if isempty(obj.cellBasedForces)
				obj.cellBasedForces = f;
			else
				obj.cellBasedForces(end + 1) = f;
			end

		end

		function AddEdgeBasedForce(obj, f)

			if isempty(obj.edgeBasedForces)
				obj.edgeBasedForces = f;
			else
				obj.edgeBasedForces(end + 1) = f;
			end

		end

		function AddSurfaceBasedForce(obj, f)

			if isempty(obj.surfaceBasedForces)
				obj.surfaceBasedForces = f;
			else
				obj.surfaceBasedForces(end + 1) = f;
			end

		end

		function AddNeighbourhoodBasedForce(obj, f)

			if isempty(obj.neighbourhoodBasedForces)
				obj.neighbourhoodBasedForces = f;
			else
				obj.neighbourhoodBasedForces(end + 1) = f;
			end

		end

		function AddTissueBasedForce(obj, f)

			if isempty(obj.tissueBasedForces)
				obj.tissueBasedForces = f;
			else
				obj.tissueBasedForces(end + 1) = f;
			end

		end

		function AddTissueLevelKiller(obj, k)

			if isempty(obj.tissueLevelKillers)
				obj.tissueLevelKillers = k;
			else
				obj.tissueLevelKillers(end + 1) = k;
			end

		end

		function AddCellKiller(obj, k)

			if isempty(obj.cellKillers)
				obj.cellKillers = k;
			else
				obj.cellKillers(end + 1) = k;
			end

		end

		function AddStoppingCondition(obj, s)

			if isempty(obj.stoppingConditions)
				obj.stoppingConditions = s;
			else
				obj.stoppingConditions(end + 1) = s;
			end

		end

		function AddSimulationModifier(obj, m)

			if isempty(obj.simulationModifiers)
				obj.simulationModifiers = m;
			else
				obj.simulationModifiers(end + 1) = m;
			end

		end

		function AddDataStore(obj, d)

			if isempty(obj.dataStores)
				obj.dataStores = d;
			else
				obj.dataStores(end + 1) = d;
			end

		end

		function AddDataWriter(obj, w)

			if isempty(obj.dataWriters)
				obj.dataWriters = w;
			else
				obj.dataWriters(end + 1) = w;
			end

		end

		function AddSimulationData(obj, d)

			% Add the simulation data calculator to the map
			% this will necessarily allow only one instance
			% of a given type of SimulationData, since the 
			% names are immutable

			% This is calculate-on-demand, so it does not have
			% an associated 'use' method here
			obj.simData(d.name) = d;

		end

		function StoreData(obj)

			for i = 1:length(obj.dataStores)
				obj.dataStores(i).StoreData(obj);
			end

		end

		function WriteData(obj)

			for i = 1:length(obj.dataWriters)
				obj.dataWriters(i).WriteData(obj);
			end

		end

		function ModifySimulationState(obj)

			for i = 1:length(obj.simulationModifiers)
				obj.simulationModifiers(i).ModifySimulation(obj);
			end

		end


		%**********************************************************
		% Not checked in this prototype so far
		%**********************************************************
		function KillCells(obj)

			% Loop through the cell killers

			% Currently the two killer types work differently.
			% This is for backward compatibility with a hack job that
			% I still need to work right now
			% Note to self, after all your work with DynamicLayer is done
			% take some time to fix this up for SquareCellJoined

			for i = 1:length(obj.tissueLevelKillers)
				obj.tissueLevelKillers(i).KillCells(obj);
			end

			killList = AbstractCell.empty();
			for i = 1:length(obj.cellKillers)
				killList = [killList, obj.cellKillers(i).MakeKillList(obj.cellList)];
			end

			obj.ProcessCellsToRemove(killList);

		end

		%**********************************************************
		% Not checked in this prototype so far
		%**********************************************************
		function ProcessCellsToRemove(obj, killList)

			% Loop from the end because we're removing cells from a list
			for i = length(killList):-1:1

				c = killList(i);

				% Do all the clean up

				% Need a separate way to handle joined cells

				% Clean up edges

				for j = 1:length(c.edgeList)

					obj.edgeList(obj.edgeList == c.edgeList(j)) = [];
					if obj.usingPartition
						obj.partition.RemoveEdgeFromPartition(c.edgeList(j));
					end

					c.edgeList(j).delete;

				end

				for j = 1:length(c.nodeList)

					obj.nodeList(obj.nodeList == c.nodeList(j)) = [];
					if obj.usingPartition
						obj.partition.RemoveNodeFromPartition(c.nodeList(j));
					end
					c.nodeList(j).delete;

				end

				% Clean up Cell

				% Since the cell List for the tissue is heterogeneous, we can't use
				% obj.cellList(obj.cellList == c) = []; to delete the cell because 
				% "one or more inputs of class 'AbstractCell' are heterogeneous
				% and 'eq' is not sealed". I have no idea what this means, but
				% it is a quirk of matlab OOP we have to work around
				for j = 1:length(obj.cellList)
					oc = obj.cellList(j);

					if oc == c
						obj.cellList(j) = [];
						break;
					end

				end
				

				c.delete;

			end

		end

		function stopped = IsStoppingConditionMet(obj)

			stopped = false;
			for i = 1:length(obj.stoppingConditions)
				if obj.stoppingConditions(i).CheckStoppingCondition(obj)
					stopped = true;
					break;
				end
			end

		end

		function Visualise(obj, varargin)


			h = figure();
			hold on

			% Intitialise the vector
			patchObjects(1) = fill([1,1],[2,2],'r');

			for i = 1:length(obj.surfList)
				s = obj.surfList(i);

				x = [s.nodeList.x];
				y = [s.nodeList.y];
				z = [s.nodeList.z];	

				patchObjects(i) = patch(x,y,z,[.5,.5,.5]);
			end

			axis equal

			drawnow

			view(3)

		end

		function Animate(obj, nsteps, sm, varargin)
			% Since we aren't storing data at this point, the only way to animate is to
			% calculate then plot

			
			xyrange = [];
			if ~isempty(varargin)
				xyrange = varargin{1};
			end

			h = figure();
			hold on
			axis equal
			view(3)

			if ~isempty(xyrange)
				xlim(xyrange(1:2));
				ylim(xyrange(3:4));
				zlim(xyrange(5:6));
			end

			patchObjects(1) = fill([1,1],[2,2],'r');
			lineObjects(1)  = line(nan,nan,nan,'Marker', 'o','MarkerSize',12);

			for i = 1:length(obj.surfList)
				s = obj.surfList(i);

				x = [s.nodeList.x];
				y = [s.nodeList.y];
				z = [s.nodeList.z];	

				patchObjects(i) = patch(x,y,z,[.5,.5,.5]);
			end

			for i = 1:length(obj.nodeList)
				n = obj.nodeList(i);

				lineObjects(i) = line(n.x,n.y,n.z,'Marker', 'o','MarkerSize',12);
			end

			totalSteps = 0;
			while totalSteps < nsteps

				obj.NTimeSteps(sm);
				totalSteps = totalSteps + sm;

				for j = 1:length(obj.surfList)
					s = obj.surfList(j);

					x = [s.nodeList.x];
					y = [s.nodeList.y];
					z = [s.nodeList.z];	

					if j > length(patchObjects)
						patchObjects(j) = patch(x,y,z,[.5,.5,.5]);
					else
						patchObjects(j).XData = x;
						patchObjects(j).YData = y;
						patchObjects(j).ZData = z;
						% patchObjects(j).FaceColor = c.GetColour();
					end
				end

				for j = 1:length(obj.nodeList)
					n = obj.nodeList(j);	

					if j > length(lineObjects)
						lineObjects(j) = line(n.x,n.y,n.z,'Marker', 'o','MarkerSize',12);
					else
						lineObjects(j).XData = n.x;
						lineObjects(j).YData = n.y;
						lineObjects(j).ZData = n.z;
						% patchObjects(j).FaceColor = c.GetColour();
					end
				end

				% Delete the line objects when there are too many
				for j = length(patchObjects):-1:length(obj.surfList)+1
					patchObjects(j).delete;
					patchObjects(j) = [];
				end

				for j = length(lineObjects):-1:length(obj.nodeList)+1
					lineObjects(j).delete;
					lineObjects(j) = [];
				end

				drawnow
				title(sprintf('t=%g',obj.t),'Interpreter', 'latex');

				pause(0.1);

			end

		end

		function AddNodesToSimulation(obj, listOfNodes)
			
			for i = 1:length(listOfNodes)
				% If any of the nodes are already in the list, don't add them
				n = listOfNodes(i);
				if sum(ismember(n, obj.nodeList)) == 0
					n.id = obj.GetNextNodeId();
					obj.nodeList = [obj.nodeList, n];
				end

			end

		end

		function AddEdgesToSimulation(obj, listOfEdges)
			
			for i = 1:length(listOfEdges)
				% If any of the Edges are already in the list, don't add them
				e = listOfEdges(i);
				if sum(ismember(listOfEdges(i), obj.edgeList)) == 0
					e.id = obj.GetNextEdgeId();
					obj.edgeList = [obj.edgeList, e];
				end

			end

		end

		function AddSurfacesToSimulation(obj, listOfSurf)
			
			for i = 1:length(listOfSurf)
				% If any of the Surf are already in the list, don't add them
				s = listOfSurf(i);
				if sum(ismember(listOfSurf(i), obj.surfList)) == 0
					s.id = obj.GetNextSurfId();
					obj.surfList = [obj.surfList, s];
				end

			end

		end

		function AddCellsToSimulation(obj, listOfCells)
			
			for i = 1:length(listOfCells)
				% If any of the Cells are already in the list, don't add them
				c = listOfCells(i);
				if sum(ismember(listOfCells(i), obj.cellList)) == 0
					c.id = obj.GetNextCellId();
					obj.cellList = [obj.cellList, c];
				end

			end

		end

		
	end

	methods (Access = protected)
		
		function id = GetNextNodeId(obj)
			
			id = obj.nextNodeId;
			obj.nextNodeId = obj.nextNodeId + 1;

		end

		function id = GetNextEdgeId(obj)
			
			id = obj.nextEdgeId;
			obj.nextEdgeId = obj.nextEdgeId + 1;

		end

		function id = GetNextSurfId(obj)
			
			id = obj.nextSurfId;
			obj.nextSurfId = obj.nextSurfId + 1;

		end

		function id = GetNextCellId(obj)
			
			id = obj.nextCellId;
			obj.nextCellId = obj.nextCellId + 1;

		end

	end

end