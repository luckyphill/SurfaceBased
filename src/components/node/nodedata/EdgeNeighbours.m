classdef EdgeNeighbours < AbstractNodeData
	% Calculates the wiggle ratio

	properties 

		name = 'edgeNeighbours'
		data = []
		radius
		steps

	end

	methods

		function obj = EdgeNeighbours(r, steps)

			% radius - the distance from the node to store neighbours
			% steps - the number of timesteps to store the current neighbour list
			
			obj.radius = r;
			obj.steps = steps;
		end

		function CalculateData(obj, n, t)
			
			% Since this is wrapped in GetData, this function will only execute one per time step
			% hence we don't need to worry about doing this repeatedly on the calculate step

			% We take the total step count modulo the storage timesteps
			% This requires t.steps to be 0 until the very end of the first
			% step, otherwise the neighbours won't be calculated at the beginning of the simulation

			if mod(t.step, obj.steps) == 0

				% Find all the neighbours within the radius
				neighbours = Edge.empty;

				for i = 1:length(t.edgeList)
					eCheck = t.edgeList(i);

					if ~ismember(n, eCheck.nodeList)

						% To check the distance between a node and an edge in 3d
						% it is sufficient to take the dot product between the
						% edge unit vector and the vector from the node to one of
						% the end points
						nToE = n.pos - eCheck.Node1.pos;
						eUnit = eCheck.GetUnitVector1to2();
						u = dot(nToE, eUnit);

						v = sqrt(eUnit.GetLength()^2 - u^2);

						if v < obj.radius
							neighbours = [neighbours, eCheck];
						end
					end

				end

				obj.data = neighbours;

			end


		end
		
	end

end