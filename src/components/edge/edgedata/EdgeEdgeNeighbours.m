classdef EdgeEdgeNeighbours < AbstractNodeData
	% Finds the edge neighbours to an edge

	% Only updates at specified intervals

	properties 

		name = 'edgeEdgeNeighbours'
		data = []
		radius
		steps

	end

	methods

		function obj = EdgeEdgeNeighbours(r, steps)

			% radius - the distance from the edge to store neighbours
			% steps - the number of timesteps to store the current neighbour list
			
			obj.radius = r;
			obj.steps = steps;
		end

		function CalculateData(obj, e1, t)
			
			% Since this is wrapped in GetData, this function will only execute one per time step
			% hence we don't need to worry about doing this repeatedly on the calculate step

			% We take the total step count modulo the storage timesteps
			% This requires t.steps to be 0 until the very end of the first
			% step, otherwise the neighbours won't be calculated at the beginning of the simulation

			if mod(t.step, obj.steps) == 0

				% Find all the neighbours within the radius
				neighbours = Edge.empty;

				for i = 1:length(t.edgeList)
					e2 = t.edgeList(i);

					% Compare each pair of edge, and decide
					% if any part of them is close enough to consider
					% in future calculations

					% To find the distance between two edges, we can
					% use the formula from https://en.wikipedia.org/wiki/Skew_lines#Distance_between_two_skew_lines

					v1 = e1.GetUnitVector1to2();
					v2 = e2.GetUnitVector1to2();

					n = cross(v1,v2);

					n = n / norm(n);

					d = norm(dot(n, (e1.Node1.pos - e2.Node1.pos)));

					if d < obj.radius
						 neighbours(end+1) = e2;
					end
				end

				obj.data = neighbours;

			end


		end
		
	end

end