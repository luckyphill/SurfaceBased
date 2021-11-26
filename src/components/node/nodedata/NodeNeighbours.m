classdef NodeNeighbours < AbstractNodeData
	% Calculates the wiggle ratio

	properties 

		name = 'nodeNeighbours'
		data = []
		radius
		steps

	end

	methods

		function obj = NodeNeighbours(r, steps)

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
				neighbours = Node.empty;

				for i = 1:length(t.nodeList)

					nCheck = t.nodeList(i);

					d = norm(n.pos - nCheck.pos);

					if d < obj.radius && d > 0
						neighbours = [neighbours, nCheck];
					end

				end

				obj.data = neighbours;

			end


		end
		
	end

end