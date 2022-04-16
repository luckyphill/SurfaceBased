classdef FaceNeighbours < AbstractNodeData
	% Finds all the faces near to a node

	properties 

		name = 'faceNeighbours'
		data = []
		radius
		steps

	end

	methods

		function obj = FaceNeighbours(r, steps)

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
				neighbours = Face.empty;


				for i = 1:length(t.faceList)
					fCheck = t.faceList(i);

					% To check the distance between a node and a face in 3d
					% we need the unit normal of the face, and the distance
					% between a node and some point on the surface
					% The dot product between these gives the distance

					if ~ismember(n, fCheck.nodeList)

						nToF = n.pos - fCheck.Node1.pos;
						fUnit = fCheck.GetUnitNormal();
						
						u = dot(nToF, fUnit);

						if abs(u) < obj.radius
							neighbours = [neighbours, fCheck];
						end
					end

				end

				obj.data = neighbours;

			end


		end
		
	end

end