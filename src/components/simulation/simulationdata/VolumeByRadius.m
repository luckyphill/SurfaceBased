classdef VolumeByRadius < AbstractSimulationData

	% Collects the volume of non-growing every cell in the
	% simulation and stores its distance from the tumour centre
	properties 

		name = 'volumeByRadius'
		data = {};

	end

	methods

		function obj = VolumeByRadius
			% No special initialisation
			
		end

		function CalculateData(obj, t)

			% This is intended for a Tumour3D simulation
			% At each time step it finds the centre of the tumour
			% It then finds the volume and distance from the centre
			% for each NodeCell.
			% This is simple for each non-growing NodeCell,
			% but when it is growing, it is actually made up of two spearate
			% NodeCells. Radius is relatively simple to define for this case
			% (just use the centre point of the two cells instead of the 
			% node position), but it is far less clear what the volume will be.
			% For the sake of getting it done, growing cells will be ignored
			% This will leave a lack of data, particularly at the early stages
			% so averaging over many simulations will be vital

			volData = [];

			% A cheeky hack.
			% In a Tumour3D simulation, the first cell is actually the membrane
			% and the first nodes in nodeList are the membrane nodes.
			% If we start the list here after the membrane nodes then the centre calculation
			% only accounts for NodeCell nodes
			n = length(t.cellList(1).nodeList) + 1;
			centre = mean(reshape([t.nodeList(n:end).pos],3,[])');

			for i = 1:length(t.cellList)

				c = t.cellList(i);
				if isa(c,'NodeCell')
					colour = c.CellCycleModel.colour;
					if colour ~= 2 % Not the best way to do it, but will skip this section if colour matches GROW phase
						vol = c.GetVolume();
						pos = c.nodeList.pos;
						radius = norm(pos-centre);
						volData(end+1,:) = [c.id, vol, radius];
					end
				end

			end

			obj.data = volData;

		end
		
	end

end