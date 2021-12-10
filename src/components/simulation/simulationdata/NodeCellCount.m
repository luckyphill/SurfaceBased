classdef NodeCellCount < AbstractSimulationData
	% Gets the number of cells

	properties 

		name = 'cellCount'
		data = []
	end

	methods

		function obj = NodeCellCount
			% No special initialisation
		end

		function CalculateData(obj, t)

			count = 0;
			for i = 1:length(t.cellList)
				c = t.cellList(i);
				if isa(c, 'NodeCell')
					count = count + 1;
				end
			end

			obj.data = count;


		end
		
	end


end