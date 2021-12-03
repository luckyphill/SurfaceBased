classdef ConfluentStoppingCondition < AbstractStoppingCondition
	% If all of the cells are in a state of contact inhibition
	% stop the simulation

	properties

		name = 'Confluent'

		minimumTime = 0

	end

	methods

		function obj = ConfluentStoppingCondition(varargin)

			% No special initialisation
			if ~isempty(varargin)
				obj.minimumTime = varargin{1};
			end

		end

		function stopped = HasStoppingConditionBeenMet(obj, t)

			stopped = true;


			if t.t > obj.minimumTime
				% If the simulation has run for a minimum specified time,
				% then it can only continue if all cells are not stopped

				% Loop through the node cell population, and if any single cell is
				% not in the stopped condition, then the simulation continues.
				for i = 1:length(t.cellList)
					c = t.cellList(i);
					if isa(c, 'NodeCell')

						if (c.CellCycleModel.colour ~= c.CellCycleModel.colourSet.GetNumber('STOPPED'))
							stopped = false;
							break;
						end

					end

				end

			else
				% Otherwise we carry on until the minimum time is reached
				stopped = false;
			end

		end

	end



end