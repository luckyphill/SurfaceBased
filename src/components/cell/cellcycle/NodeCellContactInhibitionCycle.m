classdef NodeCellContactInhibitionCycle < AbstractCellCycleModel
	% A cell cycle for node cells

	properties

		cycleLength

		fraction

		freeVolume = (4/3) * pi * 0.5^3

		pauseColour
		growthColour
		inhibitedColour
	end

	methods

		function obj = NodeCellContactInhibitionCycle(cL, f)

			obj.SetAge(cL*rand);
			obj.cycleLength = cL;
			obj.fraction = f;

			obj.colour = obj.colourSet.GetNumber('GROW');
			obj.pauseColour = obj.colourSet.GetNumber('PAUSE');
			obj.growthColour = obj.colourSet.GetNumber('GROW');
			obj.inhibitedColour = obj.colourSet.GetNumber('STOPPED');

		end

		function newCCM = Duplicate(obj)

			newCCM = NodeCellContactInhibitionCycle(obj.cycleLength, obj.fraction);
			obj.SetAge(0);
			newCCM.SetAge(0);

		end

		function ready = IsReadyToDivide(obj)

			ready = false;
			obj.colour = obj.growthColour;
			if obj.age > obj.cycleLength
				c = obj.containingCell;
				vol = c.cellData('cellVolume').GetData(c);
				if vol > obj.fraction * obj.freeVolume
					ready = true;
				else
					obj.colour = obj.inhibitedColour;
				end
			end
		end

		function fraction = GetGrowthPhaseFraction(obj)

			fraction = 1;
		end

	end


end