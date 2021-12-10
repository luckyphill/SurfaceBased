classdef WriteCellCount < AbstractDataWriter

	properties

		% No special properties
		fileNames = {'CellCount'};

		subdirectoryStructure = ''
		
	end

	methods

		function obj = WriteCellCount(sm, simName)

			obj.subdirectoryStructure = simName;
			obj.samplingMultiple = sm;
			obj.multipleFiles = false;
			obj.timeStampNeeded = true;
			obj.data = {};

		end

		function GatherData(obj, t)

			% Using the cellCount simdata object because
			% at this point we only want to count NodeCells
			% not the SphereShell

			obj.data = {t.simData('cellCount').GetData(t)};
			% obj.data = {length(t.cellList)};

		end
		
	end

end