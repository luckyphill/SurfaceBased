classdef WriteVolumeByRadius < AbstractDataWriter

	properties

		% No special properties
		fileNames = {'VolumeByRadius'};

		subdirectoryStructure = ''
		
	end

	methods

		function obj = WriteVolumeByRadius(sm, simName)

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

			obj.data = {t.simData('volumeByRadius').GetData(t)};
			% obj.data = {length(t.cellList)};

		end
		
	end

end