classdef WriteMembraneVolume < AbstractDataWriter

	properties

		% No special properties
		fileNames = {'MembraneVolume'};

		subdirectoryStructure = ''
		
	end

	methods

		function obj = WriteMembraneVolume(sm, simName)

			obj.subdirectoryStructure = simName;
			obj.samplingMultiple = sm;
			obj.multipleFiles = false;
			obj.timeStampNeeded = true;
			obj.data = {};

		end

		function GatherData(obj, t)

			% The simulation t must have a simulation data object
			% collating the complete spatial state

			obj.data = {t.simData('membraneVolume').GetData(t)};

		end
		
	end

end