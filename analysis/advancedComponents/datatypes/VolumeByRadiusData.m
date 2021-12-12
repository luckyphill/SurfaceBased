classdef VolumeByRadiusData < dataType
	% This grabs the membrane data

	properties (Constant = true)
		name = 'volumeByRadiusData';

		fileNames = 'VolumeByRadius'
	end

	methods

		function correct = verifyData(obj, data, sp)
			% All the check we're interested in to make sure the data is correct
			% Perhaps, check that there are sufficient time steps taken?

			% Check nothing yet
			correct = true;


		end

		function found = exists(obj, sp)
			% Checks if the file exists
			found = exist(obj.getFullFileName(sp), 'file');

		end
	end

	methods (Access = protected)

		function file = getFullFileName(obj,sp)
			
			file = [sp.saveLocation, obj.fileNames, '.csv'];

		end

		function data = retrieveData(obj, sp)
			% Loads the data from the file and puts it in the expected format
			% Need to read this is a format with nan padding

			% Setting the opts for readmatrix might not be necessary, but
			% doing it anyway because this thing is buggy as hell...
			% it does pad with nans though which is more useful than csvread
			opts = detectImportOptions(obj.getFullFileName(sp));
			opts.DataLines = [1 Inf];
			if strcmp(opts.VariableTypes{1}, 'char')
				opts = setvartype(opts, opts.VariableNames{1}, 'double');
			end
			data = readmatrix(obj.getFullFileName(sp),opts);

		end

		function processOutput(obj, sp)
			
			% Do nothing, simulation already puts it in the right spot

		end

	end

end