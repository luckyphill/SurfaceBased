classdef ManageTumour3D < MatlabSimulation
	% Manages the running and data handling of the Tumour3D simulation
	% for the various analyses around buckling


	properties (SetAccess = immutable)

		matlabTest = 'Tumour3D'

	end

	properties (SetAccess = private)

		% Number of cells
		radius 			double {mustBeNonnegative}

		% Cell cycle phase lengths
		p 			double {mustBeNonnegative}
		g 			double {mustBeNonnegative}

		% growthTriggerFraction
		f			double {mustBeNonnegative}

		% Force parameters
		nns			double {mustBeNonnegative}
		nms			double {mustBeNonnegative}
		mep			double {mustBeNonnegative}

		% These are solver parameters
		t 			double {mustBeNonnegative}

		% This is the RNG seed parameter
		rngSeed 	double {mustBeNumeric}

		outputLocation

		% Should really treat these as input parameters, but they won't
		% likely change, and they won't affect the simulation, maybe
		% just the amount of data produced (increase it)
		% They are just needed for the simulation command
		minimumTime
		maximumTime


	end

	properties

		% Nothing yet

	end

	methods
		function obj = ManageTumour3D(radius, p, g, f, nns, nms, mep, seed)

			obj.radius 	= radius;
			obj.p 		= p;
			obj.g 		= g;
			obj.f 		= f;
			obj.nns 	= nns;
			obj.nms		= nms;
			obj.mep		= mep;
			obj.rngSeed = seed;

			obj.minimumTime = 50;
			obj.maximumTime = 200;

			
			obj.simObj = Tumour3D(radius, p, g, f, nns, nms, mep, seed);

			obj.GenerateSaveLocation();

			% Remove the default spatial state output
			remove(obj.simObj.simData,'spatialState');
			obj.simObj.dataWriters = AbstractDataWriter.empty();

			obj.outputTypes = {CellCountData, MembraneVolumeData};

			obj.GenerateSaveLocation();

			% Normally get a warning if the data doesn't exist, so we'll
			% block it for now to save blowing out the command window
			% when there are lots of data locations to check
			warning('off','sim:LoadFailure')
			obj.LoadSimulationData();
			warning('on','sim:LoadFailure')

			% Only add the data types that are missing
			if isnan(obj.data.cellCountData)
				obj.simObj.AddSimulationData(NodeCellCount());
				obj.simObj.AddDataWriter(WriteCellCount(20,obj.simObj.pathName));
			end

			if isnan(obj.data.membraneVolumeData)
				obj.simObj.AddSimulationData(MembraneVolume());
				obj.simObj.AddDataWriter(WriteMembraneVolume(20,obj.simObj.pathName));
			end

		end

	end


	methods (Access = protected)
		% Helper methods to build the class


		function GenerateSaveLocation(obj)
			% This generates the full path to the specific data file for the simulation
			% If the path doesn't exist it creates the missing folder structure

			obj.saveLocation = obj.simObj.simulationOutputLocation;

		end

		function SimulationCommand(obj)

			obj.simObj.RunToConfluence(obj.minimumTime, obj.maximumTime);

		end

	end

end








