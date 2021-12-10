classdef Tumour3DGrowth < Analysis

	% This analysis compares buckling due to weak membrane adhesion,
	% buckling due to weak membrane tension

	properties

		% STATIC: DO NOT CHANGE
		% IF CHANGE IS NEEDED, MAKE A NEW OBJECT
		% Tumour3D(radius, p, g, f, nns, nms, mep, seed)
		radius = 1;

		p = 3;
		g = 3;

		f = 0.75;

		nns = 1;
		nms = 4;

		mep = 1:10;
		seed = 1:10;

		minimumTime = 50;
		maximumTime = 300;

		analysisName = 'Tumour3DGrowth';

		parameterSet = []
		missingParameterSet = []

		simulationRuns = 50
		slurmTimeNeeded = 72
		simulationDriverName = 'ManageTumour3D'
		simulationInputCount = 7
		

	end

	methods

		function obj = Tumour3DGrowth()

			obj.seedIsInParameterSet = false; % The seed not given in MakeParameterSet, it is set in properties
			obj.seedHandledByScript = false; % The seed will be in the parameter file, not the job script
			obj.usingHPC = true;

		end

		function MakeParameterSet(obj)

			% n, p, g, b, f, sae, spe

			params = [];

			for radius = obj.radius
				for p = obj.p
					for g = obj.g
						for mep = obj.mep
							for f = obj.f
								for nns = obj.nns
									for nms = obj.nms

										params(end+1,:) = [radius,p,g,f,nns,nms,mep];

									end
								end
							end
						end
					end
				end
			end

			

			obj.parameterSet = params;

		end

		

		function BuildSimulation(obj)

			obj.MakeParameterSet();
			obj.ProduceSimulationFiles();
			
		end

		function AssembleData(obj)

			
			obj.result = {};

			obj.missingParameterSet = [];

			if ~isempty(obj.missingParameterSet)

				obj.ProduceMissingDataSimulationFiles();
			end
			

		end

		function PlotData(obj)

			xxx = obj.result{1};
			yyy = obj.result{2};

			h = figure;

			ylabel('yyy','Interpreter', 'latex', 'FontSize', 15);
			xlabel('xxx','Interpreter', 'latex', 'FontSize', 15);
			title(sprintf('ttt, p=%g, g=%g',obj.p,obj.g),'Interpreter', 'latex', 'FontSize', 22);
			ylim([min(obj.b)-1, max(obj.b)+1]);
			xlim([min(obj.spe)-1, max(obj.spe)+1]);
			

			SavePlot(obj, h, sprintf('xxx'));


		end

	end

end