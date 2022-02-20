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
		seed = 21:100;

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

			
			MakeParameterSet(obj);

			allcount = {};
			allmemvol = {};
			allvrpairs = {};
			alltime = {};

			for i = 1:size(obj.parameterSet,1)
				s = obj.parameterSet(i,:);

				radius = s(1);
				p = s(2);
				g = s(3);
				f = s(4);
				nns = s(5);
				nms = s(6);
				mep = s(7);


				allseedcount = [];
				allseedmemvol = [];
				allseedtime = [];
				vrpairs = {};

				for j = obj.seed

					a = ManageTumour3D(radius, p, g, f, nns, nms, mep, j);
					a.LoadSimulationData();

					count = a.data.cellCountData;
					memvol = a.data.membraneVolumeData;
					volrad = a.data.volumeByRadiusData;

					if isnan(count(1,1)) || isnan(memvol(1,1)) || isnan(volrad(1,1))

						obj.missingParameterSet(end + 1,:) =[radius, p, g, f, nns, nms, mep, j];

					end
					
					if ~isnan(count(1,1))

						% Take just the count numbers
						count = count(:,2)';
						allseedcount = Concatenate(obj, allseedcount, count, nan);

					end

					if ~isnan(memvol(1,1))

						% Just take the volume
						% we only have one membrane, but in general
						% this data type can have multiple
						memvol = memvol(:,2)';
						allseedmemvol = Concatenate(obj, allseedmemvol, memvol, nan);

					end

					if ~isnan(volrad(1,1))

						% A little more involved
						% For each time step, need to collate the volume/radius pairs
						% across each seed. Parameter sets are less important here
						% but we'll do all of them to see hwo they come out

						time = volrad(:,1);
						volrad = volrad(:,2:end);
						if length(vrpairs) < length(time)
							vrpairs{length(time)} = {};
							allseedtime = time;
						end

						for k = 1:length(time)
							vr = volrad(k,:);
							vr = reshape(vr,3,[])';

							for l = 1:size(vr,1)

								if ~isnan(vr(l,1))
									vrpairs{k} = [vrpairs{k};vr(l,2:3)];
								end

							end


						end

					end

				end

				fprintf('Completed %3.2f%%\n', 100*i/length(obj.parameterSet));


				allcount{i} = allseedcount;
				allmemvol{i} = allseedmemvol;
				allvrpairs{i} = vrpairs;
				alltime{i} = allseedtime;

			end


			obj.result = {alltime, allcount, allmemvol, allvrpairs};

			if ~isempty(obj.missingParameterSet)

				obj.ProduceMissingDataSimulationFiles();
			end
			

		end

		function PlotData(obj)

			alltime = obj.result{1};
			allcount = obj.result{2};
			allmemvol = obj.result{3};
			allvrpairs = obj.result{4};

			

			MakeParameterSet(obj);

			% Plot the cell number
			h = figure;
			ax = gca;
			hold on;
			leg = {};
			for i = 1:size(obj.parameterSet,1)

				time = alltime{i};
				count = allcount{i};

				meancount = nanmean(count);

				plot(ax, time, meancount,'LineWidth', 4);
				leg{i} = num2str(obj.parameterSet(i,7));

			end

			lgd = legend(leg, 'Location','northwest');
			ylabel('Cell Count','Interpreter', 'latex', 'FontSize', 15);
			xlabel('time (hr)','Interpreter', 'latex', 'FontSize', 15);
			title(sprintf('Cell count over time for different membrane strength'),'Interpreter', 'latex', 'FontSize', 22);
			ylim([0 600]);
			xlim([0 120]); % Chosen from the plot
				

			SavePlot(obj, h, sprintf('CountvsMembraneStrength'));


			% Plot the membrane volume
			h = figure;
			hold on;
			ax = gca;
			leg = {};
			for i = 1:size(obj.parameterSet,1)

				time = alltime{i};
				memvol = allmemvol{i};

				meanmemvol = nanmean(memvol);

				plot(ax, time, meanmemvol,'LineWidth', 4);

				leg{i} = num2str(obj.parameterSet(i,7));

			end

			
			lgd = legend(leg, 'Location','northwest');	
			ylabel('Membrane Volume','Interpreter', 'latex', 'FontSize', 15);
			xlabel('time (hr)','Interpreter', 'latex', 'FontSize', 15);
			title(sprintf('Membrane volume over time for different membrane strength'),'Interpreter', 'latex', 'FontSize', 22);
			ylim([0 30]);
			xlim([0 120]); % Chosen from the plot
				

			SavePlot(obj, h, sprintf('VolumevsMembraneStrength'));


			% For one parameter set, produce a plot of the 
			% vol vs radius for a set of time steps

			for i = 1:size(obj.parameterSet,1)

				vrpairs = allvrpairs{i};

				time = alltime{i};
				tsteps = 1000*[1:floor(length(time)/1000 -1)] % The time steps to use

				h = figure;
				hold on;
				ax = gca;

				bins = 0.1:0.1:1.7; % For the radius

				allmv = {};
				leg = {};
				dt = 0.001 * 10;
				for t = tsteps

					vr = vrpairs{t};

					v = vr(:,1);
					r = vr(:,2);

					m = [];
					for i = 1:length(bins)-1
						mv(i) = mean(v(  logical( (r>bins(i)) .* (r <= bins(i+1))  )   ) );
					end
					

					allmv{end+1} = mv;

					leg{end+1} = num2str(t*dt);

				end

				for i = 1:length(allmv)

					plot(ax, bins(2:end), allmv{i},'LineWidth', 4);

				end

				lgd = legend(leg, 'Location','northwest');	
				ylabel('Volume','Interpreter', 'latex', 'FontSize', 15);
				xlabel('Distance from centre','Interpreter', 'latex', 'FontSize', 15);
				title(sprintf('Cell volume as a function of distance from centre'),'Interpreter', 'latex', 'FontSize', 22);

				SavePlot(obj, h, sprintf('VolumevsRadius_%d', i));

			end

		end

	end

end