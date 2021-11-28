classdef Visualiser < matlab.mixin.SetGet

	properties

		pathToSpatialState
		pathToOutput
		nodes
		edges
		faces
		cells
		nodeCells
		faceCells

		timeSteps

		cs = ColourSet()
		
	end

	methods

		function obj = Visualiser(v)

			% v can be either of two types:
			% a string that gives the subdirectory structure from
			% the folder SimulationOutput to the folder SpatialState for the desired simulation
			% i.e. [path]/[to]/SimulationOutput/[v]/SpatialState
			% OR
			% the simulation object handle

			if isa(v, 'AbstractCellSimulation')
				
				obj.pathToSpatialState = v.dataWriters(1).fullPath;
				rootDir = v.dataWriters(1).rootStorageLocation;

				subDir = erase(obj.pathToSpatialState, rootDir);
				subDir = subDir(1:end-13); % remove SpatialState/ from subDir

				obj.pathToOutput = [getenv('FACEDIR'),'/Images/',subDir];


			else

				% If the input is not a simulation object, then assume its a string

				if ~strcmp(v(end),'/')
					v(end+1) = '/';
				end

				v = [v, 'SpatialState/'];

				obj.pathToSpatialState = [getenv('FACEDIR'),'/SimulationOutput/',v];

				obj.pathToOutput = [getenv('FACEDIR'),'/Images/',v];

				if ~strcmp( obj.pathToOutput(end),'/' )
					obj.pathToOutput(end+1) = '/';
				end

			end

			if exist(obj.pathToOutput,'dir')~=7
				mkdir(obj.pathToOutput);
			end



			obj.LoadData();

		end

		function LoadData(obj)

			% For some reason matlab decides to ignore some lines
			% when using readmatrix, so to stop this, we need to pass in special options
			% See https://stackoverflow.com/questions/62399666/why-does-readmatrix-in-matlab-skip-the-first-n-lines?
			
			% Despite this, readmatrix is still valuable because when it feeds the file
			% into a matrix, any empty entries are filled with nans. On the other hand
			% dlmread, or csvread fill empty entries with zeros, which unfortunately are
			% also valid spatial positions, so it is difficult to reliably distinguish
			% empty values from valid zeros. To throw a complete spanner in the works
			% readmatrix can't handle exceptionally large data files (above something like
			% 100MB or 200MB), so we have to revert to dlmread in this case and just accept
			% that in certain cases cells will disappear in the visualisation because one of their
			% coordinates gets convereted to nan. Fortunately, this should be a fleeting
			% extraordinarily rare event in general. It will however be be common if
			% boundary conditions are set on the x or y axes

			% -----------------------------------------------------------------------
			% Read in the nodes file
			% in some simulations (i.e. those with fixed boundaried) a position of 0
			% will be a valid entry, so need to switch back to readmatrix
			opts = detectImportOptions([obj.pathToSpatialState, 'nodes.csv']);
			opts.DataLines = [1 Inf];
			if strcmp(opts.VariableTypes{1}, 'char')
				opts = setvartype(opts, opts.VariableNames{1}, 'double');
			end
			% nodeData = readmatrix([obj.pathToSpatialState, 'nodes.csv'],opts);
			nodeData = dlmread([obj.pathToSpatialState, 'nodes.csv']);
			nodeData(nodeData == 0) = nan;

			% -----------------------------------------------------------------------
			% Read in the edges file
			opts = detectImportOptions([obj.pathToSpatialState, 'edges.csv']);
			opts.DataLines = [1 Inf];
			if strcmp(opts.VariableTypes{1}, 'char')
				opts = setvartype(opts, opts.VariableNames{1}, 'double');
			end
			edgeData = readmatrix([obj.pathToSpatialState, 'edges.csv'],opts);


			% -----------------------------------------------------------------------
			% Read in the faces file
			opts = detectImportOptions([obj.pathToSpatialState, 'faces.csv']);
			opts.DataLines = [1 Inf];
			if strcmp(opts.VariableTypes{1}, 'char')
				opts = setvartype(opts, opts.VariableNames{1}, 'double');
			end
			faceData = readmatrix([obj.pathToSpatialState, 'faces.csv'],opts);


			% -----------------------------------------------------------------------
			% Read in the cells file
			opts = detectImportOptions([obj.pathToSpatialState, 'cells.csv']);
			opts.DataLines = [1 Inf];
			if strcmp(opts.VariableTypes{1}, 'char')
				opts = setvartype(opts, opts.VariableNames{1}, 'double');
			end
			cellData = readmatrix([obj.pathToSpatialState, 'cells.csv'],opts);
			% cellData = csvread([obj.pathToSpatialState, 'cells.csv']);


			% -----------------------------------------------------------------------
			% Modify the structure of the data so it can be matched properly

			% First remove timestamps
			obj.timeSteps = nodeData(:,1);
			nodeData = nodeData(:,2:end);
			edgeData = edgeData(:,2:end);
			faceData = faceData(:,2:end);
			cellData = cellData(:,2:end);


			% -----------------------------------------------------------------------
			% We will first build a data structure that contains all the
			% spatial positions of the nodes, and design it so we can
			% reference it by node ID and timestep to find the position
			% of a given node at a given time.

			% The exisiting data structure has ID1,x1,y1,z1,ID2,x2,y2,z2,ID3,x3,y3,z3....
			% on each row, where each row matches a time step

			% The size of the array needs to be largest ID x # timesteps x number of spatial dimensions

			
			nDim = 3;
			[tSteps,~] = size(nodeData);
			maxID = max(max(nodeData(:,1:nDim+1:end)));
			
			

			nodes = nan(maxID,tSteps,nDim);

			for i = 1:tSteps
				nD  = nodeData(i,:);
				nD = reshape(nD,nDim+1,[])';
				% First column is ID, then x and y
				
				% For each node, use the id as the first index,
				% and the second index is the time step. In that
				% position is stored the (x,y) coords
				[mnD, ~] = size(nD);

				for j = 1:mnD
					n = nD(j,:);
					if ~isnan(n(1))
						nodes(n(1),i,:) = [n(2), n(3), n(4)];
					end

				end

			end

			% Store it in the object
			obj.nodes = nodes;


			% -----------------------------------------------------------------------
			% Now make an array of edges
			% 
			obj.edges = permute(reshape(edgeData,tSteps,2,[]),[1,3,2]);


			% -----------------------------------------------------------------------
			% Now make an array of faces
			% 
			obj.faces = permute(reshape(faceData,tSteps,3,[]),[1,3,2]);
			
			
			% -----------------------------------------------------------------------
			% Now make an array of cells

			% Each row in the matrix lists the nodes for each cell. The first number is the number
			% of nodes in the cell, call it jump, then the nodes for the cell are listed, followed by
			% the cell colour
			% [m,~] = size(cellData);
			cells = {};
			for i = 1:tSteps
				a = cellData(i,:);
				j = 1;
				counter = 1;
				while j <= length(a) && ~isnan(a(j))
					jump = a(j);
					cells{i,counter} = a(j+1:j+jump+1);
					j = j + jump + 2;
					counter = counter + 1;
				end

			end

			% This provides a data structure that holds the nodes for each cell, at each timestep
			% To display, we need to break this up by NodeCells, and cells represented by surfaces
			% It is conceiveable that in the future we may allow cells to be made up of connected
			% edges in a line, but for now we assume that if a cell has more than 2 nodes, it is
			% made of faces

			% It is far more convenient to collect the faces of "face cells" rather than their
			% nodes, so we need to convert the data structure into two separate ones

			nodeCells = {};
			faceCells = {{}};

			for i = 1:tSteps

				cellsAtTime = {cells{i,:}};
				% nodeCellsRow = {};
				% faceCellsRow = {};
				nodesCounter = 1;
				facesCounter = 1;

				for j = 1:length(cellsAtTime)

					if length(cellsAtTime{j}) > 2 % Each cell has a node and a colour
						% This is a face cell
						faceCells{i,facesCounter} = cellsAtTime{j};
						facesCounter = facesCounter + 1;
					else
						% This is a node cell
						if ~isempty(cellsAtTime{j})
							nodeCells{i,nodesCounter} = cellsAtTime{j};
							nodesCounter = nodesCounter + 1;
						end

					end
				end

				% nodeCells{i} = nodeCellsRow;
				% faceCells{i} = faceCellsRow;

			end

			% To make displaying the faces easier, we need to convert the collection of nodes
			% for the whole cell, into a collection of faces
			% We don't need to be particularly efficient here since this is done once during loading

			% For each faceCell, find the ID of each face that is in the cell by comparing
			% the nodes in each each face to the nodes in each cell

			% For now, we can cheat because in Tumour3D there is only one cell with faces


			% THIS NEEDS RETHINKING< IT DOESNT WORK AT ALL
			% CURRENTLY CHEATING AND USING obj.faces DIRECTLY IN THE ANIMATION PART
			for i = 1:length(faceCells)
				faceCells{i} = [obj.faces(1,:),11]; % The faces and the face colour
			end

			obj.cells = cells;

			obj.nodeCells = nodeCells;
			obj.faceCells = faceCells;


		end


		% Really need to abstract this so the code isn't copied over and over
		function VisualiseCells(obj, varargin)

			% varargin 
			% Arg 1: [indexStart, indexEnd] - a vector of the start and ending indices. Leave empty to run the whole simulation
			% Arg 2: plot axis range in the form [xmin,xmax,ymin,ymax]

			xyrange = [];
			indices = [];
			if ~isempty(varargin)
				indices = varargin{1};
				if length(varargin) > 1
					xyzrange = varargin{2};
				end
			end

			h = figure();
			axis equal
			hold on
			ax = gca;
			if ~isempty(xyrange)
				xlim(xyzrange(1:2));
				ylim(xyzrange(3:4));
				zlim(xyzrange(5:6));
			end

			[I,~] = size(obj.cells);



			% A radius for node cells
			r = 0.25;
			% Unit sphere coords to plot node cells
			[usX,usY,usZ] = sphere(100);
			
			% Intitialise the vectors
			patchObjects(1) = patch(ax, 1,1,1, [1,1,1], 'FaceAlpha', 0.5, 'EdgeColor', [.5,.5,.5]);
			surfObjects(1)  = surf(usX,usY,nan(size(usZ)));


			% Extract the start and end points
			startI =  1;
			endI = I;
			if ~isempty(indices)
				startI = indices(1);
				endI = indices(2);
			end

			for i = startI:endI

				% Split this by updating nodeCells
				% and faceCells


				% Loop through the node cells
				[~,J] = size(obj.nodeCells);
				j = 1;
				while j <= J && ~isempty(obj.nodeCells{i,j})

					c = obj.nodeCells{i,j};
					ids = c(1:end-1);
					colour = c(end);
					nodeCoords = squeeze(obj.nodes(ids,i,:));

					x = nodeCoords(1) + r*usX;
					y = nodeCoords(2) + r*usY;
					z = nodeCoords(3) + r*usZ;

					if j > length(surfObjects)
						surfObjects(j) = surf(ax, x, y, z, 'LineStyle', 'none', 'FaceColor', obj.cs.GetRGB(colour));
					else
						surfObjects(j).XData = x;
						surfObjects(j).YData = y;
						surfObjects(j).ZData = z;
						surfObjects(j).FaceColor = obj.cs.GetRGB(colour);
						surfObjects(j).LineStyle = 'none';
					end

					j = j + 1;

				end
				% j will always end up being 1 more than the total number of non empty cells

				for k = length(surfObjects):-1:j
					surfObjects(k).delete;
					surfObjects(k) = [];
				end



				% This is a cheat way to get animations working with Tumour3D
				% all faces are displayed the same, and the number of faces never changes
				% Loop through the face cells

				[~,J,~] = size(obj.faces);
				for j = 1:J

					nIDs = obj.faces(i,j,:);
					colour = 11;
					nodeCoords = squeeze(obj.nodes(nIDs,i,:));

					x = nodeCoords(:,1);
					y = nodeCoords(:,2);
					z = nodeCoords(:,3);

					if j > length(patchObjects)
						patchObjects(j) = patch(ax, x, y, z, obj.cs.GetRGB(colour), 'FaceAlpha', 0.5, 'EdgeColor', [.5,.5,.5]);
					else
						patchObjects(j).XData = x;
						patchObjects(j).YData = y;
						patchObjects(j).ZData = z;
						patchObjects(j).FaceColor = obj.cs.GetRGB(colour);
					end

					j = j + 1;

				end
				% j will always end up being 1 more than the total number of non empty cells

				for k = length(patchObjects):-1:j
					patchObjects(k).delete;
					patchObjects(k) = [];
				end

				drawnow
				title(sprintf('t = %g',obj.timeSteps(i)),'Interpreter', 'latex');
				pause(0.1);

			end

		end

		function ProduceMovie(obj, varargin)

			% varargin 
			% Arg 1: [indexStart, indexEnd] - a vector of the start and ending indices
			% Arg 2: plot axis range in the form [xmin,xmax,ymin,ymax]
			
			xyzrange = [];
			indices = [];
			videoFormat = 'MPEG-4';
			if ~isempty(varargin)
				indices = varargin{1};
				if length(varargin) > 1
					xyzrange = varargin{2};
					if length(varargin) > 2
						videoFormat = varargin{3};
					end
				end
			end



			% Currently same as run visualiser, but saves the movie

			h = figure();
			axis equal
			axis off
			hold on
			ax = gca;

			F = getframe(h);
			
			if ~isempty(xyzrange)
				xlim(xyzrange(1:2));
				ylim(xyzrange(3:4));
				zlim(xyzrange(5:6));
			end

			% A radius for node cells
			r = 0.25;
			% Unit sphere coords to plot node cells
			[usX,usY,usZ] = sphere(100);
			
			% Intitialise the vectors
			patchObjects(1) = patch(ax, 1,1,1, [1,1,1], 'FaceAlpha', 0.5, 'EdgeColor', [.5,.5,.5]);
			surfObjects(1)  = surf(usX,usY,nan(size(usZ)));


			% Extract the start and end points
			[I,~] = size(obj.cells);
			startI =  1;
			endI = I;
			if ~isempty(indices)
				startI = indices(1);
				endI = indices(2);
			end

			for i = startI:endI

				% Split this by updating nodeCells
				% and faceCells


				% Loop through the node cells
				[~,J] = size(obj.nodeCells);
				j = 1;
				while j <= J && ~isempty(obj.nodeCells{i,j})

					c = obj.nodeCells{i,j};
					ids = c(1:end-1);
					colour = c(end);
					nodeCoords = squeeze(obj.nodes(ids,i,:));

					x = nodeCoords(1) + r*usX;
					y = nodeCoords(2) + r*usY;
					z = nodeCoords(3) + r*usZ;

					if j > length(surfObjects)
						surfObjects(j) = surf(ax, x, y, z, 'LineStyle', 'none', 'FaceColor', obj.cs.GetRGB(colour));
					else
						surfObjects(j).XData = x;
						surfObjects(j).YData = y;
						surfObjects(j).ZData = z;
						surfObjects(j).FaceColor = obj.cs.GetRGB(colour);
						surfObjects(j).LineStyle = 'none';
					end

					j = j + 1;

				end
				% j will always end up being 1 more than the total number of non empty cells

				for k = length(surfObjects):-1:j
					surfObjects(k).delete;
					surfObjects(k) = [];
				end



				% This is a cheat way to get animations working with Tumour3D
				% all faces are displayed the same, and the number of faces never changes
				% Loop through the face cells

				[~,J,~] = size(obj.faces);
				for j = 1:J

					nIDs = obj.faces(i,j,:);
					colour = 11;
					nodeCoords = squeeze(obj.nodes(nIDs,i,:));

					x = nodeCoords(:,1);
					y = nodeCoords(:,2);
					z = nodeCoords(:,3);

					if j > length(patchObjects)
						patchObjects(j) = patch(ax, x, y, z, obj.cs.GetRGB(colour), 'FaceAlpha', 0.5, 'EdgeColor', [.5,.5,.5]);
					else
						patchObjects(j).XData = x;
						patchObjects(j).YData = y;
						patchObjects(j).ZData = z;
						patchObjects(j).FaceColor = obj.cs.GetRGB(colour);
					end

					j = j + 1;

				end
				% j will always end up being 1 more than the total number of non empty cells

				for k = length(patchObjects):-1:j
					patchObjects(k).delete;
					patchObjects(k) = [];
				end

				title(sprintf('t = %g',obj.timeSteps(i)),'Interpreter', 'latex');
				F(end+1) = getframe(h);

			end

			fileName = [obj.pathToOutput,'animation'];

			if ~isempty(indices)
				ts = obj.timeSteps(tIdxStart);
				if tIdxStart == 1
					ts = 0; % A little hack to make the numbers look nice, technically its lying
				end
				te = obj.timeSteps(tIdxEnd);
				fileName = sprintf('%s_%gto%g',fileName, ts, te );
			else
				fileName = sprintf('%s_Full',fileName);
			end

			writerObj = VideoWriter(fileName,videoFormat);
			writerObj.FrameRate = 10;

			% open the video writer
			open(writerObj);
			% write the frames to the video
			for i=2:length(F)
				% convert the image to a frame
				frame = F(i) ;    
				writeVideo(writerObj, frame);
			end
			% close the writer object
			close(writerObj);

		end

		function PlotTimeStep(obj, timeStep, varargin)

			% Plots a single given timestep
			% the number timeStep must be an integer matching the row number of the
			% saved data. Usually this will be 10xt but not always

			% varargin has inputs
			% 1: plot axis range in the form [xmin,xmax,ymin,ymax]
			% 2: plot title. can include latex

			% if you want to ignore a particular input, use []

			xyrange = [];
			plotTitle = '';

			if ~isempty(varargin)
				xyrange = varargin{1};
				if length(varargin)>1
					plotTitle = varargin{2};
				end
			end

			h = figure();
			axis equal
			hold on

			i = timeStep;


			% Initialise the array with anything
			fillObjects(1) = fill([1,1],[2,2],'r');


			[~,J] = size(obj.cells);
			j = 1;
			while j <= J && ~isempty(obj.cells{i,j})

				c = obj.cells{i,j};
				ids = c(1:end-1);
				colour = c(end);
				nodeCoords = squeeze(obj.nodes(ids,i,:));

				x = nodeCoords(:,1);
				y = nodeCoords(:,2);
				
				fillObjects(j) = fill(x,y,obj.cs.GetRGB(colour));


				j = j + 1;

			end

			if ~isempty(xyrange)
				xlim(xyrange(1:2));
				ylim(xyrange(3:4));
			end
			
			% j will always end up being 1 more than the total number of non empty cells
			axis off
			drawnow
			if ~isempty(plotTitle)
				title(plotTitle,'Interpreter', 'latex', 'FontSize', 34);
			else
				title(sprintf('t = %g',obj.timeSteps(i)),'Interpreter', 'latex', 'FontSize', 34);
			end

			set(h,'Units','Inches');
			pos = get(h,'Position');
			set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
			
			fileName = sprintf('ImageAtTime_%g',obj.timeSteps(timeStep));
			fileName = strrep(fileName,'.','_'); % If any time has decimals, change the point to underscore
			fileName = sprintf('%s%s', obj.pathToOutput, fileName);
			print(fileName,'-dpdf')

		end

		function VisualiseRods(obj, r, varargin)

			% varargin 
			% Arg 1: [indexStart, indexEnd] - a vector of the start and ending indices. Leave empty to run the whole simulation
			% Arg 2: plot axis range in the form [xmin,xmax,ymin,ymax]

			xyrange = [];
			indices = [];
			if ~isempty(varargin)
				indices = varargin{1};
				if length(varargin) > 1
					xyrange = varargin{2};
				end
			end

			h = figure();
			set(gca,'Color','k');
			axis equal
			hold on

			if ~isempty(xyrange)
				xlim(xyrange(1:2));
				ylim(xyrange(3:4));
			end

			[I,~] = size(obj.cells);


			% Initialise the array with anything
			patchObjects(1) = patch([1,1],[2,2],obj.cs.GetRGB(6), 'LineWidth', 0.5);

			startI =  1;
			endI = I;
			if ~isempty(indices)
				startI = indices(1);
				endI = indices(2);
			end

			for i = startI:I
				% i is the time steps
				[~,J] = size(obj.cells);
				j = 1;
				while j <= J && ~isempty(obj.cells{i,j})

					c = obj.cells{i,j};
					ids = c(1:end-1);
					colour = c(end);
					nodeCoords = squeeze(obj.nodes(ids,i,:));

					a = nodeCoords(1,:);
					b = nodeCoords(2,:);

					if j > length(patchObjects)
						[pillX,pillY] = obj.DrawPill(a,b,r);
						patchObjects(j) = patch(pillX,pillY,obj.cs.GetRGB(colour), 'LineWidth', .5);
					else
						[pillX,pillY] = obj.DrawPill(a,b,r);
						patchObjects(j).XData = pillX;
						patchObjects(j).YData = pillY;
						patchObjects(j).FaceColor = obj.cs.GetRGB(colour);
					end

					j = j + 1;

				end
				% j will always end up being 1 more than the total number of non empty cells

				for k = length(patchObjects):-1:j
					patchObjects(k).delete;
					patchObjects(k) = [];
				end

				drawnow
				title(sprintf('t = %g',obj.timeSteps(i)),'Interpreter', 'latex');
				pause(0.1);

			end

		end

		function ProduceRodMovie(obj, r, varargin)


			% varargin 
			% Arg 1: [indexStart, indexEnd] - a vector of the start and ending indices
			% Arg 2: plot axis range in the form [xmin,xmax,ymin,ymax]

			xyrange = [];
			indices = [];
			videoFormat = 'MPEG-4';
			if ~isempty(varargin)
				indices = varargin{1};
				if length(varargin) > 1
					xyrange = varargin{2};
					if length(varargin) > 2
						videoFormat = varargin{3};
					end
				end
			end

			% Currently same as run visualiser, but saves the movie

			h = figure();
			axis equal
			% axis off
			hold on
			set(h, 'InvertHardcopy', 'off')
			set(h,'color','w');
			set(gca,'Color','k');

			if ~isempty(xyrange)
				xlim(xyrange(1:2));
				ylim(xyrange(3:4));
			end

			if ~isempty(indices)
				tIdxStart = indices(1);
				tIdxEnd = indices(2);
			else
				tIdxStart = 1;
				tIdxEnd = length(obj.timeSteps);
			end

			F = getframe(gca); % Initialise the array

			% Initialise the array with anything
			patchObjects(1) = patch([1,1],[2,2],obj.cs.GetRGB(6), 'LineWidth', 0.5);

			for i = tIdxStart:tIdxEnd
				% i is the time steps
				[~,J] = size(obj.cells);
				j = 1;
				while j <= J && ~isempty(obj.cells{i,j})

					c = obj.cells{i,j};
					ids = c(1:end-1);
					colour = c(end);
					nodeCoords = squeeze(obj.nodes(ids,i,:));

					a = nodeCoords(1,:);
					b = nodeCoords(2,:);

					if j > length(patchObjects)
						[pillX,pillY] = obj.DrawPill(a,b,r);
						patchObjects(j) = patch(pillX,pillY,obj.cs.GetRGB(colour), 'LineWidth', 0.5);
					else
						[pillX,pillY] = obj.DrawPill(a,b,r);
						patchObjects(j).XData = pillX;
						patchObjects(j).YData = pillY;
						patchObjects(j).FaceColor = obj.cs.GetRGB(colour);
					end

					j = j + 1;

				end
				% j will always end up being 1 more than the total number of non empty cells

				for k = length(patchObjects):-1:j
					patchObjects(k).delete;
					patchObjects(k) = [];
				end

				drawnow

				title(sprintf('t = %g',obj.timeSteps(i)),'Interpreter', 'latex');
				F(end+1) = getframe(gca);

			end

			fileName = [obj.pathToOutput,'animation'];

			if ~isempty(indices)
				ts = obj.timeSteps(tIdxStart);
				if tIdxStart == 1
					ts = 0; % A little hack to make the numbers look nice, technically its lying
				end
				te = obj.timeSteps(tIdxEnd);
				fileName = sprintf('%s_%gto%g',fileName, ts, te );
			else
				fileName = sprintf('%s_Full',fileName);
			end

			writerObj = VideoWriter(fileName,videoFormat);
			writerObj.FrameRate = 10;

			% open the video writer
			open(writerObj);
			% write the frames to the video
			for i=2:length(F)
				% convert the image to a frame
				frame = F(i) ;    
				writeVideo(writerObj, frame);
			end
			% close the writer object
			close(writerObj);

		end

		function ProduceRodAngleMovie(obj, r, varargin)


			% varargin 
			% Arg 1: [indexStart, indexEnd] - a vector of the start and ending indices
			% Arg 2: plot axis range in the form [xmin,xmax,ymin,ymax]

			xyrange = [];
			indices = [];
			videoFormat = 'MPEG-4';
			if ~isempty(varargin)
				indices = varargin{1};
				if length(varargin) > 1
					xyrange = varargin{2};
					if length(varargin) > 2
						videoFormat = varargin{3};
					end
				end
			end

			lineWidth = 0.5;

			% Currently same as run visualiser, but saves the movie

			h = figure();
			axis equal
			% axis off
			hold on
			set(h, 'InvertHardcopy', 'off')
			set(h,'color','w');
			set(gca,'Color','k');

			if ~isempty(xyrange)
				xlim(xyrange(1:2));
				ylim(xyrange(3:4));
			end

			if ~isempty(indices)
				tIdxStart = indices(1);
				tIdxEnd = indices(2);
			else
				tIdxStart = 1;
				tIdxEnd = length(obj.timeSteps);
			end

			F = getframe(gca); % Initialise the array

			% Initialise the array with anything
			patchObjects(1) = patch([1,1],[2,2],obj.cs.GetRGB(6), 'LineWidth', lineWidth);

			for i = tIdxStart:tIdxEnd
				% i is the time steps
				[~,J] = size(obj.cells);
				j = 1;
				while j <= J && ~isempty(obj.cells{i,j})

					c = obj.cells{i,j};
					ids = c(1:end-1);
					colour = c(end);
					nodeCoords = squeeze(obj.nodes(ids,i,:));

					a = nodeCoords(1,:);
					b = nodeCoords(2,:);

					x = nodeCoords(:,1);
					y = nodeCoords(:,2);

					angColour = 2 * abs( atan( (x(1)-x(2)) / (y(1)-y(2))) ) / pi;

					[pillX,pillY] = obj.DrawPill(a,b,r);
					

					if j > length(patchObjects)
						patchObjects(j) = patch(pillX, pillY, angColour, 'LineWidth', lineWidth);
					else
						patchObjects(j).XData = pillX;
						patchObjects(j).YData = pillY;
						patchObjects(j).FaceVertexCData = angColour;
					end

					j = j + 1;

				end
				% j will always end up being 1 more than the total number of non empty cells

				for k = length(patchObjects):-1:j
					patchObjects(k).delete;
					patchObjects(k) = [];
				end

				drawnow

				title(sprintf('t = %g',obj.timeSteps(i)),'Interpreter', 'latex');
				F(end+1) = getframe(gca);

			end

			fileName = [obj.pathToOutput,'animation_angle'];

			if ~isempty(indices)
				ts = obj.timeSteps(tIdxStart);
				if tIdxStart == 1
					ts = 0; % A little hack to make the numbers look nice, technically its lying
				end
				te = obj.timeSteps(tIdxEnd);
				fileName = sprintf('%s_%gto%g',fileName, ts, te );
			else
				fileName = sprintf('%s_Full',fileName);
			end

			writerObj = VideoWriter(fileName,videoFormat);
			writerObj.FrameRate = 10;

			% open the video writer
			open(writerObj);
			% write the frames to the video
			for i=2:length(F)
				% convert the image to a frame
				frame = F(i) ;    
				writeVideo(writerObj, frame);
			end
			% close the writer object
			close(writerObj);

		end

		function PlotRodTimeStep(obj, r, timeStep)

			% Plots a single given timestep

			h = figure();
			axis equal
			hold on
			set(h, 'InvertHardcopy', 'off')
			set(h,'color','w');
			set(gca,'Color','k');

			i = timeStep;

			lineWidth = 0.5;
			% Initialise the array with anything
			patchObjects(1) = patch([1,1],[2,2],obj.cs.GetRGB(6), 'LineWidth', lineWidth);


			[~,J] = size(obj.cells);
			j = 1;
			while j <= J && ~isempty(obj.cells{i,j})

				c = obj.cells{i,j};
				ids = c(1:end-1);
				colour = c(end);
				nodeCoords = squeeze(obj.nodes(ids,i,:));

				a = nodeCoords(1,:);
				b = nodeCoords(2,:);

				[pillX,pillY] = obj.DrawPill(a,b,r);
				patchObjects(j) = patch(pillX,pillY,obj.cs.GetRGB(colour), 'LineWidth', lineWidth);

				j = j + 1;

			end
				% j will always end up being 1 more than the total number of non empty cells
			% axis off
			drawnow
			title(sprintf('t = %g',obj.timeSteps(i)),'Interpreter', 'latex');


			set(h,'Units','Inches');
			pos = get(h,'Position');
			set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
			

			fileName = sprintf('ImageAtTime_%g',obj.timeSteps(timeStep));
			fileName = strrep(fileName,'.','_'); % If any time has decimals, change the point to underscore
			fileName = sprintf('%s%s', obj.pathToOutput, fileName);
			print(fileName,'-dpdf')

		end

		function PlotRodAngles(obj, r, timeStep)

			% Plots a single given timestep

			h = figure();
			axis equal
			hold on
			set(gca,'Color','k');
			set(h, 'InvertHardcopy', 'off')
			set(h,'color','w');

			i = timeStep;

			lineWidth = 0.5;
			% Initialise the array with anything
			patchObjects(1) = patch([1,1],[2,2],obj.cs.GetRGB(6), 'LineWidth', lineWidth);


			[~,J] = size(obj.cells);
			j = 1;
			while j <= J && ~isempty(obj.cells{i,j})

				c = obj.cells{i,j};
				ids = c(1:end-1);
				colour = c(end);
				nodeCoords = squeeze(obj.nodes(ids,i,:));

				a = nodeCoords(1,:);
				b = nodeCoords(2,:);

				x = nodeCoords(:,1);
				y = nodeCoords(:,2);

				angColour = 2 * abs( atan( (x(1)-x(2)) / (y(1)-y(2))) ) / pi;

				[pillX,pillY] = obj.DrawPill(a,b,r);
				patchObjects(j) = patch(pillX,pillY,angColour, 'LineWidth', lineWidth);

				j = j + 1;

			end
				% j will always end up being 1 more than the total number of non empty cells
			% axis off
			drawnow
			title(sprintf('t = %g',obj.timeSteps(i)),'Interpreter', 'latex');


			set(h,'Units','Inches');
			pos = get(h,'Position');
			set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
			

			fileName = sprintf('AnglesAtTime_%g',obj.timeSteps(timeStep));
			fileName = strrep(fileName,'.','_'); % If any time has decimals, change the point to underscore
			fileName = sprintf('%s%s', obj.pathToOutput, fileName);
			print(fileName,'-dpdf')

		end

		function VisualiseNodesAndEdges(obj, r, varargin)

			% r is the radius of the node
			% varargin 
			% Arg 1: [indexStart, indexEnd] - a vector of the start and ending indices. Leave empty to run the whole simulation
			% Arg 2: plot axis range in the form [xmin,xmax,ymin,ymax]

			xyrange = [];
			indices = [];
			if ~isempty(varargin)
				indices = varargin{1};
				if length(varargin) > 1
					xyrange = varargin{2};
				end
			end

			h = figure();
			axis equal
			hold on

			if ~isempty(xyrange)
				xlim(xyrange(1:2));
				ylim(xyrange(3:4));
			end

			[I,~] = size(obj.cells);


			% Initialise the array with anything
			patchObjects(1) = patch([1,1],[2,2],obj.cs.GetRGB(6), 'LineWidth', 0.5);
			lineObjects(1)  = line([1,1],[2,2],'Color', 'k', 'LineWidth', 4);

			startI =  1;
			endI = I;
			if ~isempty(indices)
				startI = indices(1);
				endI = indices(2);
			end

			for i = startI:endI
				
				% First draw node cells

				[~,J] = size(obj.cells);
				j = 1; % loops node cells
				jN = 0; % tracks the number of node cells
				jM = 0; % tracks the number of membrane 'cells'
				while j <= J && ~isempty(obj.cells{i,j})

					c = obj.cells{i,j};
					ids = c(1:end-1);
					colour = c(end);

					% This should only for the node cells
					if length(ids) < 2
						jN = jN + 1;
						a = squeeze(obj.nodes(ids,i,:))';

						if jN > length(patchObjects)
							[pillX,pillY] = obj.DrawPill(a,a,r);
							patchObjects(jN) = patch(pillX,pillY,obj.cs.GetRGB(colour), 'LineWidth', .5);
						else
							[pillX,pillY] = obj.DrawPill(a,a,r);
							patchObjects(jN).XData = pillX;
							patchObjects(jN).YData = pillY;
							patchObjects(jN).FaceColor = obj.cs.GetRGB(colour);
						end

					else

						% Doesn't quite work for closed loops
						% so need a little hack for now
						jM = jM + 1;

						nodeCoords = squeeze(obj.nodes(ids,i,:));
						x = nodeCoords(:,1);
						y = nodeCoords(:,2);

						x(end+1) = x(1);
						y(end+1) = y(1);

						if jM > length(lineObjects)
							lineObjects(jM) = line(x,y, 'Color', 'k', 'LineWidth', 4);
						else
							lineObjects(jM).XData = x;
							lineObjects(jM).YData = y;
						end

					end


					j = j + 1;

				end
				% j will always end up being 1 more than the total number of non empty cells

				for k = length(patchObjects):-1:jN+1
					patchObjects(k).delete;
					patchObjects(k) = [];
				end

				for k = length(lineObjects):-1:jM+1
					lineObjects(k).delete;
					lineObjects(k) = [];
				end

				drawnow
				title(sprintf('t = %g',obj.timeSteps(i)),'Interpreter', 'latex');
				pause(0.1);

			end

		end

		function ProduceNodesAndEdgesMovie(obj, r, varargin)


			% varargin 
			% Arg 1: [indexStart, indexEnd] - a vector of the start and ending indices
			% Arg 2: plot axis range in the form [xmin,xmax,ymin,ymax]

			xyrange = [];
			indices = [];
			videoFormat = 'MPEG-4';
			if ~isempty(varargin)
				indices = varargin{1};
				if length(varargin) > 1
					xyrange = varargin{2};
					if length(varargin) > 2
						videoFormat = varargin{3};
					end
				end
			end

			% Currently same as run visualiser, but saves the movie

			h = figure();
			axis equal
			axis off
			hold on

			if ~isempty(xyrange)
				xlim(xyrange(1:2));
				ylim(xyrange(3:4));
			end

			if ~isempty(indices)
				startI = indices(1);
				endI = indices(2);
			else
				startI = 1;
				endI = length(obj.timeSteps);
			end

			F = getframe(gca); % Initialise the array

			% Initialise the array with anything
			patchObjects(1) = patch([1,1],[2,2],obj.cs.GetRGB(6), 'LineWidth', 0.5);
			lineObjects(1)  = line([1,1],[2,2],'Color', 'k', 'LineWidth', 4);


			for i = startI:endI
				
				% First draw node cells

				[~,J] = size(obj.cells);
				j = 1; % loops node cells
				jN = 0; % tracks the number of node cells
				jM = 0; % tracks the number of membrane 'cells'
				while j <= J && ~isempty(obj.cells{i,j})

					c = obj.cells{i,j};
					ids = c(1:end-1);
					colour = c(end);

					% This should only for the node cells
					if length(ids) < 2
						jN = jN + 1;
						a = squeeze(obj.nodes(ids,i,:))';

						if jN > length(patchObjects)
							[pillX,pillY] = obj.DrawPill(a,a,r);
							patchObjects(jN) = patch(pillX,pillY,obj.cs.GetRGB(colour), 'LineWidth', .5);
						else
							[pillX,pillY] = obj.DrawPill(a,a,r);
							patchObjects(jN).XData = pillX;
							patchObjects(jN).YData = pillY;
							patchObjects(jN).FaceColor = obj.cs.GetRGB(colour);
						end

					else

						% Doesn't quite work for closed loops
						% so need a little hack for now
						jM = jM + 1;

						nodeCoords = squeeze(obj.nodes(ids,i,:));
						x = nodeCoords(:,1);
						y = nodeCoords(:,2);

						x(end+1) = x(1);
						y(end+1) = y(1);

						if jM > length(lineObjects)
							lineObjects(jM) = line(x,y, 'Color', 'k', 'LineWidth', 4);
						else
							lineObjects(jM).XData = x;
							lineObjects(jM).YData = y;
						end

					end


					j = j + 1;

				end
				% j will always end up being 1 more than the total number of non empty cells

				for k = length(patchObjects):-1:jN+1
					patchObjects(k).delete;
					patchObjects(k) = [];
				end

				for k = length(lineObjects):-1:jM+1
					lineObjects(k).delete;
					lineObjects(k) = [];
				end

				drawnow
				title(sprintf('t = %g',obj.timeSteps(i)),'Interpreter', 'latex');
				F(end+1) = getframe(gca);

			end


			fileName = [obj.pathToOutput,'animation'];

			if ~isempty(indices)
				ts = obj.timeSteps(tIdxStart);
				if tIdxStart == 1
					ts = 0; % A little hack to make the numbers look nice, technically its lying
				end
				te = obj.timeSteps(tIdxEnd);
				fileName = sprintf('%s_%gto%g',fileName, ts, te );
			else
				fileName = sprintf('%s_Full',fileName);
			end

			writerObj = VideoWriter(fileName,videoFormat);
			writerObj.FrameRate = 10;

			% open the video writer
			open(writerObj);
			% write the frames to the video
			for i=2:length(F)
				% convert the image to a frame
				frame = F(i) ;    
				writeVideo(writerObj, frame);
			end
			% close the writer object
			close(writerObj);

		end


		function PlotNodesAndEdgesTimeStep(obj, r, timeStep, varargin)


			xyrange = [];
			plotTitle = '';

			if ~isempty(varargin)
				xyrange = varargin{1};
				if length(varargin)>1
					plotTitle = varargin{2};
				end
			end


			h = figure();
			axis equal
			hold on

			i = timeStep;

			lineWidth = 0.5;
			% Initialise the array with anything
			patchObjects(1) = patch([1,1],[2,2],obj.cs.GetRGB(6), 'LineWidth', lineWidth);
			lineObjects(1)  = line([1,1],[2,2],'Color', 'k', 'LineWidth', 4);



			[~,J] = size(obj.cells);
			j = 1; % loops node cells

			while j <= J && ~isempty(obj.cells{i,j})

				c = obj.cells{i,j};
				ids = c(1:end-1);
				colour = c(end);

				% This should only for the node cells
				if length(ids) < 2

					a = squeeze(obj.nodes(ids,i,:))';

					[pillX,pillY] = obj.DrawPill(a,a,r);
					patchObjects(j) = patch(pillX,pillY,obj.cs.GetRGB(colour), 'LineWidth', .5);

				else

					nodeCoords = squeeze(obj.nodes(ids,i,:));
					x = nodeCoords(:,1);
					y = nodeCoords(:,2);

					x(end+1) = x(1);
					y(end+1) = y(1);

					% This assumes only a single membrane exists
					% to handle multiple membranes, change 1 to j
					% but this will break the "bring to front command" uistack(lineObjects(1),'top')
					lineObjects(1) = line(x,y, 'Color', 'k', 'LineWidth', 4);

				end

				j = j + 1;

			end

			uistack(lineObjects(1),'top');


			if ~isempty(xyrange)
				xlim(xyrange(1:2));
				ylim(xyrange(3:4));
			end
			
			% j will always end up being 1 more than the total number of non empty cells
			axis off
			drawnow
			if ~isempty(plotTitle)
				title(plotTitle,'Interpreter', 'latex', 'FontSize', 34);
			else
				title(sprintf('t = %g',obj.timeSteps(i)),'Interpreter', 'latex', 'FontSize', 34);
			end

			set(h,'Units','Inches');
			pos = get(h,'Position');
			set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
			
			fileName = sprintf('ImageAtTime_%g',obj.timeSteps(timeStep));
			fileName = strrep(fileName,'.','_'); % If any time has decimals, change the point to underscore
			fileName = sprintf('%s%s', obj.pathToOutput, fileName);
			print(fileName,'-dpdf')


		end


		function [pillX,pillY] = DrawPill(obj,a,b,r)

			% Draws a pill shape where the centre of the circles are at
			% a and b and the radius is r

			AtoB = b - a;
			 
			normAtoB = [-AtoB(2), AtoB(1)];
			 
			normAtoB = normAtoB / norm(normAtoB);
			if isnan(normAtoB)
				normAtoB = [1,0];
			end
			 
			R = r*normAtoB;
			% Make n equally spaced points around a circle starting from R
			
			n = 10;
			apoints = [];
			bpoints = [];
			
			rot = @(theta) [cos(theta), -sin(theta); sin(theta), cos(theta)];
			
			for i=1:n-1
				
				theta = i*pi/n;
				apoints(i,:) = rot(theta)*R' + a';
				bpoints(i,:) = -rot(theta)*R' + b';
				
			end
			pill = [ a + R; apoints; a - R; b - R; bpoints;  b + R];
			
			pillX = pill(:,1);
			pillY = pill(:,2);

		end


	end


end