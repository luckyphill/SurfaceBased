classdef SpatialState < AbstractSimulationData

	properties 

		name = 'spatialState'
		data = {};

	end

	methods

		function obj = SpatialState
			% No special initialisation
			
		end

		function CalculateData(obj, t)

			% In this case, the data is a structure containing all the node
			% positions, and a list of cells containing the nodes that make it up

			% At some point I want to add in edges as well, primarily for
			% when membrane is modelled, but I'll have to have a separate way
			% of handling a membrane object

			% This will only work when every cell has the same number of nodes
			% so it won't work with the stromal situation.
			% I will implement a way to introduce the stromal layer so it is separate
			% from the cell List

			nodeData = [[t.nodeList.id]',[t.nodeList.x]', [t.nodeList.y]', [t.nodeList.z]'];

			edgeData = [];
			for i = 1:length(t.edgeList())
				nL = t.edgeList(i).nodeList;
				edgeData(i,:) = [nL.id];
			end

			faceData = [];
			for i = 1:length(t.faceList())
				nL = t.faceList(i).nodeList;
				faceData(i,:) = [nL.id];
			end

			cellData = {};

			for i = 1:length(t.cellList())
				c = t.cellList(i);
				nL = c.nodeList;

				% A cell can have any number of nodes, but it's usually 4
				l = length(nL);
				cellData{i} = [l, nL.id, c.CellCycleModel.colour];
			end

			obj.data = {nodeData, edgeData, faceData, cellData};

		end
		
	end

end