classdef AbstractCell < handle & matlab.mixin.Heterogeneous
	% A class specifying the details about nodes

	properties
		% Essential properties of a node
		id

		age = 0

		nodeList = Node.empty()
		edgeList = Edge.empty()
		faceList = Face.empty()

		CellCycleModel

		% Determines if we are using a free or joined cell model
		freeCell = false
		newFreeCellSeparation = 0.1

		% A cell divides in 2, this will store the sister
		% cell after division
		sisterCell = AbstractCell.empty();

		% Stores the id of the cell that was in the 
		% initial configuration. Only can store the id
		% because the cell can be deleted from the simulation
		ancestorId

		% A collection objects for calculating data about the cell
		% stored in a map container so each type of data can be given a
		% meaingful name
		cellData

		% By default, the type is 1, matching a general epithelial cell
		cellType = 1
		
	end

	methods (Abstract)

		% Must implement a divide method that returns a new cell
		% and a structure containing the new components and the
		% components that need to be removed
		[newCellList, newNodeList, newEdgeList] = Divide(obj)
		
		% inside = IsPointInsideCell(obj, point) % I don't think this needs to be abstract

	end

	methods

		function delete(obj)

			clear obj;

		end

		function set.CellCycleModel( obj, v )
			% This is to validate the object given to outputType in the constructor
			if isa(v, 'AbstractCellCycleModel')
            	validateattributes(v, {'AbstractCellCycleModel'}, {});
            	obj.CellCycleModel = v;
            	v.containingCell = obj;
            else
            	error('C:NotValidCCM','Not a valid cell cycle');
            end

        end

        function currentVolume = GetVolume(obj)
			% This and the following 3 functions could be replaced by accessing the cellData
			% but they're kept here for backwards compatibility, and because
			% these types of data are fundamental enough to designate a function

			currentVolume = obj.cellData('cellVolume').GetData(obj);

		end

		% function targetArea = GetCellTargetArea(obj)
		% 	% This is so the target area can be a function of cell age

		% 	targetArea = obj.cellData('targetArea').GetData(obj);

		% end

		% function currentPerimeter = GetCellPerimeter(obj)

		% 	currentPerimeter = obj.cellData('cellPerimeter').GetData(obj);

		% end

		% function centre = GetCellCentre(obj)

		% 	centre = obj.cellData('cellCentre').GetData(obj);

		% end

		% function targetPerimeter = GetCellTargetPerimeter(obj)
		% 	% This is so the target Perimeter can be a function of cell age
		% 	targetPerimeter = obj.cellData('targetPerimeter').GetData(obj);

		% end

		function ready = IsReadyToDivide(obj)

			ready = obj.CellCycleModel.IsReadyToDivide();

		end

		function AddCellData(obj, d)

			% Need to explicitly create a map object or matlab
			% will only point to one map object for the
			% entire list of Cells...
			if isempty(obj.cellData)
				cD = containers.Map;
				for i = 1:length(d)
					cD(d(i).name) = d(i);
				end
				obj.cellData = cD;
			else
				for i = 1:length(d)
					obj.cellData(d(i).name) = d(i);
				end
			end

		end

		function AgeCell(obj, dt)

			% This will be done at the end of the time step
			obj.age = obj.age + dt;
			obj.CellCycleModel.AgeCellCycle(dt);

		end

		function age = GetAge(obj)

			age = obj.CellCycleModel.GetAge();
			
		end

		function colour = GetColour(obj)
			% Used for animating/plotting only

			colour = obj.CellCycleModel.GetColour();
		
		end

		%**********************************************************
		% Not checked in this prototype so far
		%**********************************************************
		function inside = IsNodeInsideCell(obj, n)

			% Assemble vertices in the correct order to produce a quadrilateral

			x = [obj.nodeList.x];
			y = [obj.nodeList.y];

			[inside, on] = inpolygon(n.x, n.y, x ,y);

			if inside && on
				inside = false;
			end

		end

	end

end