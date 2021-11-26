classdef Node < matlab.mixin.SetGet
	% A class specifying the details about nodes

	properties
		% Essential properties of a node
		x
		y
		z

		pos

		prevPos

		id = 0 % ids are handled by the simulation object

		force = [0, 0, 0]

		prevForce = [0, 0, 0]

		% This will be circular - each edge will have two nodes
		% each node can be part of multiple edges, similarly for cells
		edgeList = []
		faceList = []
		cellList = []

		% Each node stores it's local drag coefficient, so we can distinguish
		% between different regions in a tissue more easily
		eta = 1

		nodeAdjusted = false;
		preAdjustedPosition = [];

		nodeData

	end

	methods

		function obj = Node(x,y,z)
			% Initialise the node

			obj.x 	= x;
			obj.y 	= y;
			obj.z 	= z;

			obj.pos = [x,y,z];

			obj.prevPos = [x,y,z];

			obj.AddNodeData(NodeNeighbours(1.5, 10));
			obj.AddNodeData(FaceNeighbours(0.6, 10));

		end

		function delete(obj)

			clear obj;

		end

		function AddNodeData(obj, d)

			% Need to explicitly create a map object or matlab
			% will only point to one map object for the
			% entire list of Nodes...
			if isempty(obj.nodeData)
				nD = containers.Map;
				for i = 1:length(d)
					nD(d(i).name) = d(i);
				end
				obj.nodeData = nD;
			else
				for i = 1:length(d)
					obj.nodeData(d(i).name) = d(i);
				end
			end

		end

		function data = GetData(obj, name, t)

			data = obj.nodeData(name).GetData(obj, t);

		end

		function ApplyForce(obj, force)
			
			if sum(isnan(force)) || sum(isinf(force))
				error('N:ApplyForce:InfNaN', 'Force is inf or NaN');
			end
			obj.force = obj.force + force;

		end

		function MoveNode(obj, pos)
			% This function is used to move the position due to time stepping
			% so the force must be reset here
			% This is only to be used by the numerical integration

			obj.NewPosition(pos);
			% Reset the force for next time step
			obj.prevForce = obj.force;
			obj.force = [0,0,0];

		end

		function AdjustPosition(obj, pos)
			% Used when modifying the position manually
			% Doesn't affect previous position, or reset the force
			% But it will require fixing up the space partition

			obj.preAdjustedPosition = obj.pos;

			obj.pos = pos;

			obj.x = pos(1);
			obj.y = pos(2);
			obj.z = pos(3);

			obj.nodeAdjusted = true;

		end

		function AddEdge(obj, e)

			% e can be a vector
			if sum( ismember(e,obj.edgeList)) ~=0
				warning('N:AddEdge:EdgeAlreadyHere', 'Adding at least one edge that already appears in edgeList for Node %d. This has not been added.', obj.id);
				e(ismember(e,obj.edgeList)) = [];
			end
			obj.edgeList = [obj.edgeList , e];
			
		end

		function RemoveEdge(obj, e)
			
			% Remove the edge from the list
			if sum(obj.edgeList == e) == 0
				warning('N:RemoveEdge:EdgeNotHere', 'Edge %d does not appear in edgeList for Node %d', e.id, obj.id);
			else
				obj.edgeList(obj.edgeList == e) = [];
			end

		end

		function ReplaceEdgeList(obj, eList)

			obj.edgeList = eList;

		end

		function AddFace(obj, f)

			% c can be a vector
			if sum( ismember(f,obj.faceList)) ~=0
				warning('N:AddFace:FaceAlreadyHere', 'Adding at least one face that already appears in faceList for Node %d. This has not been added.', obj.id);
				f(ismember(f,obj.faceList)) = [];
			end
			obj.faceList = [obj.faceList , f];

		end

		function RemoveFace(obj, f)

			% Remove the face from the list
			if sum(obj.faceList == f) == 0
				warning('N:RemoveFace:FaceNotHere', 'At least one face does not appear in nodeList for Node %d', obj.id);
			else
				obj.faceList(obj.faceList == f) = [];
			end

		end

		function ReplaceFaceList(obj, fList)

			% Used for FaceFree to overwrite the existing face
			% Does not modify any links in the face, it assumes
			% they are handled in the division or creation process

			obj.faceList = fList;

		end


		function AddCell(obj, c)

			% c can be a vector
			if sum( ismember(c,obj.cellList)) ~=0
				warning('N:AddCell:CellAlreadyHere', 'Adding at least one cell that already appears in cellList for Node %d. This has not been added.', obj.id);
				c(ismember(c,obj.cellList)) = [];
			end
			obj.cellList = [obj.cellList , c];

		end

		function RemoveCell(obj, c)

			% Remove the cell from the list
			if sum(obj.cellList == c) == 0
				warning('N:RemoveCell:CellNotHere', 'At least one cell does not appear in nodeList for Node %d', obj.id);
			else
				obj.cellList(obj.cellList == c) = [];
			end

		end

		function ReplaceCellList(obj, cList)

			% Used for CellFree to overwrite the existing cell
			% Does not modify any links in the cell, it assumes
			% they are handled in the division or creation process

			obj.cellList = cList;

		end

		function SetDragCoefficient(obj, eta)

			% Use this to change the drag coefficient
			% so that the associated edges have their
			% properties updated
			obj.eta = eta;

			for i = 1:length(obj.edgeList)

				obj.edgeList(i).UpdateTotalDrag();

			end

		end

	end

	methods (Access = private)
		
		function NewPosition(obj, pos)

			% Should not be used directly, only as part of MoveNode
			obj.prevPos = obj.pos;
			obj.pos = pos;

			obj.x = pos(1);
			obj.y = pos(2);
			obj.z = pos(3);

		end

	end


end
