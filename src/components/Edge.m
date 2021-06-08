classdef Edge < matlab.mixin.SetGet
	% A class specifying the details about nodes

	properties
		% Essential porperties of an edge
		id = 0; % Handled by simulation object

		% The 'total drag' for the edge at the centre of drag
		% This will be constant for an edge, unless the 
		% coeficient of drag is changed for a node
		etaD 

		% This will be circular - each edge will have two nodes
		% each node can be part of multiple edges
		Node1
		Node2

		% Used only when an edge is modified during cell
		% division. Keep track of the old node to help with
		% adjusting the edge boxes
		oldNode1
		oldNode2

		modifiedInDivision = false

		nodeList = []
		surfList = []
		cellList = []

		torque = [0,0,0]
		previousTorque = [0,0,0]

		internal = false

		% Detemines if an edge is part of a membrane, hence
		% should be subject to membrane related forces
		isMembrane = false
		
	end

	methods

		function obj = Edge(Node1, Node2)
			% All the initilising
			% An edge will always have a pair of nodes

			% This ordering is important, because it defines the
			% orientation of the edge. If the edge is 
			% completely free in space, then the orientation is
			% arbitrary, but as soon as it is part of a cell or membrane
			% then 1 to 2 must be chosen
			% so that if u = (Node2.position - Node1.position) / length,
			% then the unit vector v = [u(2), -u(1)] points out of the tissue.
			% Then, v will point to the right if forward is 1 to 2.
			% Due to the choice of v, travelling 1 to 2 will go
			% anticlockwise around the perimeter of the cell 

			% Note: v should be perpendicular, but I can't do that in ASCII art
			%        2 o
			%         /
			%        /
			%       /\
			%      /  \v
			%    u/    \
			%    /      
			% 1 o

			obj.Node1 = Node1;
			obj.Node2 = Node2;

			obj.nodeList = [Node1, Node2];

			obj.Node1.AddEdge(obj);
			obj.Node2.AddEdge(obj);

			obj.UpdateTotalDrag();
			
		end

		function delete(obj)

			clear obj;

		end

		function UpdateTotalDrag(obj)
			% Of the three important physical quantities
			% total drag will only change if the drag coefficients
			% are explicitly changed

			obj.etaD = obj.Node1.eta + obj.Node2.eta;

		end

		function ID = GetMomentOfDrag(obj)
			% The length of the edge will change at every time
			% step, so ID needs to be calculated every time
			r1 = obj.Node1.position;
			r2 = obj.Node2.position;

			% Centre of drag
			rD = (obj.Node1.eta * r1 + obj.Node2.eta * r2) / obj.etaD;

			rDto1 = r1 - rD;
			rDto2 = r2 - rD;

			ID = obj.Node1.eta * norm(rDto1)^2 + obj.Node2.eta * norm(rDto2)^2;

		end

		function len = GetLength(obj)
			
			len = norm(obj.Node1.position - obj.Node2.position);

		end

		function direction1to2 = GetUnitVector1to2(obj)
			
			direction1to2 = obj.Node2.position - obj.Node1.position;
			direction1to2 = direction1to2 / norm(direction1to2);

		end

		function outward = GetOutwardUnitNormal(obj)
			% See constructor for discussion about why this is the outward
			% normal
			u = obj.GetVector1to2();
			outward = [u(2), -u(1)];

		end

		function midPoint = GetMidPoint(obj)

			direction1to2 = obj.Node2.position - obj.Node1.position;
			midPoint = obj.Node1.position + 0.5 * direction1to2;

		end

		function ApplyTorque(obj, torque)
			
			% Gather the torques applied to the edge.
			% This will be converted to a force on the nodes
			if sum(isnan(torque)) || sum(isinf(torque))
				error('N:ApplyTorque:InfNaN', 'Torque is inf or NaN');
			end
			obj.torque = obj.torque + torque;

		end

		function otherNode = GetOtherNode(obj, node)

			if node == obj.Node1

				otherNode = obj.Node2;

			else

				if node == obj.Node2
					otherNode = obj.Node1;
				else
					error('E:GetOtherNode:NodeNotHere', 'Node %d is not in Edge %d', node.id, obj.id);
				end

			end

		end

		function SwapNodes(obj)

			% Used when the nodes are not anticlockwise 1 -> 2
			% An edge has no concept of orientation by itself
			% so this must be controlled from outside the edge only

			a = obj.Node1;

			obj.Node1 = obj.Node2;
			obj.Node2 = a;

			obj.nodeList = [obj.Node1, obj.Node2];

		end

		function ReplaceNode(obj, oldNode, newNode)
			% Removes the old node from the edge, and replaces
			% it with a new node. This is used in cell division primarily

			% To do this properly, we need fix all the links to nodes and cells
			if oldNode == newNode
				warning('e:sameNode', 'The old node is the same as the new node. This is probably not what you wanted to do')
			else
			
				switch oldNode
					case obj.Node1
						% Remove link back to this edge
						obj.Node1.RemoveEdge(obj);

						obj.Node1 = newNode;
						newNode.AddEdge(obj);

						obj.modifiedInDivision = true;
						obj.oldNode1 = oldNode;

						obj.nodeList = [obj.Node1, obj.Node2];

					case obj.Node2
						% Remove link back to this edge
						obj.Node2.RemoveEdge(obj);

						obj.Node2 = newNode;
						newNode.AddEdge(obj);

						obj.modifiedInDivision = true;
						obj.oldNode2 = oldNode;

						obj.nodeList = [obj.Node1, obj.Node2];
						
					otherwise
						error('e:nodeNotFound','Node not in this edge')
				end
			end

		end

		function AddCell(obj, c)

			% c can be a vector
			if sum( ismember(c,obj.cellList)) ~=0
				warning('E:AddCell:CellAlreadyHere', 'Adding at least one cell that already appears in cellList for Edge %d. This has not been added.', obj.id);
				c(ismember(c,obj.cellList)) = [];
			end
			obj.cellList = [obj.cellList , c];

		end

		function ReplaceCell(obj, oldC, newC)

			% Currently the cell list has at most two entries
			if oldC == newC
				warning('E:ReplaceCell:SameCell','Both cells are the same, nothing will happen')
			else

				obj.AddCell(newC);
				obj.RemoveCell(oldC);

			end

		end

		function ReplaceCellList(obj, cellList)

			obj.cellList = cellList;
			
		end

		function RemoveCell(obj, c)

			% Remove the cell from the list
			if sum(obj.cellList == c) == 0
				warning('E:RemoveCell:CellNotHere', 'At least one cell does not appear in edgeList for Edge %d', obj.id);
			else
				obj.cellList(obj.cellList == c) = [];
			end

		end


		function AddSurface(obj, s)

			% s can be a vector
			if sum( ismember(s,obj.surfList)) ~=0
				warning('E:AddSurf:SurfAlreadyHere', 'Adding at least one surf that already appears in surfList for Edge %d. This has not been added.', obj.id);
				s(ismember(s,obj.surfList)) = [];
			end
			obj.surfList = [obj.surfList , s];

		end

		function ReplaceSurf(obj, oldS, newS)

			% Currently the surf list has at most two entries
			if oldS == newS
				warning('E:ReplaceSurf:SameSurf','Both surfaces are the same, nothing will happen')
			else

				obj.AddSurf(newS);
				obj.RemoveSurf(oldS);

			end

		end

		function ReplaceSurfList(obj, surfList)

			obj.surfList = surfList;
			
		end

		function RemoveSurf(obj, s)

			% Remove the surf from the list
			if sum(obj.surfList == s) == 0
				warning('E:RemoveSurf:SurfNotHere', 'At least one surface does not appear in edgeList for Edge %d', obj.id);
			else
				obj.surfList(obj.surfList == s) = [];
			end

		end

		function internal = IsEdgeInternal(obj)

			internal = obj.internal;

		end

	end

end