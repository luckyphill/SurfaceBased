classdef Surface < matlab.mixin.SetGet
	% A class specifying the details about nodes

	properties
		% Essential porperties of an surface
		id = 0; % Handled by simulation object

		% The 'total drag' for the surface at the centre of drag
		% This will be constant for an surface, unless the 
		% coeficient of drag is changed for a node
		etaD 

		% This will be circular - each surface will have three nodes
		% each node can be part of multiple surfaces
		Node1
		Node2
		Node3

		% Used only when an surface is modified during cell
		% division. Keep track of the old node to help with
		% adjusting the surface boxes
		oldNode1
		oldNode2
		oldNode3

		% As with nodes
		Edge1
		Edge2
		Edge3

		oldEdge1
		oldEdge2
		oldEdge3

		modifiedInDivision = false

		nodeList = []
		edgeList = []
		cellList = []

		torque = [0,0,0]
		prevTorque = [0,0,0]

		internal = false

		% Detemines if a surface is part of a membrane, hence
		% should be subject to membrane related forces
		isMembrane = false
		
	end

	methods

		function obj = Surface(Edge1,Edge2,Edge3)

			% A surface is defined by three edges
			% These edges must be connected in an anticlockwise loop
			% according to the order they are supplied

			% Comprehensive error checking must be done here to ensure
			% surfaces are constructed properly.

			% The nodes and edges must be given their local indices

			obj.AddEdgesToSurface(Edge1,Edge2,Edge3);
			
		end

		function delete(obj)

			clear obj;

		end

		function etaD = GetTotalDrag(obj)
			
			% Total drag of the surface assuming pooint drag at the nodes

			etaD = obj.Node1.eta + obj.Node2.eta + obj.Node3.eta;

		end

		function rD = GetCentreOfDrag(obj)

			% Centre of drag assumin point drag at the nodes

			rD = (obj.Node1.eta * obj.Node1.pos + obj.Node2.eta * obj.Node2.pos + obj.Node3.eta * obj.Node3.pos)/obj.GetTotalDrag();

		end

		function [rD1, rD2, rD3] = GetNodeCoordinatesRelative(obj)

			% The position of the nodes relative to the centre
			% of drag in the global axis directions

			rD 		= GetCentreOfDrag(obj);

			% Global centre to global node
			rD1 	= obj.Node1.pos - rD;
			rD2 	= obj.Node2.pos - rD;
			rD3 	= obj.Node3.pos - rD;

		end

		function ID = GetDragMatrix(obj)
			
			% Construct the drag matrix based on the current
			% configuration of the surface, with the values
			% taken in an axis centred at rD, but with unit vectors
			% identical to the fixed, global system

			% In this instance, there is no special reduction of the
			% moments due to choice of axis, so every component needs
			% to be calculated

			% The values of the products of drag are calculated without
			% the negative sign, but the sign is added to the matrix
			% This makes the calculation of angular velocities a simple
			% backslash operation

			% ID = [IDX, -IDXY, -IDXZ;
			%		-IDXY, IDY, -IDYZ;
			%		- IDXZ, IDYZ,IDZ];

			[rD1, rD2, rD3] = GetNodeCoordinatesRelative(obj);

			IDX 	= obj.Node1.eta * (rD1(2)^2 + rD1(3)^2) + obj.Node2.eta * (rD2(2)^2 + rD2(3)^2) + obj.Node3.eta * (rD3(2)^2 + rD3(3)^2);
			IDY 	= obj.Node1.eta * (rD1(1)^2 + rD1(3)^2) + obj.Node2.eta * (rD2(1)^2 + rD2(3)^2) + obj.Node3.eta * (rD3(1)^2 + rD3(3)^2);
			IDZ 	= obj.Node1.eta * (rD1(1)^2 + rD1(2)^2) + obj.Node2.eta * (rD2(1)^2 + rD2(2)^2) + obj.Node3.eta * (rD3(1)^2 + rD3(2)^2);
			IDXY 	= obj.Node1.eta * rD1(1) * rD1(2) + obj.Node2.eta * rD2(1) * rD2(2) + obj.Node3.eta * rD3(1) * rD3(2);
			IDXZ	= obj.Node1.eta * rD1(1) * rD1(3) + obj.Node2.eta * rD2(1) * rD2(3) + obj.Node3.eta * rD3(1) * rD3(3); 
			IDYZ	= obj.Node1.eta * rD1(2) * rD1(3) + obj.Node2.eta * rD2(2) * rD2(3) + obj.Node3.eta * rD3(2) * rD3(3);

			ID = [IDX, -IDXY, -IDXZ;
					-IDXY, IDY, -IDYZ;
					- IDXZ, IDYZ, IDZ];

		end

		function outward = GetUnitNormal(obj)
			
			% The unit normal to the surface pointing in the upward direction
			% as defined using the right-hand rule traversing the nodes in order 1 2 3 
			ve1 = obj.Node2.pos - obj.Node1.pos;
			ve2 = obj.Node3.pos - obj.Node2.pos; 

			% Take the cross product to get a vector normal to the surface
			u = cross(ve1,ve2);

			outward = u / norm(u);

		end

		function [u,v,w] = GetCoordinateSystem(obj)

			% Returns the coordinate system attached to the
			% surface with origin at rD.
			% The z axis is the unit normal
			% The y axis is from rD to Node2
			% The x axis is determined from these by the 
			% cross product and right hand rule

			% The vectors define the coordinate system axes
			% but are expressed in the global cordinate system

			w = obj.GetUnitNormal();
			v = obj.Node2.pos - obj.GetCentreOfDrag();
			v = v/norm(v);
			u = cross(v,w);
			u = u/norm(u);

		end

		function [r1, r2, r3] = GetNodeCoordinatesLocal(obj)

			% The position of the nodes in the local system of coordinates

			[rD1, rD2, rD3] = GetNodeCoordinatesRelative(obj);

			% Convert these to local coordinates
			[u,v,w] = GetCoordinateSystem(obj);

			r1 		= [dot(u,rD1), dot(v,rD1), 0];
			r2 		= [0, norm(rD2), 0];
			r3 		= [dot(u,rD3), dot(v,rD3), 0];

		end

		function [IDx,IDy,IDxy,IDz] = GetDragMatrixLocal(obj)
			
			% Construct the drag matrix based on the current
			% configuration of the surface, with the values
			% taken in the body (local) system of coordinates

			% In the local system of coordinates, rD is the origin
			% the z axis is the surface normal, and the y axis
			% goes from rD to Node2. The x axis follows from
			% these via a cross product and the right hand rule

			% In this situation, the drag matrix is given by 4
			% quantities: IDx, IDy, IDz, IDxy

			% If all goes to plan, IDz is not needed because
			% there should be no component of rotation around
			% the z axis, but this might be added later.
			% Either way, in this configuration, IDz = IDx + IDy
			% so we virtually get it for free

			% At this point, we are assuming the drag is concentrated
			% at the nodes only

			[r1, r2, r3] = GetNodeCoordinatesLocal(obj);
			

			IDx 	= obj.Node1.eta * r1(2)^2 + obj.Node2.eta * r2(2)^2 + obj.Node3.eta * r3(2)^2;
			IDy 	= obj.Node1.eta * r1(1)^2 + obj.Node3.eta * r3(1)^2;
			IDxy 	= obj.Node1.eta * r1(1) * r1(2) + obj.Node3.eta * r3(1) * r3(2);
			IDz		= IDx + IDy;

		end

		function midPoint = GetMidPoint(obj)

			midPoint = (obj.Node1.pos + obj.Node2.pos + obj.Node3.pos)/3;

		end

		function surfArea = GetArea(obj)

			% Returns the area of the surface
			% This can be found by taking the cross product
			% of two vectors that form adjacent sides of the triangle
			% Keep the same starting point for consistency

			v1 = obj.Node2.pos - obj.Node1.pos;
			v2 = obj.Node3.pos - obj.Node1.pos;

			surfArea = norm(cross(v1,v2))/2;

		end

		function ApplyForce(obj, F, A)

			if sum(isnan(F)) || sum(isinf(F))
				error('N:ApplyForce:InfNaN', 'Force is inf or NaN');
			end

			% Verify that A is on the surface

			% Applying a force to the surface at point A
			% This results in a linear movement, and a rotation

			% Linear movement
			obj.Node1.ApplyForce(F);
			obj.Node2.ApplyForce(F);
			obj.Node3.ApplyForce(F);

			% Rotation
			rD = obj.GetCentreOfDrag();

			rDA = A - rD;

			torque = cross(rDA, F);

			obj.ApplyTorque(torque);

		end

		function ApplyTorque(obj, torque)
			
			% Gather the torques applied to the surface.
			% This will be converted to a force on the nodes in the simulation
			obj.torque = obj.torque + torque;

		end

		function AddCell(obj, c)

			% c can be a vector
			if sum( ismember(c,obj.cellList)) ~=0
				warning('E:AddCell:CellAlreadyHere', 'Adding at least one cell that already appears in cellList for Surface %d. This has not been added.', obj.id);
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
				warning('E:RemoveCell:CellNotHere', 'At least one cell does not appear in surfaceList for Surface %d', obj.id);
			else
				obj.cellList(obj.cellList == c) = [];
			end

		end

		function internal = IsSurfaceInternal(obj)

			internal = obj.internal;

		end

		function AddEdgesToSurface(obj, e1, e2, e3)

			% Need to make sure the edges are connected and form a
			% closed loop

			% It is also important that the surface normal points in
			% a consistant direction when the surface is joined to
			% multiple other surfaces in a mesh, but we don't
			% consider that here.

			nodes = [e1.nodeList, e2.nodeList, e3.nodeList];

			% There must be exactly three unique nodes in this list,
			% any difference means the surface cannot be formed

			uNodes = unique(nodes);

			if length(uNodes) ~= 3
				error('E:ValidateEdges:Invalid', 'Edges are incorrectly formed');
			end

			% If it passes this test, then I'm almost certain the edges are
			% correctly formed for a surface. Then we just need to collate
			% the edges and nodes correctly in the local index system

			obj.Edge1 = e1;
			obj.Node1 = obj.Edge1.Node1;
			obj.Node2 = obj.Edge1.Node2;

			if obj.Node2 == e2.Node1
				% The next edge is e2
				obj.Edge2 = e2;
				obj.Edge3 = e3;
				obj.Node3 = e2.Node2;
			end

			if obj.Node2 == e2.Node2
				% The next edge is e2
				obj.Edge2 = e2;
				obj.Edge3 = e3;
				obj.Node3 = e2.Node1;
			end

			% Check if its in the
			if obj.Node2 == e3.Node1
				% The next edge is e2
				obj.Edge2 = e3;
				obj.Edge3 = e2;
				obj.Node3 = e3.Node2;
			end

			if obj.Node2 == e3.Node2
				% The next edge is e2
				obj.Edge2 = e3;
				obj.Edge3 = e2;
				obj.Node3 = e3.Node1;
			end

			obj.edgeList = [obj.Edge1, obj.Edge2, obj.Edge3];

			obj.nodeList = [obj.Node1,obj.Node2,obj.Node3];

			obj.Edge1.AddSurface(obj);
			obj.Edge2.AddSurface(obj);
			obj.Edge3.AddSurface(obj);

			obj.Node1.AddSurface(obj);
			obj.Node2.AddSurface(obj);
			obj.Node3.AddSurface(obj);

		end

	end

end