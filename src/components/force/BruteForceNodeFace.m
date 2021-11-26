classdef BruteForceNodeFace < AbstractTissueBasedForce
	% This adds a constant force to each NodeCell in the simulation
	% that points radially form a point.
	% force is a scalar magnitude - positive means it pushes away from
	% the centre, negative means it pushes towards the centre


	properties

		springRate
		dAsym
		dSep
		dLim

		c = 5

	end

	methods


		function obj = BruteForceNodeFace(springRate, dAsym, dSep, dLim)

			obj.springRate 	= springRate; 	% Spring force parameter
			obj.dAsym 		= dAsym;		% Repulsion asymptotes to inf
			obj.dSep 		= dSep;			% Preferred spearation
			obj.dLim 		= dLim;			% interaction limit
			

		end

		function AddTissueBasedForces(obj, t)

			for i = 1:length(t.cellList)
				
				c = t.cellList(i);

				if isa(c, 'NodeCell')
					
					n = c.nodeList;

					fNeighbours = n.GetData('faceNeighbours', t);
					
					for j = 1:length(fNeighbours)

						f = fNeighbours(j);

						obj.ApplyForce(n,f);

					end
				end
			end

		end

		function ApplyForce(obj, n, s)


			% For them to interact node to surf, the node must be in a region
			% normal to the surface
			% If it is not, then it may interact with an edge
			% If neither a surface or an edge, finally check node to node
			% For now, just focus on node surf interactions

			% First step is to see if the node is close enough to interact
			n1ton = n.pos - s.Node1.pos;

			% Project it onto the surface normal
			un = s.GetUnitNormal();

			% This force calculator is being used for the Tumour3D model
			% at this point, and since the SphereShell
			% has the face normals oriented outwards, the signed position
			% should be negative unless the node is outside the cell

			% This means we need to calculated the force based on -ve position
			% being the correct place. The easiest way to do this is to flip the normal

			un = -un;

			d = dot(n1ton, un);

			% This will tell us the distance of the node from the surface, signed
			% according to which side it is on, where +ve is now inside

			if d < obj.dLim

				% If it is close enough to the plane of the surface,
				% then we need to decide if it is within the boundaries of
				% the surface.

				% We can do this by finding the point of contact, and
				% testing if it is in a polygon in the local coordinate
				% system on the surface

				% The point where the force is applied
				A = n.pos - d*un;
				rD = s.GetCentreOfDrag();

				% Centre of drag to point of contact in global axes
				rDA = A - rD;

				% Transform to local coordinate system
				% The z coordinate must be zero, if not, then there
				% has been an error in making the point of contact
				[u,v,w] = s.GetCoordinateSystem();
				rA = [dot(rDA,u),dot(rDA,v),0];

				if abs(dot(rDA,w)) > 1e-6
					error('The local z coordinate of the contact point should be zero');
				end

				[r1, r2, r3] = s.GetNodeCoordinatesLocal();
				p = [r1;r2;r3];

				in = inpolygon(rA(1), rA(2), p(:,1), p(:,2));

				if in

					% To prevent attraction between nodes and
					% faces of the same cell, only do the next bit
					% if d is less than dSep

					if s.cellList ~= n.cellList || d < obj.dSep
						% Calculate the force, and apply torques
						f = obj.ForceLaw(d);

						force = f * un * sign(d);

						% Apply the force to the surface at point A
						s.ApplyForce(-force, A);

						% Apply the force to the interacting node

						n.ApplyForce(force);

					end

				end

			end

		end

		function Fa = ForceLaw(obj, x)

			% Basic non-linear spring law
			% Returns negative values if its in the repulsion region
			% Returns positive values in the attraction region

			Fa = 0;

			if (obj.dAsym < x) && ( x < obj.dSep)

				Fa = obj.springRate * log(  ( obj.dSep - obj.dAsym ) / ( x - obj.dAsym )  );
			end

			if (obj.dSep <= x ) && ( x < obj.dLim )

				Fa = obj.springRate * (  ( obj.dSep - x ) / ( obj.dSep - obj.dAsym )  ) * exp(obj.c*(obj.dSep - x)/obj.dSep );

			end

		end

	end

end