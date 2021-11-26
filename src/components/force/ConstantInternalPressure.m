classdef ConstantInternalPressure < AbstractCellBasedForce
	% A constant internal pressure for a cell made up of faces
	% There must be a clearly defined inside for this to work


	properties

		pressure

	end

	methods


		function obj = ConstantInternalPressure(p)

			obj.pressure = p;
			
		end

		function AddCellBasedForces(obj, cellList)

			for i = 1:length(cellList)
				
				c = cellList(i);
				obj.ApplyForce(c);

			end

		end

		function ApplyForce(obj, c)


			% For each face of the cell, calculate its area
			% to obtain the force on the face. Then apply
			% the force to each node in the direction of the
			% outward face normal

			% It is critical that we know which is the inside
			% so the face normals must be properly organised
			% before this point

			for i = 1:length(c.faceList)

				f = c.faceList(i);

				u = f.GetUnitNormal();

				a = f.GetArea();

				F = a * obj.pressure * u;

				f.Node1.ApplyForce(F);
				f.Node2.ApplyForce(F);
				f.Node3.ApplyForce(F);

			end


		end


	end

end