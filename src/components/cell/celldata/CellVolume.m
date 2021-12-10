classdef CellVolume < AbstractCellData
	% Calculates the volume of a cell in 3D
	% can also be used for membranes
	% This will implicitly assume that the points provided
	% form a closed surface

	properties 

		name = 'cellVolume'
		data = []

	end

	methods

		function obj = CellVolume
			% No special initialisation
			
		end

		function CalculateData(obj, c)
			% This is not ideal, but it should do for now
			% We use the built in function "alphShape" to calculate
			% the volume of the cell/membrane
			% alphaShape takes a set of points and produces a shape
			% from which we calculate the volume

			% alphaShape handles non-convex shapes, which absolutley can
			% happen for our sims, but for the time being with Tumour3D
			% as the only model, we will be dealing with convex shapes
			% the vast majority of the time

			% ************************************************************
			% WARNING!! WARNING!! WARNING!! WARNING!! WARNING!! WARNING!!
			% ************************************************************
			% We are implicitly assuming that alphaShape produces
			% the same shape given by the faces of the cell
			% but this is by no means guranteed
			% Basic testing showed that it can handle the initial shape
			% of SphereShell, as long as the alpha value is set correctly
			% But there has been no thorough testing, particularly with
			% non-convex shapes. Longer term, a custom code is needed
			% to account for the fact that we already know the faces

			% When running alphaShape without specifying alpha, it
			% defaults to a value that produces a bad, non closed surface
			% To avoid this, we use the function criticalAlpha to find the
			% minimum alpha required to produce something, then double it in
			% the hopes that it captures everything properly
			
			shp = alphaShape(reshape([c.nodeList.pos],3,[])');
			a = criticalAlpha(shp,'one-region');
			shp = alphaShape(reshape([c.nodeList.pos],3,[])',2*a);
			obj.data = volume(shp);

		end
		
	end

end