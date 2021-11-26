classdef AbstractEdgeBasedForce < matlab.mixin.Heterogeneous
	% This class gives the details for how a force will be applied
	% to each Edge (as opposed to each cell, or the whole population)


	properties


	end

	methods (Abstract)

		AddEdgeBasedForces(obj, edgeList)

	end



end