classdef MembraneVolume < AbstractSimulationData

	% Collects the volume of every membrnae in the simulation
	% assuming it is a closed surface
	properties 

		name = 'membraneVolume'
		data = {};

	end

	methods

		function obj = MembraneVolume
			% No special initialisation
			
		end

		function CalculateData(obj, t)

			% This is intended for a simulation that uses SphereShells
			% to represent a membrane. At this stage, only Tumour3D exists

			volData = [];

			for i = 1:length(t.cellList)

				c = t.cellList(i);
				if isa(c,'SphereShell')
					vol = c.GetVolume();
					volData = [vol, c.id];
				end

			end

			obj.data = {volData};

		end
		
	end

end