classdef NodeCellVolume < AbstractCellData
	% Calculates the area of the cell using 
	% neighbourhood searching

	properties 

		name = 'cellVolume'
		data = []

		simulation

		freeRadius		% The preferred seapration distance to a membrane

	end

	methods

		function obj = NodeCellVolume(simulation, freeRadius)
			
			% Need a pointer to the space simulation
			obj.simulation = simulation;

			obj.freeRadius = freeRadius;

		end

		function CalculateData(obj, c)

			if ~isa(c, 'NodeCell')
				error('NCA:NotANode','Can only use NodeCellVolume for NodeCells' );
			end

			n = c.nodeList;

			t = obj.simulation;

			nList = n.GetData('nodeNeighbours', t);
			fList = n.GetData('faceNeighbours', t);

			pos = c.nodeList.pos;
			sepN = [];

			for i = 1:length(nList)

				% Don't use the nodes on an edge in the calculation
				if isa(nList(i).cellList, 'NodeCell')
					posN = nList(i).pos;
					% Divide by two because this is centre to centre distance
					% meaning it accounts for the radius of two cells
					d = norm(pos - posN)/2;
					if d < obj.freeRadius
						sepN(end + 1) = d;
					end
				end
			end

			sepF = [];
			for i = 1:length(fList)

				f = fList(i);
				
				posN1 = f.Node1.pos;
				NtoN1 = pos - posN1;

				v = f.GetUnitNormal();
				d = abs(dot(NtoN1, v));
				if d < obj.freeRadius
					sepF(end + 1) = d;
				end

			end

			% We assume that it takes a minimum of 12 spheres to completely pack around a cell
			% This won't be true if there are growing cells neigbouring the cell of interest
			% but that gets really complicated.
			% Neighbouring nodes count for one sphere each
			% Faces are a bit different.

			% We trial three different ways of accountig for the membrane:
			% 1. Ignore it, and only caluclate volume based on cell neighbours
			% 2. If in contact with one or more membrane faces, there is a flat rate contribution
			% 3. Account for the number of membrane faces in an approximate way
			% (4. A detailed calculation taking into account the specific configuration - not tried yet) 
			
			% Type 1 is self explanitory
			% nSpheres = length(sepN);
			% Type 2 assumes that 3 sphere spots are taken up by the closest membrane face
			nSpheres = length(sepN) + 3 * (~isempty(sepF));
			% Type 3 assumes the first face takes 3 and subsequent faces take 2 sphere spots
			% nSpheres = length(sepN) + 3 * (~isempty(sepF)) + 2 * (length(sepF) - 1) * (length(sepF) > 1);

			freeSpheres = 0;
			if nSpheres < 12
				freeSpheres = 12 - nSpheres;
			end
			% [sepN, sepE, mean(sepE), obj.freeRadius*ones(1,freeSpheres) ]

			% To determine the approximate radius, we average out the values depending on the Type approach
			% Type 1: nothing special
			% radius = nanmean(  [sepN, obj.freeRadius*ones(1,freeSpheres) ]  );
			% Type 2: the smallest distance count for three spheres
			radius = nanmean(  [sepN, nanmean(sepF) * ones(1, 3), obj.freeRadius*ones(1,freeSpheres) ]  );
			% Type 3: each separation counts once for one sphere, then the remaining spheres are from the average of the separations
			%radius = nanmean(  [sepN, sepF, mean(sepF) * ones(1, 2 * length(sepF) + 1 ) , obj.freeRadius*ones(1,freeSpheres) ]  );
			
			% [radius, pi*radius^2]
			obj.data = (4*pi*radius^3)/3;

		end
		
	end

end