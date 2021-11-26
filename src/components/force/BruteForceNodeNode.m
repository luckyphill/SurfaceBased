classdef BruteForceNodeNode < AbstractTissueBasedForce
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


		function obj = BruteForceNodeNode(springRate, dAsym, dSep, dLim)

			obj.springRate 	= springRate; 	% Spring force parameter
			obj.dAsym 		= dAsym;		% Repulsion asymptotes to inf
			obj.dSep 		= dSep;			% Preferred spearation
			obj.dLim 		= dLim;			% interaction limit
			

		end

		function AddTissueBasedForces(obj, t)

			% This is only for node cells at this point
			for i = 1:length(t.cellList)
				
				c = t.cellList(i);

				if isa(c, 'NodeCell')
					
					n1 = c.nodeList;

					nNeighours = n1.GetData('nodeNeighbours', t);

					for j = 1:length(nNeighours)

						n2 = nNeighours(j);

						if isa(n2.cellList, 'NodeCell')
							obj.ApplyForce(n1, n2);
						end

					end

				end

			end

		end

		function ApplyForce(obj, n1, n2)


			% For them to interact node to surf, the node must be in a region
			% normal to the surface
			% If it is not, then it may interact with an edge
			% If neither a surface or an edge, finally check node to node
			% For now, just focus on node surf interactions

			% First step is to see if the node is close enough to interact
			n1ton2 = n2.pos - n1.pos;

			d = norm(n1ton2);

			v = n1ton2/d;

			% This will tell us the distance of the node from the surface, signed
			% according to which side it is on
			sra = obj.springRate;
			srr = obj.springRate;
			da = obj.dAsym;
			ds = obj.dSep;
			dl = obj.dLim;

			if abs(d) < obj.dLim

				if (n1.cellList.sisterCell == n2.cellList)
					% If the nodes are sister cells, then if they are
					% in the same phase, they are actually representing one
					% cell, so we have to handle the ds slightly differently

					if (n1.cellList.GetAge() < n1.cellList.CellCycleModel.growthPhaseDuration)
						% ... but only if they are in the growth phase
						dds = n1.cellList.newFreeCellSeparation;
						fr = n1.cellList.CellCycleModel.GetGrowthPhaseFraction();
						ds = dds + fr * (ds - dds);

					end


				end

				% mag = ForceLaw(obj, d);
				mag = obj.ForceLaw(d, sra, srr, da, ds, dl);

				force = mag * v;

				n1.ApplyForce(-force);
				n2.ApplyForce(force);

			end

		end

		% function Fa = ForceLaw(obj, x)

		% 	% Basic non-linear spring law
		% 	% Returns negative values if its in the repulsion region
		% 	% Returns positive values in the attraction region

		% 	Fa = 0;

		% 	if (obj.dAsym < x) && ( x < obj.dSep)

		% 		Fa = obj.springRate * log(  ( obj.dSep - obj.dAsym ) / ( x - obj.dAsym )  );
		% 	end

		% 	if (obj.dSep <= x ) && ( x < obj.dLim )

		% 		Fa = obj.springRate * (  ( obj.dSep - x ) / ( obj.dSep - obj.dAsym )  ) * exp(obj.c*(obj.dSep - x)/obj.dSep );

		% 	end

		% end

		function Fa = ForceLaw(obj, x, sra, srr, da, ds, dl)

			% This calculates the scalar force for the given separation
			% x and the controlling parameters as outline in the preamble

			% It calculates the force when it is internal to a cell

			Fa = 0;

			if (da < x) && ( x < ds)

				Fa = srr * log(  ( ds - da ) / ( x - da )  );
			end

			if (ds <= x ) && ( x < dl)

				Fa = sra * (  ( ds - x ) / ( ds - da )  ) * exp(obj.c*(ds - x)/ds );

			end

		end

	end

end