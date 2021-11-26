classdef EdgeLinearSpringForce < AbstractEdgeBasedForce
	% Treats the edge as a linear spring with a natural
	% length specified by the edges internal properties


	properties

		k

	end

	methods


		function obj = EdgeLinearSpringForce(k)

			obj.k = k;
			
		end

		function AddEdgeBasedForces(obj, edgeList)

			for i = 1:length(edgeList)
				
				e = edgeList(i);
				obj.ApplyForce(e);

			end

		end

		function ApplyForce(obj, e)


			% The edge will experience a force shrinking
			% for extending its length, pushing it towards
			% predefined natural length

			l = e.GetLength();
			n = e.naturalLength;
			u = e.GetUnitVector1to2();

			F = obj.k * (l - n) * u; % positive (l - n) means too long

			e.Node1.ApplyForce(F);
			e.Node2.ApplyForce(-F);


		end


	end

end