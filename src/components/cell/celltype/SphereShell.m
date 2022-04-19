classdef SphereShell < AbstractCell
	% A cell that is a thin shell unit sphere
	% The mesh points are hard coded
	% There are 58 nodes, 168 edges and  112 faces
	% N-E+F = 2 According to Euler's formula, which holds here

	properties


		n=[-0.4782   -0.3422   -0.8088
		    0.3695   -0.0525   -0.9277
		   -0.1244    0.0363   -0.9916
		    0.2483    0.3971   -0.8835
		   -0.6048   -0.6678   -0.4340
		   -0.1854   -0.7407   -0.6457
		    0.3010   -0.5668   -0.7669
		    0.6206   -0.6027   -0.5016
		   -0.8333   -0.1983   -0.5160
		   -0.0276   -0.3742   -0.9269
		    0.7725   -0.1747   -0.6105
		   -0.5761    0.1187   -0.8087
		    0.6905    0.2669   -0.6723
		   -0.5370    0.5949   -0.5981
		   -0.2131    0.4826   -0.8495
		    0.0301    0.7912   -0.6108
		    0.4819    0.6964   -0.5318
		   -0.2234   -0.9516   -0.2110
		    0.2348   -0.9041   -0.3571
		   -0.8733   -0.4746   -0.1099
		    0.6074   -0.7908   -0.0755
		   -0.9991   -0.0393    0.0186
		    0.9174   -0.3560   -0.1779
		   -0.9015    0.2510   -0.3526
		    0.9668    0.1001   -0.2353
		   -0.7124    0.6693   -0.2109
		    0.8011    0.5384   -0.2616
		   -0.3074    0.9158   -0.2585
		    0.1991    0.9624   -0.1848
		   -0.2198   -0.9422    0.2530
		    0.1935   -0.9736    0.1213
		   -0.6228   -0.7772    0.0902
		    0.8174   -0.5289    0.2283
		   -0.8506   -0.3850    0.3583
		    0.9686   -0.0689    0.2389
		   -0.8857    0.4292    0.1769
		    0.9073    0.3767    0.1866
		   -0.5606    0.8047    0.1953
		    0.6190    0.7801    0.0907
		   -0.1089    0.9833    0.1456
		    0.2715    0.8991    0.3433
		   -0.5045   -0.6599    0.5568
		   -0.1577   -0.5089    0.8462
		    0.1071   -0.7893    0.6046
		    0.5144   -0.7399    0.4335
		   -0.5850   -0.2523    0.7708
		    0.6909   -0.3058    0.6551
		   -0.8675    0.0810    0.4908
		   -0.1650    0.4304    0.8874
		    0.7623    0.1400    0.6319
		   -0.6005    0.5468    0.5834
		   -0.1897    0.7979    0.5721
		    0.2263    0.6067    0.7620
		    0.6137    0.5642    0.5523
		   -0.1144   -0.0370    0.9927
		    0.2857   -0.3172    0.9043
		   -0.5162    0.1674    0.8400
		    0.3115    0.1811    0.9328];

		  e=[20     5
		     5     9
		     9    20
		    23    21
		    21    33
		    33    23
		    41    54
		    54    53
		    53    41
		    52    41
		    53    52
		    24    12
		    12    14
		    14    24
		    32    20
		    20    34
		    34    32
		    45    47
		    47    33
		    33    45
		     6    10
		    10     1
		     1     6
		     8    21
		    23     8
		    23    11
		    11     8
		    48    46
		    46    34
		    34    48
		    39    54
		    41    39
		    31    44
		    44    45
		    45    31
		    31    30
		    30    44
		    33    35
		    35    47
		    14    15
		    15    12
		    15     3
		     3    12
		     4     2
		     2     3
		     3     4
		    10     2
		     2     7
		     7    10
		    52    38
		    38    51
		    51    52
		    41    40
		    40    52
		    12     9
		     9    24
		     9    22
		    22    24
		    28    38
		    38    40
		    40    28
		    32     5
		     5    18
		    18    32
		    34    42
		    42    46
		     1    12
		     3     1
		     1     5
		     9     1
		     7     6
		     6    19
		    19    18
		    18     6
		     8    19
		    19     7
		     7     8
		    50    37
		    37    54
		    54    50
		    50    47
		    35    50
		    25    37
		    37    35
		    35    25
		    56    45
		    47    56
		    16     4
		     4    15
		    15    16
		    16    17
		    17     4
		    30    32
		    18    30
		    44    42
		    42    30
		    46    57
		    57    48
		    48    36
		    36    22
		    22    48
		    57    51
		    51    48
		    39    37
		    37    27
		    27    39
		    43    46
		    42    43
		    43    56
		    56    55
		    55    43
		    18    31
		    31    19
		    19    21
		    21    31
		     2    11
		    11     7
		     2    13
		    13    11
		    36    38
		    38    26
		    26    36
		    36    51
		    24    36
		    26    24
		    17    13
		    13     4
		    46    55
		    55    57
		    58    50
		    54    58
		    53    58
		    58    47
		    58    56
		    58    55
		    25    23
		    11    25
		    27    25
		    25    13
		    13    27
		    17    27
		    41    29
		    29    40
		    39    29
		    39    17
		    17    29
		    16    29
		    28    14
		    14    26
		    26    28
		    28    16
		    16    14
		    28    29
		    49    58
		    53    49
		    49    55
		    52    49
		    49    57
		    51    49
		    20    22
		    34    22
		    10     3
		    23    35
		    32    42
		    45    21
		     5     6
		    56    44
		    44    43];

		f = [1     2     3
		     4     5     6
		     7     8     9
		    10     9    11
		    12    13    14
		    15    16    17
		    18    19    20
		    21    22    23
		    24     4    25
		    25    26    27
		    28    29    30
		    31     7    32
		    33    34    35
		    36    37    33
		    19    38    39
		    13    40    41
		    41    42    43
		    44    45    46
		    47    48    49
		    50    51    52
		    10    53    54
		    12    55    56
		    56    57    58
		    59    60    61
		    62    63    64
		    15     1    62
		    29    65    66
		    67    43    68
		    69     2    70
		    21    49    71
		    72    73    74
		    75    76    77
		    78    79    80
		    81    39    82
		    83    84    85
		    86    18    87
		    88    89    90
		    91    92    88
		    93    64    94
		    37    95    96
		    28    97    98
		    99   100   101
		    98   102   103
		   104   105   106
		   107    66   108
		   109   110   111
		    36    94   112
		   113   114   115
		   116   117    48
		   116   118   119
		   120   121   122
		   120   123    51
		   124   122   125
		   124    58   100
		    92   126   127
		    44   127   118
		    97   128   129
		   130    80   131
		   132   131     8
		   130   133    81
		   134    87   133
		   134   135   110
		   136    26   137
		   138   139   140
		   137   119   139
		   141   140   126
		   142   143    53
		   144   142    32
		   144   145   146
		   145   106   141
		   147   146    91
		   148   149   150
		   151   152   148
		   147   151   153
		   153    61   143
		   152    90    40
		   154   132   155
		   154   156   135
		   157   155    11
		   158   129   156
		   157    52   159
		   158   159   102
		     3   160    57
		    89    46    42
		    16   161   160
		    47   162    45
		     6   163    38
		    50    54    60
		    14   125   149
		    59   150   121
		    17   164    65
		    20   165     5
		    22    68   162
		    67    70    55
		    72    71    76
		   166    74    63
		    23   166    69
		    27    77   117
		    75    24   114
		    78    82    84
		    83   138   105
		   136    85   163
		    86   167    34
		    93    96   164
		    30   101   161
		    99   103   123
		   104    31    79
		   107   111   128
		   168   108    95
		   109   168   167
		   113   112    73
		    35   115   165];


		centre

	end

	methods
		
		function obj = SphereShell(pos, varargin)

			% pos: the location of the shell centre
			% varargin contains a vector [radius, ...]
			% radius: the radius of the shell

			r = 1;
			if ~isempty(varargin)
				r = varargin{1};
			end

			obj.centre = pos;

			obj.CellCycleModel = NoCellCycle();
			obj.CellCycleModel.colour = obj.CellCycleModel.colourSet.GetNumber('MEMBRANE');

			obj.n(:,1) = r * obj.n(:,1) + pos(1);
			obj.n(:,2) = r * obj.n(:,2) + pos(2);
			obj.n(:,3) = r * obj.n(:,3) + pos(3);

			for i = 1:length(obj.n)
				n = Node(obj.n(i,1),obj.n(i,2),obj.n(i,3));
				n.AddCell(obj);
				obj.nodeList(end + 1) = n;
			end

			for i = 1:length(obj.e)
				ns = obj.nodeList(obj.e(i,:));
				e = Edge(ns(1), ns(2));
				e.AddCell(obj);
				obj.edgeList(end + 1) = e;
			end

			for i = 1:length(obj.f)
				es = obj.edgeList(obj.f(i,:));
				f = Face(es(1), es(2), es(3));
				f.AddCell(obj);
				obj.faceList(end + 1) = f;
			end

			obj.OrientFaceNormals();

			cellDataArray = [CellVolume()];

			obj.AddCellData(cellDataArray);

		end

		function OrientFaceNormals(obj)

			% Need to make sure the face normals all point
			% outwards from the shell. To do this, use the centre
			% of the shell and point a vector from centre to a point
			% on the face. If the dot product is negative, swap the order
			% of Node2 and Node3 in the face, as these are used to
			% calculate the surface normal on the fly


			for i = 1:length(obj.faceList)

				f = obj.faceList(i);

				u = f.GetUnitNormal();

				ray = f.Node1.pos - obj.centre;

				if dot(ray,u) < 0
					% Unit normal is oriented the wrong way
					% Need to swap the local order of the nodes
					% to invert the normal
					n1 = f.Node1;
					n2 = f.Node3;
					n3 = f.Node2;

					f.Node2 = n2;
					f.Node3 = n3;

					f.nodeList = [n1,n2,n3];

				end

			end


		end


		function [newCell, newNodes, newEdges] = Divide(obj)

			newCell = AbstractCell.empty();
			newNodes = Node.empty();
			newEdges = Edge.empty();

		end

	end

end