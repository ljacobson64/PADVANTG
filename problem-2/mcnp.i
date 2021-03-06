BNCT modeling
c
c  Author:       Lucas Jacobson
c  Created on:   March 9, 2021
c  Last updated: October 1, 2021
c
c  Copyright 2021, Phoenix, LLC. All rights reserved.
c
c ******************************************************************************
c
c                                  CELL CARDS
c
c ******************************************************************************
c
c ------------------------------------------------------------------------------
c     Tally cells
c ------------------------------------------------------------------------------
c
  1001   4  -0.001205  -1001                                                    $ Collimator opening tally cell
c
  1101   4  -0.001205  -1101                                                    $ Tumor (replaced with air)
  1102   4  -0.001205  -1102 1101                                               $ Brain (replaced with air)
  1103   4  -0.001205  -1103 1102                                               $ Skull (replaced with air)
  1104   4  -0.001205  -1104 1103                                               $ Skin  (replaced with air)
c
c ------------------------------------------------------------------------------
c     Beam shaping assembly
c ------------------------------------------------------------------------------
c
  2001 705  -2.967686  -2001 2203                                               $ 70% AlF3 30% Al mixture
  2004  55  -8.65      -2004                                                    $ Cadmium thermal shield
  2005  28  -9.747     -2005                                                    $ Bismuth gamma shield
  2006   4  -0.001205  -2006                                                    $ Conical air beam path
  2007  55  -8.65      -2007 2000 2203                                          $ Cadmium liner
c
  2011 189 -11.35      -2011 2000 2006 2007 2203                                $ Lead reflector
c
  2021   4  -0.001205  -2021 1001 1104                                          $ Patient area
  2022 189 -11.35      -2022 2021                                               $ Lead shielding for patient area
c
c ------------------------------------------------------------------------------
c     Target and beamline
c ------------------------------------------------------------------------------
c
  2201   0             -2201                                                    $ Vacuum in target chamber
  2202   0             -2202                                                    $ Vacuum in beamline
  2203  13  -2.7       -2203 2201 2202                                          $ Beamline
c
c ------------------------------------------------------------------------------
c     Bunker
c ------------------------------------------------------------------------------
c
  5000   4  -0.001205  -5000 5050                                               $ Surrounding air
c
  5051  99  -4.5       -5051                                                    $ Heavy concrete West
  5052  99  -4.5       -5052 2203                                               $ Heavy concrete South
  5053  99  -4.5       -5053 2203                                               $ Heavy concrete North
  5054  99  -4.5       -5054                                                    $ Heavy concrete below
  5055  99  -4.5       -5055                                                    $ Heavy concrete above
c
c ------------------------------------------------------------------------------
c     Dummy cells
c ------------------------------------------------------------------------------
c
  9000   4  -0.001205  -9000 9001 9002 9003 9004 9005 9006 9007 9008 9009 9010
                             9011 9012 9013 9014 9015 9016 9017 9018 9019 9020
                             9021 9022 9023 9024 9025 9026 9027 9028 9029 9030
                             9031 9032 9033 9034 9035 9036 9037 9038            $ Bounding cell
c
  9001   4  -0.001205  -9001                                                    $ Air (dry, near sea level)
  9002   6  -2.6989    -9002                                                    $ Aluminum (Al)
  9003   7  -3.97      -9003                                                    $ Aluminum Oxide (Al2O3)
  9004  13  -2.7       -9004                                                    $ Aluminum, alloy 6061-O
  9005  25  -1.848     -9005                                                    $ Beryllium (Be)
  9006  27  -3.01      -9006                                                    $ Beryllium Oxide (BeO)
  9007  28  -9.747     -9007                                                    $ Bismuth (Bi)
  9008  55  -8.65      -9008                                                    $ Cadmium (Cd)
  9009  61  -3.18      -9009                                                    $ Calcium Fluoride (CaF2)
  9010  62  -3.3       -9010                                                    $ Calcium Oxide (CaO)
  9011  69  -1.7       -9011                                                    $ Carbon, Graphite (reactor grade)
  9012  99  -4.5       -9012                                                    $ Concrete, M-1
  9013 109  -2.3       -9013                                                    $ Concrete, Regulatory Concrete (developed for U.S. NRC)
  9014 116  -1.52      -9014                                                    $ Earth, U.S. Average
  9015 189 -11.35      -9015                                                    $ Lead (Pb)
  9016 192  -0.534     -9016                                                    $ Lithium (Li)
  9017 194  -2.635     -9017                                                    $ Lithium Fluoride (LiF)
  9018 199  -2.013     -9018                                                    $ Lithium Oxide (Li2O)
  9019 207  -1.74      -9019                                                    $ Magnesium (Mg)
  9020 208  -3.58      -9020                                                    $ Magnesium Oxide (MgO)
  9021 227  -8.902     -9021                                                    $ Nickel (Ni)
  9022 275  -0.93      -9022                                                    $ Polyethylene, Non-borated (C2H4)
  9023 362  -4.54      -9023                                                    $ Titanium (Ti)
  9024 363  -4.26      -9024                                                    $ Titanium Dioxide (TiO2)
  9025 391  -1.1044    -9025                                                    $ Water, Heavy (H2O)
  9026 392  -0.997     -9026                                                    $ Water, Liquid (H2O)
  9027 417  -1.54      -9027                                                    $ Calcium (Ca)
  9028 419  -6.11      -9028                                                    $ Vanadium (V)
  9029 505  -8.902     -9029                                                    $ Nickel-60
  9030 511  -3.148     -9030                                                    $ Magnesium Fluoride (MgF2)
  9031 512  -3.1       -9031                                                    $ Aluminum Fluoride (AlF3)
  9032 513  -3.4       -9032                                                    $ Titanium(III) Fluoride (TiF3)
  9033 521  -2.962681  -9033                                                    $ Fluental
  9034 601  -0.95      -9034                                                    $ 5% Borated Polyethylene (SWX-201)
  9035 602  -1.19      -9035                                                    $ 30% Boron Polyethylene (SWX-210)
  9036 603  -1.06      -9036                                                    $ 7.5% Lithium Polyethylene (SWX-215)
  9037 604  -2.92      -9037                                                    $ Poly-Biz Gamma Shield (SWX-217)
  9038 705  -2.967686  -9038                                                    $ 70% AlF3, 30% Al
c
c ------------------------------------------------------------------------------
c     Rest of universe
c ------------------------------------------------------------------------------
c
  9999   0              5000 9000
c
c ------------------------------------------------------------------------------

c ******************************************************************************
c
c                                SURFACE CARDS
c
c ******************************************************************************
c
c ------------------------------------------------------------------------------
c     Tally cells
c ------------------------------------------------------------------------------
c
  1001 10 rcc    95.3       0 0       0.1       0 0       7.5
c
  1101 11 sx     -4.8          1.0
  1102 11 sq      0.0236686391 0.0123456790 0.0277777778 0 0 0 -1 -1 0 0
  1103 11 sq      0.0145158949 0.0104123282 0.0216262976 0 0 0 -1  0 0 0
  1104 11 sq      0.0129132231 0.0094259591 0.0187652468 0 0 0 -1  0 0 0
c
c ------------------------------------------------------------------------------
c     Beam shaping assembly
c ------------------------------------------------------------------------------
c
  2000 20 rcc   -40.0       0 0      85.3       0 0      35.0
  2001 20 rcc   -40.0       0 0      80.0       0 0      35.0
  2004 20 rcc    40.0       0 0       0.3       0 0      35.0
  2005 20 rcc    40.3       0 0       5.0       0 0      35.0
  2006 20 trc    45.3       0 0      50.0       0 0      35.0       7.5
  2007 20 rcc   -40.0       0 0      80.3       0 0      35.1
c
  2011 20 rpp   -65.0      95.3     -60.0      60.0     -60.0      60.0
c
  2021 20 rcc    95.3       0 0      79.7       0 0      40.0
  2022 20 rpp    95.3     175.0     -60.0      60.0     -60.0      60.0
c
c ------------------------------------------------------------------------------
c     Target and beamline
c ------------------------------------------------------------------------------
c
  2201 22 rcc     0.0     -55.0       0 0     100.0       0.0       5.0
  2202 22 rcc     0.0    -120.0       0 0      65.0       0.0       5.0
  2203 22 rcc     0.0    -120.0       0 0     165.3       0.0       5.3
c
c ------------------------------------------------------------------------------
c     Bunker
c ------------------------------------------------------------------------------
c
  5000 50 rpp -30000.0  30000.0  -30000.0   30000.0  -30000.0   30000.0
c
  5050 50 rpp  -125.0     175.0    -120.0     120.0    -120.0     120.0
  5051 50 rpp  -125.0     -65.0     -60.0      60.0     -60.0      60.0
  5052 50 rpp  -125.0     175.0    -120.0     -60.0     -60.0      60.0
  5053 50 rpp  -125.0     175.0      60.0     120.0     -60.0      60.0
  5054 50 rpp  -125.0     175.0    -120.0     120.0    -120.0     -60.0
  5055 50 rpp  -125.0     175.0    -120.0     120.0      60.0     120.0
c
c ------------------------------------------------------------------------------
c     Dummy cells
c ------------------------------------------------------------------------------
c
   9000   rpp  9000.0    9039.0      -1.0       1.0      -1.0       1.0
c
   9001   sx   9001.0       0.4
   9002   sx   9002.0       0.4
   9003   sx   9003.0       0.4
   9004   sx   9004.0       0.4
   9005   sx   9005.0       0.4
   9006   sx   9006.0       0.4
   9007   sx   9007.0       0.4
   9008   sx   9008.0       0.4
   9009   sx   9009.0       0.4
   9010   sx   9010.0       0.4
   9011   sx   9011.0       0.4
   9012   sx   9012.0       0.4
   9013   sx   9013.0       0.4
   9014   sx   9014.0       0.4
   9015   sx   9015.0       0.4
   9016   sx   9016.0       0.4
   9017   sx   9017.0       0.4
   9018   sx   9018.0       0.4
   9019   sx   9019.0       0.4
   9020   sx   9020.0       0.4
   9021   sx   9021.0       0.4
   9022   sx   9022.0       0.4
   9023   sx   9023.0       0.4
   9024   sx   9024.0       0.4
   9025   sx   9025.0       0.4
   9026   sx   9026.0       0.4
   9027   sx   9027.0       0.4
   9028   sx   9028.0       0.4
   9029   sx   9029.0       0.4
   9030   sx   9030.0       0.4
   9031   sx   9031.0       0.4
   9032   sx   9032.0       0.4
   9033   sx   9033.0       0.4
   9034   sx   9034.0       0.4
   9035   sx   9035.0       0.4
   9036   sx   9036.0       0.4
   9037   sx   9037.0       0.4
   9038   sx   9038.0       0.4
c
c ------------------------------------------------------------------------------

c ******************************************************************************
c
c                                  DATA CARDS
c
c ******************************************************************************
c
c ------------------------------------------------------------------------------
c     Simulation parameters
c ------------------------------------------------------------------------------
c
  nps     1e6
  prdmp   j 1e5 1 2 j
c
  mode    n p
  dbcn    77j 1
  imp:n   1 60r 0
  imp:p   1 60r 0
c
c ------------------------------------------------------------------------------
c     Coordinate transforms
c ------------------------------------------------------------------------------
c
  *tr10    0.0        0.0        0.0                                            $ Tally cell
  *tr11  105.0        0.0        0.0                                            $ Brain phantom
  *tr20    0.0        0.0        0.0                                            $ Moderator
  *tr22    0.0        0.0        0.0                                            $ Target
  *tr50    0.0        0.0        0.0                                            $ Bunker
c
c ------------------------------------------------------------------------------
c     Read files
c ------------------------------------------------------------------------------
c
  read file = mcnp_tallies.txt
  read file = mcnp_materials.txt
  read file = mcnp_source.txt
c
c ------------------------------------------------------------------------------
c     Added by ADVANTG
c ------------------------------------------------------------------------------
c
