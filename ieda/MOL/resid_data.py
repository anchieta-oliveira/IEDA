#!/usr/bin/env python3

# Description
###############################################################################

""" Module to define the atoms of the amino acids.

This module contains the atoms of the amino acids, which are used to build the
molecular structure of the proteins. The atoms are defined as instances of the
Atom class, which is used to store and manipulate data from atoms.

"""

# Imports
###############################################################################
from ieda.MOL.atom import Atom

# License
###############################################################################
'''
IEDA
Authors: José de Anchieta de Oliveira Filho
[The Federal University of Rio de Janeiro]
Carlos Chagas Filho Institute of Biophysics
Laboratory for Molecular Modeling and Dynamics
Av. Carlos Chagas Filho 373 - CCS - bloco G1-19,
Cidade Universitária - Rio de Janeiro, RJ, CEP: 21941-902
E-mail address: anchieta.oliveira@biof.ufrj.br
This project is licensed under Creative Commons license (CC-BY-4.0) (Ver qual)
'''

# Classes
###############################################################################

# Functions
###############################################################################
## Private ##

## Public ##

# Aliases
###############################################################################
N = Atom(name="N", element= "N", coordinates=(-12.058, 34.376, 59.305))
HN = Atom(name="HN", element="H", coordinates=(-11.727, 35.343, 59.575))
HA = Atom(name="HA", element="H", coordinates=(-13.209, 33.188, 58.161))
CA = Atom(name="CA", element="C", coordinates=(-12.975, 34.150, 58.200 ))
C = Atom(name="C", element="C", coordinates=( -12.224, 34.524, 56.929 ))
O = Atom(name="O", element="O", coordinates=(-12.015, 35.710, 56.665))
CB = Atom(name = "CB", element="C", coordinates=(-14.237, 34.992, 58.377))

backbone = [N, HN, HA, CA, C, O]

resids_data = {
                #    N     HB1
                #    |    / 
                # HA-CA--CB-HB2
                #    |    \ 
                #  O=C     HB3                                                                                                       
            "ALA": [N, CA, C, O, HN, HA, #backbone                                                                    
                    CB, Atom(name="HB1", element="H", coordinates=(-14.895, 34.757, 57.663)),                         
                    Atom(name="HB2", element="H", coordinates=(-14.636, 34.807, 59.274)),                             
                    Atom(name="HB3", element="H", coordinates=(-14.002, 35.960, 58.309))], #SideChain                 


                #                           HH11    
                #                            | 
                #    N                      NH1-H12
                #    |                    //(+)
                # HA-CA--CB--CG--CD--NE--CZ
                #    |                     \
                #  O=C                      NH2-HH22
                #                             |
                #                            HH21
            "ARG": [Atom(name="N", element= "N", coordinates=(-17.231, 22.781, 72.137)), Atom(name="HN", element="H", coordinates=(-18.095, 22.267, 72.126)), 
                    Atom(name="HA", element="H", coordinates=(-16.424, 24.331, 73.176)), Atom(name="CA", element="C", coordinates=(-16.585, 23.340, 73.329)), 
                    Atom(name="C", element="C", coordinates=(-15.235, 22.681, 73.667)), Atom(name="O", element="O", coordinates=(-14.422, 23.314, 74.356)), #backbone                                                                    
                    Atom(name = "CB", element="C", coordinates=(-17.548, 23.296, 74.529)), Atom(name="HB1", element="H", coordinates=(-17.158, 23.845, 75.257)), 
                    Atom(name="HB2", element="H", coordinates=(-18.421, 23.666, 74.238)), Atom(name="CG", element="C", coordinates=(-17.784, 21.889, 75.067)), 
                    Atom(name="HG1", element="H", coordinates=(-18.204, 21.329, 74.359)), Atom(name="HG2", element="H", coordinates=(-16.915, 21.491, 75.348)),       
                    Atom(name="CD", element="C", coordinates=(-18.713, 21.897, 76.278)), Atom(name="HD1", element="H", coordinates=(-18.275, 22.409, 77.000)), 
                    Atom(name="HD2", element="H", coordinates=(-19.568, 22.312, 76.007)), Atom(name="NE", element="N", coordinates=(-18.962, 20.524, 76.737)), 
                    Atom(name="HE", element="H", coordinates=(-18.619, 19.771, 76.176)), Atom(name="CZ", element="C", coordinates=(-19.618, 20.233, 77.857)),                                                                     
                    Atom(name="NH1", element="N", coordinates=(-20.043, 21.238, 78.629)), Atom(name="HH11", element="H", coordinates=(-19.902, 22.186, 78.383)), 
                    Atom(name="HH12", element="H", coordinates=(-20.525, 21.018, 79.492)), Atom(name="NH2", element="N", coordinates=(-19.765, 18.971, 78.279)), 
                    Atom(name="HH21", element="H", coordinates=(-19.377, 18.238, 77.707)), Atom(name="HH22", element="H", coordinates=(-20.238, 18.772, 79.119))],   



                #    N   HB1        HD21
                #    |   |         /
                # HA-CA--CB--CG--ND2
                #    |   |         \
                #  O=C   HB2        HD22                                                                                                      
            "ASN": [Atom(name="N", element= "N", coordinates=(-12.973, 24.050, 63.274)), Atom(name="HN", element="H", coordinates=(-13.168, 23.073, 63.507)), 
                    Atom(name="HA", element="H", coordinates=(-13.296, 25.693, 62.179)), Atom(name="CA", element="C", coordinates=(-13.846, 24.903, 62.477)), 
                    Atom(name="C", element="C", coordinates=(-14.989, 25.440, 63.311)), Atom(name="O", element="O", coordinates=(-15.427, 24.830, 64.292)), #backbone                                                                    
                    Atom(name = "CB", element="C", coordinates=(-14.446, 24.110, 61.298)), Atom(name="HB1", element="H", coordinates=(-15.055, 23.420, 61.670)), 
                    Atom(name="HB2", element="H", coordinates=(-14.946, 24.746, 60.723)), Atom(name="CG", element="C", coordinates=(-13.379, 23.431, 60.471)), 
                    Atom(name="OD1", element="O", coordinates=(-12.668, 24.059, 59.682)), Atom(name="ND2", element="N", coordinates=(-13.224, 22.128, 60.681)), 
                    Atom(name="HD21", element="H", coordinates=(-12.529, 21.619, 60.168)), Atom(name="HD22", element="H", coordinates=(-13.796, 21.651, 61.347))],                
 

                #    N   HB1   OD1
                #    |   |    //
                # HA-CA--CB--CG
                #    |   |     \
                #  O=C   HB2    OD2(-)                                                                                                                         
            "ASP":[Atom(name="N", element= "N", coordinates=(-25.616, -1.916, 66.160)), Atom(name="HN", element="H", coordinates=(-24.916, -2.666, 65.947)), 
                   Atom(name="HA", element="H", coordinates=(-26.924, -0.529, 65.559)), Atom(name="CA", element="C", coordinates=(-26.341, -1.200, 65.111)),
                   Atom(name="C", element="C", coordinates=(-27.256, -2.198, 64.415)),  Atom(name="O", element="O", coordinates=( -27.144, -3.410, 64.609)), #backbone                                                                     
                   Atom(name="CB", element="C", coordinates=(-25.392, -0.579, 64.094)), Atom(name="HB1", element="H", coordinates=(-24.571, -0.282, 64.574)), 
                   Atom(name="HB2", element="H", coordinates=(-25.160, -1.273, 63.419)), Atom(name="CG", element='C', coordinates=(-25.984, 0.623, 63.365)),                                                                      
                   Atom(name="OD1", element="O", coordinates=(-26.927, 1.294, 63.875)), Atom(name="OD2", element="O", coordinates=(-25.456, 0.949, 62.297))],                                                                      


                #    N   HB1   OD1
                #    |   |    //
                # HA-CA--CB--CG
                #    |   |     \
                #  O=C   HB2    OD2-HD2                                                                                                                         
            "ASPP":[Atom(name="N", element= "N", coordinates=(-25.616, -1.916, 66.160)), Atom(name="HN", element="H", coordinates=(-24.916, -2.666, 65.947)), 
                   Atom(name="HA", element="H", coordinates=(-26.924, -0.529, 65.559)), Atom(name="CA", element="C", coordinates=(-26.341, -1.200, 65.111)),
                   Atom(name="C", element="C", coordinates=(-27.256, -2.198, 64.415)),  Atom(name="O", element="O", coordinates=( -27.144, -3.410, 64.609)), #backbone                                                                     
                   Atom(name="CB", element="C", coordinates=(-25.392, -0.579, 64.094)), Atom(name="HB1", element="H", coordinates=(-24.571, -0.282, 64.574)), 
                   Atom(name="HB2", element="H", coordinates=(-25.160, -1.273, 63.419)), Atom(name="CG", element='C', coordinates=(-25.984, 0.623, 63.365)),                                                                      
                   Atom(name="OD1", element="O", coordinates=(-26.927, 1.294, 63.875)), Atom(name="OD2", element="O", coordinates=(-25.456, 0.949, 62.297)),
                   Atom(name="HD2", element="H", coordinates=(-26.456, 0.949, 62.297))
                   ], 

                #    N
                #    |
                # HA-CA--CB--SG
                #    |         \
                #  O=C          HG1
            "CYS": [Atom(name="N", element= "N", coordinates=(-20.376, 30.060, 55.420)), Atom(name="HN", element="H", coordinates=(-20.447, 29.154, 54.899)), 
                    Atom(name="HA", element="H", coordinates=(-20.589, 31.138, 57.111)), Atom(name="CA", element="C", coordinates=(-20.498, 30.157, 56.864)), 
                    Atom(name="C", element="C", coordinates=(-19.287, 29.600, 57.606)), Atom(name="O", element="O", coordinates=(-19.229, 29.742, 58.828)), #backbone                                                                    
                    Atom(name="CB", element="C", coordinates=(-21.784, 29.430, 57.357)), Atom(name="HB1", element="H", coordinates=(-21.706, 28.482, 57.090)), 
                    Atom(name="HB2", element="H", coordinates=(-21.815, 29.524, 58.341)), Atom(name="SG", element="S", coordinates=(-23.272, 30.124, 56.650)),                                                                     
                    Atom(name="HG1", element="H", coordinates=(-23.436, 29.731, 55.741))],                                                                   
                                                                                                                      
 
            "CYX": [], 
 

                #    N   HB1 HG1 OE1    HE21
                #    |   |   |   ||    /
                # HA-CA--CB--CG--CD--NE2
                #    |   |   |         \
                #  O=C   HB2 HG2       HE22            
            "GLN": [Atom(name="N", element= "N", coordinates=(-12.596, 11.907, 77.489)), Atom(name="HN", element="H", coordinates=(-12.218, 12.607, 78.154)),
                    Atom(name="HA", element= "H", coordinates=(-12.099, 10.510, 76.132)), Atom(name="CA", element= "C", coordinates=(-11.876, 10.702, 77.090)), 
                    Atom(name="C", element= "C", coordinates=(-12.350, 9.519, 77.926)), Atom(name="O", element= "O", coordinates=(-12.161, 9.488, 79.147)), #backbone                                                                    
                    Atom(name="CB", element="C", coordinates=(-10.371, 10.891, 77.266)), Atom(name="HB1", element="H", coordinates=(-10.186, 11.114, 78.222)), 
                    Atom(name="HB2", element="H", coordinates=(-10.073, 11.642, 76.681)), Atom(name="CG", element="C", coordinates=(-9.582, 9.658, 76.902)), 
                    Atom(name="HG1", element="H", coordinates=(-9.892, 8.896, 77.464)), Atom(name="HG2", element="H", coordinates=(-8.615, 9.832, 77.062)),       
                    Atom(name="CD", element="C", coordinates=(-9.772, 9.295, 75.447)), Atom(name="OE1", element="O", coordinates=(-9.763, 10.173, 74.583)),                                      
                    Atom(name="NE2", element="N", coordinates=(-9.962, 8.003, 75.165)), Atom(name="HE21", element="H", coordinates=(-10.110, 7.692, 74.215)), 
                    Atom(name="HE22",element="H", coordinates=(-9.954, 7.333, 75.934))],    
                

                #    N   HB1 HG1   OE1
                #    |   |   |    //
                # HA-CA--CB--CG--CD
                #    |   |   |    \
                #  O=C   HB2 HG2   OE2(-)
            "GLU": [Atom(name="N", element="N", coordinates=(-20.660, 26.728, 52.697)),Atom(name="HN", element="H", coordinates=(-21.199, 26.270, 51.970)),
                    Atom(name="HA", element="H", coordinates=(-18.666, 26.641, 53.075)),Atom(name="CA", element="C", coordinates=(-19.300, 27.230, 52.552)),
                    Atom(name="C", element="C", coordinates=(-19.179, 28.660, 53.034)),Atom(name="O", element="O", coordinates=(-18.196, 29.015, 53.689)), #backbone                                                                    
                    Atom(name="CB", element="C", coordinates=(-18.867, 27.143, 51.087)), Atom(name="HB1", element="H", coordinates=(-17.984, 27.567, 51.008)), 
                    Atom(name="HB2", element="H", coordinates=(-19.549, 27.598, 50.544)), Atom(name="CG", element="C", coordinates=(-18.739, 25.715, 50.561)), 
                    Atom(name="HG1", element="H", coordinates=(-19.565, 25.236, 50.752)), Atom(name="HG2", element="H", coordinates=(-17.946, 25.309, 50.954)),       
                    Atom(name="CD", element="C", coordinates=(-18.548, 25.697, 49.030)), Atom(name="OE1", element="O", coordinates=(-17.661, 26.388, 48.475)), 
                    Atom(name="OE2", element="O", coordinates=(-19.316, 24.949, 48.394))],      
                                                                                                                      

            "GLY": [Atom(name="N", element="N", coordinates=(-20.165, 29.497, 52.721)), Atom(name="HN", element="H", coordinates=(-20.944, 29.103, 52.200)), 
                    Atom(name="CA", element="C", coordinates=(-20.054, 30.884, 53.140)), Atom(name="C", element="C", coordinates=(-20.163, 31.119, 54.639)), 
                    Atom(name="O", element="O", coordinates=(-20.006, 32.274, 55.069)),  Atom(name="HA1", element="H", coordinates=(-20.776, 31.436, 52.665)),  
                    Atom(name="HA", element="H", coordinates=(-19.160, 31.263, 52.809))], #backbone 


                #               HD1  HE1 
                #    N   HB1    ND1--CE1 
                #    |   |     /     ||
                # HA-CA--CB--CG      ||
                #    |   |     \\    || 
                #  O=C   HB2    CD2--NE2
                #               HD2              neutral HIS, proton on ND1                                                                                            
            "HSD":[Atom(name="N", element="N", coordinates=(-10.390, 37.499, 74.194)),Atom(name="HN", element="H", coordinates=(-10.506, 37.955, 75.102)), 
                   Atom(name="HA", element="H", coordinates=(-9.409, 37.811, 72.452)), Atom(name="CA", element="C", coordinates=(-9.164, 37.532, 73.399)), 
                   Atom(name="C", element="C", coordinates=(-8.447, 36.191, 73.328)), Atom(name="O", element="O", coordinates=(-7.538, 36.041, 72.503)), #backbone                                                                                       
                   Atom(name="CB", element="C", coordinates=(-8.207, 38.606, 73.946 )), Atom(name="HB1", element="H", coordinates=(-8.725, 39.455, 74.090)), 
                   Atom(name="HB2", element="H", coordinates=(-7.490, 38.778, 73.263)), Atom(name="CG", element="C", coordinates=(-7.552, 38.226, 75.236)),                                                                      
                   Atom(name="ND1", element="N", coordinates=(-6.396, 37.476, 75.292)), Atom(name="HD1", element="H", coordinates=(-5.901, 37.120, 74.500)),                                      
                   Atom(name="CE1", element="C", coordinates=(-6.046, 37.303, 76.557)), Atom(name="HE1", element="H", coordinates=(-5.247, 36.804, 76.877)),                                      
                   Atom(name="CD2", element="C", coordinates=(-7.887, 38.498, 76.522)), Atom(name="NE2", element="N", coordinates=(-6.940, 37.905, 77.323)), 
                   Atom(name="HD2", element="H", coordinates=(-8.672, 39.027, 76.829))],      
                

                #                  __HE1 
                #    N   HB1    ND1--CE1    
                #    |   |     /     ||
                # HA-CA--CB--CG      ||
                #    |   |     \\    || 
                #  O=C   HB2    CD2--NE2
                #               HD2  HE2              neutral His, proton on NE2                                                                                           
            "HSE":[Atom(name="N", element="N", coordinates=(-17.701, 16.227, 66.437)), Atom(name="HN", element="H", coordinates=(-18.240, 15.391, 66.181)), 
                   Atom(name="HA", element="H", coordinates=(-16.404, 17.245, 67.626)), Atom(name="CA", element="C", coordinates=(-16.861, 16.330, 67.643)), 
                   Atom(name="C", element="C", coordinates=(-15.781, 15.262, 67.779)), Atom(name="O", element="O", coordinates=(-14.704, 15.583, 68.307)), #backbone                                                                                    
                   Atom(name="CB", element="C", coordinates=(-17.712, 16.270, 68.949)), Atom(name="HB1", element="H", coordinates=(-17.082, 16.177, 69.727)), 
                   Atom(name="HB2", element="H", coordinates=(-18.292, 15.451, 68.906)), Atom(name="CG", element="C", coordinates=(-18.575, 17.465, 69.165)),                                                                      
                   Atom(name="ND1", element="N", coordinates=(-19.569, 17.852, 68.264)), Atom(name="HE2", element="H", coordinates=(-19.763, 20.117, 70.406)),                                      
                   Atom(name="CE1", element="C", coordinates=(-20.146, 18.952, 68.736)), Atom(name="HE1", element="H", coordinates=(-20.904, 19.438, 68.307)),                                      
                   Atom(name="CD2", element="C", coordinates=(-18.598, 18.361, 70.185)), Atom(name="NE2", element="N", coordinates=(-19.550, 19.298, 69.867)), 
                   Atom(name="HD2", element="H", coordinates=( -18.034, 18.342, 71.002))],      
                  
                                                                                                                      #               HD1__HE1 
            "HSP":[N, HN, HA, CA, C, O, #backbone     Protonated His                                                  #    N   HB1    ND1--CE1                                              
                   CB, Atom(name="HB1", element="H"), Atom(name="HB2", element="H"),                                  #    |   |     /     ||
                   Atom(name="CG", element="C"),                                                                      # HA-CA--CB--CG      ||
                   Atom(name="ND1", element="N"), Atom(name="HE2", element="H"), Atom(name="HD1", element="H"),       #    |   |     \\    || 
                   Atom(name="CE1", element="C"), Atom(name="HE1", element="H"),                                      #  O=C   HB2    CD2--NE2(+)
                   Atom(name="CD2", element="C"), Atom(name="NE2", element="N"), Atom(name="HD2", element="H")],      #               HD2  HE2
                    
            "HIS":[],


                #             / HG21  
                #    N     CG2--HG22
                #    |    /   \ HG23
                # HA-CA--CB-HB
                #    |    \       / HD1 
                #  O=C     CG1--CD--HD2
                #      HG11/ \HG12\ HD3                                                                                                                      
            "ILE": [Atom(name="N", element="N", coordinates=(-11.598, 24.288, 74.051)), Atom(name="HN", element="H", coordinates=(-12.546, 23.907, 74.057)), 
                    Atom(name="HA", element="H", coordinates=(-10.250, 25.541, 74.830)),Atom(name="CA", element="C", coordinates=(-10.949, 24.919, 75.210)), 
                    Atom(name="C", element="C", coordinates=(-10.218, 23.917, 76.087)), Atom(name="O", element="O", coordinates=(-9.620, 24.311, 77.105)), #backbone                                                                    
                    Atom(name="CB", element="C", coordinates=(-11.920, 25.740, 76.082)), Atom(name="CG2", element="C", coordinates=(-12.674, 26.775, 75.236)), 
                    Atom(name="HB", element="H", coordinates=(-11.379, 26.264, 76.740)), Atom(name="HG21", element="H", coordinates=(-13.217, 27.356, 75.837 )), 
                    Atom(name="HG22", element="H", coordinates=(-12.016, 27.330, 74.735)), Atom(name="HG23", element="H", coordinates=(-13.272, 26.301, 74.595)),   
                    Atom(name="CG1", element="C", coordinates=(-12.893, 24.844, 76.886)), Atom(name="HG11", element="H", coordinates=(-12.352, 24.273, 77.492)), 
                    Atom(name="HG12", element="H", coordinates=(-13.390, 24.274, 76.245)), Atom(name="CD", element="C", coordinates=(-13.883, 25.645, 77.715)), 
                    Atom(name="HD1", element="H", coordinates=(-14.427, 26.222, 77.109)), Atom(name="HD2", element="H", coordinates=(-14.481, 25.017, 78.209)),       
                    Atom(name="HD3", element="H", coordinates=(-13.382, 26.214, 78.364))],                                                                   


                #                  / HD12 
                #    N   HB1    CD1--HD12
                #    |   |     /   \ HD13
                # HA-CA--CB--CG-HG
                #        |     \   / HD21
                #        HB2    CD2--HD22
                #                  \ HD23  
            "LEU": [Atom(name="N", element="N", coordinates=(-9.535, 13.777, 80.440)), Atom(name="HN", element="H", coordinates=(-8.699, 14.256, 80.828)), 
                    Atom(name="HA", element="H", coordinates=(-10.451, 11.980, 80.402)), Atom(name="CA", element="C", coordinates=(-10.203, 12.664, 81.098)), 
                    Atom(name="C", element="C", coordinates=(-11.477, 13.128, 81.790)), Atom(name="O", element="O", coordinates=( -12.561, 12.580, 81.554)), #backbone                                                                    
                    Atom(name="CB", element="C", coordinates=(-9.250, 12.006, 82.098)), Atom(name="HB1", element="H", coordinates=(-9.265, 12.558, 82.941)),
                    Atom(name="HB2", element="H", coordinates=(-8.320, 12.052, 81.708)), Atom(name="CG", element="C", coordinates=(-9.492, 10.552, 82.508 )), 
                    Atom(name="HG", element="H", coordinates=(-10.071, 10.114, 81.824)), Atom(name="CD1", element="C", coordinates=(-8.168, 9.810, 82.601)),        
                    Atom(name="HD11", element="H", coordinates=(-7.584, 10.251, 83.284)), Atom(name="HD12", element="H", coordinates=(-8.334, 8.861, 82.869)), 
                    Atom(name="HD13", element="H", coordinates=(-7.711, 9.828, 81.711)), Atom(name="CD2", element="C", coordinates=(-10.229, 10.489, 83.841)),
                    Atom(name="HD21", element="H", coordinates=(-11.111, 10.954, 83.757)), Atom(name="HD22", element="H", coordinates=(-10.379, 9.533, 84.094)), 
                    Atom(name="HD23", element="H", coordinates=(-9.681, 10.939, 84.547))],                                  


                # HN-N   HB1 HG1 HD1 HE1    HZ1
                #    |   |   |   |   |     /
                # HA-CA--CB--CG--CD--CE--NZ--HZ2
                #    |   |   |   |   |     \
                #  O=C   HB2 HG2 HD2 HE2     HZ3
                #
            "LYS": [Atom(name="N", element="N", coordinates=(-11.402, 31.649, 80.666)),Atom(name="HN", element="H", coordinates=(-12.329, 31.218, 80.930)), 
                    Atom(name="HA", element="H", coordinates=(-10.163, 32.206, 79.167)), Atom(name="CA", element="C", coordinates=(-11.155, 32.219, 79.348)), 
                    Atom(name="C", element="C", coordinates=(-11.645, 33.660, 79.273)), Atom(name="O", element="O", coordinates=(-12.557, 34.067, 79.999 )), #backbone                                                                    #    
                    Atom(name="CB", element="C", coordinates=(-11.851, 31.408, 78.256 )), Atom(name="HB1", element="H", coordinates=(-11.681, 31.850, 77.382)), 
                    Atom(name="HB2", element="H", coordinates=(-12.825, 31.389, 78.452)), Atom(name="CG", element="C", coordinates=(-11.351, 29.973, 78.173)), 
                    Atom(name="HG1", element="H", coordinates=(-11.891, 29.471, 77.515)), Atom(name="HG2", element="H", coordinates=(-11.393, 29.555, 79.068)),       
                    Atom(name="CD", element="C", coordinates=(-9.881, 29.971, 77.702)), Atom(name="HD1", element="H", coordinates=(-9.320, 30.356, 78.434)), 
                    Atom(name="HD2", element="H", coordinates=( -9.816, 30.549, 76.890 )), Atom(name="CE", element="C", coordinates=(-9.373, 28.587, 77.366)), 
                    Atom(name="HE1", element="H", coordinates=(-9.626, 28.369, 76.418)), Atom(name="HE2", element="H", coordinates=(-9.813, 27.922, 77.977)),       
                    Atom(name="NZ", element="N", coordinates=(-7.911, 28.499, 77.516)), Atom(name="HZ1", element="H", coordinates=(-7.678, 27.743, 78.131)), 
                    Atom(name="HZ2", element="H", coordinates=(-7.487, 28.338, 76.622)), Atom(name="HZ3", element="H", coordinates=(-7.557, 29.355, 77.898))],                                                                   


                # HN-N   HB1 HG1     HE1
                #    |   |   |       |
                # HA-CA--CB--CG--SD--CE--HE3
                #    |   |   |       |
                #  O=C   HB2 HG2     HE2
            "MET": [Atom(name="N", element="N", coordinates=(-21.453, 29.117, 76.768)), Atom(name="HN", element="H", coordinates=(-21.714, 29.535, 77.682)),
                    Atom(name="HA", element="H", coordinates=(-21.146, 27.525, 75.583)), Atom(name="CA", element="C", coordinates=(-21.596, 27.704, 76.462)), 
                    Atom(name="C", element="C", coordinates=(-23.062, 27.357, 76.322)), Atom(name="O", element="O", coordinates=(-23.887, 27.799, 77.135)), #backbone                                                                    #
                    Atom(name="CB", element="C", coordinates=(-20.955, 26.895, 77.575)), Atom(name="HB1", element="H", coordinates=(-21.569, 26.918, 78.378)), 
                    Atom(name="HB2", element="H", coordinates=(-20.089, 27.346, 77.837)), Atom(name="CG", element="C", coordinates=(-20.696, 25.510, 77.188)), 
                    Atom(name="HG1", element="H", coordinates=(-20.341, 25.482, 76.224)), Atom(name="HG2", element="H", coordinates=(-21.581, 24.988, 77.165)),       
                    Atom(name="SD", element="S", coordinates=(-19.557, 24.766, 78.298)), Atom(name="CE", element="C", coordinates=(-19.952, 25.474, 79.911)), 
                    Atom(name="HE1", element="H", coordinates=(-20.857, 25.161, 80.195)), Atom(name="HE2", element="H", coordinates=(-19.269, 25.177, 80.576)),       
                    Atom(name="HE3", element="H", coordinates=(-19.946, 26.470, 79.843))],                                                                   


                #              HD1  HE1
                # HN-N   HB1   CD1--CE1
                #    |   |    //      \\
                # HA-CA--CB--CG        CZ--HZ
                #    |   |    \   __   /  
                #  O=C   HB2   CD2--CE2
                #              HD2  HE2
            "PHE": [Atom(name="N", element="N", coordinates=(-32.076, 25.289, 63.000)), Atom(name="HN", element="H", coordinates=(-32.721, 25.959, 63.455)), 
                    Atom(name="HA", element="H", coordinates=(-30.704, 23.836, 63.092)), Atom(name="CA", element="C", coordinates=(-31.380, 24.213, 63.717)), 
                    Atom(name="C", element="C", coordinates=(-32.390, 23.091, 64.002)), Atom(name="O", element="O", coordinates=(-33.331, 23.309, 64.792)), #backbone                                                                    
                    Atom(name="CB", element="C", coordinates=(-30.811, 24.715, 65.044)), Atom(name="HB1", element="H", coordinates=(-30.498, 23.928, 65.565)), 
                    Atom(name="HB2", element="H", coordinates=(-31.535, 25.190, 65.532)), Atom(name="CG", element="C", coordinates=(-29.630, 25.683, 64.935)), 
                    Atom(name="CD1", element="C", coordinates=(-28.956, 25.891, 63.735 )), Atom(name="HD1", element="H", coordinates=(-29.239, 25.430, 62.910)),       
                    Atom(name="CD2", element="C", coordinates=(-29.157, 26.302, 66.065)), Atom(name="HD2", element="H", coordinates=(-29.610, 26.132, 66.953 )),                                     
                    Atom(name="CE1", element="C", coordinates=(-27.854, 26.772, 63.688)), Atom(name="HE1", element="H", coordinates=(-27.409, 26.942, 62.799)),                                     
                    Atom(name="CE2", element="C", coordinates=(-28.064, 27.167, 66.034)), Atom(name="HE2", element="H", coordinates=(-27.767, 27.626, 66.885)),                                     
                    Atom(name="CZ", element="C", coordinates=(-27.392, 27.389, 64.841 )), Atom(name="HZ", element="H", coordinates=(-26.594, 27.979, 64.809))],                                      
                                                                                                                      

                #    N---CD -HD1/HD2   
                #    |     \
                #    |      CG -HG1
                #    |     /    
                # HA-CA--CB
                #    |   |  \    
                #  O=C   HB1 HB2                                                                                                                                                                                                                                     
            "PRO": [Atom(name="N", element="N", coordinates=(-17.312, 23.341, 68.608 )), Atom(name="HA", element="HA", coordinates=(-18.229, 22.071, 69.945)), 
                    Atom(name="CA", element="C", coordinates=(-17.292, 22.379, 69.734)), Atom(name="C", element="C", coordinates=(-16.620, 22.970, 70.961)), 
                    Atom(name="O", element="O", coordinates=(-15.551, 23.596, 70.857)), #backbone                                                                        
                    Atom(name="CB", element="C", coordinates=(-16.457, 21.218, 69.177)), Atom(name="HB1", element="H", coordinates=(-16.772, 20.343, 69.521)), 
                    Atom(name="HB2", element="H", coordinates=(-15.490, 21.341, 69.365)), Atom(name="CG", element="C", coordinates=(-16.698, 21.287, 67.639)), 
                    Atom(name="HG1", element="H", coordinates=(-17.583, 20.926, 67.434)), Atom(name="HG2", element="H", coordinates=(-15.974, 20.821, 67.176)),       
                    Atom(name="CD", element="C", coordinates=(-16.636, 22.789, 67.410)), Atom(name="HD1", element="H", coordinates=(-17.126, 23.048, 66.581)), 
                    Atom(name="HD2", element="H", coordinates=( -15.693, 23.108, 67.362))],     
                                                                                                                      

                # HN-N   HB1
                #    |   | 
                # HA-CA--CB--OG     
                #    |   |     \
                #  O=C   HB2    HG1           
            "SER": [Atom(name="N", element="N", coordinates=(-4.446, 19.119, 76.589)), Atom(name="HN", element="H", coordinates=(-3.842, 19.844, 76.173)), 
                    Atom(name="HA", element="H", coordinates=(-4.818, 17.191, 76.134)), Atom(name="CA", element="C", coordinates=(-4.105, 17.702, 76.639)), 
                    Atom(name="C", element="C", coordinates=(-4.067, 17.164, 78.064)), Atom(name="O", element="O", coordinates=(-4.243, 15.957, 78.264)), #backbone                                                                    
                    Atom(name="CB", element="C", coordinates=(-2.757, 17.475, 75.958)), Atom(name="HB1", element="H", coordinates=(-2.339, 16.659, 76.337)),                                                                
                    Atom(name="HB2", element="H", coordinates=(-2.181, 18.268, 76.112)), Atom(name="OG", element="O", coordinates=(-2.939, 17.297, 74.565)),                                                                     
                    Atom(name="HG1", element="H", coordinates=(-3.612, 17.955, 74.230))],                                                                   


                # HN-N     OG1--HG1         
                #    |    /
                # HA-CA--CB-HB
                #    |    \
                #  O=C     CG2--(HG21, HG22, HG21)
            "THR": [Atom(name="N", element="N", coordinates=(-19.234, 41.389, 63.915)), Atom(name="HN", element="H", coordinates=(-18.996, 40.461, 64.356)), 
                    Atom(name="HA", element="H", coordinates=(-20.120, 42.391, 62.411)), Atom(name="CA", element="C", coordinates=(-19.673, 41.508, 62.531)), 
                    Atom(name="C", element="C", coordinates=(-20.691, 40.411, 62.256)), Atom(name="O", element="O", coordinates=(-20.431, 39.234, 62.536)), #backbone                                                                    
                    Atom(name="CB", element="C", coordinates=(-18.495, 41.355, 61.559)), Atom(name="HB", element="H", coordinates=(-18.101, 40.421, 61.660)),                                                                 
                    Atom(name="OG1", element="O", coordinates=(-17.475, 42.308, 61.900)), Atom(name="HG1", element="H", coordinates=(-17.213, 42.206, 62.865 )),                                     
                    Atom(name="CG2", element="C", coordinates=(-18.938, 41.552, 60.087)), Atom(name="HG21", element="H", coordinates=(-19.261, 42.489, 59.966)),                                    
                    Atom(name="HG23", element="H", coordinates=(-19.674, 40.910, 59.879)), Atom(name="HG22", element="H", coordinates=(-18.159, 41.381, 59.487))],                                  


                #                       HE3
                # HN-N   HB1            CE3
                #    |   |             /   \\
                # HA-CA--CB---CG-----CD2    CZ3-HZ3
                #    |   |    ||     ||      |
                #  O=C   HB2  CD1    CE2    CH2-HH2
                #            /   \   /  \   /
                #         HD1     NE1    CZ2
                #                 HE1    HZ2
            "TRP": [Atom(name="N", element="N", coordinates=(-20.426, 30.342, 61.290)), Atom(name="HN", element="H", coordinates=(-20.100, 29.963, 60.362)), 
                    Atom(name="HA", element="H", coordinates=(-20.750, 31.919, 62.489)), Atom(name="CA", element="C", coordinates=(-20.696, 31.762, 61.500)), 
                    Atom(name="C", element="C", coordinates=(-22.040, 32.130, 60.899)), Atom(name="O", element="O", coordinates=(-22.208, 32.099, 59.673)), #backbone                                                                    
                    Atom(name="CB", element="C", coordinates=(-19.587, 32.587, 60.856)), Atom(name="HB1", element="H", coordinates=(-19.563, 32.383, 59.860)), 
                    Atom(name="HB2", element="H", coordinates=(-18.693, 32.290, 61.240)), Atom(name="CG", element="C", coordinates=(-19.619, 34.096, 60.978)), 
                    Atom(name="CD1", element="C", coordinates=(-20.224, 34.849, 61.965)), Atom(name="HD1", element="H", coordinates=(-20.728, 34.472, 62.736)),       
                    Atom(name="NE1", element="N", coordinates=(-20.028, 36.191, 61.726)), Atom(name="HE1", element="H", coordinates=(-20.383, 36.939, 62.280)),                                     
                    Atom(name="CE2", element="C", coordinates=(-19.259, 36.331, 60.594)), Atom(name="CZ2", element="C", coordinates=(-18.766, 37.499, 59.975)), 
                    Atom(name="HZ2", element="H", coordinates=(-18.967, 38.394, 60.369)), Atom(name="CH2", element="C", coordinates=(-18.027, 37.329, 58.846)), 
                    Atom(name="HH2", element="H", coordinates=(-17.663, 38.119, 58.350)), Atom(name="CZ3", element="C", coordinates=(-17.739, 36.033, 58.323)), 
                    Atom(name="HZ3", element="H", coordinates=(-17.194, 35.962, 57.496)), Atom(name="CE3", element="C", coordinates=(-18.203, 34.895, 58.944)), 
                    Atom(name="HE3", element="H", coordinates=(-17.989, 33.993, 58.577)), Atom(name="CD2", element="C", coordinates=(-18.979, 35.033, 60.100))],                                    


                # HN-N   HB1  CD1--CE1 
                #    |   |   //      \\ 
                # HA-CA--CB--CG      CZ--OH 
                #    |   |    \  __  /     \
                #  O=C   HB2  CD2--CE2      HH
                #             HD2  HE2
            "TYR": [Atom(name="N", element="N", coordinates=(-20.395, 29.269, 66.129)), Atom(name="HN", element="H", coordinates=(-20.533, 29.313, 65.084)),
                    Atom(name="HA", element="H", coordinates=(-20.694, 28.300, 67.868 )), Atom(name="CA", element="C", coordinates=(-20.677, 28.058, 66.904)), 
                    Atom(name="C", element="C", coordinates=(-19.531, 27.082, 66.673)), Atom(name="O", element="O", coordinates=(-19.125, 26.868, 65.521)), #backbone                                                                    #             HD1  HE2
                    Atom(name="CB", element="C", coordinates=(-21.989, 27.393, 66.433)), Atom(name="HB1", element="H", coordinates=(-22.146, 26.590, 66.995)), 
                    Atom(name="HB2", element="H", coordinates=(-21.880, 27.138, 65.480)), Atom(name="CG", element="C", coordinates=(-23.200, 28.302, 66.550)), 
                    Atom(name="CD1", element="C", coordinates=(-23.506, 29.212, 65.562)), Atom(name="HD1", element="H", coordinates=(-22.949, 29.263, 64.731)),       
                    Atom(name="CD2", element="C", coordinates=(-24.063, 28.182, 67.658)), Atom(name="HD2", element="H", coordinates=(-23.878, 27.497, 68.345)),                                     
                    Atom(name="CE1", element="C", coordinates=(-24.596, 30.083, 65.694)), Atom(name="HE1", element="H", coordinates=(-24.759, 30.786, 65.012)),                                     
                    Atom(name="CE2", element="C", coordinates=(-25.155, 29.019, 67.771)), Atom(name="HE2", element="H", coordinates=(-25.753, 28.957, 68.578)),                                     
                    Atom(name="CZ", element="C", coordinates=(-25.428, 29.938, 66.784)), Atom(name="OH", element="O", coordinates=(-26.516, 30.764, 66.922)), 
                    Atom(name="HH", element="H", coordinates=(-26.212, 31.715, 66.899))],        


                #            / HG11
                # HN-N     CG1-HG12
                #    |    /  \ HG13
                # HA-CA--CB-HB
                #    |    \   / HG21
                #  O=C     CG2--HG22 HG23 
            "VAL": [Atom(name="N", element="N", coordinates=(-12.728, 22.683, 72.092)), Atom(name="HN", element="H", coordinates=(-13.699, 22.795, 71.788)), 
                    Atom(name="HA", element="H", coordinates=(-10.961, 23.030, 71.175)), Atom(name="CA", element="C", coordinates=(-11.631, 23.582, 71.695)), 
                    Atom(name="C", element="C", coordinates=(-10.936, 24.184, 72.886)), Atom(name="O", element="O", coordinates=(-9.754, 24.510, 72.775)), #backbone                                                                   
                    Atom(name="CB", element="C", coordinates=(-12.156, 24.661, 70.711 )), Atom(name="HB", element="H", coordinates=(-12.769, 24.199, 70.069)),                                                                
                    Atom(name="CG1", element="C", coordinates=(-12.950, 25.720, 71.418)), Atom(name="HG11", element="H", coordinates=(-13.834, 25.851, 70.961)), 
                    Atom(name="HG12", element="H", coordinates=(-13.115, 25.449, 72.370)), Atom(name="HG13", element="H", coordinates=(-12.449, 26.589, 71.410)),  
                    Atom(name="CG2", element="C", coordinates=(-10.976, 25.348, 69.951)), Atom(name="HG21", element="H", coordinates=(-11.022, 26.325, 70.126)), 
                    Atom(name="HG22", element="H", coordinates=(-10.124, 24.969, 70.294)), Atom(name="HG23", element="H", coordinates=(-11.080, 25.159, 68.981))],  

            "ALAD": [],  
        }

