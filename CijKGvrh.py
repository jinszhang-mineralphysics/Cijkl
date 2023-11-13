# -*- coding: utf-8 -*-
"""
Created on Sun Apr 08 18:51:39 2015
Modified on Thu Jun 04 14:20:16 2015
Modified on Mon Jul 09 23:03:12 2018 
fixed bugs of monoclinic symmetry

Cij formating, include transformaiton from 1D array into 2D Cij and into 4D Cijkl

*********************************************************************************************

#####################
#Function 1: CijVRH#
#####################

Calculate aggregate elastic properties K, G, Vp and Vs from single-crystal elasticity tensor 
Cijs in the form of 1D array and density
K and G are calculated based on the VRH averaging scheme (Hill 1952):
Voight bond: constant strain
Kv = 1.0/9.0*((c11+c22+c33) + 2.0*(c12+c23+c13))
Gv = 1.0/15.0*((c11+c22+c33) - (c12+c23+c13) + 3.0*(c44+c55+c66))
Reuss bond: constant stress
Kr = 1.0/((s11+s22+s33)+2.0*(s12+s13+s23))
Gr = 15.0/(4.0*(s11+s22+s33) - 4.0*(s12+s13+s23) + 3.0*(s44+s55+s66))
Vp = ((Kvrh+4.0/3.0*Gvrh)/dens)**0.5
Vs = (Gvrh/dens)**0.5
Error is calculated based on error propagation, assuming all elastic constants are independant,
K, G and density are also independant (realsitically this is not true!).
https://en.wikipedia.org/wiki/Propagation_of_uncertainty
 
Input: c1d and c1derr both as 1x21 1D array, density and density error

Output: Aggregate elastic properties in the form of 2x9 2D array
	[dens,Kv,Kr,Kvrh,Gv,Gr,Gvrh,Vp,Vs]
	[denserr,Kvrherrb,Kvrherrb,Kvrherr,Gvrherrb,Gvrherrb,Gvrherr,Vperr,Vserr]

*********************************************************************************************

#####################
#Function 2: cijHS#
#####################

Calculate aggregate elastic properties K, G, Vp and Vs from single-crystal elasticity tensor Cijs in 
the form of 1D array and density

K and G are calculated based on the HS averaging scheme:

Input: c1d and c1derr both as 1x21 1D array, density and density error

Output: Aggregate elastic properties in the form of 2x9 2D array
	[dens,Kv,Kr,Kvrh,Gv,Gr,Gvrh,Vp,Vs]
	[denserr,Kvrherrb,Kvrherrb,Kvrherr,Gvrherrb,Gvrherrb,Gvrherr,Vperr,Vserr]

*********************************************************************************************
Redistribution and use in source and binary forms, with or without modification, are permitted 
provided that the following conditions are met:

 * Redistributions of source code must retain the above copyright notice, this list of 
   conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice, this list of 
   conditions and the following disclaimer in the documentation and/or other materials 
   provided with the distribution.
 * Neither the name of the copyright holders nor the names of any contributors may be used 
   to endorse or promote products derived from this software without specific prior written 
   permission.

 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS 
 OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF 
 MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE
 COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, 
 EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) 
 HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN
 IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

Copyright (c) 2016-2023, @author: Jin Zhang, Department of Geology and Geophysics, Texas A&M University
All rights reserved.

@author: ZhangJin
"""
import numpy as np
from CijSij124d import Cij1t2r, Cij2t4r, C1dts2d, C1dts1d

########################################################
######################Function 1: ######################
########################################################


def CijVRH(c1d,c1derr,dens,denserr):
    
    s2r = C1dts2d (c1d)
	
    #cubic: c11=c22=c33, c44=c55=c66,c12=c13=c23, 3 cosntants
    #tetragonal: c11=c22 c23=c13 c44=c55 6 constants
    c11 = c1d[0]#c11
    c22 = c1d[1]#c22
    c33 = c1d[2]#c33
    c44 = c1d[3]#c44
    c55 = c1d[4]#c55
    c66 = c1d[5]#c66
    c12 = c1d[6]#c12
    c13 = c1d[7]#c13
    c23 = c1d[8]#c23
	#Ortho:9constants; 
    c15 = c1d[9]#c15
    c25 = c1d[10]#c25
    c35 = c1d[11]#c35
    c46 = c1d[12]#c46
    #Mono: 13 constants; 	
    c14 = c1d[13]#c14
    c16 = c1d[14]#c16
    c24 = c1d[15]#c24
    c26 = c1d[16]#c26
    c34 = c1d[17]#c34
    c36 = c1d[18]#c36
    c45 = c1d[19]#c45
    c56 = c1d[20]#c56
	
	#Cij error input
    c11err = c1derr[0]
    c22err = c1derr[1]
    c33err = c1derr[2]
    c44err = c1derr[3]
    c55err = c1derr[4]
    c66err = c1derr[5]
    c12err = c1derr[6]
    c13err = c1derr[7]
    c23err = c1derr[8]
    c15err = c1derr[9]
    c25err = c1derr[10]
    c35err = c1derr[11]
    c46err = c1derr[12]
    c14err = c1derr[13]
    c16err = c1derr[14]
    c24err = c1derr[15]
    c26err = c1derr[16]
    c34err = c1derr[17]
    c36err = c1derr[18]
    c45err = c1derr[19]
    c56err = c1derr[20]
	
    Kv = 1.0/9.0*((c11+c22+c33) + 2.0*(c12+c23+c13))
    Gv = 1.0/15.0*((c11+c22+c33) - (c12+c23+c13) + 3.0*(c44+c55+c66))
	
	#cubic 3 constants; tetragonal 6 constants
    s11 = s2r[0,0] #s11
    s22 = s2r[1,1] #s22
    s33 = s2r[2,2] #s33
    s44 = s2r[3,3] #s44
    s55 = s2r[4,4] #s55
    s66 = s2r[5,5] #s66
    s12 = s2r[0,1] #s12
    s13 = s2r[0,2] #s13
    s23 = s2r[1,2] #s23
	#Ortho:9 constants; 
    s15 = s2r[0,4] #s15
    s25 = s2r[1,4] #s25
    s35 = s2r[2,4] #s35
    s46 = s2r[3,5] #s46
    #Mono: 13 constants; 	
    s14 = s2r[0,3] #s14
    s16 = s2r[0,5] #s16
    s24 = s2r[1,3] #s24
    s26 = s2r[1,5] #s26
    s34 = s2r[2,3] #s34
    s36 = s2r[2,5] #s36
    s45 = s2r[3,4] #s45
    s56 = s2r[4,5] #s56
    #Triclinic:21constants;  
	
    Kr = 1.0/((s11+s22+s33)+2.0*(s12+s13+s23))
    Gr = 15.0/(4.0*(s11+s22+s33) - 4.0*(s12+s13+s23) + 3.0*(s44+s55+s66))
	
    Kvrh = (Kr+Kv)/2.0
    Gvrh = (Gr+Gv)/2.0
     
    Kvrherra = np.absolute(Kv-Kr)/2.0 
    Kvrherrb = ((c11err**2+c22err**2+c33err**2+4.0*(c12err**2+c13err**2+c23err**2))/81.0)**0.5
    Kvrherr = Kvrherra + (1.0/4.0*Kvrherrb**2+1.0/4.0*Kvrherrb**2)**0.5
    Gvrherra = np.absolute(Gv-Gr)/2.0 
    Gvrherrb = ((c11err**2+c22err**2+c33err**2+c12err**2+c13err**2+c23err**2+9.0*(c44err**2+c55err**2+c66err**2))/225.0)**0.5
    Gvrherr = Gvrherra + (1.0/4.0*Gvrherrb**2+1.0/4.0*Gvrherrb**2)**0.5
    
    Vp = ((Kvrh+4.0/3.0*Gvrh)/dens)**0.5
    Vperr = Vp/2.0*((Kvrherr**2+16.0/9.0*Gvrherr**2)/(Kvrh+4.0/3.0*Gvrh)**2+denserr**2/dens**2)**0.5
    Vs = (Gvrh/dens)**0.5
    Vserr = Vs/2.0*(Gvrherr**2/Gvrh**2+denserr**2/dens**2)**0.5
	
    Evrh = np.array ([[dens,Kv,Kr,Kvrh,Gv,Gr,Gvrh,Vp,Vs],[denserr,Kvrherrb,Kvrherrb,Kvrherr,Gvrherrb,Gvrherrb,Gvrherr,Vperr,Vserr]])
    
    return Evrh
	
