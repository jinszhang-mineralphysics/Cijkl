# -*- coding: utf-8 -*-
"""

*********************************************************************************************
Created on Sun Apr 08 18:51:39 2015
Modified on Mon Jul 09 23:03:12 2018 
fixed bugs of monoclinic symmetry

#######################
#Function 1: UnivAniso#
#######################

Calculate aggregate elastic properties K, G, Vp and Vs from single-crystal elasticity tensor 
Cijs in the form of 1D array and density
Kv Kr and Gv Gr are calculated based on the VRH averaging scheme (Hill 1953).
Universal anisotropy index is calculated based on Ranganathan and Ostoja-Starzewski (2008):
AU = 5.0*Gv/Gr+Kv/Kr-6.0

Error is calculated based on error propagation, assuming Kv Kr Gv Gr are independant from 
each other (realsitically this is not true!).
https://en.wikipedia.org/wiki/Propagation_of_uncertainty
 
Input: c1d and c1derr both as 1x21 1D array, density and density error
Output: Universal anisotropy index and its uncertainty

*********************************************************************************************

#####################
#Function 2: VpAniso#
#####################

Transform single-crystal elasticity tensor in the form of Cij as rank 2 tensor to rank 4 tensor Cijkl

Input: c2r as 2D array
Output: c4r as 4D array

*********************************************************************************************

#####################
#Function 3: VsAniso#
#####################

Transform single-crystal elasticity tensor in the form of Cij as rank 2 tensor to rank 4 tensor Cijkl

Input: c2r as 2D array
Output: c4r as 4D array

*********************************************************************************************

########################
#Function 4: PlaneAniso#
########################

Transform single-crystal elasticity tensor in the form of Cij as rank 2 tensor to rank 4 tensor Cijkl

Input: c2r as 2D array
Output: c4r as 4D array

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

"""
import numpy as np
from CijKGvrh import CijVRH
from CijtVel import Cijtvel

########################################################
######################Function 1: ######################
########################################################


def UnivAniso(c1d,c1derr,dens,denserr):
    
    Evrh = CijVRH (c1d,c1derr,dens,denserr)
	#Evrh format: 2X9 array
	#[dens,Kv,Kr,Kvrh,Gv,Gr,Gvrh,Vp,Vs]
	#[denserr,Kvrherrb,Kvrherrb,Kvrherr,Gvrherrb,Gvrherrb,Gvrherr,Vperr,Vserr]
	
    Kv = Evrh[0,1]
    Kverr = Evrh[1,1]
    Gv = Evrh[0,4]
    Gverr = Evrh[1,4]
    Kr = Evrh[0,2]
    Krerr = Evrh[1,2]
    Gr = Evrh[0,5]
    Grerr = Evrh[1,5]	
	
    AU = 5.0*Gv/Gr+Kv/Kr-6.0
    AUerr = (25.0*((Gv/Gr)**2*(Gverr**2/Gv**2+Grerr**2/Gr**2)) + (Kv/Kr)**2*(Kverr**2/Kv**2+Krerr**2/Kr**2))**0.5
    UnivAniso = np.array([AU,AUerr])
	
	
    return UnivAniso
	
	
def AVpVs(c1d,dens,m):
    
	#m is the grid number in 1D. e.g. if n=30, then a total of 27000 sets of HKL will be calculated.
	
    # Assume no error for density and Cij
    denserr = 0.0 
    c1derr= np.zeros(21)
    Evrh = CijVRH (c1d,c1derr,dens,denserr)
	
	#Evrh format: 2X9 array
	#[dens,Kv,Kr,Kvrh,Gv,Gr,Gvrh,Vp,Vs]
	#[denserr,Kvrherrb,Kvrherrb,Kvrherr,Gvrherrb,Gvrherrb,Gvrherr,Vperr,Vserr]
	
    Vpvrh = Evrh[0,7]
    Vsvrh = Evrh[0,8]		
	# Format a 20x20x20 cartitian gird in the 3D HKL space
    vel = np.zeros((1,3))
    pol = np.zeros((3,3))
    pol1d = np.zeros((1,9))

    for nx in range (0, m+1):
        for ny in range (0, m+1):
            for nz in range (1, m+1):
				#For minerals with symmetry equal or higher than Orthorhombic: 
				#there is no difference between the velocities along [nx,ny,nz],[nx,-ny,nz],[nx,ny,-nz] and [nx,-ny,-nz]
				#For minerals with monoclinic symmetry:
				#the velocities along [nx,ny,nz]=[nx,-ny,nz] not equal to [nx,ny,-nz] =[nx,-ny,-nz]
				#Thus, for the monoclinic symmetry, we need to consider both case
				#positive
				npos = np.array([nx,ny,nz])/1.0/m 
				veltpos, poltpos = Cijtvel(c1d, npos, dens)
				pol1dtpos = np.reshape(poltpos,9)
				vel = np.append(vel,[veltpos],axis=0)
				pol1d = np.append(pol1d,[pol1dtpos],axis=0)  
				#negative 
				nneg = np.array([nx,ny,-nz])/1.0/m            
				veltneg, poltneg = Cijtvel(c1d, nneg, dens)
				pol1dtneg = np.reshape(poltneg,9)            
				vel = np.append(vel,[veltneg],axis=0)
				pol1d = np.append(pol1d,[pol1dtneg],axis=0) 
				#print vel
				#print pol1d
    VelPol = np.transpose(np.concatenate((vel.T,pol1d.T),axis=0))
    VelPol = np.delete(VelPol, (0), axis=0)
	
	
	#######"Calculating Vp azimuthal anisotropy"#######

    Vp = VelPol[:,0]
    Vpmax = Vp.max() #the calculated velocities from cij follow the sequence Vp>Vs1>Vs2
    Vpmin = Vp.min()

    VpmaxIndex = np.where(Vp == Vpmax)[0]
    VpmaxPol = np.array([VelPol[VpmaxIndex[0],3], VelPol[VpmaxIndex[0],4], VelPol[VpmaxIndex[0],5]])

    VpminIndex = np.where(Vp == Vpmin)[0]
    VpminPol = np.array([VelPol[VpminIndex[0],3], VelPol[VpminIndex[0],4], VelPol[VpminIndex[0],5]])

	#print '\n', 'Max Vp azimuthal anisotropy is:'
	#print 'Vpmax - Index numbers - Phonon/Polarization direction:', Vpmax, VpmaxIndex, VpmaxPol 
	#print 'Vpmin - Index numbers - Phonon/Polarization direction:', Vpmin, VpminIndex, VpminPol
    AVp = (Vpmax-Vpmin)/Vpvrh
	#print 'Vpmax-Vpmin is',Vpmax-Vpmin 
	#print 'Vp azimuthal anisotropy (Vpmax-Vpmin)/Vpvrh is', AVp
    AVp =([AVp, Vpmax, VpmaxPol, Vpmin, VpminPol])

	#######"Calculating Vs anisotropy"#######

    Vs1 = VelPol[:,1]
    Vs2 = VelPol[:,2]

	###"1. azimuthal anisotropy"###

    Vsmax = Vs1.max()
    Vsmin = Vs2.min()

    VsmaxIndex = np.where(Vs1 == Vsmax)[0]
    VsmaxPhon = np.array([VelPol[VsmaxIndex[0],3], VelPol[VsmaxIndex[0],4], VelPol[VsmaxIndex[0],5]])
    VsmaxPol = np.array([VelPol[VsmaxIndex[0],6], VelPol[VsmaxIndex[0],7], VelPol[VsmaxIndex[0],8]])

    VsminIndex = np.where(Vs2 == Vsmin)[0]
    VsminPhon = np.array([VelPol[VsminIndex[0],3], VelPol[VsminIndex[0],4], VelPol[VsminIndex[0],5]])
    VsminPol = np.array([VelPol[VsminIndex[0],9], VelPol[VsminIndex[0],10], VelPol[VsminIndex[0],11]])

	#print '\n', 'Max Vs azimuthal anisotropy is:'
	#print 'Vsmax - Index numbers - Phonon - Polarization direction:', '\n', Vsmax, VsmaxIndex, VsmaxPhon, VsmaxPol 
	#print 'Vsmin - Index numbers - Phonon - Polarization direction:', '\n', Vsmin, VsminIndex, VsminPhon, VsminPol
    AVs = (Vsmax-Vsmin)/Vsvrh
	#print 'Vsmax-Vsmin is',Vsmax-Vsmin 
	#print 'Vs azimuthal anisotropy (Vsmax-Vsmin)/Vsvrh is', AVs
    AVs = ([AVs, Vsmax, VsmaxPhon, VsmaxPol, Vsmin, VsminPhon, VsminPol]) 

	###"2.Vs splitting anisotropy"###

    DelVs = Vs1-Vs2 #the calculated velocities from cij follow the sequence Vp>Vs1>Vs2

    DelVsmax = DelVs.max()

    DelVsmaxIndex = np.where(DelVs == DelVsmax)[0]
    DelVsmaxPhon = np.array([VelPol[DelVsmaxIndex[0],3], VelPol[DelVsmaxIndex[0],4], VelPol[DelVsmaxIndex[0],5]])
    DelVsmaxPol1 = np.array([VelPol[DelVsmaxIndex[0],6], VelPol[DelVsmaxIndex[0],7], VelPol[DelVsmaxIndex[0],8]])
    DelVsmaxPol2 = np.array([VelPol[DelVsmaxIndex[0],9], VelPol[DelVsmaxIndex[0],10], VelPol[DelVsmaxIndex[0],11]])

	#print '\n', 'Max Vs splitting anisotropy is along phonon direction:', DelVsmaxPhon
	#print 'Index numbers in the stored results are:', DelVsmaxIndex
	#print 'Vs1 with polarization direction:', VelPol[DelVsmaxIndex[0],1], DelVsmaxPol1
	#print 'Vs2 with polarization direction:',  VelPol[DelVsmaxIndex[0],2],DelVsmaxPol2
    DVs = DelVsmax/Vsvrh
	#print 'Delta_Vsmax is',DelVsmax
	#print 'Vs splitting anisotropy (Delta_Vsmax)/Vsvrh is', DVs
    DVs = ([DVs, VelPol[DelVsmaxIndex[0],1], DelVsmaxPhon, DelVsmaxPol1, VelPol[DelVsmaxIndex[0],2], DelVsmaxPhon, DelVsmaxPol2]) 
	# AVpVs[0]: 1x5 array [AVp, Vpmax, VpmaxPol, Vpmin, VpminPol]
	# AVpVs[1]: 1x7 array [AVs, Vsmax, VsmaxPhon, Vsmaxpol, Vsmin, VsminPhon, Vsminpol]
	# AVpVs[2]: 1x7 array [DVs, VelPol[DelVsmaxIndex[0],1], DelVsmaxPhon, DelVsmaxPol1, VelPol[DelVsmaxIndex[0],2], DelVsmaxPhon, DelVsmaxPol2]
	
	
    return AVp, AVs, DVs