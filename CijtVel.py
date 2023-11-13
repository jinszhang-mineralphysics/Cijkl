# -*- coding: utf-8 -*-
"""
Created on Sun Apr 08 18:51:39 2015
Modified on Thu Jun 04 14:20:16 2015

Calculate velocities based on known Cij model

*********************************************************************************************

#####################
#Function 1: Cijtvel#
#####################

Calculate velocities based on known Cij model, phonon direction n, and density

Input: Cijs: c1d as 1D array
	   Phonon direction: n as 1D array
	   density: dens as scalar number
Output: Velocities: vel as 1D array in the sequence of Vp-Vsv-Vsh
	    Polarization directions: pol as 3x3D arrays in the sequence of 
								 pol[0]for Vp, pol[1] for Vsv, pol[2] for Vsh 

*********************************************************************************************

####################
#Function 2: DiffVs#  ## Not done yet##
####################

Calculate |Vsh-Vsv|/Vsvrh based on known Cij model and density

Input: c2r as 2D array
Output: c4r as 4D array

*********************************************************************************************

####################
#Function 3: DiffVp# ## Not done yet##
####################

Calculate |Vsh-Vsv|/Vsvrh based on known Cij model and density

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
from numpy import linalg 
import math
from CijSij124d import Cij1t2r, Cij2t4r
#from Rotation import rotation_matrix


########################################################
######################Function 1: ######################
########################################################

def Cijtvel(c1d,n,dens):

	# Converting C1d to C4r
    c2r = Cij1t2r(c1d)
    c4r = Cij2t4r(c2r)
	# Normalize the phonon direction
    n = n/linalg.norm(n)
    # Calculate the Christofol matrix: Cijkl*nj*nl
    A = np.zeros((3,3))
    for j in range (1,4):
        for l in range (1,4):
            for i in range (1,4): 
                for k in range (1,4):
                    A[i-1,k-1] += c4r[i-1,j-1,k-1,l-1]*n[j-1]*n[l-1]
    #print 'Christfol matrix is:', '\n', A
	
    # Solve the eigenvector and eigen value of Christfol matrix A, E is the psudo-elastic-moduli E = Vel**2 * dens
    E, pol = linalg.eig(A)
    # Calculate velocity 
    vel = (E/dens)**0.5
    # Sort the velocities and obtain the sorting index i for the velocities in sequence of Vs1 < Vs2 < Vp, 
    i = np.argsort(vel)
	# Output the velocities and pol in sequence of Vp, Vs1>Vs2
    vel = np.array([vel[i[2]],vel[i[1]],vel[i[0]]])
    pol = np.array([pol[:,i[2]],pol[:,i[1]],pol[:,i[0]]])
	
    # find which shear wave is Vsh, and which is Vsv, re-order the sequence of velocity array and polarization direction array in the sequence of Vp,Vsv,Vsh
    
#	M = rotation_matrix(math.pi/2,pol[:,i[2]]) #clockwise rotation matrix with Vp's polarization direction
#    polt = np.matmul(M, pol[:,i[0]])
#    if linalg.norm(polt+pol[:,i[1]])<=1e-04: # clockwise rotating pol[:,0] is -pol[:,1]
#        vel = np.array([vel[i[2]],vel[i[1]],vel[i[0]]]) #pol[:,0] is Vsh
#        pol = np.array([pol[:,i[2]],pol[:,i[1]],pol[:,i[0]]])
#    else:
#        if np.abs(linalg.norm(polt+pol[:,i[1]])-2.0)<=1e-04: # clockwise rotating pol[:,0] = pol[:,1] = 1
#            vel = np.array([vel[i[2]],vel[i[0]],vel[i[1]]]) #pol[:,0] is Vsv
#            pol = np.array([pol[:,i[2]],pol[:,i[0]],pol[:,i[1]]]) 
#        else:
#            print 'calculation error! The norms of the polarization directions are not 1!'		
    return vel, pol