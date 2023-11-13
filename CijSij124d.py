# -*- coding: utf-8 -*-
"""
Created on Sun Apr 08 18:51:39 2015
Modified on Thu Jun 04 14:20:16 2015
Modified on Mon Jul 09 23:03:12 2018 
fixed bugs of monoclinic symmetry

Cij formating, include transformaiton from 1D array into 2D Cij and into 4D Cijkl

*********************************************************************************************

#####################
#Function 1: cij1t2r#
#####################

Transform single-crystal elasticity tensor in the form of 1D array 
into Cij as rank 2 tensor to rank 4 tensor Cijkl

Input: c1d as 1D array
Output: c2r as 2D array

*********************************************************************************************

#####################
#Function 2: cij2t4r#
#####################

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

########################################################
######################Function 1: ######################
########################################################

def Cij1t2r(c1d):
    
    c2r=np.zeros((6,6))
    #cubic: c11=c22=c33, c44=c55=c66,c12=c13=c23, 3 cosntants
    #tetragonal: c11=c22 c23=c13 c44=c55 6 constants
    c11 = c2r[0,0] = c1d[0]#c11
    c22 = c2r[1,1] = c1d[1]#c22
    c33 = c2r[2,2] = c1d[2]#c33
    c44 = c2r[3,3] = c1d[3]#c44
    c55 = c2r[4,4] = c1d[4]#c55
    c66 = c2r[5,5] = c1d[5]#c66
    c12 = c2r[0,1] = c1d[6]#c12
    c13 = c2r[0,2] = c1d[7]#c13
    c23 = c2r[1,2] = c1d[8]#c23
    #Ortho:9constants; 
    c15 = c2r[0,4] = c1d[9]#c15
    c25 = c2r[1,4] = c1d[10]#c25
    c35 = c2r[2,4] = c1d[11]#c35
    c46 = c2r[3,5] = c1d[12]#c46
    #Mono: 13 constants; 	
    c14 = c2r[0,3] = c1d[13]#c14
    c16 = c2r[0,5] = c1d[14]#c16
    c24 = c2r[1,3] = c1d[15]#c24
    c26 = c2r[1,5] = c1d[16]#c26
    c34 = c2r[2,3] = c1d[17]#c34
    c36 = c2r[2,5] = c1d[18]#c35
    c45 = c2r[3,4] = c1d[19]#c45
    c56 = c2r[4,5] = c1d[20]#c56
    #Triclinic:21constants; 
    for i in range (0,6):
        for j in range (0,6):
            if i > j:
                c2r [i,j]= c2r [j,i]
            j=j+1
        i=i+1
        continue
#    print '\n','Triclinic:21constants;  ', 'Mono: 13 constants;  ',  'Ortho:9constants;', '\n', 'tetragonal: c11=c22 c23=c13 c44=c55 6 constants','\n','cubic: c11=c22=c33, c44=c55=c66,c12=c13=c23, 3 cosntants','\n'
#    print 'c11=',c11,'\n','c22=',c22,'\n','c33=',c33,'\n','c44=',c44,'\n','c55=',c55,'\n','c66=',c66,'\n','c12=',c12,'\n','c13=',c13,'\n','c23=',c23,'\n','c16=',c16,'\n','c26=',c26,'\n','c36=',c36,'\n','c45=',c45,'\n','c14=',c14,'\n','c15=',c15,'\n','c24=',c24,'\n','c25=',c25,'\n','c34=',c34,'\n','c35=',c35,'\n','c46=',c46,'\n','c56=',c56,'\n'

#    print c2r

    return c2r

########################################################
######################Function 2: ######################
########################################################

def Cij2t4r(c2r):
    
    c4r=np.zeros((3,3,3,3)) 
    r=1
    s=1
    i=1
    j=1
    k=1
    l=1
    for r in range (1,7):
        if r==1:
            i=1
            j=1
        if r==2:
            i=2
            j=2
        if r==3:
            i=3 
            j=3
        if r==4:
            i=2
            j=3
        if r==5:
            i=1 
            j=3
        if r==6:
            i=1
            j=2
        for s in range (1,7):
            if s==1:
                k=1
                l=1
            if s==2:
                k=2
                l=2
            if s==3:
                k=3
                l=3
            if s==4:
                k=2
                l=3
            if s==5:
                k=1 
                l=3
            if s==6:
                k=1
                l=2
            c4r[i-1,j-1,k-1,l-1]=c2r[r-1,s-1]
            s=s+1
        r=r+1
    for i in range (1,4):
        for j in range (1,4):
            for k in range (1,4):
                for l in range (1,4):
                    if i>j:
                        if k>l:
                            c4r[i-1,j-1,k-1,l-1]=c4r[j-1,i-1,l-1,k-1]
                        else:
                            c4r[i-1,j-1,k-1,l-1]=c4r[j-1,i-1,k-1,l-1]
                    else:
                        if k>l:
                            c4r[i-1,j-1,k-1,l-1]=c4r[i-1,j-1,l-1,k-1]
                        else:
                            c4r[i-1,j-1,k-1,l-1]=c4r[i-1,j-1,k-1,l-1]
        
#    print '\n',c4r
    return c4r
	
########################################################
######################Function 3: ######################
########################################################

def C1dts2d(c1d):
	
    c2r = Cij1t2r(c1d)
    s2r = np.linalg.inv(c2r)  
    return s2r
	
########################################################
######################Function 4: ######################
########################################################
	
def C1dts1d(c1d):
	
    c2r = Cij1t2r(c1d)
    s2r = np.linalg.inv(c2r)
    s1d = np.zeros(21)
    s1d[0] = s2r[0,0] #s11
    s1d[1] = s2r[1,1] #s22
    s1d[2] = s2r[2,2] #s33
    s1d[3] = s2r[3,3] #s44
    s1d[4] = s2r[4,4] #s55
    s1d[5] = s2r[5,5] #s66
    s1d[6] = s2r[0,1] #s12
    s1d[7] = s2r[0,2] #s13
    s1d[8] = s2r[1,2] #s23
    #Ortho:9 constants; 
    s1d[9] = s2r[0,4] #s15
    s1d[10] = s2r[1,4] #s25
    s1d[11] = s2r[2,4] #s35
    s1d[12] = s2r[3,5] #s46
    #Mono: 13 constants; 	
    s1d[13] = s2r[0,3] #s14
    s1d[14] = s2r[0,5] #s16
    s1d[15] = s2r[1,3] #s24
    s1d[16] = s2r[1,5] #s26
    s1d[17] = s2r[2,3] #s34
    s1d[18] = s2r[2,5] #s36
    s1d[19] = s2r[3,4] #s45
    s1d[20] = s2r[4,5] #s56
    #Triclinic:21constants;    
	
    return s1d