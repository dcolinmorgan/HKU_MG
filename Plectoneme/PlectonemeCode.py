#pragma rtGlobals=3		# Use modern global access method and strict wave access.
# import matplotlib.pyplot as plt
import numpy as np
import os,sys
sys.path.insert(1, './run/')
# from scipy.stats import linregress
import glob
from tqdm import tqdm
import pandas as pd
# import seaborn as sns
import scipy
from savitzky_golay import savitzky_golay
# from scipy import stats
# from sklearn import metrics



def PlectonemeCode(Swave, Twave, Wwave, Dwave, CwaveRR, CwaveTT):
    SeqLength=len(Swave)
    rise=0.339

    # make/d/o/n=(SeqLength,4) DNApath, DNApathMajorGroove  #These 2 paths will trace the center and the major groove
    # make/d/o/n=(SeqLength,2, 2) BasepairCovariance, LocalCovariance  #covariance matrices to estimate local stoffness along the sequence
    DNApathMajorGroove=(np.empty([SeqLength,4]))# for i in range(2))
    BasepairCovariance,LocalCovariance=(np.empty([SeqLength,2,2]) for i in range(2))
    # make/d/o/n=(2, 2) BendRot, Covariance, LocalCov  #more matrices for stiffness calculation
    # DNApath=0
    # make/d/o/n=(SeqLength) CurvatureSequence
    # make/d/o/n=(SeqLength)
    CurvatureSequence,Sequence_phase,Sequence_angle_energy,Sequence_angle_exp,EndEffects= (np.empty(SeqLength) for i in range(5))
    #Sequence_angle_exp will eventually store the weight assigned to a plectoneme at each position on the DNA
    # CurvatureSequence=0
    # Sequence_phase=0

    # make/d/o/n=(4,4)
    T_n, Romega_n, Q_n, Rzplus, Rx, Rzminus, Minverse_tot, M_tot= (np.empty([4,4]) for i in range(8))
    # make/d/o/n=4 StartPos, StartPosMG
    StartPos=[0,0,0,1]
    StartPosMG=[1,0,0,1]
    # make/d/o/n=4 OldPos

    DNApath= np.reshape(np.repeat(StartPos,4),[4,4])

    # OldPos=StartPos

    Minverse_tot=np.identity(4) #[[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]]
    M_tot=Minverse_tot
    T_n=M_tot#[[1,0,0,0],[0,1,0,0],[0,0,1,-rise/2],[0,0,0,1]]
    T_n[2,2]=-rise/2

    #This loop finds the 3D path of the relaxed DNA

    Letter1=0
    for ii in tqdm(np.arange(SeqLength)):
        # ii=1
        if Swave[ii]=='A':
            Letter2=0
        elif Swave[ii]=='C':
            Letter2=1
        elif Swave[ii]=='G':
            Letter2=2
        else:
            Letter2=3

        index = 4*Letter1+Letter2			#The index defines the current dinucleotide, AA=0, AC=1,...TT=15

        Sequence_phase[ii]=Sequence_phase[ii-1]+Twave[index]					#This is used to measure how far around the DNA the major groove has rotated relative to the first base pair
        BendRot=[np.cos(Sequence_phase[ii]), np.sin(Sequence_phase[ii])],[-np.sin(Sequence_phase[ii]),np.cos(Sequence_phase[ii])] #Rotation matrix
        Covariance=[CwaveRR[index], 0],[0,CwaveTT[index]]	#The covariance matrix for the current basepair, expressed in the coordinates of the current basepair
        CovRot =  np.multiply(BendRot,Covariance,np.transpose(BendRot)) #Rotating the covariance matrix so it will line up with its neighbors
        BasepairCovariance[ii]=CovRot					#Rotated covariance matrix at position is recorded

        #these next steps are outlined in the paper I sent
        omDiv2_n=Twave[index]/2

        Romega_n=[[np.cos(omDiv2_n),np.sin(omDiv2_n),0,0],[-np.sin(omDiv2_n),np.cos(omDiv2_n),0,0],[0,0,1,0],[0,0,0,1]]
        alpha_n=Wwave[index]

        beta_n=Dwave[index]-np.pi/2

        Rzplus=[[np.cos(beta_n),np.sin(beta_n),0,0],[-np.sin(beta_n),np.cos(beta_n),0,0],[0,0,1,0],[0,0,0,1]]
        Rx=[[1,0,0,0],[0,np.cos(-alpha_n),np.sin(-alpha_n),0],[0,-np.sin(-alpha_n),np.cos(-alpha_n),0],[0,0,0,1]]
        Rzminus=[[np.cos(-beta_n),np.sin(-beta_n),0,0],[-np.sin(-beta_n),np.cos(-beta_n),0,0],[0,0,1,0],[0,0,0,1]]

        Q_n = np.multiply(Rzminus, Rx)* Rzplus
        Minverse_n = np.linalg.inv( T_n* Romega_n* Q_n* Romega_n* T_n)
        Minverse_new =  np.multiply( Minverse_n,Minverse_tot)

        Minverse_tot=Minverse_new #Updating the total tranformation matrix

        CurrentPos = np.multiply(np.transpose(Minverse_tot),StartPos)   #Calculate the coordinates of the current basepair

        CurrentPosMG =  np.multiply(np.transpose(Minverse_tot) ,StartPosMG)

        DNApath=np.append(DNApath,CurrentPos)
        DNApath=np.reshape(DNApath,(int(len(DNApath)/4),4))
        DNApathMajorGroove=np.append(DNApathMajorGroove,CurrentPosMG)
        DNApathMajorGroove=np.reshape(DNApathMajorGroove,(int(len(DNApathMajorGroove)/4),4))
        Letter1=Letter2



    print("path calculated")

    #Make the curvature calculation

    # variable CurveWindow  #must be even
    Sequence_angle_energy=0
    Sequence_angle_exp=0
    CircFrac=0.667	#We assume the plectoneme tip makes a 240� arc before joining the bulk plectoneme region
    BindLength=450   #experimentally, ~450 nt are bound to the surface at each end of the DNA
    AvePlecLength=1000
    # Variable EnergyOffset
    # Variable Z1, Z2, Z3, Z4, Z5, Z6, Z7 	#These are used to assign Boltzmann weights to bending in 8 possible directions
    # Variable Cnorm, E_base

    tanlength=10 #must be even. This is number of basepairs used to calculate the local tangent vectors.
    # Variable VectorMag, CosCurve, SinCurve


    for CurveWindow in tqdm(np.arange(40,120,8)):# (CurveWindow=40; CurveWindow<120; CurveWindow+=8)
        print(CurveWindow)
        a=BasepairCovariance[:,0,0]
        b=BasepairCovariance[:,1,1]
        aa=savitzky_golay(a,CurveWindow+1,3)		#Find covariance matrix over the curvature window
        bb=savitzky_golay(b,CurveWindow+1,3)
        LocalCovariance=BasepairCovariance
        LocalCovariance[:,0,0]=aa
        LocalCovariance[:,1,1]=bb
        TanVector, NormVector, CurveVector=(np.empty([SeqLength,3]) for i in range(3))
        CurveMag, CurvePhase, HalfCurveMag, HalfCurvePhase=(np.empty([SeqLength]) for i in range(4))
        # TanVector=0
        NormVector=DNApathMajorGroove#######-DNApath  #identifies the normal vector alligned with the major groove
        # CurveVector=0
        CurrentTan, CurrentCurve, tP, tM, CurveCross, CurrentNorm=(np.empty(3) for i in range(6))

        # find the tan vectors over tanlength
        half_tan_len=int(tanlength/2)
        for ii in np.arange(half_tan_len,SeqLength-half_tan_len):
            CurrentTan=DNApath[ii+half_tan_len]-DNApath[ii-half_tan_len]
            VectorMag=np.sqrt(np.dot(CurrentTan, CurrentTan))
            TanVector[ii]=CurrentTan/VectorMag  #Normalizes tangent vector to unit length
            print(TanVector)
        half_curve=int(CurveWindow/2)
        # find the curvature vectors and values over curvewindow
        for ii in np.arange(half_curve, SeqLength-half_curve):
            tP=TanVector[ii+half_curve]					#plus tan vector
            tM=TanVector[ii-half_curve]					#minus tan vector
            CurrentCurve =[tP[1]*tM[2]-tP[2]*tM[1],tP[2]*tM[0]-tP[0]*tM[2],tP[0]*tM[1]-tP[1]*tM[0]]   #Cross product
            CurveVector[ii]=CurrentCurve						#curvature vector is recorded at this position
        #	CurveMag[ii]=acos(np.dot(tP, tM))
            CurveMag[ii]=asin(np.sqrt(np.dot(CurrentCurve, CurrentCurve)))
            CurrentCurve=CurveMag[ii]							#normalize the curvature vector to track direction

            #Calculates the phase angle of the curvature relative to major groove at start of DNA
            CurrentTan=TanVector[ii]
            CurrentNorm=NormVector[ii]
            CosCurve=np.dot(CurrentCurve, CurrentNorm)
            tP=CurrentNorm	#Only using these variables as placeholders b/c of short names, otherwise cross product calculation makes a long line of code
            tM=CurrentCurve
            CurveCross=[tP[1]*tM[2]-tP[2]*tM[1],tP[2]*tM[0]-tP[0]*tM[2],tP[0]*tM[1]-tP[1]*tM[0]]
            SinCurve=np.dot(CurveCross, CurrentTan)
            CurvePhase[ii]=mod(atan2(SinCurve,CosCurve)+Sequence_phase[ii],2*np.pi)  #curvature phase is recorded.


        EnergyOffset=25-CurveWindow*0.334*3/4.06   #adds energy penalty from pulling in DNA ends against a force

        for ii in arange(BindLength,SeqLength-BindLength):
            Covariance=LocalCovariance[ii]
            BendRot=[[np.cos(CurvePhase[ii]), np.sin(CurvePhase[ii])],[-np.sin(CurvePhase[ii]),np.cos(CurvePhase[ii])]]
            CovRot = ( BendRot * Covariance * BendRot.T)		#local covariance matrix alligned to major groove
            BendRot=[[np.cos(np.pi/4), np.sin(np.pi/4)],[-np.sin(np.pi/4),np.cos(np.pi/4)]]
            CovRot45 = ( BendRot * Covariance * BendRot.T)		#local covariance matrix alligned 45� to major groove
            Cnorm=CurveMag[ii]/(2*np.pi*CircFrac)
    #				Cnorm=0				#Uncomment to compare to straight DNA with variable stiffness
            E_base=CircFrac**2*3000/Curvewindow
            Z1=np.exp(-E_base/CovRot[0][0]*((1-Cnorm)**2)+EnergyOffset)			#Bend in direction of curve
            Z2=np.exp(-E_base/CovRot[0][0]*((1+Cnorm)**2)+EnergyOffset)			#Bend against curve
            Z3=np.exp(-E_base/CovRot[1][1]*(1-(Cnorm)**2)+EnergyOffset)			#Bend perpendicular to curve
            Z4=np.exp(-E_base/CovRot45[0][0]*(np.sqrt(Cnorm**2/2+1)-Cnorm/np.sqrt(2))**2+EnergyOffset)		#Bend at 45�, 135�, 225�, and 315�
            Z5=np.exp(-E_base/CovRot45[0][0]*(np.sqrt(Cnorm**2/2+1)+Cnorm/np.sqrt(2))**2+EnergyOffset)
            Z6=np.exp(-E_base/CovRot45[1][1]*(np.sqrt(Cnorm**2/2+1)-Cnorm/np.sqrt(2))**2+EnergyOffset)
            Z7=np.exp(-E_base/CovRot45[1][1]*(np.sqrt(Cnorm**2/2+1)+Cnorm/np.sqrt(2))**2+EnergyOffset)

            Sequence_angle_exp[ii]+=Z1+Z2+2*Z3+Z4+Z5+Z6+Z7





    Sequence_angle_energy=-np.ln(Sequence_angle_exp)  #Backing out the implied energy landscape from the summed Boltzmann weights

    EndEffects=np.max(0,np.min(1,(p-BindLength)/AvePlecLength)*np.min(1,(SeqLength-p-BindLength)/AvePlecLength))  #takes into account the effects of the handles, including limited plectoneme growth near the attachment points
    Sequence_angle_exp=EndEffects

    Sequence_angle_exp_smth=Sequence_angle_exp
    # Smooth/E=2/F/B=64 300,
    savitzky_golay(Sequence_angle_exp_smth,64,3)  #this command applies an approximate Gaussian smooth, though in reality it is 64 sequential boxcar smoothing operations
    # wavestats/q Sequence_angle_exp_smth
    np.mean(Sequence_angle_exp_smth)

