import numpy as np
import os,sys
sys.path.insert(1, './run/oric/Plectoneme/')
import glob
from tqdm import tqdm
import pandas as pd
import scipy
from savitzky_golay import savitzky_golay
from load_measures import load_measures


def PlectonemeCode(Swave):
    Twave, Wwave, Dwave, CwaveRR, CwaveTT=load_measures()
    SeqLength=len(Swave)
    rise=0.339
    DNApath,DNApathMajorGroove=(np.array([]) for i in range(2))
    BasepairCovariance,LocalCovariance=(np.zeros([SeqLength*2,2]) for i in range(2))
    CurvatureSequence,Sequence_phase,Sequence_angle_energy,Sequence_angle_exp,EndEffects= (np.zeros(SeqLength+1) for i in range(5))
    T_n, Romega_n, Q_n, Rzplus, Rx, Rzminus, Minverse_tot, M_tot= (np.zeros([4,4]) for i in range(8))
    StartPos=[0,0,0,1]
    StartPosMG=[1,0,0,1]
    Minverse_tot=np.identity(4)
    M_tot=Minverse_tot
    T_n=M_tot
    T_n[3,2]=-rise/2
    
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
    
        index = 4*Letter1+Letter2	#The index defines the current dinucleotide, AA=0, AC=1,...TT=15
    
        Sequence_phase[ii]=Sequence_phase[ii]+Twave[index]	#This is used to measure how far around the DNA the major groove has rotated relative to the first base pair
        BendRot=np.array([[np.cos(Sequence_phase[ii]), np.sin(Sequence_phase[ii])],[-np.sin(Sequence_phase[ii]),np.cos(Sequence_phase[ii])]]) #Rotation matrix
        Covariance=np.array([[CwaveRR[index], 0],[0,CwaveTT[index]]])	#The covariance matrix for the current basepair, expressed in the coordinates of the current basepair
        CovRot =  np.multiply(BendRot,Covariance,np.transpose(BendRot)) #Rotating the covariance matrix so it will line up with its neighbors
        BasepairCovariance[ii]=np.sum(CovRot,axis=1)	#Rotated covariance matrix at position is recorded
        omDiv2_n=Twave[index]/2
        Romega_n=np.array([[np.cos(omDiv2_n),np.sin(omDiv2_n),0,0],[-np.sin(omDiv2_n),np.cos(omDiv2_n),0,0],[0,0,1,0],[0,0,0,1]])
        alpha_n=Wwave[index]
        beta_n=Dwave[index]-np.pi/2
        Rzplus=np.array([[np.cos(beta_n),np.sin(beta_n),0,0],[-np.sin(beta_n),np.cos(beta_n),0,0],[0,0,1,0],[0,0,0,1]])
        Rx=np.array([[1,0,0,0],[0,np.cos(-alpha_n),np.sin(-alpha_n),0],[0,-np.sin(-alpha_n),np.cos(-alpha_n),0],[0,0,0,1]])
        Rzminus=np.array([[np.cos(-beta_n),np.sin(-beta_n),0,0],[-np.sin(-beta_n),np.cos(-beta_n),0,0],[0,0,1,0],[0,0,0,1]])
    
        Q_n = np.multiply(Rzminus, Rx, Rzplus)
        Minverse_n = np.linalg.inv( T_n* Romega_n* Q_n* Romega_n* T_n)
        Minverse_tot =  np.multiply( Minverse_n,Minverse_tot)
        Minverse_tot=Minverse_n #Updating the total tranformation matrix
        CurrentPos = np.round(np.sum(Minverse_tot,axis=1),0)*StartPos #Calculate the coordinates of the current basepair
        CurrentPosMG =  np.round(np.sum(Minverse_tot,axis=1),0)*StartPosMG
        DNApath=np.append(DNApath,CurrentPos)
        DNApathMajorGroove=np.append(DNApathMajorGroove,CurrentPosMG)
        Letter1=Letter2
    
    DNApath=np.reshape(DNApath,(int(len(DNApath)/4),4))
    DNApathMajorGroove=np.reshape(DNApathMajorGroove,(int(len(DNApathMajorGroove)/4),4))
    print("path calculated")
    #Make the curvature calculation
    # variable CurveWindow  #must be even
    CircFrac=0.667	#We assume the plectoneme tip makes a 240� arc before joining the bulk plectoneme region
    BindLength=45   #experimentally, ~450 nt are bound to the surface at each end of the DNA
    AvePlecLength=1000
    tanlength=10 #must be even. This is number of basepairs used to calculate the local tangent vectors.
    # Variable VectorMag, CosCurve, SinCurve
    
    for CurveWindow in tqdm(np.arange(40,120,8)):
        LocalCovariance[:,0]=savitzky_golay(BasepairCovariance[:,0],CurveWindow+1,3)	#Find covariance matrix over the curvature window
        LocalCovariance[:,1]=savitzky_golay(BasepairCovariance[:,1],CurveWindow+1,3)
        TanVector, NormVector,CurveVector=(np.zeros([SeqLength,3]) for i in range(3))
        CurveMag, CurvePhase, HalfCurveMag, HalfCurvePhase=(np.zeros([SeqLength]) for i in range(4))
        NormVector=DNApathMajorGroove[:,0:3]-DNApath[:,0:3]  #identifies the normal vector alligned with the major groove
        CurrentTan, CurrentCurve, tP, tM, CurveCross, CurrentNorm=(np.zeros(3).T for i in range(6))
        # find the tan vectors over tanlength
        half_tan_len=int(tanlength/2)

        for ii in np.arange(half_tan_len,SeqLength-half_tan_len):
            CurrentTan=DNApath[ii+half_tan_len]-DNApath[ii-half_tan_len]
            VectorMag=np.sqrt(np.dot(CurrentTan, CurrentTan))
            TanVector=np.append(TanVector,CurrentTan/(VectorMag))  #Normalizes tangent vector to unit length
        TanVector=np.reshape(TanVector,(int(len(TanVector)/3),3))
        half_curve=int(CurveWindow/2)
        # find the curvature vectors and values over CurveWindow
        for ii in np.arange(half_curve, SeqLength-half_curve):
            tP=TanVector[ii+half_curve]	#plus tan vector
            tM=TanVector[ii-half_curve]	#minus tan vector
            CurrentCurve =[tP[1]*tM[2]-tP[2]*tM[1],tP[2]*tM[0]-tP[0]*tM[2],tP[0]*tM[1]-tP[1]*tM[0]]   #Cross product
            CurveVector=np.append(CurveVector,CurrentCurve)	#curvature vector is recorded at this position
            CurveMag=np.append(CurveMag,np.arccos(np.dot(tP, tM)))
            CurrentCurve=CurrentCurve/np.arccos(np.dot(tP, tM))	    #normalize the curvature vector to track direction
            #Calculates the phase angle of the curvature relative to major groove at start of DNA
            CurrentTan=TanVector[ii]
            CurrentNorm=NormVector[ii]
            CosCurve=np.dot(CurrentCurve, CurrentNorm)
            tP=CurrentNorm	#Only using these variables as placeholders b/c of short names,
            #otherwise cross product calculation makes a long line of code
            tM=CurrentCurve
            CurveCross=[tP[1]*tM[2]-tP[2]*tM[1],tP[2]*tM[0]-tP[0]*tM[2],tP[0]*tM[1]-tP[1]*tM[0]]
            SinCurve=np.dot(CurveCross, CurrentTan)
            CurvePhase=np.append(CurvePhase,np.mod(np.arctan2(SinCurve,CosCurve)+Sequence_phase[ii],2*np.pi)) #curvature phase is recorded.
        CurveMag=np.reshape(CurveMag,(int(len(CurveMag)/1),1))
        EnergyOffset=25-CurveWindow*0.334*3/4.06   #adds energy penalty from pulling in DNA ends against a force
        
        for ii in np.arange(BindLength,SeqLength-BindLength):
            Covariance=LocalCovariance[ii]
            BendRot=[[np.cos(CurvePhase[ii]), np.sin(CurvePhase[ii])],[-np.sin(CurvePhase[ii]),np.cos(CurvePhase[ii])]]
            CovRot = ( BendRot * Covariance * np.transpose(BendRot))    #local covariance matrix alligned to major groove
            BendRot=[[np.cos(np.pi/4), np.sin(np.pi/4)],[-np.sin(np.pi/4),np.cos(np.pi/4)]]
            CovRot45 = ( BendRot * Covariance * np.transpose(BendRot))  #local covariance matrix alligned 45� to major groove
            Cnorm=CurveMag[ii]/(2*np.pi*CircFrac)
    #				Cnorm=0				#Uncomment to compare to straight DNA with variable stiffness
            E_base=CircFrac**2*3000/CurveWindow
            Z1=np.exp(-E_base/CovRot[0][0]*((1-Cnorm)**2)+EnergyOffset)	#Bend in direction of curve
            Z2=np.exp(-E_base/CovRot[0][0]*((1+Cnorm)**2)+EnergyOffset)	#Bend against curve
            Z3=np.exp(-E_base/CovRot[1][1]*(1-(Cnorm)**2)+EnergyOffset)	#Bend perpendicular to curve
            Z4=np.exp(-E_base/CovRot45[0][0]*(np.sqrt(Cnorm**2/2+1)-Cnorm/np.sqrt(2))**2+EnergyOffset)
            Z5=np.exp(-E_base/CovRot45[0][0]*(np.sqrt(Cnorm**2/2+1)+Cnorm/np.sqrt(2))**2+EnergyOffset)
            Z6=np.exp(-E_base/CovRot45[1][1]*(np.sqrt(Cnorm**2/2+1)-Cnorm/np.sqrt(2))**2+EnergyOffset)
            Z7=np.exp(-E_base/CovRot45[1][1]*(np.sqrt(Cnorm**2/2+1)+Cnorm/np.sqrt(2))**2+EnergyOffset)
    
            Sequence_angle_exp[ii]=Z1+Z2+2*Z3+Z4+Z5+Z6+Z7
    
    
    Sequence_angle_energy=-np.log(Sequence_angle_exp)  #Backing out the implied energy landscape from the summed Boltzmann weights
    
    Sequence_angle_exp_smth=Sequence_angle_exp
    np.nan_to_num(Sequence_angle_exp_smth, copy=True, nan=0.0, posinf=None, neginf=None)
    for i in tqdm(np.arange(300)):
        Sequence_angle_exp_smth=savitzky_golay(Sequence_angle_exp_smth,65,3)		#Find covariance matrix over the curvature window
    return Sequence_angle_exp_smth
