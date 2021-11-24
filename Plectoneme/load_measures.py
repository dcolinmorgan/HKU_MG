import pandas as pd

def load_measures():
  WW=['Dinucleotide','Direction (ϕB, degrees)','Wedge (θ, degrees)','Twist (degrees)','Tilt-Tilt covariance','Roll-Roll covariance']
  AA=['AA',-153.938,7.197,35.606,0.686,1.135]
  AC=['AC',142.942,1.100,34.386,0.649,0.999]
  AG=['AG',1.999,8.397,27.689,0.719,1.175]
  AT=['AT',0.000,2.599,31.487,0.660,0.981]
  CA=['CA',-63.974,3.499,34.486,0.970,1.450]
  CC=['CC',-56.977,2.099,33.656,0.644,1.107]
  CG=['CG',0.000,6.697,29.788,0.960,1.744]
  CT=['CT',-1.999,8.397,27.689,0.719,1.175]
  GA=['GA',119.952,5.298,36.885,0.681,1.264]
  GC=['GC',179.927,4.998,39.984,0.674,0.970]
  GG=['GG',56.977,2.099,33.656,0.644,1.107]
  GT=['GT',-142.942,1.100,34.386,0.649,0.999]
  TA=['TA',0.000,0.900,35.986,1.089,1.962]
  TC=['TC',-119.952,5.298,36.885,0.681,1.264]
  TG=['TG',63.974,3.499,34.486,0.970,1.450]
  TT=['TT',153.938,7.197,35.606,0.686,1.135]
  wave=pd.DataFrame([WW,AA,AC,AG,AT,CA,CC,CG,CT,GA,GC,GG,GT,TA,TC,TG,TT])
  wave.columns=wave.iloc[0]
  wave=wave[1:]
  Dwave=wave['Direction (ϕB, degrees)']
  Wwave=wave['Wedge (θ, degrees)']
  Twave=wave['Twist (degrees)']
  CwaveTT=wave['Tilt-Tilt covariance']
  CwaveRR=wave['Roll-Roll covariance']
  Twave.reset_index(inplace=True, drop=True)
  Wwave.reset_index(inplace=True, drop=True)
  Dwave.reset_index(inplace=True, drop=True)
  CwaveRR.reset_index(inplace=True, drop=True)
  CwaveTT.reset_index(inplace=True, drop=True)
  return Twave, Wwave, Dwave, CwaveRR, CwaveTT