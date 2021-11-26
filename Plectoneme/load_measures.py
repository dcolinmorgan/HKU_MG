import pandas as pd

def load_measures():
    ##taken from https://elifesciences.org/download/aHR0cHM6Ly9jZG4uZWxpZmVzY2llbmNlcy5vcmcvYXJ0aWNsZXMvMzY1NTcvZWxpZmUtMzY1NTctc3VwcDEtdjIuZG9jeA--/elife-36557-supp1-v2.docx?_hash=BGfe9gVycItsz84z%2F%2BXeAdeWYHgCnJhZSGZ485sHQTI%3D
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

    ##taken from https://github.com/kahutia/Plectoneme_prediction/blob/master/DinucleotideParameters.csv
#   WW=['Dinucleotide','Direction (ϕB, degrees)','Twist (degrees)','Wedge (θ, degrees)','Tilt-Tilt covariance','Roll-Roll covariance']
#   AA=['AA',1.004067109,0.612610567,0.022756309,0.685838481,1.135281571]
#   AC=['AC',-0.35877067,0.553269373,0.029824199,0.649126932,0.999217862]
#   AG=['AG',0.262119841,0.567232007,0.074089175,0.719113664,1.175291494]
#   AT=['AT',0,0.520108117,0.017453293,0.660374181,0.980849887]
#   CA=['CA',-0.019605331,0.610865238,0.089028901,0.970149706,1.449952007]
#   CC=['CC',0.139095941,0.574213324,0.088117524,0.644423114,1.1070173]
#   CG=['CG',0,0.600393263,0.095993109,0.959680465,1.743733132]
#   CT=['CT',-0.262119841,0.567232007,0.074089175,0.719113664,1.175291494]
#   GA=['GA',0.668289419,0.619591884,0.042249948,0.680589717,1.264413124]
#   GC=['GC',0,0.588175958,0.020943951,0.673715074,0.970149706]
#   GG=['GG',-0.139095941,0.574213324,0.088117524,0.644423114,1.1070173]
#   GT=['GT',0.35877067,0.553269373,0.029824199,0.649126932,0.999217862]
#   TA=['TA',0,0.65275314,0.043633231,1.088943548,1.961699774]
#   TC=['TC',-0.668289419,0.619591884,0.042249948,0.680589717,1.264413124]
#   TG=['TG',0.019605331,0.610865238,0.089028901,0.970149706,1.449952007]
#   TT=['TT',-1.004067109,0.612610567,0.022756309,0.685838481,1.135281571]



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
  
  
  ## save from https://elifesciences.org/download/aHR0cHM6Ly9jZG4uZWxpZmVzY2llbmNlcy5vcmcvYXJ0aWNsZXMvMzY1NTcvZWxpZmUtMzY1NTctc3VwcDItdjIuZG9jeA--/elife-36557-supp2-v2.docx?_hash=LCk5SxtooaKiuyzBBR1QTeWT4xs6voPDtaecWdf2uMM%3D