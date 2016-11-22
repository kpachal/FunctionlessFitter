
print 'Case 1: scaled.'

def Anum1(w,x,y,z) :
  return (z*10.0-y*10.0)

def Adenom1() :
  return (2459.0-2403.5)

def Anum2(w,x,y,z) :
  return (y*10.0-x*10.0)

def Adenom2() :
  return (2403.5-2349.0)

def Abigdenom1() :
  return (2459.0-2349.0)

def Anum3(w,x,y,z) :
  return (y*10.0-x*10.0)

def Adenom3() :
  return (2403.5-2349.0)

def Anum4(w,x,y,z) :
  return (x*10.0-w*10.0)

def Adenom4() :
  return (2349.0-2295.5)

def Abigdenom2() :
  return (2403.5-2295.5)

def Abignum1(w,x,y,z,doprint=False) :
  if doprint :
    print Anum1(w,x,y,z),"/",Adenom1(), "-", Anum2(w,x,y,z),"/",Adenom2(), "=",Anum1(w,x,y,z)/Adenom1() - Anum2(w,x,y,z)/Adenom2()
  return Anum1(w,x,y,z)/Adenom1() - Anum2(w,x,y,z)/Adenom2()

def Abignum2(w,x,y,z,doprint=False) :
  if doprint :
    print Anum3(w,x,y,z),"/",Adenom3(), "-", Anum4(w,x,y,z),"/",Adenom4(), "=", Anum3(w,x,y,z)/Adenom3() - Anum4(w,x,y,z)/Adenom4()
  return Anum3(w,x,y,z)/Adenom3()-Anum4(w,x,y,z)/Adenom4()

def AmyAll(w,x,y,z,doprint=False) :
  val = - Abignum1(w,x,y,z,True)/Abigdenom1() + Abignum2(w,x,y,z,True)/Abigdenom2()
  if doprint :
    print - Abignum1(w,x,y,z),"/",Abigdenom1(), "+", Abignum2(w,x,y,z),"/",Abigdenom2(), " = ",- Abignum1(w,x,y,z)/Abigdenom1() + Abignum2(w,x,y,z)/Abigdenom2()
  return val


w = 1000.0
x = 100.0
y = 10.0
z = 1.0

print "Result from test code:"
AmyAll(w,x,y,z,True)
print "Result using exactly the same format:", - ((z*10.0-y*10.0)/(2459.0-2403.5)-(y*10.0-x*10.0)/(2403.5-2349.0))/(2459.0-2349.0) + ((y*10.0-x*10.0)/(2403.5-2349.0)-(x*10.0-w*10.0)/(2349.0-2295.5))/(2403.5-2295.5)
#print "literally just copied",- ((1.710009615*10.0-2.07465910503*10.0)/(2459.0-2403.5)-(2.07465910503*10.0-2.45669574977*10.0)/(2403.5-2349.0))/(2459.0-2349.0) + ((2.07465910503*10.0-2.45669574977*10.0)/(2403.5-2349.0)-(2.45669574977*10.0-2.8635012513*10.0)/(2349.0-2295.5))/(2403.5-2295.5)
print "Result from my code:",0.000148406966792


print 'Case 2: unscaled.'

def Bnum1(w,x,y,z) :
  return (z*1.0-y*1.0)

def Bdenom1() :
  return (2459.0-2403.5)

def Bnum2(w,x,y,z) :
  return (y*1.0-x*1.0)

def Bdenom2() :
  return (2403.5-2349.0)

def Bbigdenom1() :
  return (2459.0-2349.0)

def Bnum3(w,x,y,z) :
  return (y*1.0-x*1.0)

def Bdenom3() :
  return (2403.5-2349.0)

def Bnum4(w,x,y,z) :
  return (x*1.0-w*1.0)

def Bdenom4() :
  return (2349.0-2295.5)

def Bbigdenom2() :
  return (2403.5-2295.5)

def Bbignum1(w,x,y,z,doprint=False) :
  if doprint :
    print Bnum1(w,x,y,z),"/",Bdenom1(), "-", Bnum2(w,x,y,z),"/",Bdenom2(), "=", Bnum1(w,x,y,z)/Bdenom1() - Bnum2(w,x,y,z)/Bdenom2()
  return Bnum1(w,x,y,z)/Bdenom1() - Bnum2(w,x,y,z)/Bdenom2()

def Bbignum2(w,x,y,z,doprint = False) :
  if doprint :
    print Bnum3(w,x,y,z),"/",Bdenom3(), "-", Bnum4(w,x,y,z),"/",Bdenom4(), "=", Bnum3(w,x,y,z)/Bdenom3() - Bnum4(w,x,y,z)/Bdenom4()
  return Bnum3(w,x,y,z)/Bdenom3() - Bnum4(w,x,y,z)/Bdenom4()

def BmyAll(w,x,y,z,doprint=False) :
  val = - Bbignum1(w,x,y,z,True)/Bbigdenom1() + Bbignum2(w,x,y,z,True)/Bbigdenom2()
  if doprint :
    print - Bbignum1(w,x,y,z),"/",Bbigdenom1(), "+", Bbignum2(w,x,y,z),"/",Bbigdenom2(), " = ",- Bbignum1(w,x,y,z)/Bbigdenom1() + Bbignum2(w,x,y,z)/Bbigdenom2()
  return val

w = 1000.0
x = 100.0
y = 10.0
z = 1.0

print "Result from test code:"
BmyAll(w,x,y,z,True)
print "Result using exactly the same format:",- ((z*1.0-y*1.0)/(2459.0-2403.5)-(y*1.0-x*1.0)/(2403.5-2349.0))/(2459.0-2349.0) + ((y*1.0-x*1.0)/(2403.5-2349.0)-(x*1.0-w*1.0)/(2349.0-2295.5))/(2403.5-2295.5)
#print "literally just copied",- ((16.9085167465*1.0-20.2242447666*1.0)/(2459.0-2403.5)-(20.2242447666*1.0-24.0734621732*1.0)/(2403.5-2349.0))/(2459.0-2349.0) + ((20.2242447666*1.0-24.0734621732*1.0)/(2403.5-2349.0)-(24.0734621732*1.0-28.4238108445*1.0)/(2349.0-2295.5))/(2403.5-2295.5)
print "Result from my code:",-1.27935856353e-17

