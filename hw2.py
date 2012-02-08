import numpy as np
import scipy.special as special
import pf as pf

def ztop(z):
    p = 0.5 + 0.5*special.erf(z/np.sqrt(2))
    return p

def ptoz(p):
    z = np.sqrt(2)*special.erfinv(p*2-1);
    return z

def afc2_p_to_dprime(ph, pf=None):
    if pf is None:
        pf = ph[:,1]
        ph = ph[:,0]

    zH=ptoz( ph )
    zF=ptoz( pf )

    dP=(zH-zF)/np.sqrt(2)
    C=-(zH+zF)/np.sqrt(2)
    lnB=(zF**2-zH**2)/2
    pC=(ph+(1.0-pf))/2

    return dP, C, lnB, pC

def afc2_dprime_to_p(dp, c):
    zH = (dp-c)/np.sqrt(2)
    zF = -(dp+c)/np.sqrt(2)
    ph = ztop( zH)
    pf = ztop( zF)
    return ph,pf

def dprime_to_p(dp, c):
    zH = -c+dp
    zF = -c
    ph = ztop( zH)
    pf = ztop( zF)
    return ph,pf

def fake_dprime(h, f=None):
    if f is None:
        p = h
    else:
        p = (h + (1-f)) / 2.0
    z = ztop(p)
    dprime = z * np.sqrt(2)
    return dprime

# p163
demo_prins_h = numpy.array( [0.6, 0.7, 0.8] )
demo_prins_f = numpy.array( [0.2, 0.2, 0.2] )
demores = afc2_p_to_dprime( demo_prins_h, demo_prins_f )

# p189
prins_tab62_h = numpy.array([0.61,0.69,0.79,0.88,0.97,0.99]) 
prins_tab62_f = numpy.array([0.53,0.42,0.33,0.18,0.06,0.03]) 
hw2res= afc2_p_to_dprime( prins_tab62_h, prins_tab62_f )
unb = afc2_dprime_to_p( hw2res[0], zeros(6) )

# Stan paper: p1425
stanh = numpy.array( [0.95] )
stanf = numpy.array( [0.50] )

figure();
plot( prins_tab62_f, prins_tab62_h, 'o', label='data' )
plot( unb[1], unb[0], '*', label='"unbiased" data' )
plot( [0,1], [0,1], '--' )
xlim( 0,1 )
ylim( 0,1 )
ylabel( 'Hit rate')
xlabel( 'False alarm rate')
title('HW2/Problem 6.2: Biased 2AFC')

# Label the d' of each point
for n in arange( len( prins_tab62_h) ):
    text( prins_tab62_f[n], prins_tab62_h[n], "d'=%.3f" % hw2res[0][n], size=9 )

for dprime in [0.5,1,2]:
    curv=afc2_dprime_to_p( dprime, np.linspace(-10,10,100) )
    plot( curv[1], curv[0], label="d'=%.1f" % dprime )

legend(loc='best')
show()

figure()
StimLevels=arange(6)
Ntrials=ones(6)
ProbitOrLogit=2
LowerAsymptote=0.5

plt.plot(StimLevels,hw2res[3],'b*', label='$Pc$')
plt.plot(StimLevels,unb[0],'rx', label='$Pc_{max}$')

smoothrang = linspace(StimLevels[0], StimLevels[-1], 100 )

# Fit to Pc
# 2 = normal slope (Palamedes way)
# [3,2.0] = guess
paramsfit = pf.fitpf( [3,2.0], StimLevels, hw2res[3], Ntrials, LowerAsymptote, ProbitOrLogit)
# Plot a smoother fitted function
LogLikX, smoothprob=pf.ProbitLogit(paramsfit, smoothrang, hw2res[3], Ntrials, LowerAsymptote, ProbitOrLogit,0)
plt.plot(smoothrang, smoothprob, 'b--', label='Fit to $Pc$');
print "Fit to Pc: %s" % str(paramsfit)

# Fit to Pc_max
paramsfit2 = pf.fitpf( [3,2.0], StimLevels, unb[0], Ntrials, LowerAsymptote, ProbitOrLogit)
LogLikX, smoothprob=pf.ProbitLogit(paramsfit2, smoothrang, hw2res[3], Ntrials, LowerAsymptote, ProbitOrLogit,0)
plt.plot(smoothrang, smoothprob, 'r--', label='Fit to $Pc_{max}$');
print "Fit to Pc_max: %s" % str(paramsfit2)

plt.legend(loc='best')
plt.show()

#fake_dprime = fake_dprime( prins_tab62_h, prins_tab62_f )

crange = linspace( -2, 2, 100. )

hrange = linspace( 0.01, 0.99, 200.)
frange = linspace( 0.01, 0.99, 200.)

print "building dprime grid"
fake_dprime_grid = numpy.array([ [fake_dprime(h,f) for f in frange] for h in hrange])
#dprime_grid = numpy.array([ [afc2_p_to_dprime(h,f)[0] for f in frange] for h in hrange])

dprime_grid = numpy.array([ [(ptoz(h)+ptoz(1-f)) for f in frange] for h in hrange]) / np.sqrt(2)

v = linspace( -3, 3, 50 )
v  = 50
plt.figure();
plt.contourf(hrange, frange, dprime_grid, v ); colorbar()
ylabel('Hit rate'); xlabel( 'False alarm rate')
title('accurate (biased) 2AFC dprime')

plt.figure();
plt.contourf(hrange, frange, fake_dprime_grid, v ); colorbar()
ylabel('Hit rate'); xlabel( 'False alarm rate')
title('naive unbiased 2AFC dprime')

plt.figure();
plt.contourf(hrange, frange, dprime_grid/fake_dprime_grid, 50 ); colorbar()
plot( prins_tab62_f, prins_tab62_h, 'x', label='table 6.2 data' )
plot( [0.5], [0.95], '*', label='Klein, from Green/Sweets p410' )

plt.show()


"""
Dear class,
My software is now working.
The problem with the 2AFC problem is that many psychophysicists are confused about it and factors of sqrt(2). So our discussion on Thursday will hopefully clarify almost everything you wanted to know about signal detection theory and 2AFC.  I'll try to undo all your confusion.  I had thought that the Exercise 2 had asked you to calculate the 2AFC dprime the normal way assuming there is no bias. Thus the assignment I had wanted you to do for  is the data on p. 189 is do calculate d' using PAL_SDT_2AFC_PHFtoDP and also using PAL_SDT_MAFC_PCtoDP.
Hint:   The answer to the latter problem is   dprimeNoBias= erfinv(2*pC-1)*2
 
 I much prefer you to calculate dprimeNoBias rather than the requested PcMax that is complicated.  Kingdom/Prins made the problem more confusing than it had to be (as seen by how simple my Hint is).
 For the plot one would plot the two dprimes (with and without bias). 
  
  Given that this is an extra thing to do (though simpler than the original assignment) you can have extra time say until noon Wednesday. Send me what you have at that point and then if you have further questions you can join me at 4:30 in my office at 420 Minor Hall (not the 520 I said before).  
   
   Stan
"""
