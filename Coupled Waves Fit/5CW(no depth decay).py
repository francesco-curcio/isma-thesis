#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May  1 15:09:28 2022

@author: aaa
"""
"""
This module defines the vector field for 5 coupled wave equations
(without a decay, Uchelnik) and first and second harmonics in the modulation; phase: 0 or pi (sign of n2).
Fit parameters are: n1,n2, d, and wavelength; 
Fit 5/(5) orders!
!!!Data: X,order,INTENSITIES
Fit  background for second orders , first and subtract it for zero orders (background fixed)
"""
from scipy.integrate import ode
from scipy import integrate
import numpy as np
import inspect,os,time
from scipy.optimize import leastsq
from scipy.optimize import least_squares
from scipy.special import erfc
import matplotlib.pyplot as plt
import matplotlib as mpl
import socket
import shutil
from scipy.optimize import curve_fit as fit
from scipy.stats import chisquare as cs
import math


sorted_fold_path="/home/aaa/Desktop/Thesis/Script/Trial/Sorted data/" #insert folder of sorted meausements files
allmeasurements = sorted_fold_path+"All measurements/"
allrenamed = allmeasurements +"All renamed/"
allmatrixes = allmeasurements + "All matrixes/"
allpictures = allmeasurements + "All pictures/"
allrawpictures = allmeasurements + "All raw pictures/"
alldata_analysis = allmeasurements + "All Data Analysis/"
allcropped_pictures = alldata_analysis + "All Cropped Pictures/"
allcontrolplots = alldata_analysis + "All Control plots/"
allcontrolfits = alldata_analysis + "All Control Fits/"
tiltangles=[0,40,48,61,69,71,79,80,81]
foldername=[]
for i in range(len(tiltangles)):
    foldername.append(str(tiltangles[i])+"deg")
foldername.append("79-4U_77c88deg")
foldername.append("79-8U_76c76deg")
foldername.append("79-12U_75c64deg")
foldername.append("79-16U_74c52deg")
n_theta=[26,46,28,17,16,20,21,20,19,48,43,59,24]  #number of measurements files for each folder (no flat, no compromised data)
n_pixel = 16384 #number of pixels in one measurement
"""
This block fits the diffraction efficiencies
"""
for k in range(1,2):#len(foldername)):
    data_analysis = sorted_fold_path+foldername[k]+"/Data Analysis/"
    diff_eff =  np.loadtxt(data_analysis+foldername[k]+'_diff_eff.mpa',skiprows=1)
    #mpl.rc('lines', linewidth=2, color='r')
    mpl.markers.MarkerStyle(marker='*', fillstyle='top') 
    #from Tkinter import Tk
    #from tkFileDialog import askopenfilename
    ################
    pi=np.pi
    deg=np.pi/180.0
    nzero=1. #Mean refractive index of Tomita SiO2
    ###############
    #wavelength=0.0047521668368307763#wavelength
    ####
    '''FIT  with or without error bars'''
    with_err='Y'#fit without error bars: Y, N: fit with.
    ###FOR 70 Degree, measurement of 2018######
    ######################
    bcr1=3.0	#+/-	9.198043505774248e-07	Modulation1
    bcr2=-0.7	#+/-	5.989619220241743e-07	Modulation2
    d=78.	#+/-	1.5379857306604519	Thickness
    #wavelength=0.004755765008537064	#+/-	0.0002979511934439464	wavelength (mum)
    #lc0=4.75e-3
    #lmin0=3e-3
    #lmax0=8e-3
    ####FOR EMG#####
    lambda_par=900.
    mu=4.5e-3
    sigma=7e-4
    ##############
    xc_est=-0.05056382236563003	#+/-	0.010645281056088298	center
    LAMBDA=1.	#+/-	0.0223260000137	LAMBDA
    b_est=2.0*pi*nzero/(mu+1/lambda_par)#propagation constant
    ###CENTER-DATA
    c0=0.#center data, should be ZERO
    bck0_est=101.#background zero order
    bck1_est=83.#.#background first orders
    bck2_est=62.#51.#background second orders
    ###Translate
    G=2.*pi/LAMBDA
    thick=d
    ###assume 1% relative error on intensities
    relaterr=0.01
    ###alpha0=0.00163118577791/2.#2*a1
    ##########################
    zeitstrg=str(time.localtime()[1])+"-"+str(time.localtime()[2])+"/"
    folderstrg="/home/aaa/Desktop/"
    if not os.path.exists(str(folderstrg)+str(zeitstrg)):
            os.mkdir(str(folderstrg)+str(zeitstrg))
    exte=".dat"
    fi="ND120801_7200_70deg_y_25-26-4fit"
    #fi="ND11504_70deg_y_29-33-4fit"
    fullfilename=str(folderstrg)+str(zeitstrg)+str(fi)+str(exte)
    #fullfilename="/home/fallym4/MyPapers/20/NanoDiamonds/ND_11504_70deg/Eval2020/5-29/"+str(fi)+str(exte)
    folder = os.path.split(fullfilename)[0]
    filesolo = os.path.split(fullfilename)[1]
    included_extensions = ['dat']
    included_startnames = ""#'1A@62-SnO3a@70'
    fn = fullfilename
    #LOAD data, filename="..."##
    mypc=socket.gethostname()
    if socket.gethostname()=='TuX':
      filename=fullfilename
    else:
      filename=fullfilename
    mdax=np.loadtxt(filename)
    mdat=np.transpose(mdax)
    markerlist=["*","o","p","<",">","+","$\\leftrightarrow$","+","$\\updownarrow$","^"]
    clist=['k','#ce1818','#7ab742','#0063a6','#fa5b00','#c864c8','#327c32','#6464c8','#fab72a','#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']
    #['1f77b4', 'ff7f0e', '2ca02c', 'd62728', '9467bd', '8c564b', 'e377c2', '7f7f7f', 'bcbd22', '17becf']
    #['r','g','b','k','y',(0.4,0.5,0.2)]
    ct=1.15
    #####DATA######
    angular_factor=0.15#Number of files into angles in Degrees
    tfac=1#every tfac point is considered
    y=mdat[2][::tfac]#intensities
    yerr=mdat[3][::tfac]#error intensities
    xe=angular_factor*mdat[0][::tfac]-ct##external angle in Degree
    print(len(xe))
    z0=-mdat[1][::tfac]#diffraction order
    tbe=np.arcsin(z0*(mu+1/lambda_par)/2.0/LAMBDA)#external Bragg angle
    ###Remove values for which yerr==0 and y==0
    x=xe#internal angles in Degrees
    #print y
    #################
    pi=np.pi
    deg=pi/180.0
    nzero=1. #Mean refractive index of Tomita SiO2
    ###############
    #wavelength=0.0047521668368307763#wavelength
    ####
    '''FIT  with or without error bars'''
    with_err='Y'#fit without error bars: Y, N: fit with.
    ###FOR 70 Degree, measurement of 2018######
    ######################
    bcr1=3.0	#+/-	9.198043505774248e-07	Modulation1
    bcr2=-0.7	#+/-	5.989619220241743e-07	Modulation2
    d=78.	#+/-	1.5379857306604519	Thickness
    #wavelength=0.004755765008537064	#+/-	0.0002979511934439464	wavelength (mum)
    #lc0=4.75e-3
    #lmin0=3e-3
    #lmax0=8e-3
    ####FOR EMG#####
    lambda_par=900.
    mu=4.5e-3
    sigma=7e-4
    ##############
    xc_est=-0.05056382236563003	#+/-	0.010645281056088298	center
    LAMBDA=1.	#+/-	0.0223260000137	LAMBDA
    b_est=2.0*pi*nzero/(mu+1/lambda_par)#propagation constant
    ###CENTER-DATA
    c0=0.#center data, should be ZERO
    bck0_est=101.#background zero order
    bck1_est=83.#.#background first orders
    bck2_est=62.#51.#background second orders
    ###Translate
    G=2.*pi/LAMBDA
    thick=d
    ###assume 1% relative error on intensities
    relaterr=0.01
    ###alpha0=0.00163118577791/2.#2*a1
    ##########################

    markerlist=["*","o","p","<",">","+","$\\leftrightarrow$","+","$\\updownarrow$","^"]
    clist=['k','#ce1818','#7ab742','#0063a6','#fa5b00','#c864c8','#327c32','#6464c8','#fab72a','#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']
    #['1f77b4', 'ff7f0e', '2ca02c', 'd62728', '9467bd', '8c564b', 'e377c2', '7f7f7f', 'bcbd22', '17becf']
    #['r','g','b','k','y',(0.4,0.5,0.2)]
    #ct=1.15
    #####DATA######
    angular_factor=0.15#Number of files into angles in Degrees
    tfac=1#every tfac point is considered
    y=mdat[2][::tfac]#intensities
    yerr=mdat[3][::tfac]#error intensities
    xe=angular_factor*mdat[0][::tfac]-ct##external angle in Degree
    print(len(xe))
    z0=-mdat[1][::tfac]#diffraction order
    tbe=np.arcsin(z0*(mu+1/lambda_par)/2.0/LAMBDA)#external Bragg angle
    ###Remove values for which yerr==0 and y==0
    x=xe#internal angles in Degrees
    #print y
    #################
    xm3=diff_eff[:,0]
    ym3=diff_eff[:,2]*0
    yem3=diff_eff[:,3]*0
    xm2=diff_eff[:,0]
    ym2=diff_eff[:,2]
    yem2=diff_eff[:,3]
    xm1=diff_eff[:,0]
    ym1=diff_eff[:,4]
    yem1=diff_eff[:,5]
    x0=diff_eff[:,0]
    y0=diff_eff[:,6]
    ye0=diff_eff[:,7]
    xp1=diff_eff[:,0]
    yp1=diff_eff[:,8]
    yep1=diff_eff[:,9]
    xp2=diff_eff[:,0]
    yp2=diff_eff[:,10]
    yep2=diff_eff[:,11]
    xp3=xm3*0
    yp3=ym3*0
    yep3=yem3*0
    for i in range(len(z0)):
      if (int(z0[i])==-3):
        xm3.append(x[i])    
        ym3.append(y[i])
        yem3.append(yerr[i])
      elif (int(z0[i])==-2) and (yerr[i]!=0):
        xm2.append(x[i])    
        ym2.append(y[i])
        yem2.append(yerr[i])
      elif (int(z0[i])==-1) and (yerr[i]!=0):
        xm1.append(x[i])    
        ym1.append(y[i])
        yem1.append(yerr[i])
      elif (int(z0[i])==0) and (yerr[i]!=0):
        x0.append(x[i])    
        y0.append(y[i])
        ye0.append(yerr[i])
      elif (int(z0[i])==1) and (yerr[i]!=0):
        xp1.append(x[i])    
        yp1.append(y[i])
        yep1.append(yerr[i])
      elif (int(z0[i])==2) and (yerr[i]!=0):
        xp2.append(x[i])    
        yp2.append(y[i])
        yep2.append(yerr[i])
      elif (int(z0[i])==3):
        xp3.append(x[i])    
        yp3.append(y[i])
        yep3.append(yerr[i])
      else:
        pass
    totl=len(ym2)+len(ym1)+len(y0)+len(yp1)+len(yp2)#TOTAL counts
    print(str(totl)),"points to be fitted alltogether."
    #Angles in degrees (neutrons: external=internal)
    xm2=np.asarray(xm2)
    xm1=np.asarray(xm1)
    x0=np.asarray(x0)
    xp1=np.asarray(xp1)
    xp2=np.asarray(xp2)
    ################
    dym2=ym2[::tfac]
    dym1=ym1[::tfac]
    dy0=y0[::tfac]
    dyp1=yp1[::tfac]
    dyp2=yp2[::tfac]
    ###ERR
    eym2=yem2[::tfac]
    eym1=yem1[::tfac]
    ey0=ye0[::tfac]
    eyp1=yep1[::tfac]
    eyp2=yep2[::tfac]
    #####
    #####
    fitx=xm2.tolist()+xm1.tolist()+x0.tolist()+xp1.tolist()+xp2.tolist()
    fito=(np.ones(len(xm2))*4).tolist()+(np.ones(len(xm1))*3).tolist()+(np.ones(len(x0))*2).tolist()+(np.ones(len(xp1))*1).tolist()+(np.zeros(len(xp2))).tolist()
    fitd=dym2+dym1+dy0+dyp1+dyp2
    fiterr=eym2+eym1+ey0+eyp1+eyp2
    total_counts=[dym2[i]+dym1[i]+dy0[i]+dyp1[i]+dyp2[i] for i in range(len(dym2))]
    ####
    reihe=np.concatenate((dym2,dym1,dy0,dyp1,dyp2))
    totl=len(reihe)#fit only +1,0,-1
    print(str(totl)),"points to be fitted alltogether."
    #################
    #x are INTERNAL angles##
    def dqo(x,b,G,o):
        "Mismatch"
        return b*(np.cos(x)-np.sqrt(1-(np.sin(x)+o*G/b)**2))
    def qz(x,b,G,o):
        "z-component of the wavevector"
        return b*np.cos(x)-dqo(x,b,G,o)
    def kappa_plus(bcr1,b):
        "This is the coupling constant"
        return pi*bcr1/b
    def kappa0(bcr2,b):
      return pi*bcr2/b
    def vectorfield(z, w, p):
        """
        Defines the differential equations for the coupled wave equations.
    
        Arguments:
            w :  vector of the state variables:
                      w = [x1,y1,x2,y2,...]
            z :  depth
            p :  vector of the parameters:
                      p = [m1,m2,k1,k2,L1,L2,b1,b2]
        """
        ym2, ym1, yz, yp1, yp2 = w
        x, b, bcr1, bcr2, d  = p
        f = [(1.0)/(1j*qz(x,b,G,-2))*b*(kappa_plus(bcr1,b)*ym1*np.exp(1j*z*dqo(x,b,G,-1)-1j*z*dqo(x,b,G,-2))+kappa_plus(bcr2,b)*yz*np.exp(-1j*z*dqo(x,b,G,-2))),
             (1.0)/(1j*qz(x,b,G,-1))*b*(kappa_plus(bcr1,b)*yz*np.exp(1j*z*dqo(x,b,G,0)-1j*z*dqo(x,b,G,-1))+kappa_plus(bcr1,b)*ym2*np.exp(1j*z*dqo(x,b,G,-2)-1j*z*dqo(x,b,G,-1))+kappa_plus(bcr2,b)*(yp1*np.exp(-1j*z*(dqo(x,b,G,-1)-dqo(x,b,G,1))))),
             (1.0)/(1j*qz(x,b,G,0))*b*(kappa_plus(bcr1,b)*yp1*np.exp(1j*z*dqo(x,b,G,1))+kappa_plus(bcr1,b)*ym1*np.exp(1j*z*dqo(x,b,G,-1))+kappa_plus(bcr2,b)*(ym2*np.exp(1j*z*(dqo(x,b,G,-2)))+yp2*np.exp(1j*z*dqo(x,b,G,2)))),
             (1.0)/(1j*qz(x,b,G,1))*b*(kappa_plus(bcr1,b)*yp2*np.exp(1j*z*dqo(x,b,G,2)-1j*z*dqo(x,b,G,1))+kappa_plus(bcr1,b)*yz*np.exp(1j*z*dqo(x,b,G,0)-1j*z*dqo(x,b,G,1))+kappa_plus(bcr2,b)*(ym1*np.exp(-1j*z*(dqo(x,b,G,1)-dqo(x,b,G,-1))))),
             (1.0)/(1j*qz(x,b,G,2))*b*(kappa_plus(bcr1,b)*yp1*np.exp(1j*z*dqo(x,b,G,1)-1j*z*(dqo(x,b,G,2)))+kappa_plus(bcr2,b)*yz*np.exp(-1j*z*(dqo(x,b,G,2))))]
        return f
    # ODE solver parameters
    abserr = 1.0e-8
    relerr = 1.0e-6
    numpoints = 2
    
    # Initial conditions
    y_start=[0.0,0.0,1.0,0.0,0.0]
    o_zahl=int(np.floor(len(y_start)/2.))
    def sol(z,p):
        x, b, bcr1, bcr2, d = p
        """definition of function for LS fit
            x angles,
            p is an array of parameters to be varied for fit"""
        # calculate ode solution, return values for each thickness of "z"
        r = ode(vectorfield).set_integrator('zvode',method='adams',atol=relerr,order=5,with_jacobian=False).set_f_params(p)
        r.set_initial_value(y_start,1e-6)
        while r.successful() and r.t<=d:
            r.integrate(r.t+d)
        return np.abs(r.y)**2
    
    def testf(x,b,bcr1,bcr2,d,ords):
        """this is mandatory to get
            proper function for fitting"""
        return sol(d,[x,b,bcr1,bcr2,d])[int(ords)]*qz(x,b,G,int(ords-o_zahl))/qz(x,b,G,0)
    #==============================================================================
    """
    Triangular distribution
    """
    def distribution(l,lmin,lc,lmax):
        if (l>lmin) and (l<=lc):
            return 2*(l-lmin)/((lmin-lmax)*(lmin-lc))
        elif (l>lc) and (l<lmax):
            return 2*(l-lmax)/((lmin-lmax)*(lmax-lc))
        elif (l<lmin):
            return 0
        elif (l>lmax):
            return 0
        else:
            return 0
    def testfi(l,x,lmin,lc,lmax,bcr1,bcr2,d,ords):#interchange variables and replace b by inverse wavelength
        return testf(x,2*pi/l,bcr1,bcr2,d,ords)*distribution(l,lmin,lc,lmax)
    def av_testf(x,lmin,lc,lmax,bcr1,bcr2,d,ords):
        return integrate.quad(testfi,lmin,lc,args=(x,lmin,lc,lmax,bcr1,bcr2,d,ords))[0]+integrate.quad(testfi,lc,lmax,args=(x,lmin,lc,lmax,bcr1,bcr2,d,ords))[0]
    # """
    # Exponentially modified Gauss
    # """
    # low_lim=1e-3
    # upp_lim=15e-3
    # wlrange=np.arange(2e-3,1.5e-2,1e-4)
    # def func(l,A,mu,sig):
    #     return A/(2.)*np.exp(A/(2.)*(2.*mu+A*sig**2-2*l))
    # def emodgauss0(l,A,mu,sig):
    #     return func(l,A,mu,sig)*erfc((mu+A*sig**2-l)/(np.sqrt(2)*sig))
    # emodgauss=np.vectorize(emodgauss0)
    # def meanwl_gauss(l,A,mu,sig):
    #     return emodgauss(l,A,mu,sig)*l
    # def av_mean_gauss(A,mu,sig):
    #     return integrate.quad(emodgauss,low_lim,upp_lim,args=(A,mu,sig))[0]
    # if av_mean_gauss(lambda_par,mu,sigma)<=0.98:
    #     print "INT_GRENZ"
    #     tttt="ACHTUNG!INTGRENZ"+str(av_mean_gauss(lambda_par,mu,sigma))
    # else:
    #     tttt='OK'
    #     pass
    # def testfi(l,x,a,mu,sig,bcr1,bcr2,d,ords):#interchange variables and replace b by inverse (mu+lambda_par)
    #         return testf(x,2*pi/l,bcr1,bcr2,d,ords)*emodgauss(l,a,mu,sig)
    # def wl_average(a,mu,sig):
    #     "Average (mu+lambda_par)"
    #     return a+mu#mu+1/lambda
    # def av_testf(x,a,mu,sig,bcr1,bcr2,d,ords):
    #     return integrate.quad(testfi,low_lim,upp_lim,args=(x,a,mu,sig,bcr1,bcr2,d,ords))[0]
    # #####
    jetzt=time.time()
    eta=np.vectorize(av_testf)
    #print(eta(x0*deg,lambda_par,mu,sigma,bcr1,bcr2,78.,1))
    ####eta=np.vectorize(testf)
    xfit=np.arange(np.min(fitx)-0.025,np.max(fitx)+0.025,0.2)
    bck2_est=np.ones(len(ym2))*bck2_est
    bck1_est=np.ones(len(ym1))*bck1_est
    bck0_est=np.ones(len(ym1))*bck0_est
    xc=0
    try:
        xc_est=xc
    except:
        pass
    ####################
    ####1) EXPORT DATA with background subtracted!
    nenn=(total_counts-2*bck2_est-2*bck1_est-bck0_est)**(-1)
    betap2=(dyp2-bck2_est)*nenn
    betap1=(dyp1-bck1_est)*nenn
    beta0=(dy0-bck0_est)*nenn
    betam1=(dym1-bck1_est)*nenn
    betam2=(dym2-bck2_est)*nenn
    ####ERRORS############
    der_nenn=(bck0_est+2*(bck1_est+bck2_est)-total_counts)**(-2)
    beyp2=np.sqrt(dyp2*nenn**2+
              bck2_est*((bck0_est+2*(bck1_est+dyp2)-total_counts)*der_nenn)**2+
              total_counts*((dyp2-bck2_est)*der_nenn)**2+
              bck0_est*((dyp2-bck2_est)*der_nenn)**2+
        bck1_est*(2*(dyp2-bck2_est)*der_nenn)**2
        )
    beyp1=np.sqrt(dyp1*nenn**2+
              bck1_est*((bck0_est+2*(bck2_est+dyp1)-total_counts)*der_nenn)**2+
              (total_counts)*((dyp1-bck1_est)*der_nenn)**2+
              bck0_est*((dyp1-bck1_est)*der_nenn)**2+
        bck2_est*(2*(dyp1-bck1_est)*der_nenn)**2
        )
    bey0=np.sqrt(dy0*nenn**2+
              bck0_est*((dy0+2*(bck2_est+bck1_est)-total_counts)*der_nenn)**2+
              (total_counts)*((dy0-bck0_est)*der_nenn)**2+
              bck1_est*(2*(dy0-bck0_est)*der_nenn)**2+
        bck2_est*(2*(dy0-bck0_est)*der_nenn)**2
        )
    beym1=np.sqrt(dym1*nenn**2+
              bck1_est*((bck0_est+2*(bck2_est+dym1)-total_counts)*der_nenn)**2+
              (total_counts)*((dym1-bck1_est)*der_nenn)**2+
              bck0_est*((dym1-bck1_est)*der_nenn)**2+
        bck2_est*(2*(dym1-bck1_est)*der_nenn)**2
        )
    beym2=np.sqrt(dym2*nenn**2+
              bck2_est*((bck0_est+2*(bck1_est+dym2)-total_counts)*der_nenn)**2+
              total_counts*((dym2-bck2_est)*der_nenn)**2+
              bck0_est*((dym2-bck2_est)*der_nenn)**2+
        bck1_est*(2*(dym2-bck2_est)*der_nenn)**2
        )
    ##############################
    #COMMENT OUT FOR FITTING
    ############
    #fig=plt.figure("PLOT2020-pre")
    #plt.clf()
    fig, (ax0,ax1)=plt.subplots(2,1,num=1,clear=True,sharex=True)
    #plt.yticks(np.arange(0,1.25,0.25))
    #ax0.errorbar(xm2,betam2,beym2,marker=markerlist[0],label='-2')
    ax0.errorbar(xm2,betam2,beym2,marker=markerlist[0],label='-2')
    #ax0.errorbar(xm1,betam1,beym1,marker=markerlist[0],label='-1')
    ax0.errorbar(xm1,betam1,beym1,marker=markerlist[0],label='-1-new')
    #ax1.errorbar(x0,beta0,bey0,marker=markerlist[0],label='0')
    ax1.errorbar(x0,beta0,bey0,marker=markerlist[0],label='0-new')
    #ax0.errorbar(xp1,betap1,beyp1,marker=markerlist[0],label='+1',alpha=0.5)
    ax0.errorbar(xp1,betap1,beyp1,marker='*',label='+1err',alpha=0.5)
    #ax0.errorbar(xp2,betap2,beyp2,marker=markerlist[0],label='+2')
    ax0.errorbar(xp2,betap2,beyp2,marker=markerlist[0],label='+2')
    ax1.legend()
    ax0.legend()
    ax1.grid()
    ax0.grid()
    ax1.plot(x0,[(dym2[i]+dym1[i]+dy0[i]+dyp1[i]+dyp2[i]) for i in range(len(dym2))])
    for i in range(5):
        if i!=2:
            ax0.plot(xfit,eta((xfit-xc_est)*deg,lambda_par,mu,sigma,bcr1,bcr2,thick,i))
        else:
            ax1.plot(xfit,eta((xfit-xc_est)*deg,lambda_par,mu,sigma,bcr1,bcr2,thick,i))
    plt.text(-1.2,0.8,"$n_1=$"+"{:.3e}".format(bcr1)+"\n$n_2=$"+"{:.2e}".format(bcr2)+"\n$d=$"+"{:.1f}".format(thick)+"$\\mu$m"+"\n$\\lambda=$"+"{:.1f}".format(mu+1/lambda_par)+"$\\mu$m")
    ax0.set_xlim(np.min(xm2),np.max(xm2))
    ax0.set_ylim(0.,np.max((dym1-bck1_est)/(total_counts-2*bck2_est-2*bck1_est-bck0_est)))
    ax1.set_ylim(np.min((dy0-bck0_est)/(total_counts-2*bck2_est-2*bck1_est-bck0_est)),1.)
    plt.xlabel("Angle in Degrees")
    plt.draw()
    #plt.pause(0.5)
    ####
    ##Scalars
    ######
    bck2=np.mean(bck2_est)#
    bck1=np.mean(bck1_est)
    bck0=np.mean([bck0_est,bck0_est])
    ######################
    #######FIT
    ###################################################################
    def residual_two_functions(pars,xm2,xm1,x0,xp1,xp2):
        mod1 = pars[0]#n1
        dora = pars[1]#thickness
        la_par=pars[2]
        wav = np.abs(pars[3])#central wavelength
        wav_max=pars[4]
        mod2 = pars[5]#n2
        cent = pars[6]#center!
    #    bck1 = pars[7]#background first orders
    #    bck2 = pars[8]#background second orders
        nenn=(total_counts-2*(bck2_est+bck1_est)-bck0_est)**(-1)
        betap2=(dyp2-bck2_est)*nenn
        betap1=(dyp1-bck1_est)*nenn
        beta0=(dy0-bck0_est)*nenn
        betam1=(dym1-bck1_est)*nenn
        betam2=(dym2-bck2_est)*nenn
    ####ERRORS############
        der_nenn=(nenn)**(2)
        beyp2=np.sqrt(dyp2*nenn**2+
                  bck2_est*((bck0_est+2*(bck1_est+dyp2)-total_counts)*der_nenn)**2+
                  total_counts*((dyp2-bck2_est)*der_nenn)**2+
                  bck0_est*((dyp2-bck2_est)*der_nenn)**2+
            bck1_est*(2*(dyp2-bck2_est)*der_nenn)**2
            )
        beyp1=np.sqrt(dyp1*nenn**2+
                  bck1_est*((bck0_est+2*(bck2_est+dyp1)-total_counts)*der_nenn)**2+
                  (total_counts)*((dyp1-bck1_est)*der_nenn)**2+
                  bck0_est*((dyp1-bck1_est)*der_nenn)**2+
            bck2_est*(2*(dyp1-bck1_est)*der_nenn)**2
            )
        bey0=np.sqrt(dy0*nenn**2+
                  bck0_est*((dy0+2*(bck2_est+bck1_est)-total_counts)*der_nenn)**2+
                  (total_counts)*((dy0-bck0_est)*der_nenn)**2+
                  bck1_est*(2*(dy0-bck0_est)*der_nenn)**2+
            bck2_est*(2*(dy0-bck0_est)*der_nenn)**2
            )
        beym1=np.sqrt(dym1*nenn**2+
                  bck1_est*((bck0_est+2*(bck2_est+dym1)-total_counts)*der_nenn)**2+
                  (total_counts)*((dym1-bck1_est)*der_nenn)**2+
                  bck0_est*((dym1-bck1_est)*der_nenn)**2+
            bck2_est*(2*(dym1-bck1_est)*der_nenn)**2
            )
        beym2=np.sqrt(dym2*nenn**2+
                  bck2_est*((bck0_est+2*(bck1_est+dym2)-total_counts)*der_nenn)**2+
                  total_counts*((dym2-bck2_est)*der_nenn)**2+
                  bck0_est*((dym2-bck2_est)*der_nenn)**2+
            bck1_est*(2*(dym2-bck2_est)*der_nenn)**2
            )
    #    print(eta(x0*deg,lmin0,lc0,lmax0,bcr1,bcr2,83.,1))
        if with_err=='Y':
            diff0 = (betap2- eta((xp2-cent)*deg,la_par,wav,wav_max,mod1,mod2,dora,0))/beyp2#eta(-2)
            diff1 = (betap1- eta((xp1-cent)*deg,la_par,wav,wav_max,mod1,mod2,dora,1))/beyp1#eta(-1)
            diff2 = (beta0-  eta((x0-cent)*deg,la_par,wav,wav_max,mod1,mod2,dora,2))/bey0#eta(0)
            diff3 = (betam1 - eta((xm1-cent)*deg,la_par,wav,wav_max,mod1,mod2,dora,3))/beym1#eta(+1)
            diff4 = (betam2 - eta((xm2-cent)*deg,la_par,wav,wav_max,mod1,mod2,dora,4))/beym2#eta(-2)
            return np.concatenate((diff0,diff1, diff2, diff3,diff4))#fit -1st AND +1st, 0th redundant
        else:
            diff0 = (betap2- eta((xp2-cent)*deg,la_par,wav,wav_max,mod1,mod2,dora,0))#eta(-2)
            diff1 = (betap1- eta((xp1-cent)*deg,la_par,wav,wav_max,mod1,mod2,dora,1))#eta(-1)
            diff2 = (beta0-  eta((x0-cent)*deg,la_par,wav,wav_max,mod1,mod2,dora,2))#eta(0)
            diff3 = (betam1 - eta((xm1-cent)*deg,la_par,wav,wav_max,mod1,mod2,dora,3))#eta(+1)
            diff4 = (betam2 - eta((xm2-cent)*deg,la_par,wav,wav_max,mod1,mod2,dora,4))#eta(-2)
            return np.concatenate((diff0,diff1, diff2, diff3,diff4))#fit -1st AND +1st, 0th redundant
    ########
    bck1=np.min(dyp1+dym1)-2#background first orders
    bck2=np.min(dyp2+dym2)-2#background second orders
    #bck0_est=bck0
    par_init = np.array([bcr1,thick,lambda_par,mu,sigma,bcr2,xc_est])#estimate for the parameters
    chi2 = np.sum((residual_two_functions(par_init,xm2,xm1,x0,xp1,xp2))**2) 
    dof=(totl-len(par_init))
    from scipy.special import gammaincc
    print ("CHi-SQuare",dof,chi2,gammaincc(dof/2.,chi2/2.))
    print ("Time before start fitting:",str(time.ctime()))
    ##offs=thick
    #################################################################
    jetzt=time.time()
    print ('Fit N or Y')
    input=raw_input('>>')
    #input='Y'
    if input in ['Y','y','YES','yes']:
        print ("Fitting is ongoing, please be patient! Might last for some ten minutes or even hours...")
        ####FIT NOW###
        best, cov, info, message, ier = leastsq(residual_two_functions,
                                                par_init, args=(xm2,xm1,x0,xp1,xp2),
                                                full_output=1)
        print(" Best-Fit Parameters: ",  best)
        fitt=(time.time()-jetzt)/60.0
        print ("Time after fitting finished:",str(time.ctime()))
        print ("Time for fitting:" +str(fitt),str("in Minutes")+"\n")
        chi2 = np.sum((residual_two_functions(best,xm2,xm1,x0,xp1,xp2))**2)
        dof =  (totl-len(best))
        print (str('DoF, ChiSquare, probability:'),dof,str(chi2),gammaincc(dof/2.,chi2/2.))
        [bcr1f,doraf,wavminf,wavf,wavmaxf,bcr2f,centf] = best
        textbest=["bcr1","d","$\\lambda_{par}$","$\\mu_(\\lambda)$","$\\sigma$","bcr2","x0","back1","back2"]
        try:
            for i in range(len(best)):
                print (str(textbest[i]),str(best[i]),"#+/-",str(np.sqrt(np.diag(cov)[i])),str(textbest[i]))
        except:
            for i in range(len(best)):
                print( str(textbest[i]),str(best[i]))
        print ("FILE="+str(filename))
        filename1=str(folder)+"/"+str(mypc)+"fit_NanoDiam2020-"+str(time.localtime()[1])+str(time.localtime()[2])+".fit"
        outf=open(filename1,'a')
        if np.sum(cov)==0:
            	cov=np.zeros([len(best),len(best)])
        else:
            pass
        print (inspect.getfile(inspect.currentframe()))
        progfile=str(inspect.getfile(inspect.currentframe()))
        try:
            gig="back1="+str(bck1_est)+"\n"+"back2="+str(bck2_est)+"\n"+"back0="+str(bck0_est)+"\n"+"ChiSQR="+str(chi2)+"\t"+"Time="+str(time.ctime())+"\t"+"Fitting time in minutes:"+str(fitt)+"\t"+"Total number of points fitted, DoF, Probability:"+str(totl)+","+str(totl-len(par_init))+","+str(gammaincc((totl-len(par_init))/2.,chi2/2.))+"\n"+str(filename)+"\n"+str(progfile)+"\n"+"Fit with Error bars: "+str(with_err)+"\n"+"\n"
            gog="\n".join(map(str,[str(textbest[i])+"="+str(best[i])+"\t"+"#+/-"+"\t"+str(np.sqrt(np.diag(cov)[i])) for i in range(len(best))]))                        
        except:#wo error if not applicable
    #        gig="n1="+str(best[0])+"\t"+"Modulation1"+"\n"+"n2="+str(best[3])+"\t"+"Modulation2"+"\n"+"d="+str(best[1])+"\t"+"Thickness"+"\n"+"lambda="+str(best[2])+"\t"+"wavelength (mum)"+"\n"+"xc="+str(best[4])+"\t"+"center"+"\n"+"ChiSQR="+str(chi2)+"\t"+"Time="+str(time.ctime())+"\t"+"Fitting time in minutes:"+str(fitt)+"\t"+"Total number of points fitted:"+str(totl)+"\n"+str(filename)+"\n"+str(progfile)+"\n"+"Fit with Error bars: "+str(with_err)+"\n"+"\n"
            gog="\n".join(map(str,[str(textbest[i])+"="+str(best[i]) for i in range(len(best))]))
        outf.write(gog)
        outf.write(gig)
        outf.close()
        print ("The END!")
        import inspect, os
        print (inspect.getfile(inspect.currentframe())) # script filename (usually with path)
        print (os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))) # script directory
    else:
        print ("Just plotted, no fit!")
    #    [n1f,doraf,wavf,n2f,centf]=[n1,thick,wavelength,n2,c0]
    ####EXPORT#####
    outputfile1=str(folder)+"/Fit_NDiam_DEQs-2020"+str(time.localtime()[1])+str(time.localtime()[2])+".dat"
    outputfile2=str(folder)+"/Fit_NDiam_DEQs-2020"+str(time.localtime()[1])+str(time.localtime()[2])+".dat2"
    fitfile=str(folder)+"/Fit_NDiam_DEQs-2020"+str(time.localtime()[1])+str(time.localtime()[2])+".txt"
    try:    
        os.remove(outputfile1)
        os.remove(outputfile2)
        os.remove(fitfile)
    except OSError:
        pass
    bck0_est=bck0_est*np.ones(len(total_counts))
    bck1_est=bck1_est*np.ones(len(total_counts))
    bck2_est=bck2_est*np.ones(len(total_counts))
    #####1) EXPORT DATA with background subtracted!
    nenn=(total_counts-2*(bck2_est+bck1_est)-bck0_est)**(-1)
    betap2=(dyp2-bck2_est)*nenn
    betap1=(dyp1-bck1_est)*nenn
    beta0=(dy0-bck0_est)*nenn
    betam1=(dym1-bck1_est)*nenn
    betam2=(dym2-bck2_est)*nenn
    ####ERRORS############
    der_nenn=(nenn)**(2)
    beyp2=np.sqrt(dyp2*nenn**2+
              bck2_est*((bck0_est+2*(bck1_est+dyp2)-total_counts)*der_nenn)**2+
              total_counts*((dyp2-bck2_est)*der_nenn)**2+
              bck0_est*((dyp2-bck2_est)*der_nenn)**2+
        bck1_est*(2*(dyp2-bck2_est)*der_nenn)**2
        )
    beyp1=np.sqrt(dyp1*nenn**2+
              bck1_est*((bck0_est+2*(bck2_est+dyp1)-total_counts)*der_nenn)**2+
              (total_counts)*((dyp1-bck1_est)*der_nenn)**2+
              bck0_est*((dyp1-bck1_est)*der_nenn)**2+
        bck2_est*(2*(dyp1-bck1_est)*der_nenn)**2
        )
    bey0=np.sqrt(dy0*nenn**2+
              bck0_est*((dy0+2*(bck2_est+bck1_est)-total_counts)*der_nenn)**2+
              (total_counts)*((dy0-bck0_est)*der_nenn)**2+
              bck1_est*(2*(dy0-bck0_est)*der_nenn)**2+
        bck2_est*(2*(dy0-bck0_est)*der_nenn)**2
        )
    beym1=np.sqrt(dym1*nenn**2+
              bck1_est*((bck0_est+2*(bck2_est+dym1)-total_counts)*der_nenn)**2+
              (total_counts)*((dym1-bck1_est)*der_nenn)**2+
              bck0_est*((dym1-bck1_est)*der_nenn)**2+
        bck2_est*(2*(dym1-bck1_est)*der_nenn)**2
        )
    beym2=np.sqrt(dym2*nenn**2+
              bck2_est*((bck0_est+2*(bck1_est+dym2)-total_counts)*der_nenn)**2+
              total_counts*((dym2-bck2_est)*der_nenn)**2+
              bck0_est*((dym2-bck2_est)*der_nenn)**2+
        bck1_est*(2*(dym2-bck2_est)*der_nenn)**2
        )
    fitd=list(betap2)+list(betap1)+list(beta0)+list(betam1)+list(betam2)
    fiterr=list(beyp2)+list(beyp1)+list(bey0)+list(beym1)+list(beym2)
    ###DIFFRACTION EFFICIENCY DATA1
    for j in range(len(fitx)):
        wri=open(str(outputfile1),'a')
        exportline=str(fitx[j]-centf)+"\t"+str(fito[j]-o_zahl)+"\t"+str(fitd[j])+"\t"+str(fiterr[j])+"\t"+str(eta((fitx[j]-centf)*deg,wavminf,wavf,wavmaxf,bcr1f,bcr2f,doraf,fito[j]))+"\n"
        wri.write(exportline)
        wri.close()
    ###DIFFRACTION EFFICIENCY DATA2
    #for j in range(len(x0)):
    #    wri=open(str(outputfile2),'a')
    #    exportline=str(x0[j]-centf)+"\t"+"\t".join(map(str,eta((x0[j]-centf)*deg,wavminf,wavf,wavmaxf,bcr1f,bcr2f,doraf,np.arange(0,5,1))))+"\n"
    #    wri.write(exportline)
    #    wri.close()
    #DIFFRACTION EFFICIENCY FIT
    xfit=np.arange(np.min(fitx)-0.025,np.max(fitx)+0.025,0.02)
    for j in range(len(xfit)):
        wri=open(str(fitfile),'a')
        exportl=str(xfit[j]-centf)+"\t"+"\t".join(map(str,eta((xfit[j]-centf)*deg,wavminf,wavf,wavmaxf,bcr1f,bcr2f,doraf,np.arange(0,5,1))))+"\n"
        wri.write(exportl)
        wri.close()
    ###########################################
    ####NOW PLOT WITH FITTED PARAMETERS
    #############################################
    ff=np.loadtxt(fitfile)####FIT
    gg=np.loadtxt(outputfile1)###DATA
    winkel=np.transpose(gg)[0]
    ordnungen=np.transpose(gg)[1]
    dounts=np.transpose(gg)[2]
    eounts=np.transpose(gg)[3]
    ws=winkel[0:len(winkel)/5]
    mata_2=dounts[0:len(dounts)/5]
    mata_1=dounts[len(dounts)/5:2*len(dounts)/5]
    nata_0=dounts[2*len(dounts)/5:3*len(dounts)/5]
    pata_1=dounts[3*len(dounts)/5:4*len(dounts)/5]
    pata_2=dounts[4*len(dounts)/5:]
    itotal=[mata_2[i]+mata_1[i]+pata_2[i]+pata_1[i]+nata_0[i] for i in range(len(mata_1))]#corr total
    #####eta
    mF_2=(mata_2)
    mF_1=(mata_1)
    mF_0=(nata_0)
    pF_1=(pata_1)
    pF_2=(pata_2)
    ####error
    emF_2=eounts[0:len(dounts)/5]
    emF_1=eounts[len(dounts)/5:2*len(dounts)/5]
    eF_0=eounts[2*len(dounts)/5:3*len(dounts)/5]
    epF_1=eounts[3*len(dounts)/5:4*len(dounts)/5]
    epF_2=eounts[4*len(dounts)/5:]
    ###HIER STARTEN 
    mpl.rc('xtick', labelsize=24) 
    mpl.rc('ytick', labelsize=24) 
    plt.figure("NDs")
    plt.clf()
    achse=plt.subplot(111)
    for j in np.arange(1,len(np.transpose(ff))):
        achse.plot(np.transpose(ff)[0],np.transpose(ff)[j],color=clist[j-1])
    achse.errorbar(ws,mF_2,emF_2,marker='*',linestyle='',label='-2',color=clist[0],markersize=10)
    achse.errorbar(ws,mF_1,emF_1,marker='*',linestyle='',label='-1',color=clist[1],markersize=10)
    achse.errorbar(ws,mF_0,eF_0,marker='*',linestyle='',label='0',color=clist[2],markersize=10)
    achse.errorbar(ws,pF_1,epF_1,marker='*',linestyle='',label='+1',color=clist[3],markersize=10)
    achse.errorbar(ws,pF_2,epF_2,marker='*',linestyle='',label='+2',color=clist[4],markersize=10)
    ##plt.text(-1.2,0.3,"$\\overline{n}_1=$"+"{:.3e}".format(n1qf)+"\n$\\overline{n}_2=$"+"{:.2e}".format(n2qf)+"\n$L=$"+"{:.1f}".format(dc)+"$\\mu$m"+"\n$d=$"+"{:.1f}".format(d)+"$\\mu$m")
    achse.text(0.02,0.98,str(time.ctime()))
    achse.set_xlim(-1.1,1.5)
    achse.set_ylim(0.,1.)
    achse.set_xticks(np.arange(-1,1.6,0.5))
    achse.set_yticks(np.arange(0,1.2,0.2))
    achse.set_xlabel("Angle in Degrees",fontsize=36)
    achse.set_ylabel("Diffraction efficiency $\\eta$",fontsize=36)
    #plt.figtext(0.01,0.01,str(time.ctime())+"\n"+str(filename))
    achse.legend(loc=6,prop={'size': 24})
    try:
        shift = max([t.get_window_extent().width for t in achse.legend_.get_texts()])
        for t in achse.legend_.get_texts():
            t.set_ha('right') # ha is alias for horizontalalignment
            t.set_position((shift,0))
    except:
        pass
    plt.grid()
    #plt.pause(0.05)
    plt.draw()
    #print "######"
    #mi=[np.min(pF_2),np.min(pF_1),np.min(mF_1),np.min(emF_2)]
    #mimi=[np.min(pF_2),np.min(pF_1),np.min(mF_1),np.min(mF_2)]
    #if np.min(mi)<0:
    #    print mi
    #elif np.min(mimi)<0:
    #    print mimi
    #else:
    #    print "All right"