import numpy as np
def seaWaterDensity(S,T,depth,lat):

	# sw_pres
	#-------------
	# BEGIN
	#-------------

	DEG2RAD = np.pi/180
	X       = np.sin(np.abs(lat)*DEG2RAD)  # convert to radians
	C1      = 5.92E-3 + ((X**2)*5.25E-3)
	p    = ( (1 - C1) - np.sqrt(((1 - C1)**2) - (8.84E-6*depth)) )/4.42E-6


	# sw_smow
	#----------------------
	# DEFINE CONSTANTS
	#----------------------
	a0 = 999.842594
	a1 =   6.793952e-2
	a2 =  -9.095290e-3
	a3 =   1.001685e-4
	a4 =  -1.120083e-6
	a5 =   6.536332e-9
	T68 = T * 1.00024
	dens_sw_smow = a0 + (a1 + (a2 + (a3 + (a4 + a5*T68)*T68)*T68)*T68)*T68

	# sw_dens0
	#----------------------
	# DEFINE CONSTANTS
	#----------------------
	T68 = T * 1.00024

	#     UNESCO 1983 eqn(13) p17.
	b0 =  8.24493e-1
	b1 = -4.0899e-3
	b2 =  7.6438e-5
	b3 = -8.2467e-7
	b4 =  5.3875e-9
	c0 = -5.72466e-3
	c1 = +1.0227e-4
	c2 = -1.6546e-6
	d0 = 4.8314e-4
	dens_sw_dens0 = dens_sw_smow + (b0 + (b1 + (b2 + (b3 + b4*T68)*T68)*T68)*T68)*S + (c0 + (c1 + c2*T68)*T68)*S*np.sqrt(S) + d0*S**2

	# sw_seck
	#--------------------------------------------------------------------
	# COMPUTE COMPRESSION TERMS
	#--------------------------------------------------------------------
	P = p/10  #convert from db to atmospheric pressure units
	T68 = T * 1.00024
	# Pure water terms of the secant bulk modulus at atmos pressure.
	# UNESCO eqn 19 p 18
	h3 = -5.77905E-7
	h2 = +1.16092E-4
	h1 = +1.43713E-3
	h0 = +3.239908   #[-0.1194975]
	AW  = h0 + (h1 + (h2 + h3*T68)*T68)*T68
	k2 =  5.2787E-8
	k1 = -6.12293E-6
	k0 =  +8.50935E-5   #[+3.47718E-5];
	BW  = k0 + (k1 + k2*T68)*T68
	e4 = -5.155288E-5
	e3 = +1.360477E-2
	e2 = -2.327105
	e1 = +148.4206
	e0 = 19652.21    #[-1930.06];
	KW  = e0 + (e1 + (e2 + (e3 + e4*T68)*T68)*T68)*T68   # eqn 19
	#--------------------------------------------------------------------
	# SEA WATER TERMS OF SECANT BULK MODULUS AT ATMOS PRESSURE.
	#--------------------------------------------------------------------
	j0 = 1.91075E-4
	i2 = -1.6078E-6
	i1 = -1.0981E-5
	i0 =  2.2838E-3
	SR = np.sqrt(S)
	A  = AW + (i0 + (i1 + i2*T68)*T68 + j0*SR)*S
	m2 =  9.1697E-10
	m1 = +2.0816E-8
	m0 = -9.9348E-7
	B = BW + (m0 + (m1 + m2*T68)*T68)*S   # eqn 18
	f3 =  -6.1670E-5
	f2 =  +1.09987E-2
	f1 =  -0.603459
	f0 = +54.6746
	g2 = -5.3009E-4
	g1 = +1.6483E-2
	g0 = +7.944E-2
	K0 = KW + (  f0 + (f1 + (f2 + f3*T68)*T68)*T68 + (g0 + (g1 + g2*T68)*T68)*SR         )*S      # eqn 16
	K = K0 + (A + B*P)*P;  # eqn 15

	# sw_dens
	#------
	# BEGIN
	#------
	densP0 = dens_sw_dens0
	P      = p/10  # convert from db to atm pressure units
	dens   = densP0/(1-P/K)

	return dens
