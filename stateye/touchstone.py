class touchstone :
	def __init__(self) :
		self.version		='v5_2'
		print "Loading touchstone module version "+self.version

		self.filename		= ''
		self.nports 		= 4
		self.victim_tx_P 	= 1		
		self.victim_tx_N 	= 3
		self.victim_rx_P 	= 2
		self.victim_rx_N 	= 4
		self.noDC 		= 1
		self.cascade		= 1
		self.DCextrapolation	= 6
		self.resoltuion 	= 1e-12	 	# no default possible

		# set these defaults to SAS
		self.txamp		= 1.0
		self.r			= 45
		self.c			= 950e-15 
		self.pole		= 10.0e9
		self.loss		= 0

		self.frequency = []
		# currently I'm keeping these as lists, but I may change this to a comple 3-d array
		self.raw_s = []
		self.s = []
		self.t = []

	def calcStep(self, txsignal, noise) :
		from numpy import concatenate, cumsum, arange, flipud, conj, real, array, zeros
		import numpy
		from string import rsplit
		from re import split
		from pylab import find
		import pdb

		resolution = self.resolution

		# currently no padding
		t = concatenate(( flipud(conj(self.t[1:])), self.t ))
		impulse = real(numpy.fft.ifft(numpy.fft.ifftshift(t)))
		timeStep = 1.0/(2.0*max(self.frequency))
		timeAxis = arange(len(impulse)) * timeStep
		self.impulse = numpy.interp(arange(timeAxis[0],timeAxis[-1],resolution),timeAxis,impulse)
		self.impulse = self.impulse * (resolution/timeStep)

		# if we don't have a provided txsignal then assume standard infinite step input
		if txsignal == [] :
			self.impulse = concatenate(( zeros(len(self.impulse)), self.impulse, zeros(len(self.impulse)) ))
			self.step = cumsum(real(self.impulse)) * 0.5 * self.txamp
		else :

			# this 0.5 comes about, because it comes from a '0' -> '1' transition of a measurement
			# and is not based on the pulse step, i.e. 'Z' -> '1'
			txsignal = 0.5 * array(txsignal)
			txsignal = txsignal - txsignal[0]
			txtime = arange(len(txsignal)) * self.resolution
			_txsignal = numpy.interp(arange(timeAxis[0],timeAxis[-1],resolution),txtime,txsignal)
			f1 = numpy.fft.fft(_txsignal)
			f2 = numpy.fft.fft(self.impulse)
			f3 = f1 * f2
			self.step = real(numpy.fft.ifft(f3))
			i = find(min(self.step)==self.step)[0]
			self.step = self.step - self.step[i]
			self.step[:i] = 0.0

		if noise != [] :
			noise = array(noise)
			noise = noise - noise[0]
			_noise = numpy.interp(arange(timeAxis[0],timeAxis[-1],resolution),txtime,noise)
			f1 = numpy.fft.fft(_noise)
			f2 = numpy.fft.fft(self.impulse)
			f3 = f1 * f2
			self.noise = real(numpy.fft.ifft(f3))


	def map(self) :
		from numpy import array,zeros
		p = [self.victim_tx_P-1,self.victim_tx_N-1,self.victim_rx_P-1,self.victim_rx_N-1 ]
		#p = (array(p)-1).tolist()
		m = [1-1,2-1,5-1,6-1]
		for _raw in self.raw_s :
			_s = (zeros(( 8,8 ))*0.0).tolist()
			_s[3-1][7-1] = 1.0
			_s[4-1][8-1] = 1.0
			_s[7-1][3-1] = 1.0
			_s[8-1][4-1] = 1.0
			for i in range(4) :
				for j in range(4) :
					_s[m[i]][m[j]] = _raw[p[i]][p[j]]

			self.s += [_s]
		self.nports = 8


	def cascadeSimple(self) :
		from numpy import zeros, dot,pi, exp
		from touchstone import s2t, t2s
		import pdb
		
		self.rl = []
		self.h = []
		self.l = []

		r = self.r
		c = self.c
		pole = self.pole
		loss = self.loss
		txamp = self.txamp

		_sPRE = zeros(( self.nports,self.nports )) * (0.0+1j*0.0)
		_sPOST = zeros(( self.nports,self.nports )) * (0.0+1j*0.0)
		_sMID = zeros(( self.nports,self.nports )) * (0.0+1j*0.0)
		for ii in range(len(self.s)) :
			_s = self.s[ii]
			_f = self.frequency[ii]
			_t = s2t(_s)
			
			Z  = 1.0 / ( 1.0/r + _f * 1j * 2 * pi * c) 
			RL = ( Z - 50 ) / ( Z + 50 ) 
			H  = 1.0 / ( 1.0 + ( 2.0 * pi * _f * 1j ) / ( 2.0 * pi * pole ) ) 
			L  = 10.0**(loss/1.0e9*_f/20) * exp(2.0 * -1j * pi * _f / (1.0 / 1e-9))
			#L  = exp(2.0 * -1j * pi * _f / (1.0 / 1e-9))
			self.rl += [RL]
			self.h  += [H]
			self.l  += [L]
			for i in range(self.nports) :
				_sPRE[i][i] = RL 
				_sPOST[i][i] = RL 
				_sMID[i][i] = 0 
			for i in range(self.nports/2) :
				_sPRE[i][i+self.nports/2] = H
				_sPRE[i+self.nports/2][i] = H
				_sPOST[i][i+self.nports/2] = H
				_sPOST[i+self.nports/2][i] = H
				_sMID[i][i+self.nports/2] = L
				_sMID[i+self.nports/2][i] = L

			# pdb.set_trace()
			_tPRE = s2t(_sPRE)
			_tMID = s2t(_sMID)
			_tPOST = s2t(_sPOST)
	
			_tCASCADE = dot( dot( dot( _tPRE, _t ), _tMID), _tPOST )

			_sCASCADE = t2s(_tCASCADE)

			self.s[ii] = _sCASCADE
					


	def extractTransfer(self) :
		from numpy import array, dot, zeros
		# this is currently hard coded and represents where the touchstone
		# matrix indexes were mapped into the 4x4 matrix
		pi = 1
		ni = 2
		pj = 5
		nj = 6
		T = array( [	[+1.0,+1.0,+0.0,+0.0],\
				[+1.0,-1.0,+0.0,+0.0],\
				[+0.0,+0.0,+1.0,+1.0],\
				[+0.0,+0.0,+1.0,-1.0]])		
		Tp = array([	[+0.5,+0.5,+0.0,+0.0],\
				[+0.5,-0.5,+0.0,+0.0],\
				[+0.0,+0.0,+0.5,+0.5],\
				[+0.0,+0.0,+0.5,-0.5]])
		map = [pi-1,ni-1,pj-1,nj-1]

		for _s in self.s :
			s = zeros((4,4)) * (0.0+0.0j)
			for i in range(4) :
				for j in range(4) :
					s[i][j] = _s[map[i]][map[j]]
			_t = dot(dot(T,s),Tp)
			# clearly this needs extending to include the other modes
			self.t += [_t[1,3]]

	def get(self,i,j) :
		r = []
		for _s in self.s :
			r += [_s[i-1][j-1]]
		return r

	# going to assume here that we only have one missing point!!! Needs extending for the final release
	def addDC(self) :
		from numpy import zeros, log10, absolute, polyfit, arange, diff, poly1d, angle, exp, transpose, flipud, floor, unwrap, pi
		#from numpy.lib import linspace
		from scipy import interp
		import pdb
		import math
		from pylab import plot,figure,hold

		if self.noDC == 0 :
			print 'Seem like we already have a DC point identified'
			return 

		n = self.DCextrapolation

		frequency = arange(0.0,self.frequency[0],diff(self.frequency[:2])[0])
		
		t = zeros(len(frequency)) * (0.0+0.0j)
		mag =[]
		for k in range(n) :
			mag += [log10((absolute(self.t[k])))]
		p = poly1d(polyfit(self.frequency[:n],mag,1))
		mag = p(frequency)
		# calculate r, when minimum frequency is too high and a 2pi wrap is necessary
		angleStep = min(diff(unwrap(angle(self.t[0:10]))))
		phase = arange(len(frequency)) * angleStep
		# pdb.set_trace()

		t = 10.0**mag * exp(1j*phase)
		figure()
		plot(self.frequency[0:20],angle(self.t[0:20]))		
		hold(1)
		plot(self.frequency[0:20],angle(self.t[0:20]),'o')		
		#print 'we are here'
		self.t = t.tolist() + self.t
		self.frequency = frequency.tolist() + self.frequency
	
		plot(self.frequency[0:20],angle(self.t[0:20]))		
		plot(self.frequency[0:20],angle(self.t[0:20]),'x')		


	def loadFile(self) :
		from string import upper,lstrip,rstrip
		from re import split, match
		from numpy import zeros
		from math import sin,cos,pi
		import sys

		i = -1
		j = -1
		part = 0

		nports 		= self.nports
		filename 	= self.filename	

		for line in open(filename) :
			line = lstrip(rstrip(line))
			
			#print 'found line <%s>'%line
			#if ((len(line)>0) and (not line.startswith('!'))) :
			if line.startswith('!') or len(line)==0:
				pass # print "Ignoring comment"
			else:
				col = split('\s*', line)
				if col[0]=='#':
					# line format is: <frequency unit> <parameter> <format> R <n>
					# specifically: [HZ/KHZ/MHZ/GHZ] [S/Y/Z/G/H] [MA/DB/RI] [R n]
					if match('HZ', upper(col[1])) :
						Kfrequency = 1.0
					if match('KHZ', upper(col[1])) :
						Kfrequency = 10.0**3
					if match('MHZ', upper(col[1])) :
						Kfrequency = 10.0**6
					if match('GHZ', upper(col[1])) :
						Kfrequency = 10.0**9
					# RI for real-imaginary, MA for magnitude-angle, DB for dB-angle
					Ktype = upper(col[3])
				else :
					#print 'found data <%s>'%col
					for _col in col :
						if (i==-1) and (j==-1) :
							i = 0
							j = 0
							self.frequency += [eval(col[0]) * Kfrequency]
							s = zeros((nports,nports)) * (0.0+0.0j)
							#print 'found frequency %e'%self.frequency[-1]
						else :
							if part==0 :
								oldcol = _col
								part = 1
							else :
								if Ktype=='RI' :
									s[i][j] = eval(oldcol) + 1.0j*eval(_col)
								elif Ktype=='MA' :
									s[i][j] = eval(oldcol) * cos(eval(_col) * pi / 180.0) \
									+ 1.0j * eval(oldcol) * sin(eval(_col) * pi / 180.0)
									#print 'added to %d,%d %e'%(i,j,s[i][j])
								elif Ktype=='DB' :
									magnitude = 10**(eval(oldcol)/20)
									#print "Magnitude %f"%magnitude
									s[i][j] = magnitude * cos(eval(_col) * pi / 180.0) \
									+ 1.0j * magnitude * sin(eval (_col) * pi / 180.0)
								else :
									print 'Unknown format%s'%Ktype
									sys.exit()
								part = 0
								i += 1
								if i==nports :
									i=0
									j+=1
									if j==nports :
										j=-1
										i=-1
										self.raw_s += [s]
										#print 'end of array'
						
								
						
def s2t(s) :
	from numpy import dot,concatenate,transpose,linalg,identity,pi

	s_a = transpose(transpose(s[0:4])[0:4])
	s_b = transpose(transpose(s[0:4])[4:8])
	s_g = transpose(transpose(s[4:8])[0:4])
	s_t = transpose(transpose(s[4:8])[4:8])
	
	si = linalg.inv(s)
	
	si_a = transpose(transpose(si[0:4])[0:4])
	si_b = transpose(transpose(si[0:4])[4:8])
	si_g = transpose(transpose(si[4:8])[0:4])
	si_t = transpose(transpose(si[4:8])[4:8])

	# print (identity(4) - dot(s_a, si_a)) 
	# print linalg.inv(identity(4) - dot(s_a, si_a)) 

	t_a = dot(dot( linalg.inv(identity(4) - dot(s_a, si_a)) , s_a), si_b)
	t_b = 		     dot( linalg.inv(identity(4) - dot(s_a, si_a)) , s_b)
	t_g = 		     dot( linalg.inv(identity(4) - dot(si_a, s_a)) , si_b)
	t_t = dot(dot( linalg.inv(identity(4) - dot(si_a, s_a)) , si_a), s_b)
	
	t = concatenate((transpose(concatenate((transpose(t_b),transpose(t_a)))), transpose(concatenate((transpose(t_t),transpose(t_g)))) )) 

	return(t)


def t2s(t) :
	from numpy import dot,concatenate,transpose,linalg,identity,pi

	t_b = transpose(transpose(t[0:4])[0:4])
	t_a = transpose(transpose(t[0:4])[4:8])
	t_t = transpose(transpose(t[4:8])[0:4])
	t_g = transpose(transpose(t[4:8])[4:8])

	t = concatenate((transpose(concatenate((transpose(t_a),transpose(t_b)))), transpose(concatenate((transpose(t_g),transpose(t_t)))) ))

	ti = linalg.inv(t)

	ti_a = transpose(transpose(ti[0:4])[0:4])
	ti_b = transpose(transpose(ti[0:4])[4:8])
	ti_g = transpose(transpose(ti[4:8])[0:4])
	ti_t = transpose(transpose(ti[4:8])[4:8])

	s_a = dot(dot( linalg.inv(identity(4) - dot(t_a,  ti_a)) , t_a),  ti_b)
	s_b = 		     dot( linalg.inv(identity(4) - dot(t_a,  ti_a)) , t_b)
	s_g = 		     dot( linalg.inv(identity(4) - dot(ti_a, t_a)) ,  ti_b)
	s_t = dot(dot( linalg.inv(identity(4) - dot(ti_a, t_a)) ,  ti_a), t_b)

	s = concatenate((transpose(concatenate((transpose(s_a),transpose(s_b)))), transpose(concatenate((transpose(s_g),transpose(s_t)))) ))

	return(s)
						
def test_s2t() :
	from numpy import array,dot
	from touchstone import s2t,t2s

	left = array([
	[-0.4555-0.0072j,-0.1140+0.1402j,-0.0000-0.0000j,-0.0000-0.0000j,+0.0002-0.0010j,-0.0004+0.0008j,-0.0000-0.0000j,-0.0000-0.0000j],
	[-0.1026+0.1374j,-0.0524-0.3005j,-0.0000-0.0000j,-0.0000-0.0000j,-0.0006+0.0008j,+0.0008-0.0001j,-0.0000-0.0000j,-0.0000-0.0000j],
	[-0.0000-0.0000j,-0.0000-0.0000j,-0.4289-0.0203j,-0.1083+0.0955j,-0.0009+0.0012j,+0.0005-0.0011j,+1.0000+0.0000j,-0.0000-0.0000j],
	[-0.0000-0.0000j,-0.0000-0.0000j,-0.0997+0.0949j,+0.0323-0.3276j,+0.0009-0.0007j,-0.0004+0.0004j,-0.0000-0.0000j,+1.0000+0.0000j],
	[+0.0002-0.0010j,-0.0006+0.0007j,-0.0008+0.0012j,+0.0009-0.0007j,+0.0532-0.3280j,+0.0015-0.0362j,-0.0008+0.0012j,+0.0009-0.0007j],
	[-0.0003+0.0008j,+0.0008-0.0001j,+0.0004-0.0011j,-0.0004+0.0004j,+0.0028-0.0374j,+0.0756+0.0069j,+0.0004-0.0011j,-0.0004+0.0004j],
	[-0.0000-0.0000j,-0.0000-0.0000j,+1.0000+0.0000j,-0.0000-0.0000j,-0.0009+0.0012j,+0.0005-0.0011j,-0.4289-0.0203j,-0.1083+0.0955j],
	[-0.0000-0.0000j,-0.0000-0.0000j,-0.0000-0.0000j,+1.0000+0.0000j,+0.0009-0.0007j,-0.0004+0.0004j,-0.0997+0.0949j,+0.0323-0.3276j]])
	
	right = array([
	[0,0,0,0,1,0,0,0],
	[0,0,0,0,0,1,0,0],
	[0,0,0,0,0,0,1,0],
	[0,0,0,0,0,0,0,1],
	[1,0,0,0,0,0,0,0],
	[0,1,0,0,0,0,0,0],
	[0,0,1,0,0,0,0,0],
	[0,0,0,1,0,0,0,0]])
	
	left_t = s2t(left)
	right_t = s2t(right)
	
	answer_t = dot(left_t,right_t)
	answer = t2s(answer_t)
	
	return(answer)
	

def printM(m) :
	for i in range(m) :
		for j in range(m[0]) :
			print '%8.3,'%(real(m[i][j])),
		print '\n'

