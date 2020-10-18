from numpy import *
from pylab import *
from scipy.special import erfinv
import pdb

def cdr_version() :
	return('v5_22')
	
def cdr (edges,k,m,name) :

	print "Running cdr version " + cdr_version()
	
	period0 = 1.0

	period = [period0]
	
	phase = [edges[0]] 
	phaseError = []
	nperiod = []
	dataperiod = []

	phaseInOld = edges[0]-period[-1]
	
	dataperiod = []

	for phaseIn in edges :
		nperiod += [ floor( (phaseIn+0.5*period[-1] - phaseInOld) / period[-1]) ]
		phaseInOld = phaseIn
		_phaseError = phaseIn - phase[-1] + period[-1]/2 
		phaseError += [ mod( _phaseError , period[-1]) - period[-1]/2 ]
		period += [period[-1] + phaseError[-1] * k]
		phase += [phase[-1] + phaseError[-1] * m + nperiod[-1]*period[-1]]

		dataperiod += [ period[-1] * nperiod[-1] ]
	

	figure()
	subplot(3,1,1)
	hold(0)
	plot(phase/mean(period))
	hold(1)
	plot(edges/mean(period))
	grid(1)
	xlabel('time [UI]')
	ylabel('Absolute Phase[UI]')
	title(name)
	
	subplot(3,1,2)
	plot(diff(edges))
	hold(1)
	plot(array(period))
	grid(1)
	xlabel('time [UI]')
	ylabel('Period\nDeviation [%mean]')
	
	subplot(3,1,3)
	plot(array(phaseError))
	grid(1)
	xlabel('time [UI]')
	ylabel('Phase Error[UI]')

	#savefig('cdrExtraction.png')
	
	return([phaseError, dataperiod])
	

def calcJitter(pe, RJ) :

	[pdf,t]=histogram(array(pe[300:]),100)

	pdf[find(pdf<5)] = 0 
	pdf = pdf*1.0 / sum(pdf)
	mid = min(find(t>0))
	left = max(find(pdf[:mid] == 0)) + 1
	rght = min(find(pdf[mid:] == 0)) - 1 + mid 
	
	leftcdf = cumsum(pdf[left:mid])
	leftcdf[find(leftcdf==1)] = 1-1e-15
	rghtcdf = flipud(cumsum(flipud(pdf[mid:rght])))
	rghtcdf[find(rghtcdf==1)] = 1-1e-15

	leftt = t[left:mid]
	rghtt = t[mid:rght]
	
	Qleft = -sqrt(2) * erfinv( 2.0 * (1 - leftcdf) -1 )
	Qrght = -sqrt(2) * erfinv( 2.0 * (1 - rghtcdf) -1 )
	
	npoints = 10
	Pleft = polyfit(leftt[0:npoints],Qleft[0:npoints],1)
	Prght = polyfit(rghtt[-npoints:],Qrght[-npoints:],1)
	
	#_Qleft = concatenate(( [-7] , Qleft ))
	#_Qrght = concatenate(( Qrght , [-7] ))
	
	if RJ == [] :
		RJ = ( 1.0 / abs(Pleft[0]) + 1.0 / abs(Prght[0]) ) / 2.0

	DJ = polyfit(Qrght[-npoints:], rghtt[-npoints:], 1)[1] - polyfit(Qleft[0:npoints], leftt[0:npoints], 1)[1]

	# pdb.set_trace()

	#_leftt = concatenate(( [ leftt[0] - (Qleft[0]+7)*RJ ] , leftt ))
	#_rghtt = concatenate(( rghtt , [ rghtt[-1] + (Qrght[-1]+7)*RJ ] ))

	print 'Extracted RJ = %0.4f, DJ = %0.4f'%(RJ, DJ)

	figure()
	hold(0)
	plot(leftt,Qleft)
	hold(1)
	plot(rghtt,Qrght)
	grid(1)
	xlabel('Time [UI]')
	ylabel('Q')
	title('Extracted Transmit Jitter, RJ = %0.4f, DJ = %0.4f'%(RJ, DJ) )
	#savefig('ExtractedJitter.png')
	
	return([RJ,DJ])	
	
