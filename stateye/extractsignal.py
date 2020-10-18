from numpy import histogram, arange, array, diff, cumsum, zeros, dot, flipud, concatenate
from pylab import figure, plot, find, grid , title, axis, hold, savefig
from string import rsplit, rstrip
from scipy import linalg, interp
from cdr import cdr, calcJitter
from re import *
import time 
import pdb

class extractsignal :

	def __init__(self) :
		self.version		= 'v5_2'
		self.filename 		= ''
		self.timecol  		= 0
		self.sigcol		= 1
		self.timestep		= 0.0
		self.samplePerBit	= 0		# should be set from stateye
		self.resolution 	= 0.0		# can be set as we don't know baudrate
		self.taplength		= 100 
		self.alpha		= 60.0
		self.beta		= 0.1
		self.yscale		= 1.0
		self.length		= 100000
		self.RJ			= []
		self.k			= 1e-3

		print "Load extractsignal module version "+self.version
		

	# from a measurement scan 
	def findtaps(self) :
		filename 	= self.filename 			
		timecol  	= self.timecol  			
		sigcol		= self.sigcol			
		timestep	= self.timestep			
		samplePerBit	= self.samplePerBit		
		resolution 	= self.resolution 		
		taplength	= self.taplength * self.samplePerBit
		alpha		= self.alpha			
		beta		= self.beta			
		length 		= self.length	 			# length of sample used to extract signal

		if (length-taplength)<1000 :
			print "warning in extract. Looks like the length of extract record vs. taplength is short."
	
		# read time and signal from file
		original_time = []
		original_signal = []
		foundBeginning = 0
		for line in open(filename) :
			if not foundBeginning :
				if search('^[0-9]',line) :
					foundBeginning = 1
				if match('Y Units:MicroVolt',line) :						
					self.yscale = 1.0e-6
			if foundBeginning :
				col = split('[,\s]',line)
				if timecol != [] :
					original_time += [eval(col[timecol])]
				original_signal += [eval(col[sigcol]) * self.yscale]
		
		# if no time column defined then use the timestep parameter to create the time axis
		if timecol==[] :
			original_time = arange(len(original_signal)) * timestep 
		
		# re-interpolate the time and signal according to the target resolution
		time = arange(original_time[0], original_time[-1], resolution)
		signal = interp(time,original_time,original_signal)
		
		# enhance the signal, to extract the underlying data 
		# this will probably not work when the channel is extremely reflective
		signal_hpf = filter(signal,alpha,beta) 
		source = array(signal_hpf)>0
		source = source * 2.0 - 1.0
	
		# extract the edges of this raw data stream and pass through CDR
		# afterwards, use the extracted clock to recreate a reasonable clock source
		# for the data
		j = (extractApproxEdge(source) ) * resolution * 6.0e9
		[pe, per] = cdr(j, 0.0125,0.0125,'test cdr')
		edges = ( cumsum(per) / resolution / 6.0e9 )
		source_sampled = zeros( int(round(edges[-1])) )
		polarity = -1.0
		for i in range(1,len(edges)) :
			source_sampled[ int(round(edges[i-1])) : int(round(edges[i])) ] = polarity
			polarity = -polarity
		source_sampled = source_sampled[ int(round(edges[0])): ]
		signal = signal[:len(source_sampled)]
	
		# use a small record set in order to extract the channel characteristic tap settings
		if len(source_sampled) < length :
			print 'record too short, would recommend at least %d samples'%length
		source_small = source_sampled[-length:]
		signal_small = signal[-length:]
		signal_small = signal_small[taplength*3/4:]
		
		# initialise taps and run interative weighted solution
		myf = figure()
		plot(signal_small)
		hold(1)
		taps = zeros(taplength) * 0.0
		for interative in range(1,21) :
			print 'Iteration %d'%interative
			tapcorrection = zeros(taplength) * 0.0
			newsignal = []
			for index in range(len(source_small)-taplength) :
				_newsignal = dot(source_small[index:index+taplength] , taps)
				signalerror =  (signal_small[index] - _newsignal )*self.k
				tapcorrection += signalerror * source_small[index:index+taplength]
				newsignal += [_newsignal]
			taps = taps + tapcorrection / len(source_small)
			if 1 : 
				plot(newsignal)

		axis([0,16000,axis()[2],axis()[3]])
		savefig('x.png')

		newsignal = []
		for index in range(len(source_small)-taplength) :
			newsignal += [ sum(source_small[index:index+taplength]*taps) ]
	
		#figure()
		#plot(newsignal)
		#plot(signal_small)
	
		# pdb.set_trace()
	
		# generate complete record
		# should start at least 500UI from beginning
	
		length2 = int(( len(source_sampled) / samplePerBit - 500 ) * samplePerBit)
	
		# pdb.set_trace()
	
		source_large = source_sampled[-length2:]
	
		signal_large = signal[-length2:]
		signal_large = signal_large[taplength*3/4:]
	
		signal_hpf_large = signal_hpf[-length2:]
		signal_hpf_large = signal_hpf_large[taplength*3/4:]
	
		print 'creating complete signal'
		newsignal = []
		for index in range(len(source_large)-taplength) :
			newsignal += [ dot(source_large[index:index+taplength],taps) ]
		
		print 'creating hpf'
		newsignal_hpf 	= filter(newsignal,alpha,beta) 
		signal_hpf_large = filter(signal_large,alpha,beta)
	
		print 'extracting jitter'
		# calculate the edge jitter 
		jitter_source		= array( extractApproxEdge(source_large) )
		jitter_signal 		= array( extractAccurateEdge(signal_hpf_large) )
		jitter_newsignal 	= array( extractAccurateEdge(newsignal_hpf) )
	
		finish = min([ len(jitter_source),len(jitter_newsignal), len(jitter_signal) ])
	
		# setting this up to optimally extract the jitter
		# this is not optimal for noise, as start must be defined to allow the initial CDR calculation
		# of the source to have settled
		jitter = jitter_source[:finish] - (jitter_newsignal[:finish] - jitter_signal[:finish])
		noise  = array(signal_large[:len(newsignal)]) - array(newsignal)
		jitter = jitter * resolution * 6.0e9
	
		#pdb.set_trace()
	
		print 'performing cdr'
		[pe,per]=cdr(jitter,0.125,0.125,'CDR Jitter Extraction')
		#pdb.set_trace()
	
		if 1 :
			figure()
			plot(newsignal)
			plot(signal_large)
			plot(noise)
			ll = len(newsignal)
			grid(1)
			axis([ll-1000,ll,-1,1])
		if 1:
			figure()
			time = arange(min([len(newsignal_hpf),len(signal_hpf_large)])) / samplePerBit
			ll = len(time)
			plot(signal_hpf_large[:ll])
			plot(newsignal_hpf[:ll])
			title('Time domain used to extract jitter')
			grid(1)
			axis([ll-1000,ll,-1,1])
	
		[self.RJ,self.DJ]	= calcJitter(pe, self.RJ)
		self.step 	= cumsum(flipud(concatenate(( zeros(len(taps)), taps, zeros(len(taps)) )) ))

	
# extract only the index values for the edges
def extractApproxEdge(x) :
	y = find( abs(diff( (array(x)>0)*1.0 )) == 1.0 )
	return(y)

# perform an accurate interpolation of the edge position
def extractAccurateEdge(x) :
	from scipy import interp
	y = []
	for i in range(len(x)-1) :
		if ( ((x[i] < 0.0) and (x[i+1] > 0.0)) or ((x[i] > 0.0) and (x[i+1] < 0.0)) ) :
			_y = interp([0],[x[i],x[i+1]],[i,i+1])[0]
			y += [_y]
	return(y)

def filter(x,alpha,beta) :
	hpf = [x[0]]
	lpf = [x[0]]
	for i in range(1,len(x)) :
		hpf += [ x[i-1] + alpha*(x[i]-x[i-1]) ]
		lpf += [ hpf[-1] * beta + (1-beta) * lpf[-1] ]

	y = lpf
	return(y)


