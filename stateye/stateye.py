from numpy import *
import pdb
import time
from matplotlib import *

############################################################
# stateye class

class stateye :

	############################################################
	# constructor
	
	def __init__(self) :
	
		self.version		= 'v5_2'
		self.nomUI		= 50
		self.pws		= 0
		self.sweepdelta		= 0.01
		self.startCursor	= -8
		self.lastCursor		= 30
		self.dj			= 0.15
		self.rj			= 0.15/14.0

		print "Load state module version "+self.version
	
		# debug file
		self.debug = open('stateye.debug','w')

		# transitionState[currentState] = [<possible next state>]
		self.transitionState	= []

		# edge[currentState] = [<edge index corresponding to transitionState>]
		self.edge 		= []

		# step[<edge index>] = [<time vs. amplitude>]
		# where @t=0;a=0, @t=inf;a=final
		self.step		= []

		# length of each step must be the same and equal to rxLength
		self.UImax		= 0

		# number of states. 
		self.nStates		= 0

		# bins construction
		# bins is a Markov pdf; bin[<state>][<amplitude index>]
		self.bins		= []
		# see binIndex and binValue for explanation of bin coefficients
		self.noBins		= 4001
		self.midBin		= 2001
		self.kBin		= 2000
		self.binMax		= 2.0
		self.binaxis		= (arange(self.noBins) - self.midBin) * self.binMax / self.kBin

		# parameters for cursor to step time conversion
		# are loaded when generating step responses
		# see cursor2index for details
		self.nomOffset		= 0		# simple offset for peak centring
		# amplitude of pulse width shrinkage
		self.pws		= 0		# in UI

	def __del__(self) :
		self.debug.flush()
		self.debug.close()

	############################################################
	# simple routine to index into pdf bins given a value
	def binIndex(self, x) :
		return( round(x / self.binMax * self.kBin) + self.midBin ) 

	############################################################
	# simple route to find value for pdf bin given a index
	def binValue(self, i):
		return( 1.0 * (i-self.midBin) / self.kBin * self.binMax )

	############################################################
	# simple routine to convert cursor index and sweep offset into time index of step response
 	def cursor2index(self, cursor, sweep) :
		return( int( round( (cursor+sweep) * self.nomUI + self.nomOffset ) ) )

	############################################################
	# simple routine to convert index to cursor
	def index2cursor(self, index) :
		return( (index - self.nomOffset) * 1.0 / self.nomUI )

	############################################################
	# shift a pdf bin description as convolution
	def binShift(self, bin, x) :
		i = self.midBin - self.binIndex(x) 
		if i < 0 :
			return( concatenate(( bin[-i:], zeros(-i) )) )
		if i > 0 :
			return( concatenate(( zeros(i), bin[0:-i] )) )
		if i==0 :
			return(bin)

	# perform stateye algorithm for single sample phase
	# clearly includes pws but not mid band jitter
	def calcpdf(self) :
                from pylab import find
		import pdb

		sweepdelta 	= self.sweepdelta
		startCursor	= self.startCursor
		lastCursor	= self.lastCursor
		dj		= self.dj
		rj 		= self.rj

		self.cdf=[]
		self.sweep = arange(-1.5,1.5,sweepdelta)
		#self.sweep = [0.0]
		self.pdf=zeros(( len(self.sweep), len(self.binaxis) ))

		# generate a pre-indexed step response for acceleration of the pdf calculation
		# should consider this for other variable as well, e.g. dfecoef
		self.stepK = []
		for _step in self.step :
			_stepK = []
			for __step in _step :
				_stepK += [self.midBin - self.binIndex(__step)]
			self.stepK += [ _stepK ]

		print 'folding %d steps'%len(self.sweep)
		binstore = [ 1.0*zeros((self.nStates,self.noBins)), 1.0*zeros((self.nStates,self.noBins)) ]
		for i_sweep in range(len(self.sweep)) :
			#self.debug.writelines('at sweep %d from %d\n'%(i_sweep,len(self.sweep)))
			#delta = time.time() - tag
			print 'folding %d / %d'%(i_sweep+1,len(self.sweep))
			#tag = time.time()
			# scan from last cursor to first cursor
			# i.e. back cursor tracing 
			# for debug 
			bintag = 0
			_bins = binstore[bintag]	
			bins = binstore[1-bintag]
			bins = self.start_bins
			_bins[:] = 0.0 
			# the value of the step for the given cursor and sweep
			_sweep		= self.sweep[i_sweep]
			for cursor in flipud(range(startCursor,lastCursor)) :
				# self.debug.writelines('at cursor %d\n'%cursor)
				#print 'at cursor %d'%cursor

				# where am I on the time axis	
				currentIndex = self.cursor2index(cursor,_sweep)	
				
				# next Markov pdf contents
				# perform Markov convolution
				# scan through each state
				for state in range(self.nStates) :

					# enable for tracking speed of exection
					# print 'at state %d'%state

					# scan each possible transition from this state
					for transition in range(len(self.transitionState[state])) :
						# sweep the pws assuming a simple dirac distribution 
						# for debug 
						# for pws in [-self.pws,self.pws] :
						for pws in [0] :
							# calculate the next state for the given state
							_nextState 	= self.transitionState[state][transition]

							# the edge needed to get to this state
							_edge		= self.edge[state][transition]

							# when pws is enabled then currentIndex clearly needs to be correctly modulated
							_delta 		= -self.stepK[ _edge ][ currentIndex + pws]
							
							# perform a convoltion using a simple shift and addition for the state given
							if _delta==0 :
								_bins[_nextState]	+=  bins[state]
							else :
								if (_delta<0) :
									_t = bins[state][-_delta:]
									_bins[_nextState][:len(_t)] +=  _t
								else :
									_t = bins[state][:-_delta]
									_bins[_nextState][-len(_t):] += _t			
								

				# enable for dumping the pdfs as they are built up
				#for _state in range(self.nStates) :
				#	if len(pylab.find(bins[_state]>0) > 0) :
				#		print '%s = %s'%(self.states[_state],array2string(self.binValue(pylab.find(bins[_state]>0))))
						
				bintag = 1-bintag
				_bins = binstore[bintag]	
				bins = binstore[1-bintag]
				_bins[:] = 0.0

				# dfe condition
				# this is also taking a second
				if cursor==2 :
					for dfeCoef in self.dfeCoef :
						for state in range(self.nStates) :

							# going to make a big assumption here!!! That the threshold for greater than and less than is the same
							# also going to make a bug assumption that the gt and lt results index are also inverse
							if 1:
								_threshold = self.binIndex(self.gt_h0[state])

								_shift 		= self.binShift( \
											concatenate(( \
											zeros(_threshold), bins[state][_threshold:] )), \
											-self.gt_true[state] * dfeCoef )
								_bins[state] = add(_bins[state], _shift)
							
								_shift 		= self.binShift( \
											concatenate(( \
											bins[state][:_threshold],  zeros(self.noBins - _threshold) )), \
											-self.gt_false[state] * dfeCoef )
								_bins[state] = add(_bins[state], _shift)
							else :								
								for binIndex in range(self.noBins) :
									if self.binValue(binIndex) > self.gt_h0[state] :
										_binValue = self.binValue(binIndex) + self.gt_true[state] * dfeCoef
									else :
										_binValue = self.binValue(binIndex) + self.gt_false[state] * dfeCoef
									_bins[state][self.binIndex(_binValue)] += bins[state][binIndex]

						bintag = 1-bintag
						_bins = binstore[bintag]	
						bins = binstore[1-bintag]
						_bins[:] = 0.0

					# enable for dumping the pdfs as they are built up
					#for _state in range(self.nStates) :
					#	if len(pylab.find(bins[_state]>0) > 0) :
					#		print '%s = %s'%(self.states[_state],array2string(self.binValue(pylab.find(bins[_state]>0))))
				
			_pdf = zeros(self.noBins)

			# typical good place to break for debugging
			# pdb.set_trace()
			# firstly what if sum = 0; secondly the sum for difference states may be different???
			# this scaling of the pdf is still not quite workig correctly
			bins = bins / bins.sum()
			for state in range(self.nStates) :
				_pdf = add(_pdf, bins[state])
				# this additional of the fliped array is mainly for 8b10b support. As we only run one set of the codes
				# we need to add in the other half. I believe this is correct, but am checking it again
				_pdf = add(_pdf, flipud(bins[state]))
			_pdf = _pdf / _pdf.sum()
			self.pdf[i_sweep] = _pdf
			self.debug.writelines('final = %s\n'%(array2string(self.binValue(find(_pdf>0)))))


		print 'folding jiiter'
		# final pdf containing the jittered version
		if 1:
        		self.p = []
        		sigma = rj;
        		mean  = dj/2;
        		for _sweep in self.sweep :
        			p = 1/(sigma*sqrt(2*pi)) * exp(-((_sweep-mean)**2)/(2*sigma**2)) + \
        				1/(sigma*sqrt(2*pi)) * exp(-((_sweep+mean)**2)/(2*sigma**2)) + \
        				1/(sigma*sqrt(2*pi)) * exp(-((_sweep)**2)/(2*sigma**2));
	        		if p>1.0e-12 :
        				self.p += [p]
        			else :
        				self.p += [0.0]

        		self.p = self.p / sum(self.p);
        
        		self.pdf_pj=zeros(( len(self.sweep), len(self.binaxis) ))
        		_jmid = len(self.sweep)/2
        		for _i in range(len(self.sweep)) :
				#print 'at sweep %d'%_i
                            	for _j in range(len(self.sweep)) :
                                	_k = _i + _j - _jmid
                                	if (_k > 0) and (_k<len(self.sweep)) :
                				self.pdf_pj[_i] += self.p[_j] * self.pdf[_k]                    
				self.pdf_pj[_i] = self.pdf_pj[_i] / self.pdf_pj[_i].sum()

	def convolveNoise(self, noise) :
		print '\nfolding noise'
		# convolve the noise into
		self.pdf_n = zeros(( len(self.sweep), len(self.binaxis) ))
		for i_sweep in range(len(self.sweep)) :
			for _noise in range(len(noise_x)) :
				if noise_y[_noise] >0 :
					_shift = self.binShift(self.pdf[i_sweep], noise_x[_noise] ) * noise_y[_noise]
					self.pdf_n[i_sweep]  = add(self.pdf_n[i_sweep], _shift)
			self.pdf_n[i_sweep] = self.pdf_n[i_sweep] / (self.pdf_n[i_sweep]).sum()


	def loadStep(self, inputStep) :
		from pylab import find

		_ui 	= self.nomUI
		_pws 	= self.pws

		# define pulse response
		self.inputStep 	= inputStep
		self.converge	= max(inputStep)
		self.pulse 	= add(-inputStep[:-_ui], +inputStep[_ui:])
		self.nomOffset 	= find(self.pulse==max(self.pulse))[0]

		# start and last cursor must allow for some margin
		#self.startCursor 	= (range(self.nomOffset,0,-_ui)[-2]-self.nomOffset)/_ui		
		#self.lastCursor		= (range(self.nomOffset,len(inputStep),_ui)[-2]-self.nomOffset)/_ui
		self.xindex		= arange(self.cursor2index(self.startCursor,0),self.cursor2index(self.lastCursor,0))

		#print 'start cursor %d, finish cursor %d'%(self.startCursor,self.lastCursor)

	def create2TapFIR(self, c, noDFEtaps) :
		from pylab import find
		# input step is the fundimental step response of the system to a 1V step, and is assumed to be 0@t=0
		# c is a 1x2 array containing the FIR coefficients
		# this function defines the transitionStates, edge transitions and generates the necessary steps

		self.nStates = 4

		# for each state x
		#                             	0      		1      		2      		3 
		x 			= [ 	[0,0], 		[0,1], 		[1,0], 		[1,1] ]
		x = array(x)*2.0 - 1.0
		self.states		= [	'00',		'01',		'10',		'11'	]
		# load the possible state transition
		self.transitionState 	= [ 	[0,1], 		[2,3], 		[0,1], 		[2,3] ]
		# define the edge used to move from state to state 
		self.edge		= [ 	[0,1], 		[2,3], 		[4,5], 		[6,7] ]
		# load empty steps
		self.step 		= [     [],[],        	[],[],		[],[],		[],[] ]

		# preload the Markov pdf with the converged values for the two stable states,
		# leave the other state empty
		self.start_bins               = zeros((self.nStates,self.noBins))
		self.start_bins[0][self.binIndex( sum(x[0]*c)*self.converge )] = 1
		self.start_bins[3][self.binIndex( sum(x[3]*c)*self.converge )] = 1

		# this could be more efficient in storage of the indexes, but for a simple example it doesn't matter
		# for each state and transitions an edge is defined
		# copy the inputStep into each step array, given the correct factor needed to move states
		for _state in range(self.nStates) :
			for _transition in range(len(self.edge[_state])) :
				k = (sum(x[self.transitionState[_state][_transition]]*flipud(c)) - sum(x[_state]*flipud(c)))
				self.debug.writelines('in state %d, transisitioning to state %d, using %0.3f\n'\
						%(_state,self.transitionState[_state][_transition],k))
				self.step[self.edge[_state][_transition]] = self.inputStep * k

		# this is the post equalised step response
		self.pulse = add( add ( \
				self.step[ self.edge[0][1] ][self.nomUI*2:] , \
				self.step[ self.edge[1][0] ][self.nomUI:-self.nomUI] ) , \
				self.step[ self.edge[2][0] ][:-self.nomUI*2] ) 
		self.pulse = self.pulse/2.0
		self.nomOffset 	= find(self.pulse==max(self.pulse))[0]

		self.dfeCoef = []
		h0 = self.pulse[self.cursor2index(0 , 0)]
		print 'found h0=%0.3f'%h0
		# 			00	01	10	11
		self.gt_h0 	= [	-h0, 	+h0,	-h0,	+h0]
		self.gt_true 	= [	-1.0, 	-1.0,	-1.0,	-1.0]
		self.gt_false 	= [	+1.0, 	+1.0,	+1.0,	+1.0]
		self.lt_h0 	= [	-h0,	+h0,	-h0, 	+h0]
		self.lt_true 	= [	+1.0,	+1.0,	+1.0, 	+1.0]
		self.lt_false 	= [	-1.0,	-1.0,	-1.0, 	-1.0]

		# clearly we need to include here the proper algorithm for finding the optimum sampling point!!
		for cursor in range(noDFEtaps) :
			self.dfeCoef += [ abs( self.pulse[self.cursor2index(cursor+1, 0)] ) ]
			print 'Extracting cursor %d, found %0.3f'%(cursor+1, self.dfeCoef[-1])


	def create8b10b_2TapFIR(self,c,noDFEtaps) :
		from pylab import find
		import pdb

		words8b10b = def8b10b()
		words = sort( words8b10b )
		states = ['x','x']

		# scan through all possible 8b10b codes, truncating to a given length l
		# collect all possible codes
		for l in range(1,11) :
			for _words in words :
				short = _words[:l]
				if not(any(array(states)==short)) :
					states += [short]

		# initialise the transition state matrix
		transitionState = []
		for i in range(len(states)) :
			transitionState += [[]]
		
		# fill transition state matrx
		for i in range(len(states)) :
			# as we search for where this code could have come from, we only start searching when the
			# code would be a minimum of 2 characters long. e.g. if the code word is 1001, we search for
			# 100 as the source of this code word
			if len(states[i])>1 :
				# find the index into the state matrix, for the source of the current word 
				source = find(states[i][:-1]==array(states))[0]
				# add this code word index to the transition matrix entry for the source of this code wor
				# clearly we will only find one single source per code word
				transitionState[source] += [i]
			# if we are at the final word, then also add the transitions to this entry in the transition
			# matrix for getting back to 0 & 1. However, as we are implementing a 2 tap FIR, we must maintain 
			# also the second entry, hence the starting states are 00,01,10 & 11
			if len(states[i])==10 :
				if states[i][-1]=='0' :
					transitionState[i] += [0]
					transitionState[i] += [1]
				if states[i][-1]=='1' :
					transitionState[i] += [2]
		 			transitionState[i] += [3]

		# as stated above we must over write the first four states to be correct
		states[0:4] = ['00','01','10','11']
		transitionState[0] = transitionState[2]
		transitionState[1] = transitionState[3]

		# generate all possible transitions
		k = []
		transitionLookUp = []
		step  = []
		for _i in range(2**3) :
			if _i > 0 :
				transitionLookUp += [binary_repr(_i)]
			else :
				transitionLookUp += ['']
			while(len(transitionLookUp[-1])<3) :
				transitionLookUp[-1] = '0' + transitionLookUp[-1]
			_k = 0
			for _j in range(2) :
				# the polarity here needs to be checked
				_k -= (eval(transitionLookUp[-1][_j])*2.0-1.0) * flipud(c)[_j] - (eval(transitionLookUp[-1][_j+1])*2.0-1.0) * flipud(c)[_j]
			k += [_k]
			step += [self.inputStep * _k]
			self.debug.writelines('edge %s/%d is %0.3f\n'%(transitionLookUp[-1],_i,_k))
		
		# scan through actual transitions and enter edge index into array
		edgeText = [] 
		edge = [] 
		for _state in range(len(states)) :
			_edgeText=[]
			_edge=[]
			for _transition in range(len(transitionState[_state])) :
				_nextstate = transitionState[_state][_transition]
				__edgeText = states[_state][-2:] + states[_nextstate][-1]
				_edgeText += [__edgeText]
				__edge = find( __edgeText == array(transitionLookUp) )
				_edge += [__edge]
				self.debug.writelines('from %s to %s using %s/%d\n'%(states[_state],states[_nextstate],__edgeText,__edge))
			edge += [_edge]
			edgeText += [_edgeText]

		self.states		= states
		self.transitionState 	= transitionState
		self.edge		= edge
		self.step 		= step

		self.nStates 			= len(self.states)
		self.start_bins               	= zeros((self.nStates,self.noBins))
		# this is the current initialisation matrix which need extending
		# see the commented conditional statements below
		if 0 :
			self.start_bins[0][self.binIndex( self.converge * sum(array([-1,-1])*c) )] = 1
			self.start_bins[3][self.binIndex( self.converge * sum(array([+1,+1])*c) )] = 1
		else :
			for _states in range(len(states)) :
				if (states[_states][-2:]=='00') :
					self.start_bins[_states][self.binIndex( self.converge * sum(array([-1,-1])*c) )] = 1
				#if (states[_states][-2:]=='01') :
				#	self.start_bins[_states][self.binIndex( self.converge * sum(array([-1,+1])*c) )] = 1
				#if (states[_states][-2:]=='10') :
				#	self.start_bins[_states][self.binIndex( self.converge * sum(array([+1,-1])*c) )] = 1
				if (states[_states][-2:]=='11') :
					self.start_bins[_states][self.binIndex( self.converge * sum(array([+1,+1])*c) )] = 1
	
		# this is the post equalised step response
		self.pulse = add( add ( \
				self.step[ 1 ][self.nomUI*2:] , \
				self.step[ 2 ][self.nomUI:-self.nomUI] ) , \
				self.step[ 4 ][:-self.nomUI*2] ) 
		self.pulse = self.pulse/2.0
		self.nomOffset 	= find(self.pulse==max(self.pulse))[0]
		self.dfeCoef = []
		h0 = self.pulse[self.cursor2index(0 , 0)]

		# clearly we need to include here the proper algorithm for finding the optimum sampling point!!
		for cursor in range(noDFEtaps) :
			self.dfeCoef += [ abs( self.pulse[self.cursor2index(cursor+1, 0)] ) ]
			print 'Extracting cursor %d, found %0.3f'%(cursor+1, self.dfeCoef[-1])

		# setup the DFE correction matrix
		self.gt_h0 	= []
		self.gt_true 	= []
		self.gt_false	= []
		self.lt_h0 	= []
		self.lt_true 	= []
		self.lt_false	= []
		#pdb.set_trace()
		for _states in states :
			if (_states[-2:]=='00') or (_states[-2:]=='10') :
				self.gt_h0 	+= [-h0]
				self.gt_true	+= [-1.0]
				self.gt_false	+= [+1.0]
				self.lt_h0 	+= [-h0]
				self.lt_true	+= [+1.0]
				self.lt_false	+= [-1.0]
			if (_states[-2:]=='01') or (_states[-2:]=='11') :
				self.gt_h0 	+= [+h0]
				self.gt_true	+= [-1.0]
				self.gt_false	+= [+1.0]
				self.lt_h0 	+= [+h0]
				self.lt_true	+= [+1.0]
				self.lt_false	+= [-1.0]



	# simple example based on step.py in steptheory
	# probably doesn't work anymore since extending the code to support more features
	def bist(self) :	
		# states are 
		# 0 = 0 0
		# 1 = 0 1
		# 2 = 1 0
		# 3 = 1 1
		self.transitionState 	= [[0,1],[2,3],[0,1],[2,3]]
		self.edge            	= [[0,1],[2,3],[4,5],[6,7]]
		self.step		= [[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], \
				[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.04, 1.3600000000000001, 1.5200000000000002, 1.6000000000000001, 1.6000000000000001, 1.6000000000000001], \
				[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.1700000000000002, -1.53, -1.7100000000000002, -1.8, -1.8, -1.8], [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.12999999999999998, -0.16999999999999996, -0.18999999999999997, -0.19999999999999996, -0.19999999999999996, -0.19999999999999996], \
				[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.12999999999999998, 0.16999999999999996, 0.18999999999999997, 0.19999999999999996, 0.19999999999999996, 0.19999999999999996], \
				[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.1700000000000002, 1.53, 1.7100000000000002, 1.8, 1.8, 1.8], [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.04, -1.3600000000000001, -1.5200000000000002, -1.6000000000000001, -1.6000000000000001, -1.6000000000000001], \
				[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]]
		self.stepLength		= 15
		self.nStates		= 4
		self.start_bins		= 1.0*zeros((self.nStates,self.noBins))
		self.start_bins[0][self.binIndex(-0.7)] = 1.0
		self.start_bins[3][self.binIndex(+0.7)] = 1.0
		self.calcpdf()


############################################################
# simple functions to load all the possible 8b10b words 
############################################################
def def8b10b() :
	word8b=[
	"00000000",
	"00000001",
	"00000010",
	"00000011",
	"00000100",
	"00000101",
	"00000110",
	"00000111",
	"00001000",
	"00001001",
	"00001010",
	"00001011",
	"00001100",
	"00001101",
	"00001110",
	"00001111",
	"00010000",
	"00010001",
	"00010010",
	"00010011",
	"00010100",
	"00010101",
	"00010110",
	"00010111",
	"00011000",
	"00011001",
	"00011010",
	"00011011",
	"00011100",
	"00011101",
	"00011110",
	"00011111",
	"00100000",
	"00100001",
	"00100010",
	"00100011",
	"00100100",
	"00100101",
	"00100110",
	"00100111",
	"00101000",
	"00101001",
	"00101010",
	"00101011",
	"00101100",
	"00101101",
	"00101110",
	"00101111",
	"00110000",
	"00110001",
	"00110010",
	"00110011",
	"00110100",
	"00110101",
	"00110110",
	"00110111",
	"00111000",
	"00111001",
	"00111010",
	"00111011",
	"00111100",
	"00111101",
	"00111110",
	"00111111",
	"01000000",
	"01000001",
	"01000010",
	"01000011",
	"01000100",
	"01000101",
	"01000110",
	"01000111",
	"01001000",
	"01001001",
	"01001010",
	"01001011",
	"01001100",
	"01001101",
	"01001110",
	"01001111",
	"01010000",
	"01010001",
	"01010010",
	"01010011",
	"01010100",
	"01010101",
	"01010110",
	"01010111",
	"01011000",
	"01011001",
	"01011010",
	"01011011",
	"01011100",
	"01011101",
	"01011110",
	"01011111",
	"01100000",
	"01100001",
	"01100010",
	"01100011",
	"01100100",
	"01100101",
	"01100110",
	"01100111",
	"01101000",
	"01101001",
	"01101010",
	"01101011",
	"01101100",
	"01101101",
	"01101110",
	"01101111",
	"01110000",
	"01110001",
	"01110010",
	"01110011",
	"01110100",
	"01110101",
	"01110110",
	"01110111",
	"01111000",
	"01111001",
	"01111010",
	"01111011",
	"01111100",
	"01111101",
	"01111110",
	"01111111",
	"10000000",
	"10000001",
	"10000010",
	"10000011",
	"10000100",
	"10000101",
	"10000110",
	"10000111",
	"10001000",
	"10001001",
	"10001010",
	"10001011",
	"10001100",
	"10001101",
	"10001110",
	"10001111",
	"10010000",
	"10010001",
	"10010010",
	"10010011",
	"10010100",
	"10010101",
	"10010110",
	"10010111",
	"10011000",
	"10011001",
	"10011010",
	"10011011",
	"10011100",
	"10011101",
	"10011110",
	"10011111",
	"10100000",
	"10100001",
	"10100010",
	"10100011",
	"10100100",
	"10100101",
	"10100110",
	"10100111",
	"10101000",
	"10101001",
	"10101010",
	"10101011",
	"10101100",
	"10101101",
	"10101110",
	"10101111",
	"10110000",
	"10110001",
	"10110010",
	"10110011",
	"10110100",
	"10110101",
	"10110110",
	"10110111",
	"10111000",
	"10111001",
	"10111010",
	"10111011",
	"10111100",
	"10111101",
	"10111110",
	"10111111",
	"11000000",
	"11000001",
	"11000010",
	"11000011",
	"11000100",
	"11000101",
	"11000110",
	"11000111",
	"11001000",
	"11001001",
	"11001010",
	"11001011",
	"11001100",
	"11001101",
	"11001110",
	"11001111",
	"11010000",
	"11010001",
	"11010010",
	"11010011",
	"11010100",
	"11010101",
	"11010110",
	"11010111",
	"11011000",
	"11011001",
	"11011010",
	"11011011",
	"11011100",
	"11011101",
	"11011110",
	"11011111",
	"11100000",
	"11100001",
	"11100010",
	"11100011",
	"11100100",
	"11100101",
	"11100110",
	"11100111",
	"11101000",
	"11101001",
	"11101010",
	"11101011",
	"11101100",
	"11101101",
	"11101110",
	"11101111",
	"11110000",
	"11110001",
	"11110010",
	"11110011",
	"11110100",
	"11110101",
	"11110110",
	"11110111",
	"11111000",
	"11111001",
	"11111010",
	"11111011",
	"11111100",
	"11111101",
	"11111110",
	"11111111"
	]
	
	word10b_p=[
	"1001110100",
	"0111010100",
	"1011010100",
	"1100011011",
	"1101010100",
	"1010011011",
	"0110011011",
	"1110001011",
	"1110010100",
	"1001011011",
	"0101011011",
	"1101001011",
	"0011011011",
	"1011001011",
	"0111001011",
	"0101110100",
	"0110110100",
	"1000111011",
	"0100111011",
	"1100101011",
	"0010111011",
	"1010101011",
	"0110101011",
	"1110100100",
	"1100110100",
	"1001101011",
	"0101101011",
	"1101100100",
	"0011101011",
	"1011100100",
	"0111100100",
	"1010110100",
	"1001111001",
	"0111011001",
	"1011011001",
	"1100011001",
	"1101011001",
	"1010011001",
	"0110011001",
	"1110001001",
	"1110011001",
	"1001011001",
	"0101011001",
	"1101001001",
	"0011011001",
	"1011001001",
	"0111001001",
	"0101111001",
	"0110111001",
	"1000111001",
	"0100111001",
	"1100101001",
	"0010111001",
	"1010101001",
	"0110101001",
	"1110101001",
	"1100111001",
	"1001101001",
	"0101101001",
	"1101101001",
	"0011101001",
	"1011101001",
	"0111101001",
	"1010111001",
	"1001110101",
	"0111010101",
	"1011010101",
	"1100010101",
	"1101010101",
	"1010010101",
	"0110010101",
	"1110000101",
	"1110010101",
	"1001010101",
	"0101010101",
	"1101000101",
	"0011010101",
	"1011000101",
	"0111000101",
	"0101110101",
	"0110110101",
	"1000110101",
	"0100110101",
	"1100100101",
	"0010110101",
	"1010100101",
	"0110100101",
	"1110100101",
	"1100110101",
	"1001100101",
	"0101100101",
	"1101100101",
	"0011100101",
	"1011100101",
	"0111100101",
	"1010110101",
	"1001110011",
	"0111010011",
	"1011010011",
	"1100011100",
	"1101010011",
	"1010011100",
	"0110011100",
	"1110001100",
	"1110010011",
	"1001011100",
	"0101011100",
	"1101001100",
	"0011011100",
	"1011001100",
	"0111001100",
	"0101110011",
	"0110110011",
	"1000111100",
	"0100111100",
	"1100101100",
	"0010111100",
	"1010101100",
	"0110101100",
	"1110100011",
	"1100110011",
	"1001101100",
	"0101101100",
	"1101100011",
	"0011101100",
	"1011100011",
	"0111100011",
	"1010110011",
	"1001110010",
	"0111010010",
	"1011010010",
	"1100011101",
	"1101010010",
	"1010011101",
	"0110011101",
	"1110001101",
	"1110010010",
	"1001011101",
	"0101011101",
	"1101001101",
	"0011011101",
	"1011001101",
	"0111001101",
	"0101110010",
	"0110110010",
	"1000111101",
	"0100111101",
	"1100101101",
	"0010111101",
	"1010101101",
	"0110101101",
	"1110100010",
	"1100110010",
	"1001101101",
	"0101101101",
	"1101100010",
	"0011101101",
	"1011100010",
	"0111100010",
	"1010110010",
	"1001111010",
	"0111011010",
	"1011011010",
	"1100011010",
	"1101011010",
	"1010011010",
	"0110011010",
	"1110001010",
	"1110011010",
	"1001011010",
	"0101011010",
	"1101001010",
	"0011011010",
	"1011001010",
	"0111001010",
	"0101111010",
	"0110111010",
	"1000111010",
	"0100111010",
	"1100101010",
	"0010111010",
	"1010101010",
	"0110101010",
	"1110101010",
	"1100111010",
	"1001101010",
	"0101101010",
	"1101101010",
	"0011101010",
	"1011101010",
	"0111101010",
	"1010111010",
	"1001110110",
	"0111010110",
	"1011010110",
	"1100010110",
	"1101010110",
	"1010010110",
	"0110010110",
	"1110000110",
	"1110010110",
	"1001010110",
	"0101010110",
	"1101000110",
	"0011010110",
	"1011000110",
	"0111000110",
	"0101110110",
	"0110110110",
	"1000110110",
	"0100110110",
	"1100100110",
	"0010110110",
	"1010100110",
	"0110100110",
	"1110100110",
	"1100110110",
	"1001100110",
	"0101100110",
	"1101100110",
	"0011100110",
	"1011100110",
	"0111100110",
	"1010110110",
	"1001110001",
	"0111010001",
	"1011010001",
	"1100011110",
	"1101010001",
	"1010011110",
	"0110011110",
	"1110001110",
	"1110010001",
	"1001011110",
	"0101011110",
	"1101001110",
	"0011011110",
	"1011001110",
	"0111001110",
	"0101110001",
	"0110110001",
	"1000110111",
	"0100110111",
	"1100101110",
	"0010110111",
	"1010101110",
	"0110101110",
	"1110100001",
	"1100110001",
	"1001101110",
	"0101101110",
	"1101100001",
	"0011101110",
	"1011100001",
	"0111100001",
	"1010110001"]
	
	word10b_n=[
	"0110001011",
	"1000101011",
	"0100101011",
	"1100010100",
	"0010101011",
	"1010010100",
	"0110010100",
	"0001110100",
	"0001101011",
	"1001010100",
	"0101010100",
	"1101000100",
	"0011010100",
	"1011000100",
	"0111000100",
	"1010001011",
	"1001001011",
	"1000110100",
	"0100110100",
	"1100100100",
	"0010110100",
	"1010100100",
	"0110100100",
	"0001011011",
	"0011001011",
	"1001100100",
	"0101100100",
	"0010011011",
	"0011100100",
	"0100011011",
	"1000011011",
	"0101001011",
	"0110001001",
	"1000101001",
	"0100101001",
	"1100011001",
	"0010101001",
	"1010011001",
	"0110011001",
	"0001111001",
	"0001101001",
	"1001011001",
	"0101011001",
	"1101001001",
	"0011011001",
	"1011001001",
	"0111001001",
	"1010001001",
	"1001001001",
	"1000111001",
	"0100111001",
	"1100101001",
	"0010111001",
	"1010101001",
	"0110101001",
	"0001011001",
	"0011001001",
	"1001101001",
	"0101101001",
	"0010011001",
	"0011101001",
	"0100011001",
	"1000011001",
	"0101001001",
	"0110000101",
	"1000100101",
	"0100100101",
	"1100010101",
	"0010100101",
	"1010010101",
	"0110010101",
	"0001110101",
	"0001100101",
	"1001010101",
	"0101010101",
	"1101000101",
	"0011010101",
	"1011000101",
	"0111000101",
	"1010000101",
	"1001000101",
	"1000110101",
	"0100110101",
	"1100100101",
	"0010110101",
	"1010100101",
	"0110100101",
	"0001010101",
	"0011000101",
	"1001100101",
	"0101100101",
	"0010010101",
	"0011100101",
	"0100010101",
	"1000010101",
	"0101000101",
	"0110001100",
	"1000101100",
	"0100101100",
	"1100010011",
	"0010101100",
	"1010010011",
	"0110010011",
	"0001110011",
	"0001101100",
	"1001010011",
	"0101010011",
	"1101000011",
	"0011010011",
	"1011000011",
	"0111000011",
	"1010001100",
	"1001001100",
	"1000110011",
	"0100110011",
	"1100100011",
	"0010110011",
	"1010100011",
	"0110100011",
	"0001011100",
	"0011001100",
	"1001100011",
	"0101100011",
	"0010011100",
	"0011100011",
	"0100011100",
	"1000011100",
	"0101001100",
	"0110001101",
	"1000101101",
	"0100101101",
	"1100010010",
	"0010101101",
	"1010010010",
	"0110010010",
	"0001110010",
	"0001101101",
	"1001010010",
	"0101010010",
	"1101000010",
	"0011010010",
	"1011000010",
	"0111000010",
	"1010001101",
	"1001001101",
	"1000110010",
	"0100110010",
	"1100100010",
	"0010110010",
	"1010100010",
	"0110100010",
	"0001011101",
	"0011001101",
	"1001100010",
	"0101100010",
	"0010011101",
	"0011100010",
	"0100011101",
	"1000011101",
	"0101001101",
	"0110001010",
	"1000101010",
	"0100101010",
	"1100011010",
	"0010101010",
	"1010011010",
	"0110011010",
	"0001111010",
	"0001101010",
	"1001011010",
	"0101011010",
	"1101001010",
	"0011011010",
	"1011001010",
	"0111001010",
	"1010001010",
	"1001001010",
	"1000111010",
	"0100111010",
	"1100101010",
	"0010111010",
	"1010101010",
	"0110101010",
	"0001011010",
	"0011001010",
	"1001101010",
	"0101101010",
	"0010011010",
	"0011101010",
	"0100011010",
	"1000011010",
	"0101001010",
	"0110000110",
	"1000100110",
	"0100100110",
	"1100010110",
	"0010100110",
	"1010010110",
	"0110010110",
	"0001110110",
	"0001100110",
	"1001010110",
	"0101010110",
	"1101000110",
	"0011010110",
	"1011000110",
	"0111000110",
	"1010000110",
	"1001000110",
	"1000110110",
	"0100110110",
	"1100100110",
	"0010110110",
	"1010100110",
	"0110100110",
	"0001010110",
	"0011000110",
	"1001100110",
	"0101100110",
	"0010010110",
	"0011100110",
	"0100010110",
	"1000010110",
	"0101000110",
	"0110001110",
	"1000101110",
	"0100101110",
	"1100010001",
	"0010101110",
	"1010010001",
	"0110010001",
	"0001110001",
	"0001101110",
	"1001010001",
	"0101010001",
	"1101001000",
	"0011010001",
	"1011001000",
	"0111001000",
	"1010001110",
	"1001001110",
	"1000110001",
	"0100110001",
	"1100100001",
	"0010110001",
	"1010100001",
	"0110100001",
	"0001011110",
	"0011001110",
	"1001100001",
	"0101100001",
	"0010011110",
	"0011100001",
	"0100011110",
	"1000011110",
	"0101001110"]

	codewords=[	
	"0011110100",
	"1100001011",
	"0011111001",
	"1100000110",
	"0011110101",
	"1100001010",
	"0011110011",
	"1100001100",
	"0011110010",
	"1100001101",
	"0011111010",
	"1100000101",
	"0011110110",
	"1100001001",
	"0011111000",
	"1100000111",
	"1110101000",
	"0001010111",
	"1101101000",
	"0010010111",
	"1011101000",
	"0100010111",
	"0111101000",
	"1000010111"]

	   
	return(concatenate(( word10b_p, word10b_n, codewords )))																	

