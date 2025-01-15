"""SDFA_KLx.py

implementing the approximate KL algorithm for stochastic
deterministic finite acceptors (SDFAs), from Cortes&al'08, Mohri'02.

Given SDFAs m1 and m2, klx(m1,m2,e) returns their relative entropy,
sometimes called Kullback-Liebler (KL) divergence, up to
error threshold epsilon e.

To run a simple example, type:
  > python SDFA_KLx.py
Other examples, displays of the relevant machines, etc can
be uncommented. Set VERBOSE=true for lots of info.

The Cortes&al calculations use weights in the "entropy semiring" K,
sometimes also called the "expectation semiring".  K weights are pairs
of floating point numbers representing (probability,entropy). 
K also includes +infty,-infty, but here we just let infinite values
raise exceptions.  The semiring zero is (0.,0.), and one is (1.,0.).

This approximate algorithm uses Bellman-Ford breadth-first queue
updating to a threshold epsilon. Here Python deques implement queues:
enqueue e with appendleft(e) and dequeue from the right with pop().

Cortes&al 2008.   https://cs.nyu.edu/~mohri/pub/nkl.pdf
                  https://doi.org/10.1142/s0129054108005644
Mohri 2002.       https://cs.nyu.edu/~mohri/pub/jalc.pdf
                  https://doi.org/10.25596/jalc-2002-321
Allauzen&al 2007. https://doi.org/10.1007/978-3-540-76336-9_3

############

Each weighted DFA is represented by a 3-tuple:

( InitialState,
  {FinalState -> Weight},
  {(SourceState,Symbol) -> (TargetState,Weight)} )

where symbols are strings; states are strings or ints; weights are
    sometimes probabilities (i.e. floats 0<=p<=1),
    sometimes log2(probabilities) (i.e. floats),
    sometimes entropic weights (i.e. pairs of floats)

Since this structure provides exactly one initial state, since
transitions are given by a function (i.e. a dict), and since every
transition consumes a symbol (i.e. no epsilon transitions), this
notation defines only DFAs.

Cortes&al'08 treat the larger class of unambiguous finite
acceptors, but here we treat the slightly simpler class of DFAs.

Functions with names beginning with 'k' are defined specifically for
entropic weights and the DFAs that use them.

Additional requirements:
  * probabilities in the input SDFA should be normalized -- i.e., for
    every state, the sum of the probabilities of outgoing arcs and
    final weight should sum to one
  * any cycle should have probability strictly less than one
  * the two input SDFAs m1 and m2 should have the same alphabet.
    In fact, if there is a string that only one of the input SDFAs
    maps to a positive probability, the relative entropy of m1,m2 is
    infinite, and so the algorithm defined here will throw an
    exception.

"""
from numpy import log2, zeros, array, ndarray
from itertools import product # cartesian product of lists
from collections import deque # to implement queue
from SDFA_display import * # optional, for graphical displays
from SDFA_YuExamples import mIDS, mADS

E = 0.0001 # default epsilon for testing klx

VERBOSE = True
if VERBOSE: print('VERBOSE is set to True')

""" EXAMPLE SDFAs based on the machine in Figure 1 of Mohri 2002
       (initial state, final weights, transitions)
    where transitions is a dict mapping
       (state,input) -> (state,probability)
"""
fig1a = (0, {1:0.1}, { (0,'a'):(1,1.), (1,'b'):(1,0.9) })
fig1b = (0, {1:0.9}, { (0,'a'):(1,1.), (1,'b'):(1,0.1) })

""" Basic functions on DFAs
"""
def states(m):
  """ return the states of m """
  (init, finals, transitions) = m
  qs = set(finals.keys())
  qs.add(init)
  for ((qi,_), (qj,_)) in transitions.items():
    qs.add(qi)
    qs.add(qj)
  return list(qs)

def vocabulary(m):
  """ return the vocabulary of m """
  (init, finals, transitions) = m
  vocabulary = set([])
  for ((_,ii), (_,_)) in transitions.items():
    vocabulary.add(ii)
  return list(vocabulary)

def finalize(m):
  """ Add new final symbol and state to guarantee unique final with
      weight 1.0, to simplify entropy calculations
  """
  if '>' in vocabulary(m):
    raise RuntimeError('finalize: final symbol name ">" conflict')
  if '$' in states(m):
    raise RuntimeError('finalize: final state name "$" conflict')
  transitions = dict(m[2]) # copy the transition function of m
  for q in m[1].keys():    # add transitions to new final state
    if m[1][q] != 1:
      transitions[(q,'>')] = ('$', m[1][q])
    else:
      transitions[(q,'>')] = ('$', 1.)
  return(m[0], {'$':1.0}, transitions)

def logDFA(m):
  """ return DFA like m but with log2(probability) weights """
  finals = {}
  for q in m[1].keys():
    if m[1][q] == 0.: finals[q] = 0.
    else: finals[q] = log2(m[1][q])
  transitions = {}
  for si in m[2].keys():
    if m[2][si][1] == 0.:
      transitions[si] = (m[2][si][0],0.)
    else:
      transitions[si] = (m[2][si][0],log2(m[2][si][1]))
  return (m[0], finals, transitions)

def phi1(m):
  """ return k-DFA like m but with entropic weights (w,0.) """
  return( (m[0],
    dict([(q, (m[1][q],0.)) for q in m[1].keys()]),
    dict([(si, (m[2][si][0], (m[2][si][1],0.))) for si in m[2].keys()]) ) )

def phi2(m):
  """ return k-DFA like m but with entropic weights (1.,w) """
  return( (m[0],
    dict([(q, (1.,m[1][q])) for q in m[1].keys()]),
    dict([(si, (m[2][si][0], (1.,m[2][si][1]))) for si in m[2].keys()]) ) )

""" Operations on the entropy semiring
"""

def kplus(w1,w2):
  """ entropy weight sum: sum is coordinate-wise """
  (x1,y1) = w1
  (x2,y2) = w2
  return ( x1+x2, y1+y2 )

def ktimes(w1,w2):
  """ entropy weight product """
  (x1,y1) = w1
  (x2,y2) = w2
  return ( x1*x2, (x1*y2) + (x2*y1) )

""" Operations on K semiring DFAs
"""

def kcomplete(m):
  """ complete the transitions of entropic m *in place*
      by adding transitions to a new sink state 'x' if necessary
  """
  qs = states(m)
  v = vocabulary(m)
  newTransitions = False
  for q in qs:
    defined = [ii for ((qi,ii), (qj,pi)) in m[2].items() if qi==q]
    for a in v:
      if not(a in defined):
        m[2][(q,a)] = ('x',(0.,0.))
        newTransitions = True
  if newTransitions:
    for a in v:
        m[2][('x',a)] = ('x',(0.,0.))

def kintersect(k1,k2):
  """ return intersection of entropic DFAs k1, k2 """
  init = (k1[0], k2[0])
  finals = {}
  for q1 in k1[1].keys():
    for q2 in k2[1].keys():
      finals[(q1,q2)] = ktimes(k1[1][q1],k2[1][q2])
  transitions = {}
  for ((qi1,ii1),(qj1,w1)) in k1[2].items():
    for ((qi2,ii2),(qj2,w2)) in k2[2].items():
      if ii1 == ii2:
        transitions[((qi1,qi2),ii1)] = ((qj1,qj2), ktimes(w1,w2))
  return ( (init, finals, transitions ) )

def kbuild(m1,m2):
  """ return Cortes&al entropic intersection machines and their states """
  if VERBOSE:
    print('\n--- k1 = kcomplete(kintersect( phi1(m1), phi2(logDFA(m1)) ) )')
  k1a = phi1(finalize(m1))
  #dotDFA(k1a,states(k1a))
  k1b = phi2(logDFA(finalize(m1)))
  k1 = kintersect(k1a,k1b)
  kcomplete(k1)
  qs1 = states(k1) # matrix indices are positions in this list
  #dotDFA(k1,qs1)
  if VERBOSE:
    print('--- k2 = kcomplete(kintersect( phi1(m1), phi2(logDFA(m2)) ) )\n')
  k2b = phi2(logDFA(finalize(m2)))
  k2 = kintersect(k1a,k2b)
  kcomplete(k2)
  qs2 = states(k2) # matrix indices are positions in this list
  return (k1, qs1, k2, qs2)

def ks_approx(m, qs, e):
  """ approximate distance s to final state, up to epsilon e """
  # initialize vectors for m
  n = len(qs)
  d = zeros((n,2))
  r = zeros((n,2))
  initial,final = (qs.index(m[0]), qs.index(list(m[1].keys())[0]))
  d[initial] = (1.,0.)
  r[initial] = (1.,0.)
  # convert transitions to {sourceState -> [(target,weight)] list}
  edges = {}
  for i in range(n): edges[i] = []
  for t in m[2].items():
    ((qi,ii),(qf,w)) = t
    i,f = qs.index(qi),qs.index(qf)
    edges[i] += [(f,w)]
  if VERBOSE:
    if n < 40: print('states =',qs)
    print('initial,final state indices =',(initial,final))
    print('initialized vectors d and r identically with %d weights' % n)
    if n < 20: kshowMx([d])
  # update vectors with distances, in place
  if VERBOSE: print('updating in place, up to epsilon %.2f...' % e)
  s = deque([initial]) # to start: initial state in queue
  iteration = 0
  while s:
    q = s.pop()
    r1 = r[q].copy()
    r[q] = (0.,0.)
    for t in edges[q]:
      (f,w) = t
      newf = kplus( d[f], ktimes(r1, w) )
      changef = abs(newf[0]-d[f][0]) + abs(newf[1]-d[f][1])
      if changef > e:
        d[f] = newf
        r[f] = kplus(r[f], ktimes(r1, w))
        if not(f in s): s.appendleft(f)
    iteration += 1
  if VERBOSE:
    print('%d iterations up to epsilon %f' % (iteration,e))
    if n < 20: kshowMx([d])
    print('distance s = d[%d] = %f,%f\n' % (final,d[final][0],d[final][1]))
  return d[final][1]

def klx(m1, m2, e):
  """ approximate relative entropy of SDFAs m1,m2 up to error e """
  (k1,qs1,k2,qs2) = kbuild(m1,m2)
  #if len(qs1) < 40: dotDFA(k1,qs1)
  if VERBOSE: print('calculate approx distances in k1...')
  s1 = ks_approx(k1, qs1, e)
  if VERBOSE: print('calculate approx distances in k2...')
  s2 = ks_approx(k2, qs2, e) 
  change = s1-s2
  if VERBOSE: print('s1 - s2 = %f bits\n' % change)
  return abs(change)

""" line mode display functions
"""

def printTransition(t):
  ((qi,ii),(qf,w)) = t
  if (isinstance(w,tuple) or isinstance(w,ndarray)) and len(w) == 2:
    print('   (%s,%s) -> (%s,(%s,%s))' % (qi,ii,qf,'%.2f' % w[0],'%.2f' % w[1]))
  elif isinstance(w,float):
    print('   (%s,%s) -> (%s,%s)' % (qi,ii,qf,'%.2f' % w))
  else: raise runtimeError()

def printDFA(m):
  """ print m to console in a pretty format """
  print('initial state:', m[0])
  print('transitions:')
  for t in m[2].items(): printTransition(t)
  print('final states (with weights):')
  for f in m[1].keys(): print('   ',f,'--',m[1][f])

def kshowMx(mx):
  """ print out a matrix of pairs of floats in more readable form """
  for row in mx:
    print('\t'.join([('(%.2f,%.2f)' % (i[0],i[1])) for i in row]))

""" EXAMPLES
"""

def firstApprox():
  """ compute klx of Mohri'02 Figure 1 with cycle probabilities p1, p2 """
  p1,p2 = 0.9,0.1 # when p1,p2 = 0.9,0.1, these are examples fig1a,fig1b above
  m1 = (0, {1:p1}, { (0,'a'):(1,1.), (1,'b'):(1,1-p1) })
  #printDFA(m1)
  #dotDFA(m1,states(m1))
  m2 = (0, {1:p2}, { (0,'a'):(1,1.), (1,'b'):(1,1-p2) })
  #dotDFA(m2,states(m2))
  dot2DFA(m1,states(m1),m2,states(m2))
  print(klx(m1,m2,E))

def cycleApprox():
  """ plot klx for pairs of Mohri'02 Fig 1 with various cycle probabilities """
  p1s = [0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
  kls = []
  for x in p1s:
    p1,p2 = x, 0.5
    m1 = (0, {1:p1}, { (0,'a'):(1,1.), (1,'b'):(1,1-p1) })
    m2 = (0, {1:p2}, { (0,'a'):(1,1.), (1,'b'):(1,1-p2) })
    diff = klx(m1,m2,E)
    print('klx(fig1(%.2f),fig1(%.2f),%.2f) = %f' % (p1,p2,E,diff))
    kls.append(diff)
  kls2 = []
  for x in p1s:
    p1,p2 = 0.5, x
    m1 = (0, {1:p1}, { (0,'a'):(1,1.), (1,'b'):(1,1-p1) })
    m2 = (0, {1:p2}, { (0,'a'):(1,1.), (1,'b'):(1,1-p2) })
    diff = klx(m1,m2,E)
    print('klx(fig1(%.2f),fig1(%.2f),%.2f) = %f' % (p1,p2,E,diff))
    kls2.append(diff)
  # matplotlib plt is imported from SDFA_display
  plt.xlabel("probability p")
  plt.ylabel("approx kl in bits")
  plt.plot(p1s, kls, 'r', label="klx(fig1(%.2f), fig1(p), %.4f) for 0<p<1" % (p1,E))
  plt.plot(p1s, kls2, 'g', label="klx(fig1(p), fig1(%.2f), %.4f) for 0<p<1" % (p1,E))
  plt.legend()
  plt.title("varying the difference in loop probabilities in Mohri'02 Fig 1")
  plt.show()

def yuApprox():
  """ compare the slightly bigger machines from Yu&al'25 """
  print('--- m1 = the acceptor corresponding to tIDS')
  m1 = mIDS
  print('--- m2 = the acceptor corresponding to tADS')
  m2 = mADS
  print(klx(m1,m2,E))

if __name__ == '__main__':
  if not(VERBOSE): print('To see calculations, edit this file to set VERBOSE = True')
  firstApprox()
  #cycleApprox()
  #yuApprox()
