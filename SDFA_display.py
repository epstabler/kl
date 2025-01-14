""" DFA graphical display
"""
from graphviz import Source
import matplotlib.pyplot as plt
from numpy import ndarray

def dotDFA(m, qs):
  """ display m on screen, where states and/or weights might be pairs of numbers """
  dotfile = ''
  dotfile += 'digraph finite_state_machine {\n'
  dotfile += 'rankdir = LR;\n'
  # initial state with bold circle, all states numbered by position in qs
  dotfile += 'node [shape=circle,style=bold] %s;\n' % qs.index(m[0])
  # final states with double circle
  finalStates = ' '.join([str(qs.index(k)) for k in m[1].keys()])
  dotfile += 'node [shape=doublecircle,style=solid] %s;\n' % finalStates
  # all other states with plain circle
  dotfile += 'node [shape=circle,style=solid];\n'
  # now list all states with their labels
  for q in qs:
    if isinstance(qs.index(q),tuple) and len(qs.index(q))==2:
      dotfile += '%s [label="(%s,%s)"];\n' % (str(qs.index(q)),str(q[0]),str(q[1]))
    else:
      dotfile += '%s [label="%s"];\n' % (str(qs.index(q)),str(q))
  # show each edge with symbol,weight
  for e in m[2].items():
    ((src,input),(target,p)) = e
    if (isinstance(p,tuple) or isinstance(p,ndarray)) and len(p)==2:
      edgeLabelParts = (str(qs.index(src)),
                        str(qs.index(target)),
                        str(input),
                        ('(%.2f,%.2f)' % (p[0],p[1])))
      dotfile += '%s -> %s [label="%s,%s"];\n' % edgeLabelParts
    else:
      edgeLabelParts = (str(qs.index(src)),
                        str(qs.index(target)),
                        str(input),
                        ('%.2f' % p))
      dotfile += '%s -> %s [label="%s,%s"];\n' % edgeLabelParts
  dotfile += '}\n'
  #print('------ dot file:\n' + dotfile + '------')
  s = Source(dotfile, filename="/tmp/dfa", format="png")
  s.view()

def new(q,i):
  """ to show a second machine in same display, rename the states of m2.
  """
  if isinstance(q,int): return i+q  # we add len(qs1) to positions in qs2 """
  elif isinstance(q,str): return 2*q
  elif isinstance(q,tuple) and len(q)==2: return (new(q[0]),new(q[1]))
  else: raise runtimeError("new: unexpected state of machine 2")

def dot2DFA(m1, qs1, m2, qs2):
  """ display two DFAs """
  # rename states in m2 with simple strategy (error if there are clashes)
  m2increment = len(qs1)
  print('qs2 =',qs2)
  newqs2 = [new(q,m2increment) for q in qs2]
  print('newqs2 =',newqs2)
  if [x for x in qs1 if x in newqs2]: raise RuntimeError('State name clash')
  qs12 = qs1 + newqs2
  m20 = new(m2[0],m2increment)
  m21 = dict([(new(f,m2increment),w) for (f,w) in m2[1].items()])
  m22 = dict([((new(qi,m2increment),ii),(new(qf,m2increment),w)) for ((qi,ii),(qf,w)) in m2[2].items()])
  # write graph to file
  dotfile = ''
  dotfile += 'digraph m12 {\n'
  dotfile += 'rankdir = LR;\n'
    # initial states with bold circle, all states numbered by positions in qs1,qs2
  initialStates = [m1[0], m20]
  initialStateString = ' '.join([str(qs12.index(k)) for k in initialStates])
  dotfile += 'node [shape=circle,style=bold] %s;\n' % initialStateString
  # final states with double circle
  finalStates = list(m1[1].keys()) + list(m21.keys())
  finalStateString = ' '.join([str(qs12.index(k)) for k in finalStates])
  dotfile += 'node [shape=doublecircle,style=solid] %s;\n' % finalStateString
  # all other states with plain circle
  otherStates = [x for x in qs12 if not(x in initialStates) and not(x in finalStates)]
  dotfile += 'node [shape=circle,style=solid];\n'
  # now list all states with their labels
  for q in qs12:
    if isinstance(qs12.index(q),tuple) and len(qs12.index(q))==2:
      dotfile += '%s [label="%s,%s"];\n' % (str(qs12.index(q)),str(q[0]),str(q[1]))
    else:
      dotfile += '%s [label="%s"];\n' % (str(qs12.index(q)),str(q))
  # show each edge with symbol,weight
  for e in list(m1[2].items()) + list(m22.items()):
    ((src,input),(target,p)) = e
    if (isinstance(p,tuple) or isinstance(p,ndarray)) and len(p)==2:
      edgeLabelParts = (str(qs12.index(src)),
                        str(qs12.index(target)),
                        str(input),
                        ('(%.2f,%.2f)' % (p[0],p[1])))
      dotfile += '%s -> %s [label="%s,%s"];\n' % edgeLabelParts
    else:
      edgeLabelParts = (str(qs12.index(src)),
                        str(qs12.index(target)),
                        str(input),
                        ('%.2f' % p))
      dotfile += '%s -> %s [label="%s,%s"];\n' % edgeLabelParts
  dotfile += '}\n'
  #print('------ dot file:\n' + dotfile + '------')
  s = Source(dotfile, filename="/tmp/dfa", format="png")
  s.view()
