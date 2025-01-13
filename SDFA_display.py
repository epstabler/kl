""" DFA graphical display
"""
from graphviz import Source
import matplotlib.pyplot as plt
from numpy import ndarray

def dotDFA(m, qs):
  """ display m on screen, where states and/or weights might be pairs of floats """
  dotfile = ''
  dotfile += 'digraph finite_state_machine {\n'
  dotfile += 'rankdir = LR;\n'
  # initial state with bold circle, all states numbered by position in qs
  dotfile += 'node [shape=circle,style=bold] %s;\n' % qs.index(m[0])
  if (isinstance(qs.index(m[0]),tuple) or isinstance(qs.index(m[0]),ndarray)) and len(qs.index(m[0]))==2:
    dotfile += '%s [label="(%s,%s)"];\n' % (str(qs.index(m[0])),str(m[0][0]),str(m[0][1]))
  else:
    dotfile += '%s [label="%s"];\n' % (str(qs.index(m[0])),str(m[0]))
  # final states with double circle
  finalStates = ' '.join([str(qs.index(k)) for k in m[1].keys()])
  dotfile += 'node [shape=doublecircle,style=solid] %s;\n' % finalStates
  for q in m[1].keys():
    if (isinstance(qs.index(m[0]),tuple) or isinstance(qs.index(m[0]),ndarray)) and len(qs.index(m[0]))==2:
      dotfile += '%s [label="(%s,%s)"];\n' % (str(qs.index(q)),str(q[0]),str(q[1]))
    else:
      dotfile += '%s [label="%s"];\n' % (str(qs.index(q)),str(q))
  # all other states with plain circle
  dotfile += 'node [shape=circle,style=solid];\n'
  for q in [x for x in qs if x!=m[0] and not(x in m[1].keys())]:
    if (isinstance(q,tuple) or isinstance(q,ndarray)) and len(q)==2:
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
  #if VERBOSE: print('------ dot file:\n' + dotfile + '------')
  s = Source(dotfile, filename="/tmp/dfa", format="png")
  s.view()
