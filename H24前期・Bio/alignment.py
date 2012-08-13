#!/usr/bin/env python 
#coding=utf-8
#inspired by https://github.com/alevchuk/pairwise-alignment-in-python/
from numpy import *

class LocalAlignment(object):

  MATCH_AWARD      = 2
  MISMATCH_PENALTY = -1
  GAP_PENALTY      = -1 
  SCORE_LOWBOUND   = 0

  TRACEBACK_END  = 0
  TRACEBACK_UP  = 1
  TRACEBACK_LEFT = 2
  TRACEBACK_DIAGONAL = 3

  debug = False

  def alignment(self, sa, sb):
    m,n = len(sa), len(sb)
    score = zeros((m+1,n+1))
    pointer = zeros((m+1,n+1))
    max_score = max_i = max_j = 0

    for i in xrange(1, m+1):
      for j in xrange(1, n+1):
        # Calculate score table
        d = score[i-1,j-1] + \
            (self.MATCH_AWARD if sa[i-1] == sb[j-1] else self.MISMATCH_PENALTY)
        l = score[i-1][j] + self.GAP_PENALTY
        u = score[i][j-1] + self.GAP_PENALTY
        score[i,j] = max(d, l, u, self.SCORE_LOWBOUND)

        # Calculate trackback table
        if score[i,j] == u: pointer[i,j] = self.TRACEBACK_UP
        if score[i,j] == l: pointer[i,j] = self.TRACEBACK_LEFT
        if score[i,j] == d: pointer[i,j] = self.TRACEBACK_DIAGONAL
        if score[i,j] == 0: pointer[i,j] = self.TRACEBACK_END
        if score[i,j] >= max_score:
          # Log down the position of max score cell
          max_i, max_j, max_score = i, j, score[i,j]

    if self.debug: print score.transpose()
    # traceback
    i, j = max_i, max_j
    align_a = align_b = ""
    while pointer[i,j] != self.TRACEBACK_END:
      if pointer[i,j] == self.TRACEBACK_DIAGONAL:
        if self.debug: print "cell(%d,%d):diagonal" % (i,j)
        align_a = sa[i-1] + align_a
        align_b = sb[j-1] + align_b
        i -= 1; j -=1
      elif pointer[i,j] == self.TRACEBACK_LEFT:
        if self.debug: print "cell(%d,%d):left" % (i,j)
        align_a = sa[i-1] + align_a
        align_b = '-' + align_b
        i -= 1
      elif pointer[i,j] == self.TRACEBACK_UP:
        if self.debug: print "cell(%d,%d):up" % (i,j)
        align_a = '-' + align_a
        align_b = sb[j-1] + align_b
        j -= 1
    if self.debug: print "max score:%d" % max_score
    if self.debug: print "result:\n%s\n%s" % (align_a, align_b)
    return align_a, align_b

# test
la = LocalAlignment()
la.debug = True
a,b = la.alignment("ACACACTA","AGCACACA")
c,d = la.alignment("GCTCGTTG","AACCGTAA")
