#!/usr/bin/env python

import sys
import dendropy
from dendropy.calculate.treecompare import false_positives_and_negatives

argc = len(sys.argv)
assert (argc == 3), "Usage compare_tree.py tree1 tree2"

tax = dendropy.TaxonNamespace()
tr1 = dendropy.Tree.get(path=sys.argv[1],
                       schema='newick',
                       rooting="force-unrooted",
                       taxon_namespace=tax)
tr2 = dendropy.Tree.get(path=sys.argv[2],
                       schema='newick',
                       rooting="force-unrooted",
                       taxon_namespace=tax)
tr1_len=tr1.length()
for e in tr1.postorder_edge_iter():
    if e.length != None:
        e.length = e.length / tr1_len
tr2_len=tr2.length()
for e in tr2.postorder_edge_iter():
    if e.length != None:
        e.length = e.length / tr2_len
ed = dendropy.calculate.treecompare.euclidean_distance(tr1,tr2)

tr1.encode_bipartitions()
tr2.encode_bipartitions()

[fp, fn] = false_positives_and_negatives(tr1, tr2)
# assert (fp == fn), "Trees may not be binary or over the same namespace"

nl = len(tax)
rf = (fp + fn) / (2. * nl - 6)
len1=tr1.length()
len2=tr2.length()
sys.stdout.write('%d, %d, %d, %1.6f, %e, %e, %e\n' % (nl, fp, fn, rf, tr1_len, tr2_len, ed))
