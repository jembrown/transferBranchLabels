#! /usr/bin/env python

#################################################################################
#	transferBranchLabels.py			 											#
#									 											#
#	Copyright Jeremy M. Brown, 2010												#
#	jeremymbrown@gmail.com														#
#																				#
#  This program is free software; you can redistribute it and/or modify			#
#  it under the terms of the GNU General Public License as published by			#			
#  the Free Software Foundation; either version 3 of the License, or			#
#  (at your option) any later version.											#
#																				#
#  This program is distributed in the hope that it will be useful,				#
#  but WITHOUT ANY WARRANTY; without even the implied warranty of				#
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the				#
#  GNU General Public License for more details.									#
#																				#
#  You should have received a copy of the GNU General Public License along		#
#  with this program. If not, see <http://www.gnu.org/licenses/>.				#
#																				#
#################################################################################

"""
Program to take branch labels (e.g., support values) from an arbitrary number of labeled trees
(e.g., majority-rule consensus from bootstrapped trees or a posterior distribution) and place any 
relevant labels on another (e.g., ML) tree.  Topologies need not be identical or nested.

NOTE: Requires the previous installation of Dendropy v3
"""

import dendropy
import sys

def main(argv):

	if (argv[0] == "--help"):
		usage()
	else:
		print """
# # # # # # # # # # # # # # # # #
#                           	#
#   transferBranchLabels.py     #
#                             	#
#   Jeremy M. Brown         	#
#   jeremymbrown@gmail.com  	#
#                           	#
# # # # # # # # # # # # # # # # #
		"""
	
		transferLabels(argv)

		print """Execution complete...
		"""
	
	sys.exit(0)
		
		
		
def transferLabels(argv):
	"""
	Transferring branch labels from arbitrary number of input trees with support
	"""
	
	# Defining general taxonSet and treeList
	taxa = dendropy.TaxonSet()
	trees = dendropy.TreeList(taxon_set=taxa)

	# Reading in target (e.g., ML) tree and placing first in the treeList
	targetTree = dendropy.Tree.get_from_path(argv[0],"nexus",taxon_set=taxa)
	trees.append(targetTree)	
	
	# Reading in all consensus trees provided at the command line
	for conTree in argv[1:]:
		trees.append(dendropy.Tree.get_from_path(conTree,"nexus",taxon_set=taxa))	
	
	# Creates reference list of all leaves
	all_leaves = trees[0].leaf_nodes()

	# Creates dictionary to hold internal bitmasks as keys and Node OIDs as values
	target_bitmasks = {}
	
	# Create list of bitmasks for each internal node in target tree
	for i in trees[0].internal_nodes():
		if i is not trees[0].seed_node:
			bitmask = bit_mask([j.taxon.label for j in i.leaf_nodes()],[k.taxon.label for k in all_leaves])
			target_bitmasks["".join([str(num) for num in bitmask])] = i.oid
		i.label = []
		for conTree in argv[1:]:	# Initializes label list for internal nodes to "-" for each conTree
			i.label.append("-")
	
	# Iterate over internal nodes in each consensus tree, comparing against target bitmasks
	# 	Store support from consensus as new Node attribute when matching bitmasks found
	conTreeIndex = 0
	for conTree in trees[1:]:
		conCladeInTarget = 0
		conCladeNotInTarget = 0
		for int in conTree.internal_nodes():
			if int is not conTree.seed_node:
				bitmask = bit_mask([j.taxon.label for j in int.leaf_nodes()],[k.taxon.label for k in all_leaves])
				try:
					focalNodeOID = target_bitmasks["".join([str(num) for num in bitmask])]					
					conCladeInTarget += 1
					oidFilter = lambda x: x.oid == focalNodeOID
					focalNode = trees[0].find_node(filter_fn=oidFilter)
					if (int.label is not None):
						focalNode.label[conTreeIndex] = int.label
				except KeyError:
					conCladeNotInTarget += 1
								
		print("conTreeIndex: %d" % conTreeIndex)
		print("conCladeInTarget: %d" % conCladeInTarget)
		print("conCladeNotInTarget: %d" % conCladeNotInTarget)
		print"""
		"""
		conTreeIndex += 1
		
	# Two new tree-writing functions:
	#	(1) Support for each consensus type is written to all branches in target tree
	#	(2) Support is only written to target tree branches with at least one new consensus value

	treeString = writeWithLabels(trees[0].seed_node,"/")
	outFile = open(argv[0]+"_labeled.tre",'w')
	outFile.write(treeString)
	outFile.close()

def writeWithLabels(node,delim):
	"""
	Novel tree-writing method for multiple internal node labels stored as a list (separated by delim in output)
	"""
	if node.is_internal():
		treeString = "("
		for child in node.child_nodes():
			treeString += writeWithLabels(child,delim)
			if child is not node.child_nodes()[len(node.child_nodes())-1]: # Outputs a "," after non-terminal child nodes
				treeString += ","
		treeString += ")"
		for label in range(len(node.label)):
			if node.label[0] is not "-":
				treeString += '%.2f'% float(node.label.pop(0))
			else:
				treeString += node.label.pop(0)
			if len(node.label) > 0: # Outputs a delim character after non-terminal consensus values
				treeString += delim
		treeString += ":"
		if node.edge_length is not None:
			treeString += '%f' % node.edge_length
	if node.is_leaf():
		treeString = node.taxon.label
		treeString += ":"
		treeString += '%f' % node.edge_length
	return treeString
	
def bit_mask(focal,ref):
	"""
	Generic function to create bitmask corresponding to a particular bipartition
	"""
	bitmask = []
	for i in ref:
		if i == ref[0]:
			bitmask.append(1)
		elif (i in focal) == (ref[0] in focal):
			bitmask.append(1)
		else:
			bitmask.append(0)
	return bitmask
	

def usage():
	"""
	Reminds user of proper command-line syntax
	"""
	print """
transferBranchLabels.py targetTreeFile consensusTreeFile1 [consensusTreeFile2 [consensusTreeFile3 [...]]]
	"""

if __name__ == "__main__":
	main(sys.argv[1:])