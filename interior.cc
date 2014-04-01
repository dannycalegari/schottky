/* interior.cc 

	Algorithm to use the method of Efficient Covers to test whether
	a point z is in (the interior of) Set A.
	
	Here is a description of the algorithm, in slightly abstract terms.
	
	We have a 2-generator semigroup G:=<f,g> depending on a complex
	parameter z where
	
	f: x -> zx and g: x -> z(x-1) + 1
	
	A disk D is good if fD and gD are contained in D. If D is good, define
	
	L: = the intersection of G_n D over all n
	
	where 
		G_n means words in G of length n,
		G_n D means the union of translates wD over all w in G_n
		
	Note that L does not depend on the choice of a good disk D. If we define
	E to be the set of right-infinite words in the alphabet {f,g}, then there
	is a map p: E -> L which takes an infinite word w_infty to the limit
	w_n(1/2) where w_n is the initial subword of w of length n.
 	
	Set A is the set of complex numbers z for which L is disconnected
	(equivalently, for which L is a Cantor set). Set A is closed, and 
	contains the set of z with |z|>=1/sqrt{2}=0.7071067811 (approximately)
	So we can always just fix D to be the disk of radius 4 centered at 1/2.
	
	A word w of length n is *good for z* if there is some m>n so that for every
	word wu of length m beginning with w, the disk wuD is completely contained
	in the union of disks of the form w'D as w' ranges over G_n - w.
	
	We know that z is in Set A if there is some word w of some length
	which is good for z.
	
	The algorithm proceeds as follows: first decide on n and a bound d for m-n. 
	Then for each word w in G_n we do the following:
		Test if w is good for z as follows.
		Are wfD and wgD contained in some w'D with w' in G_n-w? To test, first
		see which of fD and gD they are contained in, then which of
		ffD, fgD, gfD, ggD (pruning the tree as we go) to find if they are in
		some w'D. If they both are, then w is good for z with m=1.
		Otherwise, if wfD is not in some w'D as above, check whether wffD and
		wfgD are in some w'D as above, and similarly for wgD. Basically,
		refine the tree of descendents of w of length at most n+d until we
		get a word which is contained in some w'D.
	If we succeed for some w in G_n	then z is in Set A.

 */
