This folder contains 4 algorithms (each in subfolders) and 2 *.pdf files
which describe the four methods.

% **** UNCONSTRAINED OPTIMIZATION ************************************

The file, "AIAA-2010-2912 Unconstrained.pdf" describes the unconstrained
multifidelity optimization method which is contained in the folder "Unconstrained".
Within that folder is another Readme file describing all of the included files and
how to run the method.


% ****** CONSTRAINED OPTIMIZATION ************************************

The file, "MDAO_2010 Constrained.pdf" describes the three constrained formulations:
	* Algorithm 1 in that paper is to optimize a high-fidelity objective function
		using multifidelity methods. This method is contained in the folder 
		"AppObjective". A readme file is within that folder describing all of the
		codes.

	* The folder "AppConstraint" performs a simplified version of Algorithm 2 in that
		paper. It optimizes a single-fidelity objective subjet to a multifidelity
		constraint and other easy constraints. There is a readme file contained
		in that directly explaining what all of the files are and how to run the
		codes.

	* The folder "AppEverything" performs the full version of Algorithm 2 in that paper.
		It optimizes a multifidelity objective function subject to multifidelity
		constraint and other easy constraints. That folder also contains a readme
		file describing all of the codes and how to run them.