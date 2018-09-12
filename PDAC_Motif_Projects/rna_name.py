# import OS module
import os

# removes the "unc.edu." part from file names. 
for filename in os.listdir("."):
	#if filename.startswith("unc.edu."):
	if filename.endswith(".rsem.genes.normalized_results"):
		os.rename(filename, filename[0:])

