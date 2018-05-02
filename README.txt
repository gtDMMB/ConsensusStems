README for operating ConsensusStems.pl

ConsensusStems code for Academic Users, Version 1.0 (April 2018)                                               

Usage: perl ConsensusStems.pl [options]... directory_of_sequences

Required: perl, Boltzmann sampling with constraints, RNA profiling, a MFE prediction program
        Boltzmann sampling with constraints is implemented by Sfold, which can be requested from http://sfold.wadsworth.org/cgi-bin/index.pl
        RNA profiling is available at https://github.com/gtfold/RNAStructProfiling
        A recommended MFE prediction program is GTfold's gtmfe option, available at https://github.com/gtfold/gtfold
                                                                                                               
Input a directory containing all sequences (in FASTA format) to be processed.                                  
        Recommended: at least 8 related sequences.                                                             
        ConsensusStems is not designed nor tested with less than 8 sequences                                   
                                                                                                               
Options:                                                                                                       
        -s <path>        Path to invoke Sfold, eg -s /usr/bin/sfold (default)                                  
        -p <path>        Path to invoke RNAprofiling, eg -p ../RNAprofiling (default is ./RNAprofile)          
        -m <path>        Path to invoke MFE program, eg -m /usr/bin/gtmfe (default)                            
        -d <name>        Path to MFE NNTM parameters eg -d Desktop/data/Turner99 (default is ./data)
	-o <path>        Path to store output files, eg -o ../output (default is ./output)                                     
        -f <name>        Family name, eg -f tRNA (default is last directory name if exists or RNA)  
	-v		 Flag to enable verbose output concerning details of the method

FILES:
Generated in the <outputname> given with the -o option (or default in ./output):

--ConsensusStems_output.txt (containing the final clusters)

--seqs (directory containing the made-up sequences derived from the search windows, used to run MFE validation)

For each <seq> in <dir_of_sequences>:
--<seq>_constraints.txt (optional; containing the last constraints file used to resample before termination)
--<seq>.out (containing the result of the initial profiling output given Sfold sample of <seq>)
--<seq>_resample (optional; directory containing the result of last Sfold sampling with constraints)
--<seq>_resample.out (optional; containing the result of the profiling output given the last Sfold resample of <seq>)
--Sfold_<seq> directories containing output of initial Sfold sampling on <seq>


