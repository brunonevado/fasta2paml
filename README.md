# fasta2paml: prepare input files for paml

  
Author: B. Nevado  

Makes phylip file from fasta file after:  
      - removing last codon if it is a stop  
      - masking any stop codons within alignment  
      - removing sequences with ambiguous bases % above threshold (-cutSeqs)  
      - removing codon positions with ambiguous proportion above threshold (-cutCodons)  
Usage: fasta2paml -infile in.fas -outfile out_prefix -names names.txt -tree tree.tre -cut n -dropBranches 1/0 -unroot 1/0 -strictNames 1/0  
Args: -infile: fasta file with sequences.  
      -outfile: prefix for output files (will write to outfile.phy and outfile.tre).  
      -names: name of sequences to write - sequences or leaves not named will be removed.  
      -tree: input tree file, newick format, with or without branch lengths.  
      -cutSeqs|cutCodons: n threshold to remove sequences/codons (e.g. with 50, codons where more than 50% individuals have missing data, or sequences with > 50% missing bases, are removed).  
      -dropBranches: whether to keep (0) or remove (1) branch lengths if present in tree file.  
      -unroot: whether to unroot the tree before printing (otherwise the tree is printed as read, after removing tips).  
      -strictNames: if set to 1, species in '-names' missing in alignments will raise error, otherwise (set to 0) will raise warning only.  
Only sequences listed in file names.txt are printed (1 sequence name per line).  
Using '-names all' will print all sequences in infile.  
Sequences removed due to missing data are also removed from tree file.  
Sequences in tree that are not present in names file (or fasta file) are also removed from tree.  
  
  
Requirements: BIO++ (bio-core and bio-phyl)  
  
Installation (Linux / Mac):    
If bpp-core and bpp-seq libraries are installed and available on PATH, this should work (from within /source):    
 g++ -O3 -std=c++11 *cpp -lbpp-core -lbpp-phyl -o fasta2paml  

  
 A static-compiled version is provided in /bin    
  

