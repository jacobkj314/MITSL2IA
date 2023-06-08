# Subregular experiments

This is a notebook set up to replicate the experiments for an implementation of De Santo and Aksënova (2012)'s MITSL 2,2 learner tested on an array of languages belonging to separate subregular classes (SL, TSL, SP, ITSL, MTSL, MITSL). We also provide a Python 3 implementation of the learner following the standards of SigmaPie (https://github.com/alenaks/SigmaPie).

The notebook in general requires SigmaPie (https://github.com/alenaks/SigmaPie).  Our own implementation of the MITSL learner is in local_sigmapie/code/mitsl_class.py. 

The general evaluation pipeline is adapted following Chapter 3 of Aksënova (2020).
All code in the first few sections of the notebook was taken verbatim from https://github.com/alenaks/subregular-experiments.

The section 'Experiment 10: Multi-tier Input-Sensitive Harmony' is adapted from De Santo & Aksënova (2021), with code taken verbatim from https://github.com/alenaks/2IMTSL.

The MITSL Experiments section is modelled after the set-up in Aksënova (2020).

Artificial samples for the subregular languages are generated within the notebook itself using SigmaPie, but some external files are needed to provide the natural language corpora. They can be found in the natural_data folder.

Sources of the real linguistic data:
* German data (`german.txt`): https://github.com/enz/german-wordlist
* Finnish data (`finnish.txt`): https://github.com/douglasbuzatto/WordLists
* Turkish data (`turkish.txt`): http://www.swarthmore.edu/SocSci/harmony/public_html/dummyresults.html