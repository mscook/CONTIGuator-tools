CONTIGuator tools
=================

We're increasingly using CONTIGuator in our work! Here are some tools/wrappers.
Eventually we should look at getting these incorporated into the CONTIGuator
core.

**At the moment this only adds the excluded contigs to 
PseudoContig.embl file**

Usage::

    $ python CONTIGuator_tools.py -h
    usage: CONTIGuator_tools.py [-h] [-o OUT] [-v] EMBL Fasta

    Add CONTIGuator excluded contigs to PseudoContig.embl file

    positional arguments:
      EMBL               Path to PseudoContig.embl
      Fasta              Path to Excluded.fsa

    optional arguments:
      -h, --help         show this help message and exit
      -o OUT, --out OUT  Output PseudoContig_Excluded.embl to this file and
                         location
      -v, --verbose      Be loud & noisy

    By: Mitchell Stanton-Cook, E: m.stantoncook@gmail.com

