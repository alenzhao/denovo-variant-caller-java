denovo-variant-caller-java |Build Status|_ |Build Coverage|_
============================================================

.. |Build Status| image:: http://img.shields.io/travis/googlegenomics/denovo-variant-caller-java.svg?style=flat
.. _Build Status: https://travis-ci.org/googlegenomics/denovo-variant-caller-java

.. |Build Coverage| image:: http://img.shields.io/coveralls/googlegenomics/denovo-variant-caller-java.svg?style=flat
.. _Build Coverage: https://coveralls.io/r/googlegenomics/denovo-variant-caller-java?branch=master


Calls de novo variants using information from a mother, father and child trio.

Uses a bayes net encoded with the inheritence relationship in the trio in order
to judge the denovo calls. 

**NOTE** : Currently under development. Usage should be considered experimental.

Documentation
-------------
Documentation for the project can be found in `denovo.pdf`_.

.. _denovo.pdf: https://raw.githubusercontent.com/googlegenomics/denovo-variant-caller-java/master/denovo.pdf

Getting started
---------------

This Java program allows you to discover denovo variants using Bayesian de novo
variant calling.

* To use, first build the client using `Apache Maven`_::

    cd denovo-variant-caller-java
    mvn package

* Then, follow the `sign up instructions`_ to generate a valid
  ``client_secrets.json`` file.

* Move the ``client_secrets.json`` file into the client-java directory.
  (Authentication will take place the first time you make an API call.)

There are three modes for Denovo calling

* **Variants** Based - Examines variant calls and filters based on mendelian inheritance rules.

* **Reads** based - Examines reads for candidate positions and filters based based on Bayesian evidence weighting. Lower false positive rate but more expensive to compute. **Note** that this step requires a  pre selected list of candidate positions ``--input_calls_file`` which can be obtained from the previous variant based step.

* **Full** - A utility mode that runs both the variant and the reads mode for you such that output of variants mode is piped to reads mode ::

    java -jar target/denovo-variant-caller-0.1.jar --caller full \
    --client_secrets_filename ${HOME}/Downloads/client_secrets.json \
    --dataset_id 3049512673186936334 \
    --dad_callset_name NA12891 \
    --mom_callset_name NA12892 \
    --child_callset_name NA12878 \
    --chromosome chr1 \
    --start_position 1 \
    --end_position 14000000 \
    --log_level debug \
    --num_threads 25 \
    --output_file NA12878_full.calls


Additional Options
------------------

To speed up execution increase the number of threads with the ``--num_threads`` 
option. 

To restrict to one or more chromosomes use the ``--chromosome`` flag.

See below for all options ::

    Usage: DenovoMain [flags...]
     --caller [VARIANT | READ | FULL]       : The caller mode
     --child_callset_name <name>            : Child's callset name e.g. NA12879
     --chromosome <name>                    : specify the chromosomes to search
                                              (specify multiple times for multiple
                                              chromsomes)
     --client_secrets_filename <file>       : Path to client_secrets.json
     --dad_callset_name <name>              : Dad's callset name e.g. NA12877
     --dataset_id <id>                      : Dataset id
     --denovo_mut_rate <rate>               : Specify the denovo mutation rate
                                              (default 1e-8)
     --end_position <position>              : end position ( usually set
                                              automatically )
     --inference_method [MAP | BAYES | LRT] : Inference method (map | bayes | lrt)
     --input_calls_file <file>              : File to read from
     --log_file <file>                      : specify the log file
     --log_level [ERROR | INFO | DEBUG]     : specify the logging level
     --lrt_threshold <sig_level>            : likelihood ratio test significance
                                              level (default 1. ;higher the
                                              stricter)
     --max_api_retries <num>                : max api retry count (default 5)
     --max_variant_results <num>            : max variants returned per request
                                              (default 10000)
     --mom_callset_name <name>              : Mom's callset name e.g. NA12878
     --num_threads <num>                    : Specify the number of threads
                                              (default 1 ; 1 to 50 suggested)
     --output_dir <dir>                     : File to write results
     --output_file <file>                   : File to write results
     --seq_err_rate <rate>                  : Specify the sequence error rate
                                              (default 1e-2)
     --start_position <position>            : start position ( usually 1 )
	
.. _Google Genomics API: https://developers.google.com/genomics
.. _Apache Maven: http://maven.apache.org/download.cgi
.. _sign up instructions: https://developers.google.com/genomics


Building Documentation
----------------------

The documentation in this repository relies on the 
`LaTeX maven plugin <http://mojo.codehaus.org/latex-maven-plugin>`_

To build the documentation you need to first have ``pdflatex`` available on your system. 
Try `MacTeX <http://www.tug.org/mactex/>`_ for Macs or
`TeX Live <http://mirror.utexas.edu/ctan/systems/texlive/Images/>`_ for Windows.

::
    
    cd denovo-variant-caller-java
    mvn latex:latex
    cp target/denovo.pdf denovo.pdf

Todos / Next Steps
------------------
* The caller currently calls SNPs and ignores indels. This feature can be added by carefully 
  treating structural variations.
* Parameters in the bayes net are fixed and not learned. Baseline mutation rates
  could be learned for the trio under study.
* Additional supervised classifiers could be added to the set of callers. It 
  should be sufficient to derive from ``DenovoCaller`` class and initialized by
  ``DenovoCallers`` static factory.
* To get a correct estimate of the precision/recall values of the caller a gold
  standard dataset with de novo mutations is needed. Unfortunately, none such 
  exists. It can be closely approximated with blood derived DNA samples from 
  multiple trios of siblings.

    
The mailing list
----------------

The `Google Genomics Discuss mailing list <https://groups.google.com/forum/#!forum/google-genomics-discuss>`_ is a good
way to sync up with other people who use genomics-tools including the core developers. You can subscribe
by sending an email to ``google-genomics-discuss+subscribe@googlegroups.com`` or just post using
the `web forum page <https://groups.google.com/forum/#!forum/google-genomics-discuss>`_.
