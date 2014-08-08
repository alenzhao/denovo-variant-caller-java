denovo-variant-caller
=====================

Calls de novo variants using information from a mother, father and child trio.

Uses a bayes net encoded with the inheritence relationship in the trio in order
to judge the denovo calls. 

**NOTE** : Currently under development. Usage should be considered experimental.

Getting started
---------------

This Java program allows you to discover denovo variants using Bayesian de novo
variant calling.

* To use, first build the client using `Apache Maven`_::

    cd api-client-java
    mvn package

* Then, follow the `sign up instructions`_ to generate a valid
  ``client_secrets.json`` file.

* Move the ``client_secrets.json`` file into the client-java directory.
  (Authentication will take place the first time you make an API call.)

There are two stages in de novo Calling:

Stage1
    Filter Candidate Variant Sites using VCF file called variants

Stage 2
    Apply Bayesian Infererence from raw reads using the list of candidate 
    variants from stage1

Run Stage 1 which outputs candidate calls to ``--output_file`` ::

    java -jar target/denovo-variant-caller-0.1.jar stage1 \
    --client_secrets_filename ${HOME}/Downloads/client_secrets.json \
    --job_name myStage1Job --output_file stage1.calls \
    --debug_level 1

Run Stage 2 which reads in candidate calls  from ``--candidates_file``::

    java -jar target/denovo-variant-caller-0.1.jar stage2 \
    --client_secrets_filename ${HOME}/Downloads/client_secrets.json \
    --job_name MyStage2Job \
    --input_file stage1.calls \
    --output_file stage2.calls \
    --debug_level 1

Additional Options
------------------

To speed up execution increase the number of threads with the ``--num_threads`` 
option. 

To restrict to one or more chromosomes use the ``--chromosome`` flag.

See below for all options ::

	Usage: GenomicsExperiment stage_id [flags...]
	 <stage id>                             : The stage of the calling pipeline ;
		                                  usually stage1 or stage2
	 --chromosome <chromosome>              : specify the chromosomes to search
		                                  (specify multiple times for multiple
		                                  chromsomes)
	 --client_secrets_filename              : Path to client_secrets.json
	 <client_secrets_filename>                 
	 --debug_level <debug_level>            : specify the debug level (0 for no
		                                  debug spew)
	 --denovo_mut_rate <denovo_mut_rate>    : Specify the denovo mutation rate
		                                  (default 1e-8)
	 --inference_method <map|bayes|lrt>     : Inference method (map | bayes | lrt)
	 --input_file <file>                    : File to read from
	 --job_name <job name>                  : Name of your job
	 --lrt_threshold <lrt_sig_level>        : likelihood ratio test significance
		                                  level (default 1. ;higher the
		                                  stricter)
	 --num_threads <num_threads>            : Specify the number of threads
		                                  (default 1 ; 1 to 50 suggested)
	 --output_file <file>                   : File to write results
	 --seq_err_rate <seq_err_rate>          : Specify the sequence error rate
		                                  (default 1e-2)

.. _Google Genomics API: https://developers.google.com/genomics
.. _Apache Maven: http://maven.apache.org/download.cgi
.. _sign up instructions: https://developers.google.com/genomics


 
