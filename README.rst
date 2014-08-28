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

    cd denovo-variant-caller
    mvn package

* Then, follow the `sign up instructions`_ to generate a valid
  ``client_secrets.json`` file.

* Move the ``client_secrets.json`` file into the client-java directory.
  (Authentication will take place the first time you make an API call.)

There are two modes for Denovo calling

* **Variant** Based - Examines variant calls and filters based on mendelian inheritance rules ::

    java -jar target/denovo-variant-caller-0.1.jar --caller variant \
    --client_secrets_filename ${HOME}/Downloads/client_secrets.json \
    --debug_level 1 \
    --chromosome chr1 \
    --output_file NA12878_stage1.calls \
    --num_threads 25 \
    --dad_callset_name NA12891 \
    --mom_callset_name NA12892 \
    --child_callset_name NA12878 \
    --dataset_id ${DATASET_ID} \
    --project_id ${PROJECT_ID}

* **Reads** based - Examines reads for candidate positions and filters based based on Bayesian evidence weighting. Lower false positive rate but more expensive to compute. Not that this step requires a  pre selected list of candidate positions ``--input_calls_file`` which can be obtained from the  previous variant based step. Default behavior is to invoke variants method if candidates are not  present. ::

    java -jar target/denovo-variant-caller-0.1.jar --caller read \
    --client_secrets_filename ${HOME}/Downloads/client_secrets.json \
    --input_calls_file  NA12878_candidates.calls
    --debug_level 1 \
    --chromosome chr1 \
    --output_file NA12878_stage1.calls \
    --num_threads 25 \
    --dad_callset_name NA12891 \
    --mom_callset_name NA12892 \
    --child_callset_name NA12878 \
    --dataset_id ${DATASET_ID} \
    --project_id ${PROJECT_ID}
    --inference_method bayes

Additional Options
------------------

To speed up execution increase the number of threads with the ``--num_threads`` 
option. 

To restrict to one or more chromosomes use the ``--chromosome`` flag.

See below for all options ::

  Usage: DenovoMain [flags...]
   --caller [VARIANT | READ]              : The caller to use (variant or read)
                                            based
   --child_callset_name <name>            : Child's callset name e.g. NA12879
   --chromosome [CHR1 | CHR2 | CHR3 |     : specify the chromosomes to search
   CHR4 | CHR5 | CHR6 | CHR7 | CHR8 |       (specify multiple times for multiple
   CHR9 | CHR10 | CHR11 | CHR12 | CHR13     chromsomes)
   | CHR14 | CHR15 | CHR16 | CHR17 |         
   CHR18 | CHR19 | CHR20 | CHR21 | CHR22     
   | CHRX | CHRY | CHRM]                     
   --client_secrets_filename <file>       : Path to client_secrets.json
   --dad_callset_name <name>              : Dad's callset name e.g. NA12877
   --dataset_id <id>                      : Dataset id
   --debug_level <level>                  : specify the debug level (0 for no
                                            debug spew)
   --denovo_mut_rate <rate>               : Specify the denovo mutation rate
                                            (default 1e-8)
   --end_position <position>              : end position ( usually set
                                            automatically )
   --inference_method [MAP | BAYES | LRT] : Inference method (map | bayes | lrt)
   --input_calls_file <file>              : File to read from
   --lrt_threshold <sig_level>            : likelihood ratio test significance
                                            level (default 1. ;higher the
                                            stricter)
   --max_api_retries <num>                : max api retry count (default 5)
   --max_variant_results <num>            : max variants returned per request
                                            (default 10000)
   --mom_callset_name <name>              : Mom's callset name e.g. NA12878
   --num_threads <num>                    : Specify the number of threads
                                            (default 1 ; 1 to 50 suggested)
   --output_file <file>                   : File to write results
   --project_id <id>                      : Project id
   --seq_err_rate <rate>                  : Specify the sequence error rate
                                            (default 1e-2)
   --start_position <position>            : start position ( usually 1 )
	
.. _Google Genomics API: https://developers.google.com/genomics
.. _Apache Maven: http://maven.apache.org/download.cgi
.. _sign up instructions: https://developers.google.com/genomics


The mailing list
----------------

The `Google Genomics Discuss mailing list <https://groups.google.com/forum/#!forum/google-genomics-discuss>`_ is a good
way to sync up with other people who use genomics-tools including the core developers. You can subscribe
by sending an email to ``google-genomics-discuss+subscribe@googlegroups.com`` or just post using
the `web forum page <https://groups.google.com/forum/#!forum/google-genomics-discuss>`_.
