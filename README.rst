denovo-variant-caller
=====================

Calls de novo variants using information from a mother, father and child trio.
Uses a bayesian inference method to make the calls. 

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
    Apply Bayesian Infererence from raw reads

* You can then run Stage 1 which outputs candidate calls to ``candidates_file`` ::

    java -jar target/denovo-v1beta.jar stage1 --candidates_file candidate.calls --client_secrets_filename /home/subhodeep/Downloads/client_secrets.json --require_all_scopes

* You can then run Stage 2 which reads in candidate calls  from ``candidates_file``::

    java -jar target/denovo-v1beta.jar stage2 --candidates_file candidate.calls --client_secrets_filename /home/subhodeep/Downloads/client_secrets.json --require_all_scopes

.. _Google Genomics API: https://developers.google.com/genomics
.. _Apache Maven: http://maven.apache.org/download.cgi
.. _sign up instructions: https://developers.google.com/genomics


 
