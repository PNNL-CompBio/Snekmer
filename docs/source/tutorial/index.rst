Tutorial
========

Running the Tutorial
--------------------

To run an example set of jobs use the following commands:

.. code-block:: bash

   source ~/snekmer_env/bin/activate
   cd resources/tutorial/demo_example
   python3 run_demo.py

This will execute the snekmer model, search, and cluster modes in succession
on a set of three input families and produce output files for each in output.


Running the Learn/Apply Tutorial  
--------------------  

To run an example set of fasta files use the following commands:  
  
.. code-block:: bash

   source ~/snekmer_env/bin/activate
   cd resources/tutorial
   jupyter notebook snekmer_learnapp_tutorial.ipynb




Evaluating Results
------------------

Apply writes a summary of results to ``output/Snekmer_Apply_Report.html``. Below is an example of a report file.

.. image:: resources/images/apply_report_example.png 

See the :ref:`Accessing Results <usage-results>` section for more details.

