******************
GIMSAN Readme
******************

################
Installation
################

1. Install required Python libraries::

    sudo pip install -r pip_requirements.txt

#. Compile GibbsMarkov and ColumnDependency::
    
    ./compile_clean.sh

################
Usage Example
################

1. Set and verify GIMSAN job configuration file ``conf_examples/test_window_sampling.cfg``. 
   Ensure your ``gimsan_home`` and ``r_path`` directory/path are pointing to the
   correct location.

#. Submit GIMSAN job using window-sampling::
    
    ./gimsan_submit.py --conf=conf_examples/test_window_sampling.cfg -v
    
#. Generate results and HTML output for GIMSAN job. If parameter ``main_output_dir=testout/``
   is set in configuration file, then run::

    ./gimsan_result.py --dir=testout 
    
#. Open HTML file ``testout/ABF1_YPD_mod/output.html``
