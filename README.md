GIMSAN
===

http://www.cs.cornell.edu/~ppn3/gimsan/

Installation
---

1. Install required Python libraries:

        sudo pip install -r pip_requirements.txt

2. Compile GibbsMarkov and ColumnDependency:

        ./compile_clean.sh

Usage Example
---

1. Set and verify GIMSAN job configuration file `conf_examples/test_window_sampling.cfg`.
   Ensure your `gimsan_home` and `r_path` directory/path are pointing to the
   correct location.

2. Submit GIMSAN job using window-sampling:

         ./gimsan_submit.py --conf=conf_examples/test_window_sampling.cfg -v

3. Generate results and HTML output for GIMSAN job. If parameter `main_output_dir=testout/`
   is set in configuration file, then run:

         ./gimsan_result.py --dir=testout

4. Open HTML file `testout/ABF1_YPD_mod/output.html`
