This script allows you to run variants of the algorithms on our case studies.
By modifying the variables, e.g. reps or output_dir, you can tune the script to your needs.
Commenting and uncommenting the configurations and models allows selecting the combinations of those that you are interested in.
Use
''' python3 run_benchmarks.py run '''
to execute every uncommented configuration on every uncommented model.
Use 
''' python3 run_benchmarks.py run '''
to read the resulting output files and compile them into a file output_dir+"/times.csv" with the runtimes and output_dir+"/values.csv" with the results.
