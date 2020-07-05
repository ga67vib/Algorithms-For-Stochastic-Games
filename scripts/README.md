This script allows you to run variants of the algorithms on our case studies.
By modifying the variables, e.g. reps or output_dir, you can tune the script to your needs.
Commenting and uncommenting the configurations and models allows selecting the combinations of those that you are interested in.
Use<br>
''' python3 run_benchmarks.py run '''<br>
to execute every uncommented configuration on every uncommented model.
Use <br>
''' python3 run_benchmarks.py read '''<br>
to read the resulting output files and compile them into a file output_dir+"/times.csv" with the runtimes and output_dir+"/values.csv" with the results.
