# SW-SpTRSV
Sparse triangular solver for Sunway architecture

1.Login to Sunway TaihuLight

- Land the homepage of Sunway TaihuLight Supercomputer: http://www.nsccwx.cn/wxcyw/

- Select one VPN service: 'Telecom', 'Unicom' or 'China Mobile' on the top of the website. Please choose the best one for better connection.

- Login to Sunway TaihuLight Supercomputer: ssh 41.0.0.188

Please email to clarencewxl@gmail.com for further questions.

2.Experiments for the Temperature Map

- Build the binary:

cd ./SpTRSV_map
makeall.sh

- Run a single benchmark

./test.sh -I benchmark_name
Note: The benchmarks' names are listed in the file: benchmarks\_list in the same directory

- Use nohup to run all the benchmarks

nohup ./total_run.sh &
Note: Use bjobs to check whether the jobs finish or not.

- The referenced performance results are listed in the file 'ref_map' in the current directory.
Note: It will take about 24 hours to run all the benchmarks.

3.Experiments for the Scatter Plot

- Build the binary:

cd ./SpTRSV_tuning
make

- Run a single benchmark

./run.sh -I benchmark_name
Note: The benchmarks' names are listed in the file: benchmarks\_list in the same directory

- Use nohup to run all the benchmarks

nohup ./total_run.sh &
Note: Use bjobs to check whether the jobs finish or not.

- The referenced performance results are listed in the file 'ref_SwSpTRSV' in the current directory.
Note: It will take about 3 hours to run all the benchmarks.
