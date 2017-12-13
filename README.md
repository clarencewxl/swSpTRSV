# swSpTRSV
Sparse triangular solve for Sunway architecture

1.Login to Sunway TaihuLight

- Land the homepage of the National Supercomputing Center in Wuxi: http://www.nsccwx.cn/wxcyw/

- Select one VPN service: 'Telecom', 'Unicom' or 'China Mobile' on the top of the website. Please choose the best one for better connection.

- Login to Sunway TaihuLight Supercomputer: ssh 41.0.0.188

Please email to clarencewxl@gmail.com for further questions.

2.Experiments for the Temperature Map

- Build the binary:

cd ./SpTRSV_map

makeall.sh

- Run a single benchmark

./test.sh -I benchmark_name

Note: The names of all the benchmarks are listed in the file: benchmarks_list in the same directory

- The referenced output can be found in the following format:

Filename: ../triangular_files/atmosmodd.cscu, PRODUCER_CONSUMER_ROWS: 8 CACHE_X_V4: 1024 Average time is 0.005880s, Average MFlops is 1499.029826

- Use nohup to run all the benchmarks

nohup ./total_run.sh &

Note: Use bjobs to check whether the jobs finish or not. Each job will produce a new file in the subdirectory 'results' in the current directory.

- The referenced performance results are listed in the file 'ref_map' in the current directory.

Note: There are totally 4x6x2057=49368 works and will cost about 36 hours. You do not need to wait for finishing all the test cases and can compare the result of any case with the referenced result as long as it is finished.

3.Experiments for the Scatter Plot

- Build the binary:

cd ./SpTRSV_tuning

make

- Run a single benchmark

./run.sh -I benchmark_name

Note: The names of all the benchmarks are listed in the file: benchmarks_list in the same directory

- The referenced output can be found in the following format:

SwSpTRSV: Filename is ../triangular_files/atmosmodd.cscu Average time is 0.005869s, Average MFlops is 1501.908293

Serial: Filename is ../triangular_files/atmosmodd.cscu Average time is 0.110587s, Average MFlops is 79.710070

level-sets: Filename is ../triangular_files/atmosmodd.cscu Average time is 0.055369s, Average MFlops is 159.201148

- Use nohup to run all the benchmarks

nohup ./total_run.sh &

Note: Use bjobs to check whether the jobs finish or not. Each job will produce a new file in the subdirectory 'results' in the current directory.

- The referenced performance results are listed in the file 'ref_SwSpTRSV', 'ref_serial' and 'ref_levelsets' in the current directory. 

Note: There are totally $2057$ works and will cost about 2 hours. You do not need to wait for finishing all the test cases and can compare the result of any case with the referenced result as long as it is finished.
