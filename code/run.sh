#!/bin/bash

if [ ! -d ./Eigen ]; then
    ln -s /Eigen/Eigen ./Eigen
fi


# check to see if all required arguments were provided
if [ $# -eq 1 ]; then
    # assign the provided argument to variables
    param_file=$1
else
    # assign the default value to a variable
    param_file="InputHelixParams.txt"
fi

make

echo "Running tendonmech_test with argument $param_file:"
./tendonmech_test $param_file | tee ../results/tendonmech_test_results.txt

echo "Running tendonmech_test with argument $param_file:"
./tendonmech $param_file | tee ../results/tendonmech_results.txt
