#!/bin/bash

ls utest*.py > utests.txt

while read test; do
  echo "Running test: "$test
  python $test
done < utests.txt

rm utests.txt
