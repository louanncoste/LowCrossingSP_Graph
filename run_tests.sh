#!/bin/bash

ALPHA=0.75

MAX_DURATION_HOURS=9

MAX_DURATION=$(( MAX_DURATION_HOURS * 3600 ))   #in seconds

declare -a test_names=( 
"BIOGRID-MV-Physical-3.5.169.edgelist"
"BIOGRID-SYSTEM-Affinity_Capture-MS-3.5.169.edgelist"
"BIOGRID-SYSTEM-Affinity_Capture-RNA-3.5.169.edgelist"
"buddha-w.txt"
"ca-HepPh.edgelist"
"CAIDA_as_20130601.edgelist"
"com-dblp.ungraph.txt"
"DIMES_201204.edgelist"
"dip20170205.edgelist"
"epinions1-d.txt"
"facebook_combined.txt"
"froz-w.txt"
"gnutella31-d.txt"
"grid300-10.txt"
"oregon2_010331.txt"
"p2p-Gnutella09.edgelist"
"twitter_combined-d.txt"
"xgrid500-10.txt"
"z-alue7065.txt"
"t.FLA-w.txt"
"notreDame-d.txt"
"t.CAL-w.txt"
)



for test in "${test_names[@]}"; do
    for (( i=1; i<=10; i++ )); do
        echo "Job started: $test (Iteration $i)"

        # Start the timer and run the command in the background
        start_time=$(date +%s)
        ./target/release/low_crossing_path data/Simple/$test ${test}_$i $ALPHA > big_tests/${test}_$i.out &
        pid=$!

        # Wait for the process to finish or timeout
        timed_out=false
        while kill -0 $pid 2> /dev/null; do
            current_time=$(date +%s)
            elapsed_time=$(( current_time - start_time ))
            if (( elapsed_time > MAX_DURATION )); then
                echo "Job $test (Iteration $i) exceeded maximum duration of $MAX_DURATION_HOURS hours. Skipping to next file."
                kill -9 $pid
                timed_out=true
                break
            fi
            sleep 1
        done

        # Check if the process was killed due to timeout
        if ! $timed_out; then
            echo "Job completed: $test (Iteration $i) in $((elapsed_time / 3600)) hours and $(((elapsed_time % 3600) / 60)) minutes."
        elif (( i == 1 )); then
            # If the iteration timed out and it was first iteration, break out of the loop for this file
            break
        fi
    done
done


