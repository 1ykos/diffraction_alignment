#!/bin/bash
echo 541 514 | ./find_directbeam2 <(gunzip -c < $1) > $(basename $1 .gz).directbeam
# # # ls -1 ./mica_5/onepanel*.gz | xargs -P 16 -n 1 ./find_directbeam2.sh 2>/dev/null &
