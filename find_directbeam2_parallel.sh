#!/bin/bash
ls -1 ./mica_5/onepanel*.gz | xargs -P $(nproc) -n 1 ./find_directbeam2.sh 2>/dev/null &
