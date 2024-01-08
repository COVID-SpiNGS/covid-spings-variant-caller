#!/bin/bash

mkdir -p ./log
mkdir -p ./tmp
mkdir -p ./output

echo "$PWD"
./run_server.sh &
#python3 -m client_server.live_server
./run_watcher.sh "$1" &
#python3 -m watcher.watcher /home/vickylara/Documents/test
