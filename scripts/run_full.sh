#!/bin/bash

mkdir -p ./log
mkdir -p ./tmp
mkdir -p ./output

echo "$PWD"
./run_server.sh &
./run_watcher.sh "$1" &
