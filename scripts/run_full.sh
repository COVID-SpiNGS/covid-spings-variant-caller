#!/bin/bash

mkdir -p ./log
mkdir -p ./tmp
mkdir -p ./output

echo "$PWD"
./scripts/run_server.sh &
./scripts/run_watcher.sh "$1" &
