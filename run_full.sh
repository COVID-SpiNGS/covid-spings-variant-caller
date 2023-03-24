#!/bin/bash

mkdir -p ./log
mkdir -p ./tmp
mkdir -p ./output

python -m client_server.live_server
python -m watcher.watcher /home/vickylara/Documents/test