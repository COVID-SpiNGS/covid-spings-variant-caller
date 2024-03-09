#!/bin/bash

mkdir -p ./log
mkdir -p ./tmp
mkdir -p ./output

python -m src.architecture.server.live_server

