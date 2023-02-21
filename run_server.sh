#!/bin/bash

mkdir -p ./log
mkdir -p ./tmp
mkdir -p ./output

python -m main_live_server

python -m client_server.live_server

