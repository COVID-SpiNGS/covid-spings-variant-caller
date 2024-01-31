#!/bin/bash

mkdir -p ./log
mkdir -p ./tmp
mkdir -p ./output

python -m src.client_server.live_server

