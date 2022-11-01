#!/usr/bin/env bash

pushd ../src
# build
make tests
make rootcat
make verify
popd
