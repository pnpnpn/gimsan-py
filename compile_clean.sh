#!/bin/bash -x

pushd gibbsmarkov/code/ && make clean && make && popd
pushd column_dependency_app/code/ && make clean && make && popd
