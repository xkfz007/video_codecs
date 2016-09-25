#!/bin/bash
for i in *
do
    test -d $i || continue
    \rm -rv $i/*
done
