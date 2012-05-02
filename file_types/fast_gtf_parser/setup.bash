#!/usr/bin/bash
gcc -O2 -fpic -shared -Wl,-soname,libgtf.so.1 -o libgtf.so gtf_parser.c
