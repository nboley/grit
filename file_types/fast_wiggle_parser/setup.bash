#!/usr/bin/bash
gcc -fpic -shared -Wl,-soname,libgtf.so.1 -o libbedgraph.so bedgraph_parser.c
