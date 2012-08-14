#!/usr/bin/bash
gcc -fpic -shared -Wl,-soname,libbedgraph.so.1 -o libbedgraph.so bedgraph_parser.c
