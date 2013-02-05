cython -a ~/grit/sparsify/sparsify_support_fns.pyx
gcc -shared -pthread -fPIC -fwrapv -O3 -Wall -fno-strict-aliasing \
    -I/usr/include/python2.7 \
    -o ~/grit/sparsify/sparsify_support_fns.so \
    ~/grit/sparsify/sparsify_support_fns.c
