cp -r grit/ grit-0.1.2/
cd grit-0.1.2/
`find ./ | grep -P "~$"`
rm -rf ./build/
rm -rf ./grit/build/
rm `find ./ | grep -P ".pyc$"`
rm ./grit/sparsify_support_fns.c
cd ..
tar -cf grit-0.1.2.tar grit-0.1.2/
rm -rf grit-0.1.2/
gzip grit-0.1.2.tar
