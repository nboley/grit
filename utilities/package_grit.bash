cp -r ./grit/ ../grit-0.1.1/
cd ../grit-0.1.1/
rm `find | grep -P '.pyc$'`
rm `find ./ | grep '~$'`
rm build/ -rf
rm build/grit/ -rf
cd ..
tar -cvf grit-0.1.1.tar ./grit-0.1.1/
gzip grit-0.1.1.tar
rm -rf ./grit-0.1.1/
