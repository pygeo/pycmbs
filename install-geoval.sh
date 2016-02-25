#!/bin/sh
# install pycmbs from source
#~ set -ex
#~ wget https://protobuf.googlecode.com/files/protobuf-2.4.1.tar.gz
#~ tar -xzvf protobuf-2.4.1.tar.gz
#~ cd protobuf-2.4.1 && ./configure --prefix=/usr && make && sudo make install

git clone https://github.com/pygeo/geoval.git
cd geoval
sudo /usr/bin/python setup.py install
sudo sh compile_extensions.sh
