./configure --prefix=$PREFIX --enable-shared --disable-java --disable-dependency-tracking --with-readline=$PREFIX --with-xml-prefix=$PREFIX CFLAGS="-I${PREFIX}/include -I${PREFIX}/include/libxml2 $CFLAGS"
export CFLAGS="-I${PREFIX}/include -I${PREFIX}/include/libxml2 $CFLAGS"
export XML_CFLAGS="-I${PREFIX}/include -I${PREFIX}/include/libxml2"
export XML_LIBS="-L${PREFIX}/lib -lxml2 -lz -llzma -lpthread -liconv -licui18n -licuuc -licudata -lm"

make MOTIF_APS="" PARTS="mdsshr mdsdcl treeshr tdishr mdstcpip mdslib" MISC_PARTS="tdi include"
make install MOTIF_APS="" PARTS="mdsshr mdsdcl treeshr tdishr mdstcpip mdslib" MISC_PARTS="tdi include"
cd mdsobjects/python
python setup.py install --single-version-externally-managed --record record.txt
