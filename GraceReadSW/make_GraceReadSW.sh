#!/bin/bash -ue

# keeping track of where I am
DIR_NOW="$(cd $(dirname $BASH_SOURCE); pwd)"

#update this as needed
VERSION="2010-03-31"
FILE_ROOT="GraceReadSW_L1_"
FILE_EXT=".tar"

#derived parameters
LOCAL_DIR="$(pwd)"
FILENAME=$FILE_ROOT$VERSION$FILE_EXT
LOCAL_FILE="$LOCAL_DIR/$FILE_ROOT$VERSION$FILE_EXT"
SW_DIR="$LOCAL_DIR/RELEASE_$VERSION"

#remove software directory 
[[ ! "${@/--clean/}" == "$@" ]] && rm -fvr "$SW_DIR"

#decompress package (if not already)
[ ! -d "$SW_DIR" ] && tar -x -C "$LOCAL_DIR" -f "$LOCAL_FILE" -v


#need to tweak the make files
find "$SW_DIR" -name Makefile | while read i
do 

	  echo "$i"
	  echo " "
	  if [ ! -z "$(grep '^CFLAGS = $' "$i")" ] || [ ! -z "$(grep '^CFLAGS =$' "$i")" ]
	  then
		CFLAGS=$(echo "
		  -Wno-return-type
		  -Wno-pointer-sign
		  -Wno-implicit-int
		  -Wno-implicit-function-declaration
		  -Wno-format
		  -Wno-pointer-bool-conversion
		  -Wno-comment
		  -Wno-incompatible-pointer-types
		" |tr '\n' ' '|tr -s ' ')
		echo "Adding appropriate compiler flags to $i"
		sed "s/CFLAGS =/CFLAGS = ${CFLAGS}/" "$i" > "$i".tmp
		mv -fv "$i".tmp "$i"
	  fi

done

#need to tweak some code
for i in gps1x2rnx.c Bin2AsciiLevel1.c
do
  i="$(find "$SW_DIR" -name "$i")"
  if [ ! -z "$(grep 'int8_t \*argv\[\]' "$i")" ]
  then
    sed 's/int8_t \*argv\[\]/char *argv[]/' "$i" > "$i".tmp
    mv -fv "$i".tmp "$i"
  fi
done

#update the location of the header files
i="$(find "$SW_DIR" -name GRACEsyspath.h)"
if [ ! -z "$(grep '/goa/local/grace/includes' "$i")" ]
then
  sed "s:/goa/local/grace/includes:${SW_DIR}:" "$i" > "$i".tmp
  mv -fv "$i".tmp "$i"
fi

#compile the software
make -C "$SW_DIR"

#create links to this dir
#ln -svf "$(basename $SW_DIR)"/*.e .
