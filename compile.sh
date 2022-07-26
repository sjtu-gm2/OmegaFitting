#!/bin/bash

ALL="NO"

TARGET=./OmegaFitting/main.cpp
while [[ $# -gt 0 ]]; do
  case $1 in
    -t|--target)
      TARGET="$2"
      shift # past argument
      shift # past value
      ;;
    -s|--searchpath)
      SEARCHPATH="$2"
      shift # past argument
      shift # past value
      ;;
    --all)
      ALL="YES"
      shift # past argument
      ;;    
    -*|--*)
      echo "Unknown option $1"
      exit 1
      ;;
    *)
      POSITIONAL_ARGS+=("$1") # save positional arg
      shift # past argument
      ;;
  esac
done

BASENAME=`basename ${TARGET}`
EXE=${BASENAME%.*}.exe

[ -d objs ] || mkdir objs
[ -d output ] || mkdir output

if [[ ${ALL} == 'YES' ]];then
  g++ -g -Wall -Wextra -O3 -funroll-loops -finline-functions -fomit-frame-pointer -fpic -I./OmegaFitting/rlib/include -c \
    -o ./objs/Random.o ./OmegaFitting/rlib/src/Random.cpp &&\
  g++ -I OmegaFitting/rlib/include -I /usr/local/opt/openssl/include OmegaFitting/Blinders.cc -std=c++11 -Wall -Wextra -Werror -pedantic-errors -fpic -c \
    -o ./objs/Blinders.o &&\
  g++ -shared -L/usr/local/opt/openssl/lib -lssl -lcrypto \
    -o ./objs/libBlinders.so ./objs/Blinders.o ./objs/Random.o
fi

g++ `root-config --cflags --glibs` -I ./OmegaFitting/rlib/include -L./objs -lBlinders -I ./OmegaFitting -I/usr/local/include -L/usr/local/lib -ljsoncpp \
  -o ${EXE} ./OmegaFitting/FitTools.cpp ${TARGET} \
&& echo -e "Generated ${EXE} from ${TARGET}"