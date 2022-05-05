[ -d objs ] || mkdir objs
[ -d output ] || mkdir output
if [[ ${1} == 'All' ]];then
	g++ -g -Wall -Wextra -O3 -funroll-loops -finline-functions -fomit-frame-pointer -fpic -I./srcs/rlib/include -c -o ./objs/Random.o ./srcs/rlib/src/Random.cpp &&\
	g++ -I srcs/rlib/include -I /usr/local/opt/openssl/include srcs/Blinders.cc -std=c++11 -Wall -Wextra -Werror -pedantic-errors -fpic -c -o ./objs/Blinders.o &&\
	g++ -shared -o ./objs/libBlinders.so ./objs/Blinders.o ./objs/Random.o -L/usr/local/opt/openssl/lib -lssl -lcrypto
fi

g++ `root-config --cflags --glibs` -I ./srcs//rlib/include -L./objs -lBlinders -I ./srcs -o main.exe ./srcs/FitTools.cpp ./srcs/main.cpp