#!/bin/sh

gcc -cpp -fPIC -shared calc_adv.c -lm -o calc_adv.so -O3
gcc -cpp -fPIC -shared calc_charge.c -lm -o calc_charge.so -O3
gcc -cpp -fPIC -shared calc_gauss.c -lm -o calc_gauss.so -O3
gcc -cpp -fPIC -shared flash.c -lm -o flash.so -O3
gcc -cpp -fPIC -shared lookup.c -lm -o lookup.so -O3

time=`date +%s`
mv nohup.out $time.log