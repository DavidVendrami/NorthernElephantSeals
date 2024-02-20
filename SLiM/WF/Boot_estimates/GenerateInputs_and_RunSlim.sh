#!/bin/bash

cat Inputs.txt | while read i;
do
XXXXX=$(echo $i | cut -d ' ' -f1)
TOT_LEN=$(echo $i | cut -d ' ' -f2)
YYYYY=$(echo $i | cut -d ' ' -f3)
LGM=$(echo $i | cut -d ' ' -f4)
WWWWW=$(echo $i | cut -d ' ' -f5)
ZZZZZ=$(echo $i | cut -d ' ' -f6)
BOTT=$(echo $i | cut -d ' ' -f7)
AAAAA=$(echo $i | cut -d ' ' -f8)
BBB=$(echo $i | cut -d ' ' -f9)
CCC=$(echo $i | cut -d ' ' -f10)
DDD=$(echo $i | cut -d ' ' -f11)
EEE=$(echo $i | cut -d ' ' -f12)
FFF=$(echo $i | cut -d ' ' -f13)
GGG=$(echo $i | cut -d ' ' -f14)
HHH=$(echo $i | cut -d ' ' -f15)
III=$(echo $i | cut -d ' ' -f16)
JJJ=$(echo $i | cut -d ' ' -f17)
KKK=$(echo $i | cut -d ' ' -f18)
LLL=$(echo $i | cut -d ' ' -f19)
MMM=$(echo $i | cut -d ' ' -f20)
NNN=$(echo $i | cut -d ' ' -f21)
OOO=$(echo $i | cut -d ' ' -f22)
PPP=$(echo $i | cut -d ' ' -f23)
QQQ=$(echo $i | cut -d ' ' -f24)
RRR=$(echo $i | cut -d ' ' -f25)
SSS=$(echo $i | cut -d ' ' -f26)
TTT=$(echo $i | cut -d ' ' -f27)
UUU=$(echo $i | cut -d ' ' -f28)
VVV=$(echo $i | cut -d ' ' -f29)
POSTG=$(echo $i | cut -d ' ' -f30)
GRBOT=$(echo $i | cut -d ' ' -f31)
N=$(echo $i | cut -d ' ' -f32)
sed "s/XXXXX/${XXXXX}/g" WF_delet_along_whole_genome_F_repN.txt > st1.txt
sed -i "s/TOT_LEN/${TOT_LEN}/g" st1.txt
sed -i "s/YYYYY/${YYYYY}/g" st1.txt
sed -i "s/LGM/${LGM}/g" st1.txt
sed -i "s/WWWWW/${WWWWW}/g" st1.txt
sed -i "s/ZZZZZ/${ZZZZZ}/g" st1.txt
sed -i "s/BOTT/${BOTT}/g" st1.txt
sed -i "s/AAAAA/${AAAAA}/g" st1.txt
sed -i "s/BBB/${BBB}/g" st1.txt
sed -i "s/CCC/${CCC}/g" st1.txt
sed -i "s/DDD/${DDD}/g" st1.txt
sed -i "s/EEE/${EEE}/g" st1.txt
sed -i "s/FFF/${FFF}/g" st1.txt
sed -i "s/GGG/${GGG}/g" st1.txt
sed -i "s/HHH/${HHH}/g" st1.txt
sed -i "s/III/${III}/g" st1.txt
sed -i "s/JJJ/${JJJ}/g" st1.txt
sed -i "s/KKK/${KKK}/g" st1.txt
sed -i "s/LLL/${LLL}/g" st1.txt
sed -i "s/MMM/${MMM}/g" st1.txt
sed -i "s/NNN/${NNN}/g" st1.txt
sed -i "s/OOO/${OOO}/g" st1.txt
sed -i "s/PPP/${PPP}/g" st1.txt
sed -i "s/QQQ/${QQQ}/g" st1.txt
sed -i "s/RRR/${RRR}/g" st1.txt
sed -i "s/SSS/${SSS}/g" st1.txt
sed -i "s/TTT/${TTT}/g" st1.txt
sed -i "s/UUU/${UUU}/g" st1.txt
sed -i "s/VVV/${VVV}/g" st1.txt
sed -i "s/POSTG/${POSTG}/g" st1.txt
sed "s/GRBOT/${GRBOT}/g" st1.txt > WF_delet_along_whole_genome_F_rep$N.txt
slim WF_delet_along_whole_genome_F_rep$N.txt > rep_N/out.txt
cp -r rep_N ./rep_$N
done





