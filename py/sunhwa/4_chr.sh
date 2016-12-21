for i in Vr01 Vr02 Vr03 Vr04 Vr05 Vr06 Vr07 Vr08 Vr09 Vr10 Vr11 scaff
do
cat ${1} | grep ${i} > ${1}.${i}.chr
done
