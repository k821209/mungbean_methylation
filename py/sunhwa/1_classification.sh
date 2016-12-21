cat ${1} | grep 'CG' > ${1}.CG
python3 1_window_count.py ${1}.CG 
cat ${1} | grep 'CHG' > ${1}.CHG
python3 1_window_count.py ${1}.CHG
cat ${1} | grep 'CHH' > ${1}.CHH
python3 1_window_count.py ${1}.CHH
