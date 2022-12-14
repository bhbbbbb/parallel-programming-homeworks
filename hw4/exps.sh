#/src/bin/bash
for i in {1..32}
do
    echo i = $i
    $1 $i >> $2
done
