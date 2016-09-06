head -n 2 configforveryshort5.cfg>tmp.cfg
source tmp.cfg
arr=$(echo $FILENAMES | tr "," "\n")
echo ${arr[*]}


if [ "$ZIPPED" = "True" ]
	then
		gunzip -c ${arr[*]} | ./secondtrie||echo "File not found. Check input file parameters."; exit

else
	for f in ${arr[*]}
    do
        echo "Processing $f"
        sed -i '' -e '$a\' $f
    done
    cat ${arr[*]} | ./secondtrie
fi
rm tmp.cfg

