head -n 2 config.cfg>tmp.cfg
source tmp.cfg
arr=$(echo $FILENAMES | tr "," "\n")
echo ${arr[*]}
command="./secondtrie "
for f in ${arr[*]}
     do
	if [ "$ZIPPED" = "True" ]
        then
	    command="$command <(gunzip -c $f) "
           
    	else

		command="$command <(cat $f) "
	fi
    done

eval $command   

rm tmp.cfg

