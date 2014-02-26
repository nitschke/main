while true
  do
    eval "python evalCSV.py &"
    pid=$!
    while ps -p $pid > /dev/null
    do
      sleep 1
    done
  done
