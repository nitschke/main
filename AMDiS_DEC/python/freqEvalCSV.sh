while true
  do
    eval "python evalCSV.py &"
    sleep $1
    kill $!
  done
