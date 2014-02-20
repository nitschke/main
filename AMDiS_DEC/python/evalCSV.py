import csv
with open('../meshStatsSphereDivBy4.csv', 'rb') as f:
    reader = csv.reader(f)
    n = 0
    for col in reader[0]:
        if (col == "MaxMaxAngle"):
          break
        n = n + 1

print n
