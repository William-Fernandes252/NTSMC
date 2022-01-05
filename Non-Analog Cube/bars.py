import matplotlib.pyplot as plt
import numpy as np
import csv

# Create the figure that will show the graph:
fig, graph = plt.subplots(figsize = (7, 5))

# Read the data from the report file:
with open('report.csv', 'r') as report:
    reader = csv.reader(report, delimiter = ';', quoting=csv.QUOTE_NONE)
    rows = [line for line in reader]
quantities = [rows[0][0], rows[0][1], rows[0][2], rows[0][3], rows[0][4]]
mean = np.array([rows[1][0], rows[1][1], rows[1][2], rows[1][3], rows[1][4]]).astype(np.float64)
std = np.array([rows[2][0], rows[2][1], rows[2][2], rows[2][3], rows[2][4]]).astype(np.float64)
pos = np.arange(len(quantities))

# Create the graph:
graph.bar(pos, mean, yerr = std, orientation = 'vertical',
          align = 'center', alpha = 0.5,
          color = ['orange', 'blue', 'green', 'yellow', 'red'], 
          ecolor = 'black', capsize = 10, log = True)

# Set the labels:
graph.set_ylabel("Occurences")
graph.set_xticks(pos)
graph.set_xticklabels(quantities)
graph.yaxis.grid(True)

# Tittle of the graph:
graph.set_title("Estimatives of the average number of occurences of each type of event")

# Show the graph:
plt.show()
