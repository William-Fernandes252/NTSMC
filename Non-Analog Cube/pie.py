import matplotlib.pyplot as plt
import numpy as np
import csv

def computePercentage(pct, data, err):
    abs = int(pct/100.*np.sum(data))
    return "{:.5f}%\n({:.5e})".format(pct, abs)

# Create the figure that will show the graph:
fig, graph = plt.subplots(figsize = (8, 6), subplot_kw = dict(aspect = "equal"))

# Read the data from the report file:
with open('report.csv', 'r') as report:
    reader = csv.reader(report, delimiter = ';', quoting=csv.QUOTE_NONE)
    rows = [line for line in reader]
quantities = [rows[0][0], rows[0][1], rows[0][2], rows[0][3], rows[0][4]]
mean = np.array([rows[1][0], rows[1][1], rows[1][2], rows[1][3], rows[1][4]]).astype(np.float64)
std = np.array([rows[2][0], rows[2][1], rows[2][2], rows[2][3], rows[2][4]]).astype(np.float64)

# Create the graph:
wedges, texts, autotexts = graph.pie(mean, autopct = lambda pct: computePercentage(pct, mean, std), 
                                     textprops = dict(color = "k"), radius = 1.1)

# Defining the subtittles box:
graph.legend(wedges, quantities, title = "Eventos", loc = "center left", bbox_to_anchor = (1, 0, 0.5, 1))

# Setting the font and the size of the text inside the graph:
plt.setp(autotexts, size = 5, weight = "bold")

# Tittle of the graph:
graph.set_title("Average number of occurences of each type of event:")

plt.show()
