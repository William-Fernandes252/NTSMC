import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

def create_graph(option, x, x_label, y, y_label, err, color):

    # Create the figure that will show the graph:
    fig, graph = plt.subplots()

    # Create the graph, set logarithm scale depending on the option parameter:
    if option == 'n':
        graph.errorbar(x, y, yerr = err, color = color, ecolor = 'black', 
                       capsize = 2)
    elif option == 's':
        graph.bar(x, y, yerr = err, align = 'center', color = color, 
                  ecolor = 'black', capsize = 1, log = True)

    # Set the labels:
    axis = plt.gca()
    axis.set_xlabel(x_label)
    axis.set_ylabel(y_label)
    axis.set_xticks(np.arange(0, 101, 5))
    graph.grid(True)

    # Tittle of the graph:
    graph.set_title("Estimatives of the " + y_label.lower() + " per " 
                    + x_label.lower())

    plt.tight_layout()
    plt.show()

def set_variables(results):
    if 'Layer' in results:
        x = results['Layer']
        x_label = "Layer"
        if 'Collisions mean' in results:
            y = results['Collisions mean']
            y_label = "Collisions"
            color = 'green'
        elif 'Escapes mean' in results:
            y = results['Escapes mean']
            y_label = "Escapes"
            color = 'blue'
    elif 'Radius' in results:
        x = results['Radius']
        x_label = "Radius (cm)"
        y = results['Total flux mean']
        y_label = "Total flux"
        color = 'orange'
    err = results['Std. deviation']
    return x, x_label, y, y_label, err, color

# Store the data in a dataframe:
results = pd.read_csv("report.csv", sep = ';')

# Set the variables for analysis:
x, x_label, y, y_label, err, color = set_variables(results)

# Ask for the type of analysis desired:
while True:
    option = input("Order of magnitude analysis? (s/n)\n")
    if(option == 's'):
        break
    elif(option == 'n'):
        break

# Created and plot the graph:
create_graph(option, x, x_label, y, y_label, err, color)
